from enum import IntEnum
from abc import ABC, abstractmethod
from mewpy.utils.constants import EAConstants, ModelConstants
from mewpy.optimization.ea import Solution
from mewpy.simulation import get_simulator
from collections import OrderedDict
import warnings
import numpy as np
import copy


class Strategy(IntEnum):
    """
    The available optimization strategies
    """
    KO = 1
    OU = 2


class KOBounder(object):
    """
    A bounder of possible indexes in an enumeration
    """

    def __init__(self, lower_bound, upper_bound):

        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.range = self.upper_bound - self.lower_bound + 1

    # make it callable
    def __call__(self, candidate, args):
        bounded_candidate = set()
        for val in candidate:
            if val > self.upper_bound or val < self.lower_bound:
                val = val % self.range + self.lower_bound
            bounded_candidate.add(val)
        return bounded_candidate


class OUBounder(object):
    """
    A bounder for (int,int) representations
    """

    def __init__(self, lower_bound, upper_bound):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.range = [self.upper_bound[i]-self.lower_bound[i] +
                      1 for i in range(len(self.lower_bound))]

    def __call__(self, candidate, args):
        bounded_candidate = set()
        for idx, lv in candidate:
            if idx > self.upper_bound[0] or idx < self.lower_bound[0]:
                idx = idx % self.range[0] + self.lower_bound[0]
            if lv > self.upper_bound[1] or idx < self.lower_bound[1]:
                lv = lv % self.range[1] + self.lower_bound[1]
            bounded_candidate.add((idx, lv))
        return bounded_candidate


class AbstractProblem(ABC):
    """
    Base class for optimization problems.

    :param model: A Metabolic model.
    :param list fevaluation: A list of callable EvaluationFunctions.
    :param **kwargs: Additional parameters dictionary
    """

    def __init__(self, model, fevaluation=None,**kwargs):
        self.model = model
        self.fevaluation = fevaluation
        self.number_of_objectives = len(self.fevaluation)
        
        # simulation context : defines the simulations environment
        self.simul_context = None
        # The target product reaction id may be specified when optimizing for a single product.
        # Only required for probabilistic modification targeting.
        self.product = kwargs.get('product', None)
        # Environmental conditions
        self.environmental_conditions = kwargs.get('envcond', None)
        # Additional persistent constraints
        self.persistent_constraints = kwargs.get('constraints', None)
        # Reference reaction fluxes
        self._reference = None
        # solution size
        self.candidate_min_size = kwargs.get(
            'candidate_min_size', EAConstants.MIN_SOLUTION_SIZE)
        self.candidate_max_size = kwargs.get(
            'candidate_max_size', EAConstants.MAX_SOLUTION_SIZE)
        # non target
        self.non_target = kwargs.get('non_target', None)
        # targets
        # If not provided, targets are build in the context of the problem.
        # Objectives are not automatically removed from the targets... this should be a user concern!?
        self._trg_list = kwargs.get('target', None)
        # the EA representation bounds
        self._bounder = None
        # scaling factor
        self.scalefactor = kwargs.get('scalefactor', None)
        # required simulations
        methods = []
        for f in self.fevaluation:
            methods.extend(f.required_simulations())
        self.methods = list(set(methods))

    @abstractmethod
    def generator(self, random, args):
        """The generator function for the problem."""
        raise NotImplementedError

    @abstractmethod
    def translate(self, candidate, reverse=False):
        """The generator function for the problem."""
        raise NotImplementedError

    @abstractmethod
    def decode(self, candidate):
        """The decoder function for the problem."""
        raise NotImplementedError

    def get_name(self):
        """The problem name."""
        return self.__class__.__name__

    def pre_process(self):
        """ Defines pre processing tasks
        """
        self.target_list
        self.reset_simulator()

    @property
    def simulator(self):
        if self.simul_context is None:
            self.simul_context = get_simulator(
                self.model, reference=self._reference)
            self.simul_context.environmental_conditions = self.environmental_conditions
            self.simul_context.constraints = self.persistent_constraints
        return self.simul_context

    def reset_simulator(self):
        self.simul_context = None

    def __str__(self):
        if self.number_of_objectives > 1:
            return '{0} ({1}  objectives)'.format(self.__class__.__name__, self.number_of_objectives)
        else:
            return '{0} '.format(self.__class__.__name__)

    def __repr__(self):
        return self.__class__.__name__

    @property
    def target_list(self):
        "list of modification targets"
        if self._trg_list is None:
            self._build_target_list()
            if not self._trg_list:
                raise RuntimeError("Could not build target list.")
        return self._trg_list

    @abstractmethod
    def _build_target_list(self):
        raise NotImplementedError

    def get_constraints(self, solution):
        """
        :returns: The constrainst enconded into an individual.
        """
        return self.decode(solution.candidate)

    def get_environmental_conditions(self):
        return self.environmental_conditions

    def get_persistent_constraints(self):
        return self.persistent_constraints

    def evaluate_solution(self, solution, decode=True):
        """
        Evaluates a single solution, a list of constraints.

        :param solution: The solution to be evaluated.
        :param decode: If the solution needs to be decoded.
        :returns: A list of fitness.
        """
        p = []
        # decoded constraints
        if decode:
            constraints = self.decode(solution)
        else:
            constraints = solution
        # pre simulation
        simulation_results = OrderedDict()
        try:
            for method in self.methods:
                simulation_result = self.simulator.simulate(
                    constraints=constraints, method=method, scalefactor=self.scalefactor)
                simulation_results[method] = simulation_result
            # apply the evaluation function(s)
            for f in self.fevaluation:
                p.append(f(simulation_results, solution,
                           scalefactor=self.scalefactor))
        except Exception as e:
            for f in self.fevaluation:
                p.append(f.worst_fitness)
            if EAConstants.DEBUG:
                warnings.warn(f"Solution couldn't be evaluated [{e}]\n {constraints}")
        return p

    @property
    def is_maximization(self):
        return all([f.maximize for f in self.fevaluation])

    def simplify(self, solution, tolerance=1e-6,):
        """
        Simplify a solution by removing the modification that do not affect the final fitness value.
        Two solutions are considered different if the maximum allowed difference between objective values is exceeded.
        Tolerance may be defined by a single float value, or per objective by setting a list of floats of size equal
        to the number of objectives.

        :param solution: the solution to be simplified
        :param float tolerance: The maximum allowed difference between objective values.
        :returns: A list of simplified solutions.

        """

        values = self.translate(solution.values, reverse=True)
        fitness = self.evaluate_solution(values)
        one_to_remove = {}
        # single removal
        for entry in values:
            simul_constraints = copy.copy(values)
            simul_constraints.remove(entry)
            fit = self.evaluate_solution(simul_constraints)
            diff = np.abs(np.array(fit)-np.array(fitness))
            is_equal = False
            if isinstance(tolerance, float):
                is_equal = np.all(diff <= tolerance)
            else:
                is_equal = np.all(diff <= np.array(tolerance))
            if is_equal:
                one_to_remove[entry] = fit

        simul_constraints = copy.copy(values)
        for entry in one_to_remove.keys():
            simul_constraints.remove(entry)

        # test all simultaneous removal
        fit = self.evaluate_solution(simul_constraints)
        diff = np.abs(np.array(fit)-np.array(fitness))
        is_equal = False
        if isinstance(tolerance, float):
            is_equal = np.all(diff <= tolerance)
        else:
            is_equal = np.all(diff <= np.array(tolerance))

        if is_equal:
            v = self.translate(simul_constraints)
            c = self.decode(simul_constraints)
            simplification = Solution(v, fitness, c)
            return [simplification]
        else:
            res = []
            for entry, fit in one_to_remove.items():
                simul_constraints = copy.copy(values)
                simul_constraints.remove(entry)
                v = self.translate(simul_constraints)
                c = self.decode(simul_constraints)
                simplification = Solution(v, fitness, c)
                res.append(simplification)
            return res


class AbstractKOProblem(AbstractProblem):
    """
    Base class for Knockout optimization problems.

    :param model: The constraint metabolic model.
    :param list fevaluation: A list of callable EvaluationFunctions.

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(AbstractKOProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        self.strategy = Strategy.KO

    @abstractmethod
    def decode(self, candidate):
        """The decode function for the problem."""
        raise NotImplementedError

    @property
    def bounder(self):
        """
        The KO list index bounder
        """
        if self._bounder is None:
            max = len(self.target_list)-1
            self._bounder = KOBounder(0, max)
        return self._bounder

    def generator(self, random, args):
        """
        Generates a solution, a random int set with length in range min_solution_size to max_solution_size.
        """
        solution = set()
        solution_size = random.uniform(
            self.candidate_min_size, self.candidate_max_size)
        while len(solution) < solution_size:
            solution.add(random.randint(0, len(self.target_list)-1))
        return solution

    def translate(self, candidate, reverse=False):
        """
        Translates a candidate solution in problem specific representation to
        an iterable of ids, or (ids, folds).

        :param candidate: The candidate representation.
        :param boolean reverse: Performs a reverse translation.

        """
        if not reverse:
            return {self.target_list[idx] for idx in candidate}
        else:
            return [self.target_list.index(k) for k in candidate]


    def solution_to_constraints(self,solution):
        """"Transforms a solution for the problem into metabolic constraints.

        :param dict solution: A dictionary of genetic modifications.
        :returns: A dictionary of metabolic constraints that may be applied to the model.
        """
        return self.decode(self.translate(solution,True))

class AbstractOUProblem(AbstractProblem):
    """ Base class for Over/Under expression optimization problems
    """

    def __init__(self, model, fevaluation=None, **kwargs):
        """
        :param model: The constraint metabolic model.
        :param list fevaluation: A list of callable EvaluationFunctions.

        Optional:

        :param dic reference: Dictionary of flux values to be used in the over/under expression values computation.
        :param list levels: Over/under expression levels (Default EAConstants.LEVELS).

        """
        super(AbstractOUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        self.strategy = Strategy.OU
        self.levels = kwargs.get('levels', EAConstants.LEVELS)
        self._reference = kwargs.get('reference', None)

    @abstractmethod
    def decode(self, candidate):
        """The decoder function for the problem. Needs to be implemented by extending classes."""
        raise NotImplementedError

    @property
    def bounder(self):
        """
        The list and levels index bounder

        :returns: a OUBounder object.
        """
        if self._bounder is None:
            max_idx = len(self.target_list)-1
            max_lv = len(self.levels)-1
            self._bounder = OUBounder([0, 0], [max_idx, max_lv])
        return self._bounder

    def generator(self, random, args):
        """
        Generates a solution, a random (int,int) set with length in range min_solution_size to max_solution_size.

        :param random: A random number generator. Needs to implement uniform and randint generators methods.
        :param dict args: A dictionary of additional parameters.
        :returns: A new solution.
        """
        solution = set()
        solution_size = random.uniform(
            self.candidate_min_size, self.candidate_max_size)
        while len(solution) < solution_size:
            idx = random.randint(0, len(self.target_list)-1)
            lv = random.randint(0, len(self.levels)-1)
            solution.add((idx, lv))
        return solution

    @property
    def reference(self):
        """
        :returns: A dictionary of reference flux values.
        """
        if not self._reference:
            self._reference = self.simulator.reference
        return self._reference

    def ou_constraint(self, level, wt):
        """ Computes the bounds for a reaction.

        :param float level: The expression level for the reaction.
        :param float wt: The reference reaction flux.
        :returns: A tupple, flux bounds for the reaction.

        """
        if level > 1:
            if wt >= 0:
                return (level*wt, ModelConstants.REACTION_UPPER_BOUND)
            else:
                return (-1*ModelConstants.REACTION_UPPER_BOUND, level*wt)
        else:
            return (0, level * wt) if wt >= 0 else (level * wt, 0)

    def reaction_constraints(self, rxn, lv):
        """
        Converts a (reaction, level) pair into a constraint
        If a reaction is reversible, the direction with no or less wild type flux
        is knocked out.

        :param rxn: The reaction identifier.
        :param lv: the wild type multiplier factor. The reaction bounds are altered accordingly.
        :returns: A dictionary of reaction constraints.
        """
        constraints = {}
        fluxe_wt = self.reference[rxn]
        rev_rxn = self.simulator.reverse_reaction(rxn)
        if lv == 0:
            # KO constraint
            constraints[rxn] = (0, 0)
        elif lv == 1:
            # No contraint is applyed
            pass
        elif rev_rxn is None or rev_rxn == rxn:
            # if there is no reverse reaction
            constraints[rxn] = self.ou_constraint(lv, fluxe_wt)
        else:
            # there's a reverse reaction...
            # one of the two reactions needs to be KO, the one with no flux in the wt.
            rev_fluxe_wt = self.reference[rev_rxn]
            if abs(fluxe_wt) >= abs(rev_fluxe_wt):
                ko_rxn, ou_rxn, fwt = rev_rxn, rxn, fluxe_wt
            else:
                rxn, rev_rxn, rev_fluxe_wt
            constraints[ko_rxn] = (0, 0)
            constraints[ou_rxn] = self.ou_constraint(lv, fwt)
        return constraints

    def translate(self, candidate, reverse=False):
        """
        Translates a candidate solution in problem specific representation to
        an iterable of ids, or (ids, folds).

        :param iterable candidate: The candidate representation.
        :param boolean reverse: Performs the reverse translation.
        :returns: If reverse, a list of index tupple (modification_target_index,level_index). The indexes are
                  problem dependent.
                  If not reverse, a dictionary of {modification_target: level}
        """
        if not reverse:
            return {self.target_list[idx]: self.levels[lv_idx]
                    for idx, lv_idx in candidate}
        else:
            return [(self.target_list.index(k), self.levels.index(lv))
                    for k, lv in candidate.items()]

    def solution_to_constraints(self,solution):
        """"Transforms a solution for the problem into metabolic constraints.

        :param dict solution: A dictionary of genetic modifications.
        :returns: A dictionary of metabolic constraints that may be applied to the model.
        """
        return self.decode(self.translate(solution,True))
