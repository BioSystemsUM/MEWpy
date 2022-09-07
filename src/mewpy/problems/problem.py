import copy
import warnings
from abc import ABC, abstractmethod
from enum import Enum
import numpy as np
from ..optimization.ea import Solution, filter_duplicates
from ..simulation import get_simulator
from ..util.constants import EAConstants, ModelConstants
from ..util.process import get_fevaluator


class Strategy(Enum):
    KO = 'KO'
    OU = 'OU'

    def __eq__(self, other):
        """Overrides equal to enable string name comparison.
        Allows to seamlessly use:
            Strategy.KO = Strategy.KO
            Strategy.KO = 'KO'.
        """
        if isinstance(other, Strategy):
            return super().__eq__(other)
        elif isinstance(other, str):
            return self.name == other
        else:
            return False

    def __hash__(self):
        return hash(self.name)


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
        self.range = [self.upper_bound[i] - self.lower_bound[i] +
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
    :param kwargs: Additional parameters dictionary
    """

    def __init__(self, model, fevaluation=None, **kwargs):
        self.model = model
        self.fevaluation = [] if fevaluation is None else fevaluation
        self.number_of_objectives = len(self.fevaluation)

        # simulation context : defines the simulations environment
        self._reset_solver = kwargs.get('reset_solver', ModelConstants.RESET_SOLVER)
        self._simul = None
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
        # problem specific EA operators
        self.operators = None
        self.operators_param = None

    @abstractmethod
    def generator(self, random, **kwargs):
        """The generator function for the problem."""
        raise NotImplementedError

    @abstractmethod
    def encode(self, candidate):
        """The generator function for the problem."""
        raise NotImplementedError

    @abstractmethod
    def decode(self, candidate):
        """The decoder function for the problem."""
        raise NotImplementedError

    @abstractmethod
    def solution_to_constraints(self, solution):
        """Converts a decoded solution to metabolict constraints."""
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
        if self._simul is None:
            self._simul = get_simulator(
                self.model, envcond=self.environmental_conditions,
                constraints=self.persistent_constraints,
                reference=self._reference, reset_solver=self._reset_solver)
        return self._simul

    def simulate(self, *args, **kwargs):
        return self._simul.simulate(*args, **kwargs)

    def reset_simulator(self):
        self._simul = None

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
        return self.solution_to_constraints(self.decode(solution.candidate))

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
        decoded = {}
        # decoded constraints
        if decode:
            decoded = self.decode(solution)
            constraints = self.solution_to_constraints(decoded)
        else:
            constraints = solution

        # pre simulation
        simulation_results = dict()
        try:
            p = []
            for method in self.methods:
                simulation_result = self.simulator.simulate(
                    constraints=constraints, method=method, scalefactor=self.scalefactor)
                simulation_results[method] = simulation_result
            # apply the evaluation function(s)
            for f in self.fevaluation:
                v = f(simulation_results, decoded, scalefactor=self.scalefactor)
                p.append(v)
        except Exception as e:
            p = []
            for f in self.fevaluation:
                p.append(f.worst_fitness)
            if EAConstants.DEBUG:
                warnings.warn(f"Solution couldn't be evaluated [{e}]\n {constraints}")
        return p

    @property
    def is_maximization(self):
        return all([f.maximize for f in self.fevaluation])

    def simplify(self, solution, tolerance=1e-6, ):
        """
        Simplify a solution by removing the modification that do not affect the final fitness value.
        Two solutions are considered different if the maximum allowed difference between objective values is exceeded.
        Tolerance may be defined by a single float value, or per objective by setting a list of floats of size equal
        to the number of objectives.

        :param solution: the solution to be simplified
        :param float tolerance: The maximum allowed difference between objective values.
        :returns: A list of simplified solutions.

        """

        enc_values = self.encode(solution.values)
        fitness = self.evaluate_solution(enc_values)
        one_to_remove = {}
        # single removal
        for entry in enc_values:
            simul_enc_values = copy.copy(enc_values)
            simul_enc_values.remove(entry)
            fit = self.evaluate_solution(simul_enc_values)
            diff = np.abs(np.array(fit) - np.array(fitness))
            is_equal = False
            if isinstance(tolerance, float):
                is_equal = np.all(diff <= tolerance)
            else:
                is_equal = np.all(diff <= np.array(tolerance))
            if is_equal:
                one_to_remove[entry] = fit

        simul_enc_values = copy.copy(enc_values)
        for entry in one_to_remove.keys():
            simul_enc_values.remove(entry)

        # test all simultaneous removal
        fit = self.evaluate_solution(simul_enc_values)
        diff = np.abs(np.array(fit) - np.array(fitness))
        is_equal = False
        if isinstance(tolerance, float):
            is_equal = np.all(diff <= tolerance)
        else:
            is_equal = np.all(diff <= np.array(tolerance))

        if is_equal:
            v = self.decode(simul_enc_values)
            c = self.solution_to_constraints(v)
            simplification = Solution(v, fitness, c)
            return [simplification]
        else:
            res = []
            for entry, fit in one_to_remove.items():
                simul_enc_values = copy.copy(enc_values)
                simul_enc_values.remove(entry)
                v = self.decode(simul_enc_values)
                c = self.solution_to_constraints(v)
                simplification = Solution(v, fitness, c)
                res.append(simplification)
            return res

    def simplify_population(self, population, n_cpu=1):
        """Simplifies a population of solutions

        Args:
            population (list): List of mewpy.optimization.ea.Solution
            n_cpu (int): Number of CPUs.
        Returns:
            list: Simplified population
        """
        pop = []
        for solution in population:
            try:
                res = self.simplify(solution)
                pop.extend(res)
            except Exception:
                pop.append(solution)
        return filter_duplicates(pop)


class AbstractKOProblem(AbstractProblem):

    def __init__(self, model, fevaluation=None, **kwargs):
        """
        Base class for Knockout optimization problems.

        :param model: The constraint metabolic model.
        :param list fevaluation: A list of callable EvaluationFunctions.

        """
        super(AbstractKOProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        self.strategy = Strategy.KO

    def decode(self, candidate):
        decoded = dict()
        for idx in candidate:
            try:
                decoded[self.target_list[idx]] = 0
            except IndexError:
                raise IndexError("Index out of range: {} from {}".format(
                    idx, len(self.target_list[idx])))
        return decoded

    def encode(self, candidate):
        """
        Translates a candidate solution in problem specific representation to
        an iterable of ids, or (ids, folds).

        :param candidate: The candidate representation.
        """
        return set([self.target_list.index(k) for k in candidate])

    def solution_to_constraints(self, decoded_candidate):
        """
        Converts a candidate, a dictionary of reactions, into a dictionary of constraints.
        This is problem specific. By default return the decoded candidate.
        """
        return decoded_candidate

    @property
    def bounder(self):
        """
        The KO list index bounder
        """
        if self._bounder is None:
            max = len(self.target_list) - 1
            self._bounder = KOBounder(0, max)
        return self._bounder

    def generator(self, random, **kwargs):
        """
        Generates a solution, a random int set with length in range min_solution_size to max_solution_size.
        """

        if self.candidate_min_size == self.candidate_max_size:
            solution_size = self.candidate_min_size
        else:
            solution_size = random.randint(
                self.candidate_min_size, self.candidate_max_size)

        solution = set(random.sample(range(len(self.target_list)), solution_size))
        return solution


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
        :param boolean twostep: If deletions should be applied before identifiying reference flux values.
        :param dict partial_solution: A partial solution to be appended to any other solution
        """
        super(AbstractOUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        self.strategy = Strategy.OU
        self.levels = kwargs.get('levels', EAConstants.LEVELS)
        self._reference = kwargs.get('reference', None)
        self.twostep = kwargs.get('twostep', False)
        self._partial_solution = kwargs.get('partial_solution', dict())

    def decode(self, candidate):
        """The decoder function for the problem. Needs to be implemented by extending classes."""
        decoded = self._partial_solution.copy()
        for idx, lv_idx in candidate:
            try:
                rxn = self.target_list[idx]
                lv = self.levels[lv_idx]
                decoded[rxn] = lv
            except IndexError:
                raise IndexError("Index out of range")
        return decoded

    def encode(self, candidate):
        """
        Translates a candidate solution in problem specific representation to
        an iterable of ids, or (ids, folds).

        :param iterable candidate: The candidate representation.
        :returns: a list of index tupple (modification_target_index,level_index). The indexes are
                  problem dependent.
        """

        return set([(self.target_list.index(k), self.levels.index(lv))
                    for k, lv in candidate.items()])

    def solution_to_constraints(self, decoded_candidate):
        """
        Decodes a candidate, a dictionary of reactions, into a dictionary of constraints.
        This is problem specific. By default return the decoded candidate.
        """
        return decoded_candidate

    @property
    def bounder(self):
        """
        The list and levels index bounder

        :returns: a OUBounder object.
        """
        if self._bounder is None:
            max_idx = len(self.target_list) - 1
            max_lv = len(self.levels) - 1
            self._bounder = OUBounder([0, 0], [max_idx, max_lv])
        return self._bounder

    def generator(self, random, **kwargs):
        """
        Generates a solution, a random (int,int) set with length in range min_solution_size to max_solution_size.

        :param random: A random number generator. Needs to implement uniform and randint generators methods.
        :param dict args: A dictionary of additional parameters.
        :returns: A new solution.
        """

        if self.candidate_min_size == self.candidate_max_size:
            solution_size = self.candidate_min_size
        else:
            solution_size = random.randint(
                self.candidate_min_size, self.candidate_max_size)

        keys = random.sample(range(len(self.target_list)), solution_size)
        values = random.choices(range(len(self.levels)), k=solution_size)
        solution = set(zip(keys, values))
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
            if wt > 0:
                return (level * wt, ModelConstants.REACTION_UPPER_BOUND)
            elif wt < 0:
                return (-1 * ModelConstants.REACTION_UPPER_BOUND, level * wt)
        else:
            if wt > 0:
                return (0, level * wt)
            elif wt < 0:
                return (level * wt, 0)

    def reaction_constraints(self, rxn, lv, reference):
        """
        Converts a (reaction, level) pair into a constraint
        If a reaction is reversible, the direction with no or less wild type flux
        is knocked out.

        :param rxn: The reaction identifier.
        :param lv: the wild type multiplier factor. The reaction bounds are altered accordingly.
        :returns: A dictionary of reaction constraints.
        """
        constraints = {}
        fluxe_wt = reference[rxn]
        rev_rxn = self.simulator.reverse_reaction(rxn)
        if lv == 0:
            # KO constraint
            constraints[rxn] = (0, 0)
        elif lv == 1 or fluxe_wt == 0:
            # No contraint is applyed
            pass
        elif rev_rxn is None or rev_rxn == rxn:
            # if there is no reverse reaction
            constraints[rxn] = self.ou_constraint(lv, fluxe_wt)
        else:
            # there's a reverse reaction...
            # one of the two reactions needs to be KO, the one with no flux in the wt.
            rev_fluxe_wt = reference[rev_rxn]
            if abs(fluxe_wt) >= abs(rev_fluxe_wt):
                ko_rxn, ou_rxn, fwt = rev_rxn, rxn, fluxe_wt
            else:
                rxn, rev_rxn, rev_fluxe_wt
            constraints[ko_rxn] = (0, 0)
            constraints[ou_rxn] = self.ou_constraint(lv, fwt)
        return constraints
