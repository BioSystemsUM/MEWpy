from enum import Enum
from abc import ABC, abstractmethod
from mewpy.utils.constants import EAConstants, ModelConstants
from mewpy.simulation import get_simulator
from collections import OrderedDict
import warnings
import numpy as np
import copy


class Strategy(Enum):
    """
    The available optimization strategies
    """
    KO = 'KO'
    OU = 'OU'


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
    Optimization Problem base class

    parameters:

    model (metabolic model)
    fevaluation (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        self.model = model
        self.fevaluation = fevaluation
        self.number_of_objectives = len(self.fevaluation)

        # simulation context : defines the simulations environment
        self.simul_context = None
        self.environmental_conditions = kwargs.get('envcond', None)
        self.persistent_constraints = kwargs.get('constraints', None)
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
    def decode(self, candidate):
        """The decoder function for the problem."""
        raise NotImplementedError

    def get_name(self):
        """The decoder function for the problem."""
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
        "list of allowed OU"
        if self._trg_list is None:
            self._build_target_list()
        return self._trg_list

    @abstractmethod
    def _build_target_list(self):
        raise NotImplementedError

    def get_constraints(self, solution):
        """
        returns the constrainst enconded into an individual
        """
        return self.decode(solution.candidate)

    def get_environmental_conditions(self):
        return self.environmental_conditions

    def get_persistent_constraints(self):
        return self.persistent_constraints

    def evaluate_solution(self, solution):
        """
        Evaluates a single solution, a list of constraints
        returns a list of fitness
        """
        p = []
        # decoded constraints
        constraints = self.decode(solution)
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
            warnings.warn(
                f"Solution couldn't be evaluated [{e}]\n {constraints}")
        return p

    @property
    def is_maximization(self):
        return all([f.maximize for f in self.fevaluation])

    def simplify(self, solution, tolerance=1e-6):
        """
        Simplify a solution by removing the modification that not affect the final fitness value.
        Args:
            solution : the solution to be simplified
            tolerance: max allowed objective difference values for two solutions to be considered diferent.
                       Tolerance may be defined by a single float value, or per objective by means of a list of floats
                       of size equal to the number of objectives.
        Returns: a list of simplified constraints
        """

        constraints = dict(copy.copy(self.get_constraints(solution)))
        constraints_to_remove = []
        if self.number_of_objectives == 1:
            fitness = [solution.fitness]
        else:
            fitness = solution.fitness
        for key in constraints.keys():
            simul_constraints = copy.copy(constraints)
            del simul_constraints[key]
            # pre simulation
            simulation_results = OrderedDict()
            for method in self.methods:
                simulation_result = self.simulator.simulate(
                    constraints=simul_constraints, method=method, scalefactor=self.scalefactor)
                simulation_results[method] = simulation_result
            # apply the evaluation function(s)
            fit = []
            for f in self.fevaluation:
                fit.append(f(simulation_results, solution,
                             scalefactor=self.scalefactor))
            diff = np.abs(np.array(fit)-np.array(fitness))
            is_equal = False
            if isinstance(tolerance, float):
                is_equal = np.all(diff <= tolerance)
            else:
                is_equal = np.all(diff <= np.array(tolerance))
            if is_equal:
                constraints_to_remove.append(key)
        for key in constraints_to_remove:
            del constraints[key]
        return constraints


class AbstractKOProblem(AbstractProblem):
    """
    Base class for Knockout optimization problems

    args:

        model (metabolic model): The constraint based metabolic model.
        fevaluation (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness


    **kwargs:

        envcond (OrderedDict): environmental conditions.
        constraints (OrderedDict): additional constraints to be applied to the model.
        candidate_min_size (int) : The candidates minimum size.
        candidate_min_size (int) : The candidates maximum size.
        target (list): List of target reactions.
        non_target (list): List of non target reactions. Not considered if a target list is provided.
        scalefactor (floaf): a scaling factor to be used in the LP formulation. 
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
        Generates a solution, a random int set with length in range min_solution_size to max_solution_size
        """
        solution = set()
        solution_size = random.uniform(
            self.candidate_min_size, self.candidate_max_size)
        while len(solution) < solution_size:
            solution.add(random.randint(0, len(self.target_list)-1))
        return solution


class AbstractOUProblem(AbstractProblem):
    """ Base class for Over/Under expression optimization problems 

    arguments:

        * model* (metabolic model): the constraint metabolic model
        * fevaluation* (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness

    **args:

        envcond (OrderedDict): environmental conditions
        constraints (OrderedDict): additional constraints to be applied to the model 
        candidate_min_size : The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
        candidate_max_size : The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
        non_target (list): list of non target reactions
        levels (list): over/under expression levels (Default EAConstants.LEVELS)
        reference (dic): dictionary of reference flux values

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(AbstractOUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        self.strategy = Strategy.OU
        self.levels = kwargs.get('levels', EAConstants.LEVELS)
        self._reference = kwargs.get('reference', None)

    @abstractmethod
    def decode(self, candidate):
        """The decoder function for the problem."""
        raise NotImplementedError

    @property
    def bounder(self):
        """
        The list and levels index bounder  
        """
        if self._bounder is None:
            max_idx = len(self.target_list)-1
            max_lv = len(self.levels)-1
            self._bounder = OUBounder([0, 0], [max_idx, max_lv])
        return self._bounder

    def generator(self, random, args):
        """
        Generates a solution, a random (int,int) set with length in range min_solution_size to max_solution_size
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
        if not self._reference:
            self._reference = self.simulator.reference
        return self._reference

    def ou_constraint(self, level, wt):
        if level > 1:
            return (level*wt, ModelConstants.REACTION_UPPER_BOUND) if wt >= 0 else (-1*ModelConstants.REACTION_UPPER_BOUND, level*wt)
        else:
            return (0, level * wt) if wt >= 0 else (level * wt, 0)

    def reaction_constraints(self, rxn, lv):
        """
        Converts a (reaction, level) pair into a constraint
        If a reaction is reversible, the direction with no or less wild type flux
        is knocked out.
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
            # one of the two reactions needs to be KO, the one with no (or lesser) flux in the wt
            rev_fluxe_wt = self.reference[rev_rxn]
            if abs(fluxe_wt) >= abs(rev_fluxe_wt):
                ko_rxn, ou_rxn, fwt = rev_rxn, rxn, fluxe_wt
            else:
                rxn, rev_rxn, rev_fluxe_wt
            constraints[ko_rxn] = (0, 0)
            constraints[ou_rxn] = self.ou_constraint(lv, fwt)
        return constraints
