import warnings
from collections import OrderedDict
import numpy as np
from .problem import AbstractKOProblem, AbstractOUProblem


class RKOProblem(AbstractKOProblem):
    """
    Reaction Knockout Optimization Problem.

    :param model: The constraint metabolic model.
    :param list fevaluation: A list of callable EvaluationFunctions.

    Optional parameters:

    :param OrderedDict envcond: Environmental conditions.
    :param OrderedDict constraints: Additional constraints to be applied to the model.
    :param int candidate_min_size: The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
    :param int candidate_max_size: The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
    :param list target: List of modification target reactions.
    :param list non_target: List of non target reactions. Not considered if a target list is provided.
    :param float scalefactor: A scaling factor to be used in the LP formulation.

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(RKOProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)

    def _build_target_list(self):
        """Default modification target builder.
        Removes drains, transport and essential reactions
        """
        reactions = set(self.simulator.reactions)
        essential = set(self.simulator.essential_reactions())
        drains = set(self.simulator.get_drains())
        transport = set(self.simulator.get_transport_reactions())
        target = reactions - essential - drains - transport
        if self.non_target is not None:
            target = target - set(self.non_target)
        target = list(target)
        self._trg_list = target


class ROUProblem(AbstractOUProblem):
    """
    Reaction Over/Under Expression Optimization Problem

    :param model: The constraint metabolic model.
    :param list fevaluation: A list of callable EvaluationFunctions.

    Optional parameters:

    :param OrderedDict envcond: Environmental conditions.
    :param OrderedDict constraints: Additional constraints to be applied to the model.
    :param int candidate_min_size: The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
    :param int candidate_max_size: The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
    :param list target: List of modification target reactions.
    :param list non_target: List of non target reactions. Not considered if a target list is provided.
    :param float scalefactor: A scaling factor to be used in the LP formulation.
    :param dic reference: Dictionary of flux values to be used in the over/under expression values computation.
    :param list levels: Over/under expression levels (Default EAConstants.LEVELS)

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(ROUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)

    def _build_target_list(self):
        reactions = set(self.simulator.reactions)
        # drains = set(self.simulator.get_drains())
        target = reactions  # - drains
        if self.non_target is not None:
            target = target - set(self.non_target)
        target = list(target)
        self._trg_list = target

    def solution_to_constraints(self, candidate):
        """
        Decodes a candidate, an dict {idx:lv} into a dictionary of constraints
        Suposes that reverseble reactions have been treated and bounded with positive flux values
        """
        constraints = dict()
        # print(type(candidate), candidate)
        for rxn, lv in candidate.items():
            rev_rxn = self.simulator.reverse_reaction(rxn)
            # skips if the reverse reaction was already processed
            if rev_rxn and rev_rxn in constraints.keys():
                continue
            elif lv < 0:
                raise ValueError("All UO levels should be positive")
            else:
                constraints.update(self.reaction_constraints(rxn, lv))
        return constraints


class MediumProblem(AbstractOUProblem):
    """
    Medium Optimization Problem. Try to find an optimized uptake configuration. 
    By default all uptake reactions are considered. Uptake reactions not included on 
    a solution candidate are KO.

    :param model: The constraint metabolic model.
    :param list fevaluation: A list of callable EvaluationFunctions.

    Optional parameters:

    :param OrderedDict envcond: Environmental conditions.
    :param OrderedDict constraints: Additional constraints to be applied to the model.
    :param int candidate_min_size: The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
    :param int candidate_max_size: The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
    :param list target: List of modification target reactions.
    :param list non_target: List of non target reactions. Not considered if a target list is provided.
    :param float scalefactor: A scaling factor to be used in the LP formulation.
    :param dic reference: Dictionary of flux values to be used in the over/under expression values computation.
    :param list levels: Over/under expression levels (Default [0,0.1,0.2,...,9.9,10.0])

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(ROUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        self.levels = kwargs.get('levels',np.linspace(0,10,101))
        self.candidate_max_size = kwargs.get(
            'candidate_max_size', len(self.target))

    def _build_target_list(self):
        target = set(self.simulator.get_uptake_reactions())
        if self.non_target is not None:
            target = target - set(self.non_target)
        target = list(target)
        self._trg_list = target

    def solution_to_constraints(self, candidate):
        """
        Decodes a candidate, a dict {idx:lv} into a dictionary of constraints
        Suposes that reversible reactions have been treated and bounded with positive flux values
        """
        constraints = dict()
        from mewpy.util.constants import ModelConstants
        for rxn in self.target_list:
            if rxn in candidate.items():
                constraint[rxn] = (-1*lv, ModelConstants.REACTION_UPPER_BOUND)
            else:
                constraint[rxn] = (0,0)
        return constraints
