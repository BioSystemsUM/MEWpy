from .problem import AbstractKOProblem, AbstractOUProblem
from collections import OrderedDict
import warnings


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
        reactions = set(self.simulator.reactions)
        essential = set(self.simulator.essential_reactions)
        drains = set(self.simulator.get_drains())
        target = reactions - essential - drains
        if self.non_target is not None:
            target = target - set(self.non_target)
        target =  list(target)
        try:
            from mewpy.utils.constants import EAConstants
            if EAConstants.PROB_TARGET and self.product:
                from mewpy.utils.graph import probabilistic_reaction_targets
                target = probabilistic_reaction_targets(self.model,self.product,target)        
        except Exception as e:
            warnings.warn(str(e))
            
        self._trg_list = target

    def decode(self, candidate):
        """
        Decodes a candidate, an integer set, into a dictionary of constraints
        """
        constraints = OrderedDict()
        for idx in candidate:
            try:
                constraints[self.target_list[idx]] = 0
            except IndexError:
                raise IndexError("Index out of range: {} from {}".format(
                    idx, len(self.target_list[idx])))
        return constraints


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
        target =  list(target)
        try:
            from mewpy.utils.constants import EAConstants
            if EAConstants.PROB_TARGET and self.product:
                from mewpy.utils.graph import probabilistic_reaction_targets
                target = probabilistic_reaction_targets(self.model,self.product,target)        
        except Exception as e:
            warnings.warn(str(e))
            
        self._trg_list = target

    def decode(self, candidate):
        """
        Decodes a candidate, an set (idx,lv) into a dictionary of constraints
        Suposes that reverseble reactions have been treated and bounded with positive flux values
        """
        constraints = OrderedDict()
        for idx, lv_idx in candidate:
            try:
                rxn = self.target_list[idx]
                lv = self.levels[lv_idx]
                rev_rxn = self.simulator.reverse_reaction(rxn)
                # skips if the reverse reaction was already processed
                if rev_rxn and rev_rxn in constraints.keys():
                    continue
                elif lv < 0:
                    raise ValueError("All UO levels should be positive")
                else:
                    constraints.update(self.reaction_constraints(rxn, lv))
            except IndexError:
                raise IndexError("Index out of range")
        return constraints
