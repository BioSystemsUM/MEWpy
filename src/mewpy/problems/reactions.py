from .problem import AbstractKOProblem, AbstractOUProblem
from collections import OrderedDict


class RKOProblem(AbstractKOProblem):
    """
    Reaction Knockout Optimization Problem.

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
        super(RKOProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)

    def _build_target_list(self):
        reactions = set(self.simulator.reactions)
        essential = set(self.simulator.essential_reactions)
        drains = set(self.simulator.get_drains())
        target = reactions - essential - drains
        if self.non_target is not None:
            target = target - set(self.non_target)
        self._trg_list = list(target)

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

    arguments:

        * model* (metabolic model): the constraint metabolic model
        * fevaluation* (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness

    **args:

        *envcond* (OrderedDict): environmental conditions
        *constraints* (OrderedDict): additional constraints to be applied to the model 
        *candidate_min_size* : The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
        *candidate_max_size* : The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
        *non_target* (list): list of non target reactions
        *levels* (list): over/under expression levels (Default EAConstants.LEVELS)
    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(ROUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)

    def _build_target_list(self):
        reactions = set(self.simulator.reactions)
        #drains = set(self.simulator.get_drains())
        target = reactions  # - drains
        if self.non_target is not None:
            target = target - set(self.non_target)
        self._trg_list = list(target)

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
