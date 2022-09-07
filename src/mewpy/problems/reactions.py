import numpy as np
from .problem import AbstractKOProblem, AbstractOUProblem
from ..simulation import SStatus


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
        print("Building modification target list.")
        reactions = set(self.simulator.reactions)
        essential = set(self.simulator.essential_reactions())
        drains = set(self.simulator.get_exchange_reactions())
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
    :param boolean twostep: If deletions should be applied before identifiying reference flux values.
    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(ROUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)

    def _build_target_list(self):
        print("Building modification target list.")
        reactions = set(self.simulator.reactions)
        # drains = set(self.simulator.get_exchange_reactions())
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
        # computes reference fluxes based on deletions
        reference = self.reference
        if self.twostep:
            try:
                deletions = {rxn: 0 for rxn, lv in candidate.items() if lv == 0}
                sr = self.simulator.simulate(constraints=deletions, method='pFBA')
                if sr.status in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL):
                    reference = sr.fluxes
            except Exception as e:
                print(e)

        for rxn, lv in candidate.items():
            rev_rxn = self.simulator.reverse_reaction(rxn)
            # skips if the reverse reaction was already processed
            if rev_rxn and rev_rxn in constraints.keys():
                continue
            elif lv < 0:
                raise ValueError("All UO levels should be positive")
            else:
                constraints.update(self.reaction_constraints(rxn, lv, reference))
        return constraints


class MediumProblem(AbstractOUProblem):
    """
    Medium Optimization Problem. Try to find an optimized uptake configuration.
    By default all uptake reactions are considered. Uptake reactions not included in
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
        super(MediumProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        self.levels = kwargs.get('levels', np.linspace(0, 10, 101))
        self.candidate_max_size = kwargs.get(
            'candidate_max_size', len(self.target_list))
        self.simulator._allow_env_changes = True

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
            if rxn in candidate.keys():
                lv = candidate[rxn]
                constraints[rxn] = (-1 * lv, ModelConstants.REACTION_UPPER_BOUND)
            else:
                constraints[rxn] = (0, 0)
        return constraints
