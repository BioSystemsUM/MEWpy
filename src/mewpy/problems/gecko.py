import warnings

from .problem import AbstractKOProblem, AbstractOUProblem
from ..util.constants import ModelConstants
from ..simulation import SStatus


class GeckoKOProblem(AbstractKOProblem):
    """Gecko KnockOut Optimization Problem

    :param model: The constraint based metabolic model.
    :param list fevaluation: a list of callable EvaluationFunctions.

    Optional:

    :param OrderedDict envcond: environmental conditions.
    :param OrderedDict constraints: additional constraints to be applied to the model.
    :param int candidate_min_size: The candidates minimum size.
    :param int candidate_max_size: The candidates maximum size.
    :param list target: List of target reactions.
    :param list non_target: List of non target reactions. Not considered if a target list is provided.
    :param float scalefactor: a scaling factor to be used in the LP formulation.
    :param str prot_prefix: the protein draw reaction prefix. Default `draw_prot_`.


    Note:

    Targets as well as non target proteins are defined using their prot id, ex  `P0351`, and not by the associated draw
    reaction id, ex `draw_prot_P0351`.

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        # if isinstance(model,GeckoModel):
        super(GeckoKOProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        # else:
        #    raise Exception("The model should be an instance of GeckoModel")
        # problem simulation context
        self.prot_prefix = kwargs.get('prot_prefix', 'draw_prot_')
        self.simulator.prot_prefix = self.prot_prefix

    def _build_target_list(self):
        """
        If not provided, targets are all non essential proteins.
        """
        print("Building modification target list.")
        proteins = set(self.simulator.proteins)
        print("Computing essential proteins.")
        essential = self.simulator.essential_proteins()
        target = proteins - set(essential)
        if self.non_target:
            target = target - set(self.non_target)
        self._trg_list = list(target)

    def decode(self, candidate):
        """
        Decodes a candidate, an integer set, into a dictionary of constraints
        """
        decoded_candidate = dict()
        for idx in candidate:
            try:
                decoded_candidate["{}{}".format(
                    self.prot_prefix, self.target_list[idx])] = 0
            except IndexError:
                raise IndexError(
                    f"Index out of range: {idx} from {len(self.target_list[idx])}")
        return decoded_candidate

    def encode(self, candidate):
        """
        Translates a candidate solution in problem specific representation to
        an iterable of ids, or (ids, folds).

        :param candidate: The candidate representation.
        """
        p_size = len(self.prot_prefix)
        return set([self.target_list.index(k[p_size:]) for k in candidate])

    def solution_to_constraints(self, candidate):
        """
        Converts a candidate, a dictionary of reactions, into a dictionary of constraints.
        This is problem specific. By default return the decoded candidate.
        """
        # check prefix
        def add_prefix(prot):
            if prot.startswith(self.prot_prefix):
                return prot
            else:
                return f"{self.prot_prefix}{prot}"

        _candidate = {add_prefix(k): v for k, v in candidate.items()}
        return _candidate

class GeckoOUProblem(AbstractOUProblem):
    """
    Gecko Under/Over expression Optimization Problem


    :param model (metabolic model): The constraint based metabolic model.
    :param fevaluation (list): a list of callable EvaluationFunctions.

    Optional:

    :param OrderedDict envcond: environmental conditions.
    :param OrderedDict constraints: additional constraints to be applied to the model.
    :param int candidate_min_size: The candidates minimum size.
    :param int candidate_max_size: The candidates maximum size.
    :param list target: List of target reactions.
    :param list non_target: List of non target reactions. Not considered if a target list is provided.
    :param float scalefactor: a scaling factor to be used in the LP formulation.
    :param dic reference: Dictionary of flux values to be used in the over/under expression values computation.
    :param str prot_prefix: the protein draw reaction prefix. Default `draw_prot_`.
    :param boolean twostep: If deletions should be applied before identifiying reference flux values.

    Note:
    Target as well as non target proteins are defined with their prot id, ex `P0351`, and with the associated reaction
    id, ex `draw_prot_P0351`.

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        # if isinstance(model,GeckoModel):
        super(GeckoOUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        # else:
        #    raise Exception("The model should be an instance of GeckoModel")
        # problem simulation context
        self.prot_rev_reactions = None
        self.prot_prefix = kwargs.get('prot_prefix', 'draw_prot_')

    def _build_target_list(self):
        """
        If not provided, targets are all non essential proteins.
        """
        proteins = set(self.simulator.proteins)
        target = proteins
        if self.non_target:
            target = target - set(self.non_target)
        if self._partial_solution:
            target = target - set(self._partial_solution.keys())
        self._trg_list = list(target)

    def decode(self, candidate):
        """
        Decodes a candidate, an integer set, into a dictionary of constraints
        """
        decoded_candidate = dict()
        for idx, lv_idx in candidate:
            try:

                decoded_candidate["{}{}".format(
                    self.prot_prefix, self.target_list[idx])] = self.levels[lv_idx]

            except IndexError:
                raise IndexError(
                    f"Index out of range: {idx} from {len(self.target_list[idx])}")
        return decoded_candidate

    def encode(self, candidate):
        """
        Translates a candidate solution in problem specific representation to
        an iterable of ids, or (ids, folds).

        :param iterable candidate: The candidate representation.
        :returns: a list of index tupple (modification_target_index,level_index). The indexes are
                  problem dependent.
        """
        p_size = len(self.prot_prefix)
        return set([(self.target_list.index(k[p_size:]), self.levels.index(lv))
                    for k, lv in candidate.items()])

    def solution_to_constraints(self, candidate):
        """
        Converts a candidate, a dict {protein:lv}, into a dictionary of constraints
        Reverseble reactions associated to proteins with over expression are KO
        according to the flux volume in the wild type.

        :param candidate: The candidate to be decoded.
        :returns: A dictionary of metabolic constraints.
        """
        constraints = dict()
        if len(candidate) == 0:
            return constraints

        # check prefix
        def add_prefix(prot):
            if prot.startswith(self.prot_prefix):
                return prot
            else:
                return f"{self.prot_prefix}{prot}"

        _candidate = {add_prefix(k): v for k, v in candidate.items()}

        reference = self.reference
        if self.twostep:
            try:
                deletions = {rxn: 0 for rxn, lv in _candidate.items() if lv == 0}
                if deletions and len(deletions) < len(_candidate):
                    sr = self.simulator.simulate(constraints=deletions, method='pFBA')
                    if sr.status in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL):
                        reference = sr.fluxes
            except Exception as e:
                print(e)

        if self.prot_rev_reactions is None:
            self.prot_rev_reactions = self.simulator.protein_rev_reactions

        for rxn, lv in _candidate.items():
            fluxe_wt = reference[rxn]
            prot = rxn[len(self.prot_prefix):]
            if lv < 0:
                raise ValueError("All UO levels should be positive")
            # a level = 0 is interpreted as KO
            elif lv == 0:
                constraints[rxn] = 0.0
            # under expression
            elif lv < 1:
                constraints[rxn] = (0.0, lv * fluxe_wt)
            # TODO: Define how a level 1 is tranlated into constraints...
            elif lv == 1:
                continue
            else:
                constraints[rxn] = (
                    lv * fluxe_wt, ModelConstants.REACTION_UPPER_BOUND)
                # Deals with reverse reactions associated with the protein.
                # Strategy: The reaction direction with no flux in the wild type (reference) is KO.
                if prot in self.prot_rev_reactions.keys():
                    reactions = self.prot_rev_reactions[prot]
                    for r, r_rev in reactions:
                        if reference[r] == 0 and reference[r_rev] == 0:
                            continue
                        elif reference[r] > 0 and reference[r_rev] == 0:
                            constraints[r_rev] = 0.0
                        elif reference[r] == 0 and reference[r_rev] > 0:
                            constraints[r] = 0.0
                        else:
                            warnings.warn(
                                f"Reactions {r} and {r_rev}, associated with the protein {prot},\
                                both have fluxes in the WT.")

        return constraints
