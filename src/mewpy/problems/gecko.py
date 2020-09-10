from .problem import AbstractKOProblem, AbstractOUProblem
from mewpy.utils.parsing import GeneEvaluator, build_tree, Boolean
from mewpy.utils.constants import ModelConstants
from collections import OrderedDict
import warnings


class GeckoKOProblem(AbstractKOProblem):
    """Gecko KnockOut Optimization Problem 

    :param model: The constraint based metabolic model.
    :param list fevaluation: a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness.


    Optional:

    :param OrderedDict envcond: environmental conditions.
    :param OrderedDict constraints: additional constraints to be applied to the model.
    :param int candidate_min_size: The candidates minimum size.
    :param int candidate_max_size: The candidates maximum size.
    :param list target: List of target reactions.
    :param list non_target: List of non target reactions. Not considered if a target list is provided.
    :param float scalefactor: a scaling factor to be used in the LP formulation. 
    :param str prot_prefix: the protein draw reaction prefix. Default `draw_prot_`  


    Note:
    
    Targets as well as non target proteins are defined using their prot id, ex  `P0351`, and not by the associated draw reaction id, ex `draw_prot_P0351`. 
    
    """

    def __init__(self, model, fevaluation=None, **kwargs):
        # if isinstance(model,GeckoModel):
        super(GeckoKOProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        # else:
        #    raise Exception("The model should be an instance of GeckoModel")
        # problem simulation context
        self.prot_prefix = kwargs.get('prot_prefix', 'draw_prot_')

    def _build_target_list(self):
        """
        If not provided, targets are all non essential proteins.
        """
        proteins = set(self.simulator.proteins)
        # as draw_prot_XXXXXX
        ess = self.simulator.essential_proteins(self.prot_prefix)
        # remove 'draw_prot_'
        n = len(self.prot_prefix)
        essential = set([p[n:] for p in ess])
        target = proteins - essential
        if self.non_target:
            target = target - set(self.non_target)
        self._trg_list = list(target)

    def decode(self, candidate):
        """
        Decodes a candidate, an integer set, into a dictionary of constraints
        """
        constraints = OrderedDict()
        for idx in candidate:
            try:
                constraints["{}{}".format(
                    self.prot_prefix, self.target_list[idx])] = 0
            except IndexError:
                raise IndexError(
                    f"Index out of range: {idx} from {len(self.target_list[idx])}")
        return constraints


class GeckoOUProblem(AbstractOUProblem):
    """
    Gecko Under/Over expression Optimization Problem 

  
    :param model (metabolic model): The constraint based metabolic model.
    :param fevaluation (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness.


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


    Note:
    Target as well as non target proteins are defined with their prot id, ex `P0351`, and with the associated reaction id, ex `draw_prot_P0351`. 
    
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
        self._trg_list = list(target)

    def decode(self, candidate):
        """
        Decodes a candidate, a set (idx,lv), into a dictionary of constraints
        Reverseble reactions associated to proteins with over expression are KO 
        according to the flux volume in the wild type.

        Note: Fluxes in Yeast7 gecko model are always non negative 
        """
        constraints = OrderedDict()

        if self.prot_rev_reactions is None:
            self.prot_rev_reactions = self.simulator.protein_rev_reactions

        for idx, lv_idx in candidate:
            try:
                prot = self.target_list[idx]
                rxn = self.prot_prefix+prot
                lv = self.levels[lv_idx]
                fluxe_wt = self.reference[rxn]
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
                        lv*fluxe_wt, ModelConstants.REACTION_UPPER_BOUND)
                    # Deals with reverse reactions associated with the protein.
                    # Strategy: The reaction direction with no flux in the wild type (reference) is KO.
                    if prot in self.prot_rev_reactions.keys():
                        reactions = self.prot_rev_reactions[prot]
                        for r, r_rev in reactions:
                            if self.reference[r] == 0 and self.reference[r_rev] == 0:
                                continue
                            elif self.reference[r] > 0 and self.reference[r_rev] == 0:
                                constraints[r_rev] = 0.0
                            elif self.reference[r] == 0 and self.reference[r_rev] > 0:
                                constraints[r] = 0.0
                            else:
                                warnings.warn(
                                    f"Reactions {r} and {r_rev}, associated with the protein {prot}, both have fluxes in the WT.")
            except IndexError:
                raise IndexError("Index out of range")

        return constraints
