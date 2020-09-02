from .problem import AbstractKOProblem, AbstractOUProblem
from mewpy.utils.parsing import GeneEvaluator, build_tree, Boolean
from collections import OrderedDict
from mewpy.utils.constants import EAConstants
import logging

logger = logging.getLogger('mewpy.problems.genes')

class GKOProblem(AbstractKOProblem):
    """
    Gene Knockout Optimization Problem

    arguments:

        * model* (metabolic model): the constraint metabolic model
        * fevaluation* (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness

    **args:

        *envcond* (OrderedDict): environmental conditions
        *constraints* (OrderedDict): additional constraints to be applied to the model 
        *candidate_min_size* : The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
        *candidate_max_size* : The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
        *non_target* (list): list of non target genes
    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(GKOProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)

    def _build_target_list(self):

        genes = set(self.simulator.genes)
        essential = set(self.simulator.essential_genes)
        target = genes - essential
        if self.non_target:
            target = target - set(self.non_target)
        self._trg_list = list(target)

    def decode(self, candidate):
        """
        Decodes a candidate, an set of genes into a dictionary of constraints
        """
        genes = [self.target_list[idx] for idx in candidate]
        active_genes = set(self.simulator.genes) - set(genes)
        active_reactions = self.simulator.evaluate_gprs(active_genes)
        inactive_reactions = set(
            self.simulator.reactions) - set(active_reactions)
        gr_constraints = {rxn: 0 for rxn in inactive_reactions}
        return gr_constraints


class GOUProblem(AbstractOUProblem):
    """
    Gene Over/Under expression Optimization Problem

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
        reference (dic): Dictionary of flux values to be used in the over/under expression values computation.
        operators (tuple of function): (and, or) operations. Default (MIN, MAX) 

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(GOUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        # operators to replace 'and'/'or'. By default min/max
        self._operators = kwargs.get('operators', None)

    def _build_target_list(self):

        target = set(self.simulator.genes)
        if self.non_target:
            target = target - set(self.non_target)
        self._trg_list = list(target)

    def decode(self, candidate):
        """
        Decodes a candidate, a set of (genes,lv) into a dictionary of reaction constraints
        """
        gr_constraints = OrderedDict()
        genes = {self.target_list[idx]: self.levels[lv_idx]
                 for idx, lv_idx in candidate}
        # evaluate gpr
        if not self._operators:
            self._operators = (lambda x, y: min(x, y), lambda x, y: max(x, y))

        evaluator = GeneEvaluator(
            genes, self._operators[0], self._operators[1])
        for rxn_id in self.simulator.reactions:
            gpr = self.simulator.get_gpr(rxn_id)
            if gpr:
                tree = build_tree(gpr, Boolean)
                # apply the operators to obtain a level for the reaction
                # if a gene as no level associated its factor is 1 (see GeneEvaluator)
                lv = tree.evaluate(evaluator.f_operand, evaluator.f_operator)
                # debugging 
                logger.debug(f"{gpr}\n{genes}\nlevel:{lv}\n")
                
                # adds the reaction constraint
                rev_rxn = self.simulator.reverse_reaction(rxn_id)
                # skips if the reverse reaction was already processed
                if rev_rxn and rev_rxn in gr_constraints.keys():
                    continue
                elif lv < 0:
                    raise ValueError("All UO levels should be positive")
                else:
                    gr_constraints.update(
                        self.reaction_constraints(rxn_id, lv))

        return gr_constraints
