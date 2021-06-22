import logging
import warnings
import itertools

from .problem import AbstractKOProblem, AbstractOUProblem
from ..util.parsing import GeneEvaluator, build_tree, Boolean

logger = logging.getLogger(__name__)


def gene_has_associated_enzyme(model, gene):
    if any([gene in x.composition for x in model.enzymes]):
        try:
            return model._get_translation_name(gene_id)
        except:
            return None
    return None


def associated_enzyme(model, gene):
    return [x for x in model.enzymes if gene in x.composition]


class ETFLGKOProblem(AbstractKOProblem):
    """
    Gene Knockout Optimization Problem.

    :param model: The constraint metabolic model.
    :param list fevaluation: A list of callable EvaluationFunctions.

    Optional:

    :param OrderedDict envcond: Environmental conditions.
    :param OrderedDict constraints: Additional constraints to be applied to the model.
    :param int candidate_min_size: The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
    :param int candidate_max_size: The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
    :param list target: List of modification target genes.
    :param list non_target: List of non target genes. Not considered if a target list is provided.
    :param float scalefactor: A scaling factor to be used in the LP formulation.
    :param boolean only_gpr: Only uses GPRs and do not alter pseudo translation reactions bounds (Default False).

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(ETFLGKOProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        self._only_gpr = kwargs.get('only_gpr', False)
        self.gene_reaction_mapping()

    def gene_reaction_mapping(self):
        enzyme_reaction = {}
        for rx in self.model.reactions:
            if hasattr(rx, 'enzymes') and rx.enzymes:
                for e in rx.enzymes:
                    enzyme_reaction.setdefault(e.id, []).append(rx.id)
        gene_reaction = {}
        self.has_enzyme = []
        for g in self.simulator.genes:
            gene_reaction[g] = []
            ee = associated_enzyme(self.model, g)
            if ee:
                self.has_enzyme.append(g)
            for e in ee:
                gene_reaction[g].extend(enzyme_reaction[e.id])
        self.gene_reaction = gene_reaction

    def _build_target_list(self):
        genes = set(self.simulator.genes)
        essential = set(self.simulator.essential_genes())
        transport = set(self.simulator.get_transport_genes())
        target = genes - essential - transport
        if self.non_target:
            target = target - set(self.non_target)
        target = list(target)
        self._trg_list = target

    def _trans_solution_to_constraints(self, candidate):
        """
        Converts a candidate, dict of genes:0 into a dictionary of constraints.
        """
        genes = list(candidate.keys())
        gr_constraints = dict()
        no_tans = []
        # Translation
        for g in genes:
            if g in self.has_enzyme:
                try:
                    rx = self.model._get_translation_name(g)
                    gr_constraints[rx] = 0
                except:
                    no_trans.append(g)
        # GPR based reaction KO
        active_genes = set(self.simulator.genes) - set(genes)
        active_reactions = self.simulator.evaluate_gprs(active_genes)
        catalyzed_reactions = set(itertools.chain.from_iterable(
            [self.gene_reaction[g] for g in genes if g not in no_trans]))
        inactive_reactions = set(self.simulator.reactions) - set(active_reactions) - catalyzed_reactions
        gr_constraints = {rxn: 0 for rxn in inactive_reactions}
        return gr_constraints

    def _gpr_solution_to_constraints(self, candidate):
        """
        Converts a candidate, dict of genes:0 into a dictionary of constraints.
        """
        genes = list(candidate.keys())
        gr_constraints = dict()
        trans = [gene_has_associated_enzyme(self.model, g) for g in genes]
        # USE TRANSLATION REACTION
        if all(trans):
            gr_constraints = {rxn: 0 for rxn in trans}
        # USE GPR
        else:
            active_genes = set(self.simulator.genes) - set(genes)
            active_reactions = self.simulator.evaluate_gprs(active_genes)
            inactive_reactions = set(
                self.simulator.reactions) - set(active_reactions)
            gr_constraints = {rxn: 0 for rxn in inactive_reactions}
        return gr_constraints

    def solution_to_constraints(self, candidate):
        if self._only_gpr:
            return self._gpr_solution_to_constraints(candidate)
        else:
            return self._trans_solution_to_constraints(candidate)


class ETFLGOUProblem(AbstractOUProblem):
    """Gene Over/Under expression Optimization Problem

    :param model: The constraint metabolic model.
    :param list fevaluation: A list of callable EvaluationFunctions.

    Optional:

    :param OrderedDict envcond: Environmental conditions.
    :param OrderedDict constraints: Additional constraints to be applied to the model.
    :param int candidate_min_size: The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
    :param int candidate_max_size: The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
    :param list target: List of modification target genes.
    :param list non_target: List of non target genes. Not considered if a target list is provided.
    :param float scalefactor: A scaling factor to be used in the LP formulation.
    :param dic reference: Dictionary of flux values to be used in the over/under expression values computation.
    :param tuple operators: (and, or) operations. Default (MIN, MAX).
    :param list levels: Over/under expression levels (Default EAConstants.LEVELS).
    :param boolean only_gpr: Only uses GPRs and do not alter pseudo translation reactions bounds (Default False).

    Note:  Operators that can not be pickled may be defined by a string e.g. 'lambda x,y: (x+y)/2'.

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(ETFLGOUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        # operators to replace 'and'/'or'. By default min/max
        self._temp_op = kwargs.get('operators', None)
        self._operators = None
        self._only_gpr = kwargs.get('only_gpr', False)
        self.gene_reaction_mapping()

    def gene_reaction_mapping(self):
        enzyme_reaction = {}
        for rx in self.model.reactions:
            if hasattr(rx, 'enzymes') and rx.enzymes:
                for e in rx.enzymes:
                    enzyme_reaction.setdefault(e.id, []).append(rx.id)
        gene_reaction = {}
        self.has_enzyme = []
        for g in self.simulator.genes:
            gene_reaction[g] = []
            ee = associated_enzyme(self.model, g)
            if ee:
                self.has_enzyme.append(g)
            for e in ee:
                gene_reaction[g].extend(enzyme_reaction[e.id])
        self.gene_reaction = gene_reaction

    def _build_target_list(self):

        genes = set(self.simulator.genes)
        transport = set(self.simulator.get_transport_genes())
        target = genes - transport
        if self.non_target:
            target = target - set(self.non_target)
        target = list(target)
        self._trg_list = target

    def __op(self):
        # set default operators as configurable options
        if self._operators:
            pass
        elif not self._temp_op or None in self._temp_op or len(self._temp_op) < 2:
            self._operators = (lambda x, y: min(x, y), lambda x, y: max(x, y))
        else:
            ops = []
            for i in [0, 1]:
                op = None
                if isinstance(self._temp_op[i], str):
                    op = eval(self._temp_op[i])
                else:
                    op = self._temp_op[i]
                if callable(op):
                    ops.append(op)
                else:
                    raise ValueError(f"The operator at index {i} is not callable.")
            self._operators = tuple(ops)

    def _gpr_solution_to_constraints(self, candidate):
        """
        Decodes a candidate, a dict of genes:lv into a dictionary of reaction constraints
        """
        gr_constraints = dict()
        # USE TRANSLATION REACTION
        trans = [gene_has_associated_enzyme(self.model, g) for g in candidate.keys()]
        if all(trans):
            for gene_id, lv in candidate.items():
                rxn_id = self.model._get_translation_name(gene_id)
                gr_constraints.update(
                    self.reaction_constraints(rxn_id, lv))
        # USE GPR
        else:
            # operators check
            self.__op()
            # evaluate gpr
            evaluator = GeneEvaluator(
                candidate, self._operators[0], self._operators[1])
            for rxn_id in self.simulator.reactions:
                gpr = self.simulator.get_gpr(rxn_id)
                if gpr:
                    tree = build_tree(gpr, Boolean)
                    # apply the operators to obtain a level for the reaction
                    # if a gene as no level associated its factor is 1 (see GeneEvaluator)
                    lv = tree.evaluate(evaluator.f_operand, evaluator.f_operator)

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

    def _trans_solution_to_constraints(self, candidate):
        """
        Decodes a candidate, a dict of genes:lv into a dictionary of reaction constraints
        """
        gr_constraints = dict()
        no_trans = []
        # translation reactions
        for gene_id, lv in candidate.items():
            if gene_id in self.has_enzyme:
                try:
                    rx = self.model._get_translation_name(gene_id)
                    gr_constraints.update(
                        self.reaction_constraints(rx, lv))
                except Exception as ex:
                    no_trans.append(gene_id)
        catalyzed_reactions = set(itertools.chain.from_iterable(
            [self.gene_reaction[g] for g in candidate if g not in no_trans]))
        # GPR based reaction
        self.__op()
        # evaluate gpr
        evaluator = GeneEvaluator(
            candidate, self._operators[0], self._operators[1])
        for rxn_id in self.simulator.reactions:
            if rxn_id not in catalyzed_reactions:
                gpr = self.simulator.get_gpr(rxn_id)
                if gpr:
                    tree = build_tree(gpr, Boolean)
                    # apply the operators to obtain a level for the reaction
                    # if a gene as no level associated its factor is 1 (see GeneEvaluator)
                    lv = tree.evaluate(evaluator.f_operand, evaluator.f_operator)
                    # debugging
                    # logger.debug(f"{gpr}\n{genes}\nlevel:{lv}\n")

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

    def solution_to_constraints(self, candidate):
        if self._only_gpr:
            return self._gpr_solution_to_constraints(candidate)
        else:
            return self._trans_solution_to_constraints(candidate)
