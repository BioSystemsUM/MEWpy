import logging
from .problem import AbstractKOProblem, AbstractOUProblem
from ..util.parsing import GeneEvaluator, build_tree, Boolean
from ..simulation import SStatus
logger = logging.getLogger(__name__)


class GKOProblem(AbstractKOProblem):
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

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(GKOProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)

    def _build_target_list(self):
        print("Building modification target list.")
        genes = set(self.simulator.genes)
        essential = set(self.simulator.essential_genes())
        transport = set(self.simulator.get_transport_genes())
        target = genes - essential - transport
        if self.non_target:
            target = target - set(self.non_target)
        target = list(target)
        self._trg_list = target

    def solution_to_constraints(self, candidate):
        """
        Converts a candidate, dict of genes:0 into a dictionary of constraints.
        """
        genes = list(candidate.keys())
        active_genes = set(self.simulator.genes) - set(genes)
        active_reactions = self.simulator.evaluate_gprs(active_genes)
        inactive_reactions = set(
            self.simulator.reactions) - set(active_reactions)
        gr_constraints = {rxn: 0 for rxn in inactive_reactions}
        return gr_constraints


class GOUProblem(AbstractOUProblem):
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
    :param boolean twostep: If deletions should be applied before identifiying reference flux values.

    Note:  Operators that can not be pickled may be defined by a string e.g. 'lambda x,y: (x+y)/2'.

    """

    def __init__(self, model, fevaluation=None, **kwargs):
        super(GOUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        # operators to replace 'and'/'or'. By default min/max
        self._temp_op = kwargs.get('operators', None)
        self._operators = None

    def _build_target_list(self):

        genes = set(self.simulator.genes)
        transport = set(self.simulator.get_transport_genes())
        target = genes - transport
        if self.non_target:
            target = target - set(self.non_target)
        if self._partial_solution:
            target = target - set(self._partial_solution.keys())
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

    def solution_to_constraints(self, candidate):
        """
        Decodes a candidate, a dict of genes:lv into a dictionary of reaction constraints
        """
        gr_constraints = dict()
        genes = candidate

        # Computes reference fluxes based on deletions
        reference = self.reference
        if self.twostep:
            try:
                deletions = [gene for gene, lv in candidate.items() if lv == 0]
                active_genes = set(self.simulator.genes) - set(deletions)
                active_reactions = self.simulator.evaluate_gprs(active_genes)
                inactive_reactions = set(self.simulator.reactions) - set(active_reactions)
                gr_constraints = {rxn: 0 for rxn in inactive_reactions}
                sr = self.simulator.simulate(constraints=gr_constraints, method='pFBA')
                if sr.status in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL):
                    reference = sr.fluxes
            except Exception as e:
                print(e)
        # operators check
        self.__op()
        # evaluate gpr
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
                        self.reaction_constraints(rxn_id, lv, reference))

        return gr_constraints
