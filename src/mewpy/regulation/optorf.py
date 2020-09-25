from mewpy.problems.problem import AbstractKOProblem


class OptOrfProblem(AbstractKOProblem):
    """
    EA OptOrf problem based on a RFBA problem.
    Currently only consideres deletions.
    """

    def __init__(self, model, fevaluation, regmodel, **kwargs):

        super(OptOrfProblem, self).__init__(model, fevaluation, **kwargs)

        self.regmodel = regmodel

        # ignore user defined target list
        self._trg_list = None

    def _build_target_list(self):

        """
        The EA target list is the combination [metabolic regulatory genes]+[regulators] for a integrated model which
        is the basis for a RFBA model
        """

        self._trg_list = self.regmodel.genes + [regulator for regulator in self.regmodel.regulators.keys()
                                                if regulator not in self.regmodel.metabolic_regulatory_reactions and
                                                regulator not in self.regmodel.metabolic_regulatory_metabolites and
                                                regulator not in self.regmodel.metabolic_regulatory_genes and
                                                regulator not in self.regmodel.regulatory_conditions]

    def decode(self, candidate):
        """

        OptOrf decode based on a RFBA model

        :param candidate:
        :return:
        """

        regulatory_candidate = {}
        metabolic_candidate = {}

        for t_idx in candidate:

            if t_idx > len(self.regmodel.genes) - 1:
                regulatory_candidate[self.target_list[t_idx]] = 0
            else:
                metabolic_candidate[self.target_list[t_idx]] = 0

        # solve first the regulatory model
        state = self.regmodel.solve_regulatory_model(regulatory_candidate, silently=True)

        # update the resulting state with the metabolic share of the candidate
        state.update(metabolic_candidate)

        # find the affected reactions
        constraints = self.regmodel.decode_metabolic_state(state)

        return constraints
