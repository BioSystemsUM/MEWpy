from typing import Union, Dict

from mewpy.mew.lp import ConstraintContainer, VariableContainer, LinearProblem
from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel
from mewpy.solvers.solution import Solution
from mewpy.solvers.solver import VarType, Solver


class FBA(LinearProblem):

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = False,
                 attach: bool = False):
        """
        Flux Balance Analysis (FBA) of a metabolic model. Regular implementation of a FBA for a metabolic model.

        For more details consult: https://dx.doi.org/10.1038%2Fnbt.1614

        :param model: a mewpy Model, MetabolicModel, RegulatoryModel or all. The model is used to retrieve
        variables and constraints to the linear problem
        :param solver: A Solver, CplexSolver, GurobiSolver or OptLangSolver instance.
        Alternatively, the name of the solver is also accepted.
        The solver interface will be used to load and solve a linear problem in a given solver.
        If none, a new solver is instantiated. An instantiated solver may be used,
        but it will be overwritten if build is true.
        :param build: Whether to build the linear problem upon instantiation. Default: False
        :param attach: Whether to attach the linear problem to the model upon instantiation. Default: False
        """
        super().__init__(model=model, solver=solver, build=build, attach=attach)

    def _build_mass_constraints(self):
        gene_state = {gene.id: max(gene.coefficients) for gene in self.model.yield_genes()}

        constraints = {metabolite.id: ConstraintContainer(name=metabolite.id, lbs=[0.0], ubs=[0.0], coefs=[{}])
                       for metabolite in self.model.yield_metabolites()}
        variables = {}

        for reaction in self.model.yield_reactions():
            if reaction.gpr.is_none:
                lb, ub = reaction.bounds

            else:
                res = reaction.gpr.evaluate(values=gene_state)
                if not res:
                    lb, ub = 0.0, 0.0
                else:
                    lb, ub = reaction.bounds

            variable = VariableContainer(name=reaction.id, sub_variables=[reaction.id],
                                         lbs=[float(lb)], ubs=[float(ub)], variables_type=[VarType.CONTINUOUS])
            variables[reaction.id] = variable

            for metabolite, stoichiometry in reaction.stoichiometry.items():
                constraints[metabolite.id].coefs[0][reaction.id] = stoichiometry

        self.add_variables(*variables.values())
        self.add_constraints(*constraints.values())
        return

    def _build(self):
        """
        It builds the linear problem from the model. The linear problem is built from the model
        variables and constraints. The linear problem is then loaded into the solver.
        :return:
        """
        if self.model.is_metabolic():
            # mass balance constraints and reactions' variables
            self._build_mass_constraints()

            self._linear_objective = {var.id: value for var, value in self.model.objective.items()}
            self._minimize = False

        return

    def _optimize(self, solver_kwargs: Dict = None, **kwargs) -> Solution:
        """
        It optimizes the linear problem. The linear problem is solved by the solver interface.
        :param solver_kwargs: A dictionary of keyword arguments to be passed to the solver.
        :return: A Solution instance.
        """
        return self.solver.solve(**solver_kwargs)
