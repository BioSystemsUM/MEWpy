from typing import Union, Dict

from mewpy.mew.analysis import FBA
from mewpy.mew.lp import ConstraintContainer, VariableContainer
from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel
from mewpy.solvers.solution import Solution, Status
from mewpy.solvers.solver import VarType, Solver


class pFBA(FBA):

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = False,
                 attach: bool = False):
        """
        Parsimonious Flux Balance Analysis (FBA) of a metabolic model.
        Regular implementation of a pFBA for a metabolic model.

        This pFBA implementation was heavily inspired by pFBA implementation of reframed python package. Take a look at
        the source: https://github.com/cdanielmachado/reframed and https://reframed.readthedocs.io/en/latest/

        For more details consult: https://doi.org/10.1038/msb.2010.47

        :param model: a mewpy Model, MetabolicModel, RegulatoryModel or all. The model is used to retrieve
        variables and constraints to the linear problem
        :param solver: A Solver, CplexSolver, GurobiSolver or OptLangSolver instance.
        Alternatively, the name of the solver is also accepted.
        The solver interface will be used to load and solve a linear problem in a given solver.
        If none, a new solver is instantiated. An instantiated solver may be used, but it will be overwritten
        if build is true.
        :param build: Whether to build the linear problem upon instantiation. Default: False
        :param attach: Whether to attach the linear problem to the model upon instantiation. Default: False
        """
        super().__init__(model=model, solver=solver, build=build, attach=attach)

    def _wt_bounds(self, fraction: float = None, solver_kwargs: Dict = None):
        """
        It builds the linear problem from the model. The linear problem is built from the model
        variables and constraints. The linear problem is then loaded into the solver.
        :return:
        """
        if not solver_kwargs:
            solver_kwargs = {}

        sol = FBA(model=self.model, build=True, attach=False).optimize(solver_kwargs=solver_kwargs, to_solver=True)
        if sol.status != Status.OPTIMAL:
            lb, ub = 0.0, 0.0

        else:
            if fraction is None:
                lb, ub = float(sol.fobj), float(sol.fobj)
            else:
                lb, ub = float(sol.fobj) * fraction, float(sol.fobj)
        return lb, ub

    def _build_pfba_constrains(self, fraction: float = None, solver_kwargs: Dict = None):
        """
        It builds the pfba constraints of the linear problem.
        :return:
        """
        if not solver_kwargs:
            solver_kwargs = {}

        lb, ub = self._wt_bounds(fraction, solver_kwargs)

        if 'linear' in solver_kwargs:
            coef = solver_kwargs['linear'].copy()
        else:
            coef = {variable.id: val for variable, val in self.model.objective.items()}

        if 'constraints' in solver_kwargs:
            constraints = solver_kwargs['constraints'].copy()
        else:
            constraints = {}

        constraint = ConstraintContainer(name='pfba_constraints', coefs=[coef], lbs=[lb], ubs=[ub])
        variable = VariableContainer(name='pfba_variables', sub_variables=[], lbs=[], ubs=[], variables_type=[])
        objective = {}
        for reaction in self.model.yield_reactions():

            if reaction.reversibility:
                rxn_forward = f'{reaction.id}_forward'
                rxn_reverse = f'{reaction.id}_reverse'

                rxn_ub = float(constraints.get(reaction.id, reaction.bounds)[1])

                variable.sub_variables.extend([rxn_forward, rxn_reverse])
                variable.lbs.extend([0.0, 0.0])
                variable.ubs.extend([rxn_ub, rxn_ub])
                variable.variables_type.extend([VarType.CONTINUOUS, VarType.CONTINUOUS])

                constraint.lbs.extend([0.0, 0.0])
                constraint.ubs.extend([rxn_ub, rxn_ub])
                constraint.coefs.extend([{reaction.id: -1, rxn_forward: 1},
                                         {reaction.id: 1, rxn_reverse: 1}])

                objective[rxn_forward] = 1
                objective[rxn_reverse] = 1

            else:
                objective[reaction.id] = 1

        self.add_variables(variable)
        self.add_constraints(constraint)
        self._linear_objective = objective
        self._minimize = True

    def _build(self):
        """
        It builds the linear problem from the model. The linear problem is built from the model
        variables and constraints. The linear problem is then loaded into the solver.
        :return:
        """
        if self.model.is_metabolic():
            # mass balance constraints and reactions' variables
            self._build_mass_constraints()

            # pFBA constraints can be added again to the linear problem during optimization
            self._build_pfba_constrains()

        return

    def _optimize(self, fraction: float = None, solver_kwargs: Dict = None, **kwargs) -> Solution:
        """
        It optimizes the linear problem. The linear problem is solved by the solver interface.
        :param solver_kwargs: A dictionary of keyword arguments to be passed to the solver.
        :return: A Solution instance.
        """
        if not solver_kwargs:
            solver_kwargs = {}

        linear = solver_kwargs.get('linear')
        constraints = solver_kwargs.get('constraints')

        # if linear and constraints are not provided, build new pfba constraints and solver
        replace_pfba_constraints = [x for x in (fraction, linear, constraints) if x is not None]

        if replace_pfba_constraints:
            self._build_pfba_constrains(solver_kwargs=solver_kwargs)
            self.build_solver()

        solution = self.solver.solve(**solver_kwargs)

        # restore the pfba constraints and solver to the previous state
        if replace_pfba_constraints:
            self._build_pfba_constrains()
            self.build_solver()

        return solution
