from typing import Union, Dict

from mewpy.mew.analysis import FBA
from mewpy.mew.lp import LinearMixIn
from mewpy.mew.solution import ModelSolution
from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel
from mewpy.solvers import Solution
from mewpy.solvers.solver import Solver


class SRFBA(FBA, LinearMixIn):

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = False,
                 attach: bool = False):
        """
        Steady-state Regulatory Flux Balance Analysis (SRFBA) of a metabolic-regulatory model.
        Implementation of a steady-state version of SRFBA for an integrated metabolic-regulatory model.

        For more details consult Shlomi et al. 2007 at https://dx.doi.org/10.1038%2Fmsb4100141

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

    def _build_algebraic_constraints(self):
        """
        It builds the algebraic constraints for SRFBA, namely the GPR constraints and the interaction constraints.
        It is called automatically upon instantiation if build is True.
        :return:
        """
        variables = []
        constraints = []

        for reaction in self.model.yield_reactions():
            gpr_variables, gpr_constraints = self.gpr_constraint(reaction)
            variables.extend(gpr_variables)
            constraints.extend(gpr_constraints)

        for interaction in self.model.yield_interactions():
            interaction_variables, interaction_constraints = self.interaction_constraint(interaction)
            variables.extend(interaction_variables)
            constraints.extend(interaction_constraints)

        self.add_variables(*variables)
        self.add_constraints(*constraints)

    def _build(self):
        """
        It builds the linear problem for SRFBA. It is called automatically upon instantiation if build is True.
        The SRFBA problem is a mixed-integer linear problem (MILP) with the following structure:
            - metabolic constraints
            - GPR constraints
            - interaction constraints

        :return:
        """
        if self.model.is_metabolic() and self.model.is_regulatory():
            self._build_mass_constraints()
            self._build_algebraic_constraints()

            self._linear_objective = {var.id: value for var, value in self.model.objective.items()}
            self._minimize = False

    def _optimize(self,
                  to_solver: bool = False,
                  solver_kwargs: Dict = None,
                  initial_state: Dict[str, float] = None,
                  **kwargs) -> Union[ModelSolution, Solution]:
        """
        It solves the linear problem. The linear problem is solved using the solver interface.

        The optimize method allows setting temporary changes to the linear problem. The changes are
        applied to the linear problem reverted to the original state afterward.
        Objective, constraints and solver parameters can be set temporarily.

        The solution is returned as a ModelSolution instance, unless to_solver is True. In this case,
        the solution is returned as a SolverSolution instance.

        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False
        :param solver_kwargs: A dictionary of solver parameters to be set temporarily. Default: None
        :param initial_state: a dictionary of variable ids and their values to set as initial state
        :return: A ModelSolution instance or a SolverSolution instance if to_solver is True.
        """
        if not initial_state:
            initial_state = {}

        if not solver_kwargs:
            solver_kwargs = {}

        if 'constraints' not in solver_kwargs:
            solver_kwargs['constraints'] = initial_state.copy()
        else:
            solver_kwargs['constraints'].update(initial_state)

        solution = self.solver.solve(**solver_kwargs)
        return solution
