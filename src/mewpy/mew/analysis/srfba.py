from typing import Union, Dict, Tuple

from mewpy.solvers import Solution
from mewpy.solvers.solver import Solver
from mewpy.mew.lp import Notification, MetabolicLinearizer, InteractionLinearizer, GPRLinearizer
from mewpy.model import Model, MetabolicModel, RegulatoryModel
from mewpy.mew.solution import ModelSolution


class SRFBA(MetabolicLinearizer, GPRLinearizer, InteractionLinearizer):

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = True,
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
        If none, a new solver is instantiated. An instantiated solver may be used but it will be overwritten
        if build is true.

        :param build: Whether to build the linear problem upon instantiation. Default: False
        :param attach: Whether to attach the linear problem to the model upon instantiation. Default: False
        """
        super().__init__(model=model,
                         solver=solver,
                         build=build,
                         attach=attach)

    def build(self):
        """
        It builds the linear problem for SRFBA. It is called automatically upon instantiation if build is True.
        The SRFBA problem is a mixed-integer linear problem (MILP) with the following structure:
            - metabolic constraints
            - GPR constraints
            - interaction constraints

        :return:
        """
        if self._variables or self._constraints:
            self.clean()

        if self.model.is_metabolic():
            self.metabolite_reaction_lookup(reactions=list(self.model.yield_reactions()),
                                            metabolites=list(self.model.yield_metabolites()))

            self.add_gprs(list(self.model.yield_reactions()))

        if self.model.is_regulatory():
            self.add_interactions(list(self.model.yield_interactions()))

        self.update()

        self.set_objective(linear={var.id: value for var, value in self.model.objective.items()},
                           minimize=False)

    def notification(self, notification: Notification):
        """
        It is called by the model when a notification is sent. It updates the linear problem accordingly.
        :param notification: A Notification instance
        :return:
        """
        if notification.content_type == 'coefficients' and notification.action == 'set':

            # if the notification is a coefficient object, the whole problem must be build again,
            # as lower and upper bounds of variables are also encoded into the constraints

            self.update_milp_constraints(notification.content.id)

            lb, ub = notification.content.bounds

            return super(SRFBA, self).set_bounds(notification.content.id, lb=lb, ub=ub)

        else:
            return super(SRFBA, self).notification(notification)

    def update_milp_constraints(self, variable_id):
        """
        It updates the MILP constraints for a given variable. It is called when a coefficient is set.
        :param variable_id: The variable id
        :return:
        """
        variable = self.model.get(variable_id)

        to_add = []

        if variable.is_regulator():

            for interaction in variable.yield_interactions():
                to_add.append(interaction)

        if variable.is_target():
            to_add.append(variable.interaction)

        self.add_interactions(to_add)

        if variable.is_reaction():
            self.add_gprs([variable])

    def _optimize_with_constraints(self,
                                   objective=None,
                                   minimize=False,
                                   constraints=None,
                                   get_values=True,
                                   shadow_prices=False,
                                   reduced_costs=False,
                                   pool_size=0,
                                   pool_gap=None):

        constraints = constraints.copy()

        old_simulators = self.model.simulators.copy()

        self.model._simulators = [self]

        with self.model:

            for variable, bounds in constraints.items():
                variable_obj = self.model.get(variable)

                # sending notification. Notification processing will update variables and constraints accordingly
                if isinstance(bounds, (tuple, list, set)):
                    variable_obj.coefficient.coefficients = bounds

                else:
                    variable_obj.coefficient.coefficients = (bounds,)

            # updating the last modifications
            self.update()

            # solving for the last modifications
            solver_solution = self._solver.solve(linear=objective,
                                                 quadratic=None,
                                                 minimize=minimize,
                                                 model=None,
                                                 constraints=None,
                                                 get_values=get_values,
                                                 shadow_prices=shadow_prices,
                                                 reduced_costs=reduced_costs,
                                                 pool_size=pool_size,
                                                 pool_gap=pool_gap)

        # context exit. restore. update
        self.update()

        self.model._simulators = old_simulators

        return solver_solution

    def optimize(self,
                 initial_state: Union[Dict[str, float], Dict[str, Tuple]] = None,
                 objective: Union[str, Dict[str, float]] = None,
                 minimize: bool = False,
                 constraints: Dict[str, Tuple[float, float]] = None,
                 to_solver: bool = False,
                 get_values: bool = True,
                 shadow_prices: bool = False,
                 reduced_costs: bool = False,
                 pool_size: int = 0,
                 pool_gap: float = None) -> Union[ModelSolution, Solution]:
        """
        It solves the linear problem. The linear problem is solved using the solver interface.

        The optimize method allows setting temporary changes to the linear problem. The changes are
        applied to the linear problem reverted to the original state afterward.
        Objective, constraints and solver parameters can be set temporarily.

        The solution is returned as a ModelSolution instance, unless to_solver is True. In this case,
        the solution is returned as a SolverSolution instance.

        :param initial_state: a dictionary of variable ids and their values to set as initial state
        :param objective: A dictionary with the objective coefficients. The keys are the variable ids and the values
        are the coefficients. Alternatively, a string with the variable id can be provided.
        In this case, the coefficient is set to 1.0. Default: None
        :param minimize: Whether to minimize the objective. Default: False
        :param constraints: A dictionary with the constraints bounds. The keys are the constraint ids and the values
        are tuples with the lower and upper bounds. Default: None
        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False
        :param get_values: Whether to retrieve the solution values. Default: True
        :param shadow_prices: Whether to retrieve the shadow prices. Default: False
        :param reduced_costs: Whether to retrieve the reduced costs. Default: False
        :param pool_size: The size of the solution pool. Default: 0
        :param pool_gap: The gap between the best solution and the worst solution in the pool. Default: None
        :return: A ModelSolution instance or a SolverSolution instance if to_solver is True.
        """

        self._set_objective(linear=objective, minimize=minimize)

        if not initial_state:
            initial_state = {}

        if not constraints:
            constraints = {}

        constraints = {**constraints, **initial_state}

        if constraints:
            solver_solution = self._optimize_with_constraints(objective=objective,
                                                              minimize=minimize,
                                                              constraints=constraints,
                                                              get_values=get_values,
                                                              shadow_prices=shadow_prices,
                                                              reduced_costs=reduced_costs,
                                                              pool_size=pool_size,
                                                              pool_gap=pool_gap)

        else:
            solver_solution = self._solver.solve(linear=objective,
                                                 quadratic=None,
                                                 minimize=minimize,
                                                 model=None,
                                                 constraints=None,
                                                 get_values=get_values,
                                                 shadow_prices=shadow_prices,
                                                 reduced_costs=reduced_costs,
                                                 pool_size=pool_size,
                                                 pool_gap=pool_gap)

        if to_solver:
            return solver_solution

        if minimize:
            sense = 'minimize'

        else:
            sense = 'maximize'

        return ModelSolution.from_solver(method=self.method, solution=solver_solution, model=self.model,
                                         objective_direction=sense)
