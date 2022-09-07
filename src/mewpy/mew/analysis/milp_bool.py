from typing import Union, Dict, Tuple

from mewpy.solvers.solver import Solver

from mewpy.model import Model, MetabolicModel, RegulatoryModel
from mewpy.solvers.solution import Solution
from mewpy.mew.lp import Notification, InteractionLinearizer, GPRLinearizer
from mewpy.mew.solution import ModelSolution


class milpBool(GPRLinearizer, InteractionLinearizer):

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = True,
                 attach: bool = False):
        """
        Mixed-Integer Boolean simulation (milpBool) of a regulatory model.
        This analysis method solves a set of interactions by linearization using mixed-integer constraints,
        so that regulators are variables of the linear problem to be solved.

        Check the mewpy.mew.lp.linearizers module for more detail regarding interactions linearization

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
        It builds the linear problem by adding variables and constraints to the solver interface.
        :return:
        """
        if self._variables or self._constraints:
            self.clean()

        if self.model.is_regulatory():
            self.add_interactions(self.model.yield_interactions())

        self.update()

    def notification(self, notification: Notification):
        """
        It processes a notification by updating the linear problem accordingly.

        :param notification: a Notification instance
        :return:
        """

        if notification.content_type == 'coefficients' and notification.action == 'set':

            # if the notification is a coefficient object, the whole problem must be build again,
            # as lower and upper bounds of variables are also encoded into the constraints

            self.update_milp_constraints(notification.content.id)

            lb, ub = notification.content.bounds

            return super(milpBool, self).set_bounds(notification.content.id, lb=lb, ub=ub)

        else:
            return super(milpBool, self).notification(notification)

    def update_milp_constraints(self, variable_id):
        """
        It builds a new problem if there is a regulator or target. It is safer and easier to build a new problem rather
        than update the existing one, as lower and upper bounds of variables are also encoded into the constraints.

        This method is called when a notification is received.
        :param variable_id: the id of the variable that has been modified
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

    def _optimize_with_state(self,
                             initial_state: Union[Dict[str, float], Dict[str, Tuple]] = None,
                             objective: Union[str, Dict[str, float]] = None,
                             minimize: bool = False,
                             get_values: bool = True,
                             shadow_prices: bool = False,
                             reduced_costs: bool = False,
                             pool_size: int = 0,
                             pool_gap: float = None):
        """
        It solves the linear problem with a given initial state.
        To use an initial state, the solver must be able to set the initial state of the problem.
        This is performed by temporarily adding the initial state to the model and build a new problem.
        The initial state is then removed from the model.

        Internal use only!!
        :param initial_state: a dictionary of variable ids and their values
        :param objective: the objective function
        :param minimize: whether to minimize or maximize the objective function
        :param get_values: whether to retrieve the values of the variables
        :param shadow_prices: whether to retrieve the shadow prices of the constraints
        :param reduced_costs: whether to retrieve the reduced costs of the variables
        :param pool_size: the size of the solution pool
        :param pool_gap: the gap between the best and the worst solution in the pool
        :return: a ModelSolution instance
        """
        initial_state = initial_state.copy()

        old_simulators = self.model.simulators.copy()

        self.model._simulators = [self]

        with self.model:

            for variable, bounds in initial_state.items():
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

        # exiting the context exit sends a signal to restore the model
        # Then, this linear problem is updated
        self.update()

        # restoring the simulators too
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
        if minimize is None:
            minimize = self._minimize

        self._set_objective(linear=objective, minimize=minimize)

        if initial_state:
            solver_solution = self._optimize_with_state(objective=objective,
                                                        minimize=minimize,
                                                        initial_state=initial_state,
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
