from typing import TYPE_CHECKING, Tuple, Union, Dict

from mewpy.solvers.solution import Status
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.solvers import Solution
    from mewpy.mew.lp import LinearProblem


def decode_solver_solution(solution: 'Solution') -> Tuple[float, str]:
    """
    It decodes the solution of the solver and returns the objective value.

    :param solution: the solution of the solver
    :return: the objective value and the status of the solution
    """
    sol_status = solution.status.value

    if solution.status == Status.OPTIMAL:
        sol_f_obj = solution.fobj

    elif solution.status == Status.UNBOUNDED or solution.status == Status.INF_OR_UNB:

        sol_f_obj = ModelConstants.REACTION_UPPER_BOUND

    else:
        sol_f_obj = None

    return sol_f_obj, sol_status


def run_method_and_decode(method: 'LinearProblem',
                          objective: Union[str, Dict[str, float]] = None,
                          constraints: Dict[str, Tuple[float, float]] = None,
                          **kwargs) -> Tuple[float, str]:
    """
    It runs a method and decodes the objective value and status returned by the solver.
    :param method: the method to be run
    :param objective: an alternative temporary objective function
    :param constraints: alternative temporary constraints
    :param kwargs: additional arguments to be passed to the method
    :return: the objective value and the status of the solution
    """
    solver_kwargs = {'get_values': False}

    if objective:
        if hasattr(objective, 'keys'):
            solver_kwargs['linear'] = objective.copy()
        else:
            solver_kwargs['linear'] = {str(objective): 1.0}

    if constraints:
        solver_kwargs['constraints'] = constraints

    if 'minimize' in kwargs:
        solver_kwargs['minimize'] = kwargs['minimize']

    solution = method.optimize(to_solver=True, solver_kwargs=solver_kwargs, **kwargs)
    objective_value, status = decode_solver_solution(solution=solution)
    return objective_value, status
