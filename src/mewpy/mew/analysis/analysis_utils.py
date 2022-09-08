from typing import TYPE_CHECKING, Tuple

from mewpy.solvers.solution import Status
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.solvers import Solution


def decode_solver_solution(solution: 'Solution', minimize: bool = True) -> Tuple[float, str]:
    """
    It decodes the solution of the solver and returns the objective value.

    :param solution: the solution of the solver
    :param minimize: whether the problem is a minimization problem
    :param status: whether the status of the solution should be returned
    :return: the objective value and the status of the solution
    """
    sol_status = solution.status.value

    if solution.status == Status.OPTIMAL:
        sol_f_obj = solution.fobj

    elif solution.status == Status.UNBOUNDED or solution.status == Status.INF_OR_UNB:

        if minimize:

            sol_f_obj = ModelConstants.REACTION_LOWER_BOUND

        else:

            sol_f_obj = ModelConstants.REACTION_UPPER_BOUND

    else:
        sol_f_obj = None

    return sol_f_obj, sol_status
