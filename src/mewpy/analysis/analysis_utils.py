from mewpy.solvers.solution import Status

from mewpy.util import SLIM_LB, SLIM_UB


def decode_solver_solution(solution, minimize=True, status=False):

    sol_status = solution.status.value

    if solution.status == Status.OPTIMAL:
        sol_f_obj = solution.fobj

    elif solution.status == Status.UNBOUNDED or solution.status == Status.INF_OR_UNB:

        if minimize:

            sol_f_obj = SLIM_LB

        else:

            sol_f_obj = SLIM_UB

    else:
        sol_f_obj = None

    if status:
        return sol_f_obj, sol_status

    return sol_f_obj
