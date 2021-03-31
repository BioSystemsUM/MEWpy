from mewpy.solvers.solution import Status
from mewpy.util.constants import ModelConstants


def decode_solver_solution(solution, minimize=True, status=False):

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

    if status:
        return sol_f_obj, sol_status

    return sol_f_obj
