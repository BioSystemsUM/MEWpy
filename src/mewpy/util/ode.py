"""ODE Solver
"""
from enum import Enum


class SolverConfigurations:
    ABSOLUTE_TOL = 1e-9
    RELATIVE_TOL = 1e-6
    N_STEPS = 10000


class SolverMethod(Enum):
    LSODA = 1
    LSODAR = 2
    LSODE = 3
    HEUN = 4
    EULER = 5
    RK4 = 6
    DORMAN_PRINCE = 7
    RKFehlberg = 8
    Dopri5 = 9
    Dop853 = 10
    Vode = 11
    Radau5 = 12
    AdamsBashforth2 = 13
    AdamsBashMoulton2 = 14


class KineticConfigurations:
    STEADY_STATE_TIME = 1e9
    SOLVER = "odespy"
    SOLVER_METHOD = SolverMethod.HEUN  # ode solver method used in the phenotype simulation
    SOLVER_TIMEOUT = 6000  # maximum time allowed by simulation


class ODEStatus(Enum):
    ERROR = 0
    OPTIMAL = 1


class ODESpySolver:
    """
    ODE solver method implemented on odespy package.
    """

    def __init__(self, solverMethod):
        self.solverMethod = solverMethod

    def get_solver(self, func):
        """
        Returns the solver method from odespy package.

        :param func: function with ODE system.
        :return: an instance of odeSolver

        """
        try:
            import odespy
        except Exception:
            raise RuntimeError("ODEspy package is required")

        if self.solverMethod is SolverMethod.LSODA:
            solver = odespy.Lsoda(func)
        elif self.solverMethod is SolverMethod.LSODAR:
            solver = odespy.Lsodar(func)
        elif self.solverMethod is SolverMethod.LSODE:
            solver = odespy.Lsode(func)
        elif self.solverMethod is SolverMethod.HEUN:
            solver = odespy.Heun(func)
        elif self.solverMethod is SolverMethod.EULER:
            solver = odespy.Euler(func)
        elif self.solverMethod is SolverMethod.RK4:
            solver = odespy.RK4(func)
        elif self.solverMethod is SolverMethod.DORMAN_PRINCE:
            solver = odespy.DormandPrince(func)
        elif self.solverMethod is SolverMethod.RKFehlberg:
            solver = odespy.RKFehlberg(func)
        elif self.solverMethod is SolverMethod.Dopri5:
            solver = odespy.Dopri5(func)
        elif self.solverMethod is SolverMethod.Dop853:
            solver = odespy.Dop853(func)
        elif self.solverMethod is SolverMethod.Vode:
            solver = odespy.Vode(func)
        elif self.solverMethod is SolverMethod.AdamsBashforth2:
            solver = odespy.AdamsBashforth2(func, method='bdf')
        elif self.solverMethod is SolverMethod.Radau5:
            solver = odespy.Radau5(func)
        elif self.solverMethod is SolverMethod.AdamsBashMoulton2:
            solver = odespy.AdamsBashMoulton2(func)

        # update default parameters
        solver.nsteps = SolverConfigurations.N_STEPS
        solver.atol = SolverConfigurations.ABSOLUTE_TOL
        solver.rtol = SolverConfigurations.RELATIVE_TOL

        return solver

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
