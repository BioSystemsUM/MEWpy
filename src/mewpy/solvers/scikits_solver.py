"""scikits ODE Solver
"""
from .ode import ODEMethod, SolverConfigurations, ODESolver
from scikits.odes.odeint import odeint
import warnings
methods = {
}


class ScikitsODESolver(ODESolver):
    """
    ODE solver method implemented on odespy package.
    """

    def __init__(self, func, method):
        self.func = func
        if method in methods.keys():
            self.method = method
        else:
            warnings.warn(f'Method {method} is unavailable.')

        self.initial_condition = None

    def set_initial_condition(self, initial_condition):
        self.initial_condition = initial_condition

    def solve(self, y0, t_span, **kwargs):
        """
        Returns the solver method from odespy package.

        :param func: function with ODE system.
        :return: an instance of odeSolver

        """
        sol = odeint(self.func, t_span, y0)
        C = [c[-1] for c in sol.y]
        return C, sol.t, sol.y
