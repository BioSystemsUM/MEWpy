"""odespy ODE Solver
"""
from .ode import ODEMethod, SolverConfigurations, ODESolver
import numpy as np
import odespy

methods = {
    ODEMethod.LSODA: odespy.Lsoda,
    ODEMethod.LSODAR: odespy.Lsodar,
    ODEMethod.LSODE: odespy.Lsode,
    ODEMethod.HEUN: odespy.Heun,
    ODEMethod.EULER: odespy.Euler,
    ODEMethod.RK4: odespy.RK4,
    ODEMethod.DORMAN_PRINCE: odespy.DormandPrince,
    ODEMethod.RKFehlberg: odespy.RKFehlberg,
    ODEMethod.Dopri5: odespy.Dopri5,
    # ODEMethod.Dop853: odespy.Dop853,
    ODEMethod.Vode: odespy.Vode,
    ODEMethod.Radau: odespy.Radau5,
    ODEMethod.AdamsBashMoulton2: odespy.AdamsBashMoulton2,
    ODEMethod.AdamsBashforth2: odespy.AdamsBashforth2
}


class ODESpySolver(ODESolver):
    """
    ODE solver method implemented on odespy package.
    """

    def __init__(self, func, method):
        self.func = func
        if method in methods.keys():
            self.method = method
        else:
            raise ValueError(f'Method {method} is unavailable.')

        self.initial_condition = None

    def set_initial_condition(self, initial_condition):
        self.initial_condition = initial_condition

    def solve(self, y0, t_points, **kwargs):
        """
        Solves the ODE
        """
        def f(u, t):
            return self.func(t, u)

        try:
            if self.method == ODEMethod.AdamsBashforth2:
                solver = methods[self.method](f, method='bdf')
            else:
                solver = methods[self.method](f)
            
            # update default parameters
            time_points = t_points
            solver.atol = SolverConfigurations.ABSOLUTE_TOL
            solver.rtol = SolverConfigurations.RELATIVE_TOL
            solver.set_initial_condition(y0)
            y, t = solver.solve(time_points)
            C = [c[-1] for c in y]
            return C, t, y
        except Exception as e:
            print(e)
            raise(Exception)
