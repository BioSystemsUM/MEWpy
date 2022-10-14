"""scipy ODE solver"""
from .ode import ODEMethod, ODESolver
from scipy.integrate._ivp import solve_ivp

methods = {
    ODEMethod.RK45: 'RK45',
    ODEMethod.RK23: 'RK23',
    ODEMethod.DOP853: 'DOP853',
    ODEMethod.Radau: 'Radau',
    ODEMethod.BDF: 'BDF',
    ODEMethod.LSODA: 'LSODA',
}


class ScipySolver(ODESolver):

    def __init__(self, func, method='LSODA'):
        self.func = func
        self.method = method
        self.initial_condition = None

    def set_initial_condition(self, initial_condition):
        self.initial_condition = initial_condition

    def solve(self, y0, t_points, **kwargs):
        t_span=[t_points[0],t_points[-1]]
        sol = solve_ivp(self.func, t_span, y0, method=methods[self.method], t_eval=t_points,**kwargs)
        C = [c[-1] for c in sol.y]
        return C, sol.t, sol.y
