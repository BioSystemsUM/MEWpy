# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
##############################################################################
Interface for scipy ODE solver

Author: Vitor Pereira
##############################################################################
"""
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
