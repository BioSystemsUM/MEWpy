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
Interface for scikits ODE Solver

Author Vitor Pereira
##############################################################################
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
