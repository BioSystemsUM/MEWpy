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
from scikits.odes import ode
import numpy as np
import warnings

methods = {
    ODEMethod.BDF : 'cvode',
    
}


class ScikitsODESolver(ODESolver):
    """
    ODE solver method implemented on odespy package.
    """

    def __init__(self, func, method):
        """
          Integrate a system of ordinary differential equations.\n,
          *odeint* is a wrapper around the ode class, as a confenience function to,
          quickly integrate a system of ode.
          Solves the initial value problem for stiff or non-stiff systems,
          of first order ode's:,
              rhs = dy/dt = fun(t, y),
          where y can be a vector, then rhsfun must be a function computing rhs with\n",
          signature:,
              rhsfun(t, y, rhs)",
          storing the computated dy/dt in the rhs array passed to the function.
          All sundials methods, except the Runge-Kutta method, are implicit.
        """
            
        self.func = func
        
        # makes the function implicit
        def rhs(t,y,out):
            res = self.func(t,y)
            for i, v in enumerate(res):
                out[i]=v
        
        if method in methods.keys():
            self.method = method
        else:
            self.method='cvode'
            
        self.initial_condition = None
        self.solver = ode(self.method,rhs)

    def set_initial_condition(self, initial_condition):
        self.initial_condition = initial_condition

    def solve(self, y0, t_points, **kwargs):
        """
        Returns the solver method from package.

        :param func: function with ODE system.
        :return: an instance of odeSolver

        """
        sol = self.solver.solve(t_points,y0)
        array = np.array(sol.values.y)
        transposed_array = array.T
        y = transposed_array.tolist()
        C = sol.values.y[-1]
        t = sol.values.t
        
        return C, t , y

        