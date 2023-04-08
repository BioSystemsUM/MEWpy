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
Author: Vítor Pereira
##############################################################################
"""

from .simulator import get_simulator, get_container
from .simulation import Simulator, SimulationMethod, SStatus, SimulationResult
from .environment import Environment
from .sglobal import __MEWPY_sim_solvers__


default_solver = None


def get_default_solver():
    """
    Returns:
        [type]: [description]
    """
    global default_solver

    if default_solver:
        return default_solver

    solver_order = ['cplex', 'gurobi', 'glpk']

    for solver in solver_order:
        if solver in __MEWPY_sim_solvers__:
            default_solver = solver
            break

    if not default_solver:
        raise RuntimeError("No solver available.")

    return default_solver


def set_default_solver(solvername):
    """ Sets default solver.
    Arguments:
        solvername : (str) solver name (currently available: 'gurobi', 'cplex')
    """

    global default_solver

    if solvername.lower() in __MEWPY_sim_solvers__:
        default_solver = solvername.lower()
    else:
        raise RuntimeError(f"Solver {solvername} not available.")

