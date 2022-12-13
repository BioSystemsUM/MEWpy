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
An interface for ODE solvers

Author: Vitor Pereira
##############################################################################
"""
from enum import Enum
from abc import ABC, abstractmethod


class ODEMethod(Enum):
    RK45 = 'RK45',
    RK23 = 'RK23',
    DOP853 = 'DOP853',
    Radau = 'Radau',
    BDF = 'BDF',
    LSODA = 'LSODA',
    LSODAR = 'LSODAR',
    LSODE = 'LSODE',
    HEUN = 'HEUN',
    EULER = 'EULER',
    RK4 = 'RK4',
    DORMAN_PRINCE = 'DORMAN_PRINCE',
    RKFehlberg = 'RKFehlberg',
    Dopri5 = 'Dopri5',
    Vode = 'Vode',
    CVode = 'cvode'
    Radau5 = 'Radau5',
    Ida ='ida',
    AdamsBashforth2 = 'AdamsBashforth2',
    AdamsBashMoulton2 = 'AdamsBashMoulton2'

    def __eq__(self, other):
        """Overrides equal to enable string name comparison"""
        if isinstance(other, ODEMethod):
            return super().__eq__(other)
        elif isinstance(other, str):
            return self.name == other
        else:
            return False

    def __hash__(self):
        return hash(self.name)


class SolverConfigurations:
    ABSOLUTE_TOL = 1e-9
    RELATIVE_TOL = 1e-6
    N_STEPS = 20


class KineticConfigurations:
    SOLVER_METHOD = ODEMethod.LSODA
    STEADY_STATE_TIME = 1e9
    SOLVER_TIMEOUT = 6000


class ODEStatus(Enum):
    ERROR = 0
    OPTIMAL = 1


class ODESolver(ABC):

    @abstractmethod
    def solve(self, y0, t_points, **kwargs):
        raise NotImplementedError
