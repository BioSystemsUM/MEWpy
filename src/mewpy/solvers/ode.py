"""An interface for ODE solvers"""
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
