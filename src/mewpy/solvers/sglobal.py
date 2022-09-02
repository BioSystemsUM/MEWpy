from mewpy.util.utilities import Singleton


class MEWPYSolvers(Singleton):

    def __init__(self):
        self._mewpy_solvers = dict()

    def build(self):
        try:
            from .gurobi_solver import GurobiSolver
            self._mewpy_solvers['gurobi'] = GurobiSolver
        except ImportError:
            pass

        try:
            from .cplex_solver import CplexSolver
            self._mewpy_solvers['cplex'] = CplexSolver
        except ImportError:
            pass

        try:
            from .optlang_solver import OptLangSolver
            self._mewpy_solvers['optlang'] = OptLangSolver
        except ImportError:
            pass

    def get_solvers(self):
        if not self._mewpy_solvers:
            self.build()
        return self._mewpy_solvers


class MEWPYODESolvers(Singleton):

    def __init__(self):
        self._mewpy_ode_solvers = dict()

    def build(self):
        try:
            from .scikits_solver import ScikitsODESolver
            self._mewpy_ode_solvers['scikits'] = ScikitsODESolver
        except ImportError:
            pass

        try:
            from .scipy_solver import ScipySolver
            self._mewpy_ode_solvers['scipy'] = ScipySolver
        except ImportError:
            pass

        try:
            from .odespy_solver import ODESpySolver
            self._mewpy_ode_solvers['odespy'] = ODESpySolver
        except ImportError:
            pass

    def get_solvers(self):
        if not self._mewpy_ode_solvers:
            self.build()
        return self._mewpy_ode_solvers


__MEWPY_solvers__ = MEWPYSolvers().get_solvers()
__MEWPY_ode_solvers__ = MEWPYODESolvers().get_solvers()




