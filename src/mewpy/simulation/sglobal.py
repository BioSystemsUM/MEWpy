from mewpy.util.utilities import Singleton


class MEWPYSimSolvers(Singleton):

    def __init__(self):
        self._mewpy_sim_solvers = []

    def build(self):
        try:
            import gurobipy
            self._mewpy_sim_solvers.append('gurobi')
        except ImportError as e:
            pass
        try:
            import cplex
            self._mewpy_sim_solvers.append('cplex')
        except ImportError as e:
            pass

        try:
            import swiglpk
            self._mewpy_sim_solvers.append('glpk')
        except ImportError:
            pass

    def get_solvers(self):
        if not self._mewpy_sim_solvers:
            self.build()
        return self._mewpy_sim_solvers


__MEWPY_sim_solvers__ = MEWPYSimSolvers().get_solvers()
