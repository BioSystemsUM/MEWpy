from mewpy.problems.problem import AbstractKOProblem, AbstractOUProblem
from mewpy.optimization.evaluation import KineticEvaluationFunction
from mewpy.simulation.kinetic import KineticSimulation
from mewpy.solvers import KineticConfigurations
from re import search

class KineticKOProblem(AbstractKOProblem):

    def __init__(self, model, fevaluation=None, **kwargs):
        super(KineticKOProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        kinetic_parameters = kwargs.get('kparam', None)
        t_points = kwargs.get('t_points', None)
        timeout = kwargs.get('timeout', KineticConfigurations.SOLVER_TIMEOUT)
        self.kinetic_sim = KineticSimulation(model, 
                                             parameters=kinetic_parameters,
                                             timeout=timeout,
                                             t_points=t_points)
        for f in self.fevaluation:
            if isinstance(f,KineticEvaluationFunction):
                f.kinetic = True
            else:
                raise ValueError(f"The optimization function {f} is not intended for kinetic optimization.")

    def _build_target_list(self):
        p = list(self.model.get_parameters(exclude_compartments=True))
        target =[]
        for k in p:
            if search(r'(?i)[rv]max',k):
                target.append(k)
        if self.non_target is not None:
            target = target - set(self.non_target)
        self._trg_list = list(target)

    def decode(self, candidate):
        factors = {self.target_list[idx]: 0 for idx in candidate}
        return factors

    def evaluate_solution(self, solution, decode=True):
        p = []
        factors = self.decode(solution) if decode else solution
        simulation_results = self.kinetic_sim.simulate(factors=factors)
        for f in self.fevaluation:
            p.append(f(simulation_results, factors))
        return p


class KineticOUProblem(AbstractOUProblem):

    def __init__(self, model, fevaluation=None, **kwargs):
        super(KineticOUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        kinetic_parameters = kwargs.get('kparam', None)
        t_points = kwargs.get('t_points', None)
        timeout = kwargs.get('timeout', KineticConfigurations.SOLVER_TIMEOUT)
        self.kinetic_sim = KineticSimulation(model, 
                                             parameters=kinetic_parameters,
                                             timeout=timeout,
                                             t_points=t_points)
        for f in self.fevaluation:
            if isinstance(f,KineticEvaluationFunction):
                f.kinetic = True
            else:
                raise ValueError(f"The optimization function {f} is not intended for kinetic optimization.")

    def _build_target_list(self):
        p = list(self.model.get_parameters(exclude_compartments=True))
        target =[]
        for k in p:
            if search(r'(?i)max',k):
                target.append(k)
        if self.non_target is not None:
            target = set(target) - set(self.non_target)
        self._trg_list = list(target)

    def decode(self, candidate):
        factors = {self.target_list[idx]: self.levels[lv_idx] for idx, lv_idx in candidate}
        return factors

    def evaluate_solution(self, solution, decode=True):
        p = []
        factors = self.decode(solution) if decode else solution
        simulation_results = self.kinetic_sim.simulate(factors=factors)
        for f in self.fevaluation:
            p.append(f(simulation_results, factors))
        return p
