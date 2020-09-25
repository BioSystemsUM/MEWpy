from .problem import AbstractKOProblem, AbstractOUProblem
from mewpy.simulation.kinetic import KineticSimulation
from collections import OrderedDict


class KineticKOProblem(AbstractKOProblem):

    def __init__(self, model, fevaluation=None, **kwargs):
        super(KineticKOProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        kinetic_parameters = kwargs.get('kparam', None)
        tSteps = kwargs.get('tSteps', None)
        self.kinetic_sim = KineticSimulation(model, parameters=kinetic_parameters, tSteps=tSteps)

    def _build_target_list(self):
        target = set(self.model.reactions.keys())
        if self.non_target is not None:
            target = target - set(self.non_target)
        self._trg_list = list(target)

    def decode(self, candidate):
        factors = OrderedDict([(self.target_list[idx], 0) for idx in candidate])
        return factors

    def evaluate_solution(self, solution, decode=True):
        """
        overrides evaluate_solution
        """
        p = []
        factors = self.decode(solution) if decode else solution
        simulation_results = self.kinetic_sim.simulate(factors=factors)
        for f in self.fevaluation:
            p.append(f(simulation_results, solution))
        return p


class KineticOUProblem(AbstractOUProblem):

    def __init__(self, model, fevaluation=None, **kwargs):
        super(KineticOUProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)
        kinetic_parameters = kwargs.get('kparam', None)
        tSteps = kwargs.get('tSteps', None)
        self.kinetic_sim = KineticSimulation(model, parameters=kinetic_parameters, tSteps=tSteps)

    def _build_target_list(self):
        target = set(self.model.reactions.keys())
        if self.non_target is not None:
            target = target - set(self.non_target)
        self._trg_list = list(target)

    def decode(self, candidate):
        factors = {self.target_list[idx]: self.levels[lv_idx] for idx, lv_idx in candidate}
        return factors

    def evaluate_solution(self, solution, decode=True):
        """
        overrides evaluate_solution
        """
        p = []
        factors = self.decode(solution) if decode else solution
        simulation_results = self.kinetic_sim.simulate(factors=factors)
        for f in self.fevaluation:
            p.append(f(simulation_results, solution))
        return p
