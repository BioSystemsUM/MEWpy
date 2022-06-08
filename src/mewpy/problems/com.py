from copy import deepcopy
from warnings import warn
from collections import OrderedDict
from .problem import AbstractKOProblem
from mewpy.model import CommunityModel
from mewpy.simulation import get_simulator
from mewpy.simulation.simulation import Simulator



class CommunityKOProblem(AbstractKOProblem):
    """
    Community Knockout Optimization Problem.

    :param models: A list of metabolic models.
    :param list fevaluation: A list of callable EvaluationFunctions.

    Optional parameters:

    :param OrderedDict envcond: Environmental conditions.
    :param OrderedDict constraints: Additional constraints to be applied to the model.
    :param int candidate_min_size: The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
    :param int candidate_max_size: The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
    :param list target: List of modification target reactions.
    :param list non_target: List of non target reactions. Not considered if a target list is provided.
    :param float scalefactor: A scaling factor to be used in the LP formulation.
    
    """

    def __init__(self, models: list, fevaluation=[], copy_models=True, **kwargs):
        super(CommunityKOProblem, self).__init__(
            None, fevaluation=fevaluation, **kwargs)
        self.organisms = OrderedDict()
        self.model_ids = list({model.id for model in models})

        if len(self.model_ids) < len(models):
            warn("Model ids are not unique, repeated models will be discarded.")

        for model in models:

            m = model if isinstance(model, Simulator) else get_simulator(model)
            self.organisms[m.id] = deepcopy(m) if copy_models else m

    def _build_target_list(self):
        """Target organims, i.e., organisms that may be removed from the community.
        """
        print("Building modification target list.")
        target = set(self.model_ids)
        if self.non_target is not None:
            target = target - set(self.non_target)
        self._trg_list = list(target)

    def solution_to_constraints(self, candidate):
        """Returns a community model that includes all non
        KO organisms.
        TODO: this is a naif approach. Insteas consider using weights
        to turn organism ON/OFF within the community.

        :param candidate: [description]
        :return: [description]
        """
        ko_organisms = list(candidate.keys())
        models = [x for x in self.model_ids if x not in ko_organisms]
        cmodel = CommunityModel(models)
        return cmodel.get_community_model()

    def evaluate_solution(self, solution, decode=True):
        """
        Evaluates a single solution, a community.

        :param solution: The solution to be evaluated (a community model).
        :param decode: If the solution needs to be decoded (convert a list of model ids to a community model).
        :returns: A list of fitness values.
        """
        decoded = {}
        # decoded constraints
        if decode:
            decoded = self.decode(solution)
            cmodel = self.solution_to_constraints(decoded)
        else:
            cmodel = solution
        try:
            p = []
            simulation_results = dict()
            for method in self.methods:
                sim = get_simulator(cmodel, self.environmental_conditions)
                simulation_result = self.simulator.simulate(method=method, scalefactor=self.scalefactor)
                simulation_results[method] = simulation_result
            # apply the evaluation function(s)
            for f in self.fevaluation:
                v = f(simulation_results, decoded,
                      scalefactor=self.scalefactor)
                p.append(v)
        except Exception as e:
            p = []
            for f in self.fevaluation:
                p.append(f.worst_fitness)
        return p


