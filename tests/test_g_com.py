import unittest

MODELS_PATH = 'tests/data/'
EC_CORE_MODEL = MODELS_PATH + 'e_coli_core.xml.gz'


class TestComReframed(unittest.TestCase):

    def setUp(self):
        """Set up"""
        from reframed.io.sbml import load_cbmodel
        from mewpy.problems import CommunityKOProblem
        model1 = load_cbmodel(EC_CORE_MODEL)
        model1.set_flux_bounds('R_ATPM', 0, 0)
        model1.id = 'm1'

        model2 = model1.copy()
        model2.id = 'm2'

        model3 = model1.copy()
        model3.id = 'm3'
        self.models = [model1, model2, model3]
        self.problem = CommunityKOProblem(self.models)
        cmodel = self.problem.merge_models()
        # Define the environmental conditions
        from mewpy.simulation.environment import Environment
        medium = Environment.from_model(model1).get_compounds()
        env = Environment.from_compounds(medium, prefix='R_')
        env = env.apply(cmodel, inplace=False)
        self.env = env

    def test_FBA(self):
        from mewpy.simulation import get_simulator
        cmodel = self.problem.merge_models()
        sim = get_simulator(cmodel, envcond=self.env)
        res = sim.simulate()
        self.assertGreater(res.objective_value, 0)

    def test_EA(self):
        from mewpy.optimization.evaluation import BPCY, TargetFlux
        PRODUCT = 'R_EX_succ_e'
        BIOMASS = 'community_growth'
        obj1 = BPCY(BIOMASS, PRODUCT)
        obj2 = TargetFlux(PRODUCT)
        from mewpy.problems import CommunityKOProblem
        problem = CommunityKOProblem(self.models, [obj1, obj2], envcond=self.env, candidate_max_size=2)
        from mewpy.optimization import EA
        ea = EA(problem, max_generations=2)
        final_pop = ea.run()


class TestComCobra(TestComReframed):

    def setUp(self):
        """Set up"""
        from cobra.io.sbml import read_sbml_model
        from mewpy.problems import CommunityKOProblem
        model1 = read_sbml_model(EC_CORE_MODEL)
        model1.reactions.get_by_id("ATPM").lower_bound = 0
        model1.reactions.get_by_id("ATPM").upper_bound = 0
        model1.id = 'm1'
        model2 = model1.copy()
        model2.id = 'm2'
        model3 = model1.copy()
        model3.id = 'm3'
        self.models = [model1, model2, model3]
        # build the problem
        self.problem = CommunityKOProblem(self.models)
        # Define the environmental conditions
        from mewpy.simulation.environment import Environment
        medium = Environment.from_model(model1).get_compounds()
        # environmental conditionc from cobrapy models that drop the 'R_' prefix
        env = Environment.from_compounds(medium, prefix='')
        cmodel = self.problem.merge_models()
        # apply environmental conditions to a reframed model uses the 'R_' prefix
        env = env.apply(cmodel, inplace=False, prefix='R_')
        self.env = env
