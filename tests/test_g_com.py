import unittest

MODELS_PATH = 'tests/data/'
EC_CORE_MODEL = MODELS_PATH + 'e_coli_core.xml.gz'


class TestCommReframed(unittest.TestCase):

    def setUp(self):
        """Set up"""
        from reframed.io.sbml import load_cbmodel
        from mewpy.model import CommunityModel
        model1 = load_cbmodel(EC_CORE_MODEL)
        model1.set_flux_bounds('R_ATPM', 0, 0)
        model1.id = 'm1'

        model2 = model1.copy()
        model2.id = 'm2'

        model3 = model1.copy()
        model3.id = 'm3'
        self.models = [model1, model2, model3]
        self.comm = CommunityModel(self.models, flavor='reframed')

    def FBA(self):
        sim = self.comm.get_community_model()
        res = sim.simulate()
        self.assertGreater(res.objective_value, 0)

    def SteadyCom(self):
        from mewpy.cobra.com.steadycom import SteadyCom
        SteadyCom(self.comm)

    def SteadyComVA(self):
        from mewpy.cobra.com.steadycom import SteadyComVA
        SteadyComVA(self.comm)



class TestCommCobra(TestCommReframed):

    def setUp(self):
        """Set up"""
        from cobra.io.sbml import read_sbml_model
        from mewpy.model import CommunityModel
        model1 = read_sbml_model(EC_CORE_MODEL)
        model1.set_flux_bounds('ATPM', 0, 0)
        model1.id = 'm1'

        model2 = model1.copy()
        model2.id = 'm2'

        model3 = model1.copy()
        model3.id = 'm3'
        self.models = [model1, model2, model3]
        self.comm = CommunityModel(self.models)
