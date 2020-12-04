import unittest

MODELS_PATH = 'tests/data/'
EC_CORE_MODEL = MODELS_PATH + 'e_coli_core.xml.gz'
MIN_GROWTH = 0.1

      

class TestReframedSimul(unittest.TestCase):

    def setUp(self):
        from reframed.io.sbml import load_cbmodel
        model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.simulation import get_simulator
        self.simul = get_simulator(model)
        self.BIOMASS_ID = model.biomass_reaction

    def test_essential_reactions(self):
        essential = self.simul.essential_reactions
        self.assertGreater(len(essential), 0)

    def test_essential_genes(self):
        essential = self.simul.essential_genes
        self.assertGreater(len(essential), 0)

    def test_fba(self):
        res=self.simul.simulate()
        self.assertGreater(res.objective_value, MIN_GROWTH)

    def test_pfba(self):
        res=self.simul.simulate(method='pFBA')
        self.assertGreater(res.fluxes[self.BIOMASS_ID], MIN_GROWTH)

    def test_moma(self):
        res=self.simul.simulate(method='MOMA')
        self.assertGreater(res.fluxes[self.BIOMASS_ID], MIN_GROWTH)

    def test_lmoma(self):
        res=self.simul.simulate(method='lMOMA')
        self.assertGreater(res.fluxes[self.BIOMASS_ID], MIN_GROWTH)

    def test_room(self):
        res=self.simul.simulate(method='ROOM')
        self.assertGreater(res.fluxes[self.BIOMASS_ID], MIN_GROWTH)


class TestCobra(TestReframedSimul):

    def setUp(self):
        from cobra.io.sbml import read_sbml_model
        model = read_sbml_model(EC_CORE_MODEL)
        from mewpy.simulation import get_simulator
        self.simul = get_simulator(model)
        k = list(self.simul.objective.keys())
        self.BIOMASS_ID = k[0]

    

class TestGeckoLoad(unittest.TestCase):

    def test_gecko(self):
        from mewpy.model.gecko import GeckoModel
        GeckoModel('single-pool')

    def test_simulator(self):
        from mewpy.model.gecko import GeckoModel
        model = GeckoModel('single-pool')
        from mewpy.simulation import get_simulator
        get_simulator(model)

        

class TestGeckoSimul(unittest.TestCase):

    def setUp(self):
        from mewpy.model.gecko import GeckoModel
        model = GeckoModel('single-pool')
        from mewpy.simulation import get_simulator
        self.simul=get_simulator(model)

    def test_essential_proteins(self):
        essential = self.simul.essential_proteins()
        self.assertGreater(len(essential), 0)



if __name__ == '__main__':
    unittest.main()
