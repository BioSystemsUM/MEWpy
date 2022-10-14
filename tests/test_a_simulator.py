import unittest
from pathlib import Path

MODELS_PATH = 'tests/data/'
EC_CORE_MODEL = MODELS_PATH + 'e_coli_core.xml.gz'
EC_CORE_MODEL2 = Path(__file__).parent.joinpath('data', 'e_coli_core.xml')
MIN_GROWTH = 0.1


class TestReframedSimul(unittest.TestCase):
    """ Tests the REFRAMED Simulator
    """

    def setUp(self):
        """Set up
        Loads a model
        """
        from reframed.io.sbml import load_cbmodel
        model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.simulation import get_simulator
        self.simul = get_simulator(model)
        self.BIOMASS_ID = model.biomass_reaction
        self.SUCC = 'R_EX_succ_e'

    def test_essential_reactions(self):
        """Tests essential reactions
        """
        essential = self.simul.essential_reactions()
        self.assertGreater(len(essential), 0)

    def test_essential_genes(self):
        """Tests essential genes
        """
        essential = self.simul.essential_genes()
        self.assertGreater(len(essential), 0)

    def test_uptake_reactions(self):
        """Tests uptake reactions
        """
        uptake_reactions = self.simul.get_uptake_reactions()
        self.assertGreater(len(uptake_reactions), MIN_GROWTH)

    def test_transport_reactions(self):
        """Tests transport reactions
        """
        transport_reactions = self.simul.get_transport_reactions()
        self.assertGreater(len(transport_reactions), MIN_GROWTH)

    def test_fba(self):
        """Tests FBA
        """
        res = self.simul.simulate()
        self.assertGreater(res.objective_value, MIN_GROWTH)

    def test_pfba(self):
        """Tests pFBA
        """
        res = self.simul.simulate(method='pFBA')
        self.assertGreater(res.fluxes[self.BIOMASS_ID], MIN_GROWTH)

    def test_moma(self):
        """Tests MOMA
        """
        res = self.simul.simulate(method='MOMA')
        self.assertGreater(res.fluxes[self.BIOMASS_ID], MIN_GROWTH)

    def test_lmoma(self):
        """Tests lMOMA
        """
        res = self.simul.simulate(method='lMOMA')
        self.assertGreater(res.fluxes[self.BIOMASS_ID], MIN_GROWTH)

    def test_room(self):
        """Tests ROOM
        """
        res = self.simul.simulate(method='ROOM')
        self.assertGreater(res.fluxes[self.BIOMASS_ID], MIN_GROWTH)

    def test_FVA(self):
        self.simul.FVA(reactions=[self.SUCC])

    def test_envelope(self):
        from mewpy.visualization.envelope import plot_flux_envelope
        plot_flux_envelope(self.simul, self.BIOMASS_ID, self.SUCC)

    def test_solver(self):
        from mewpy.solvers import solver_instance
        solver = solver_instance(self.simul)
        solver.solve()


class TestCobra(TestReframedSimul):
    """Tests COBRApy Simulator
    """

    def setUp(self):
        """Set up
        Loads a model
        """
        from cobra.io.sbml import read_sbml_model
        model = read_sbml_model(EC_CORE_MODEL)
        from mewpy.simulation import get_simulator
        self.simul = get_simulator(model)
        k = list(self.simul.objective.keys())
        self.BIOMASS_ID = k[0]
        self.SUCC = 'EX_succ_e'


class TestGERM(TestReframedSimul):
    """Tests GERM Simulator
    """

    def setUp(self):
        """Set up
        Loads a model
        """
        from mewpy.io import read_sbml
        model = read_sbml(EC_CORE_MODEL2, regulatory=False)
        from mewpy.simulation import get_simulator
        self.simul = get_simulator(model)
        k = list(self.simul.objective.keys())
        self.BIOMASS_ID = k[0]
        self.SUCC = 'EX_succ_e'

    def test_essential_reactions(self):
        """Tests essential reactions
        """
        essential = self.simul.essential_reactions()
        self.assertGreater(len(essential), 0)

    def test_essential_genes(self):
        """Tests essential genes
        """
        essential = self.simul.essential_genes()
        self.assertGreater(len(essential), 0)

    def test_uptake_reactions(self):
        """Tests uptake reactions
        """
        uptake_reactions = self.simul.get_uptake_reactions()
        self.assertGreater(len(uptake_reactions), MIN_GROWTH)

    def test_transport_reactions(self):
        """Tests transport reactions
        """
        transport_reactions = self.simul.get_transport_reactions()
        self.assertGreater(len(transport_reactions), MIN_GROWTH)

    def test_fba(self):
        """Tests FBA
        """
        res = self.simul.simulate()
        self.assertGreater(res.objective_value, MIN_GROWTH)

    def test_pfba(self):
        """Tests pFBA
        """
        res = self.simul.simulate(method='pFBA')
        self.assertGreater(res.fluxes[self.BIOMASS_ID], MIN_GROWTH)

    def test_moma(self):
        """Tests MOMA
        """
        pass

    def test_lmoma(self):
        """Tests lMOMA
        """
        pass

    def test_room(self):
        """Tests ROOM
        """
        pass

    def test_FVA(self):
        self.simul.FVA(reactions=self.simul.reactions[0:2])

    def test_envelope(self):
        pass

    def test_solver(self):
        from mewpy.solvers import solver_instance
        solver = solver_instance(self.simul)
        solver.solve()


class TestGeckoLoad(unittest.TestCase):
    """Tests GECKO simulator
    """

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
        self.simul = get_simulator(model)

    def test_essential_proteins(self):
        """
        Can not run on community CPLEX
        """
        # essential = self.simul.essential_proteins()
        # self.assertGreater(len(essential), 0)
        pass


if __name__ == '__main__':
    unittest.main()
