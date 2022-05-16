from mewpy.simulation.kinetic import KineticSimulation
import unittest

MODELS_PATH = 'tests/data/'
MODEL = MODELS_PATH + 'chassagnole2002.xml'


class TestKineticSimulation(unittest.TestCase):

    def setUp(self):
        from mewpy.io.sbml import load_ODEModel
        self.model = load_ODEModel(MODEL)

    def test_build_ode(self):
        self.model.build_ode()

    def test_simulation(self):
        from mewpy.simulation.kinetic import KineticSimulation
        sim = KineticSimulation(self.model)
        sim.simulate()
