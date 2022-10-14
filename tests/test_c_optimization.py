import unittest
from pathlib import Path

MODELS_PATH = 'tests/data/'
EC_CORE_MODEL = MODELS_PATH + 'e_coli_core.xml.gz'
EC_CORE_MODEL2 = Path(__file__).parent.joinpath('data', 'e_coli_core.xml')
BIOMASS_ID = 'R_BIOMASS_Ecoli_core_w_GAM'
SUCC_ID = 'R_EX_succ_e'
MIN_GROWTH = 0.1


class TestOptInspyred(unittest.TestCase):
    """ Unittests of Inspyred based optimizations.
    """

    def setUp(self):
        """Sets up the the model
        """
        from reframed.io.sbml import load_cbmodel
        self.model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.optimization.settings import set_default_population_size
        set_default_population_size(10)
        from mewpy.optimization import set_default_engine, get_available_engines
        if len(get_available_engines()):
            set_default_engine('inspyred')

    def test_engine(self):
        """Assert the availability of optimization engines
        """
        from mewpy.optimization import get_available_engines
        eng = get_available_engines()
        self.assertGreater(len(eng), 0)

    def test_KOProblem(self):
        """Tests KO problems
        """
        from mewpy.optimization.evaluation import BPCY, WYIELD
        f1 = BPCY(BIOMASS_ID, SUCC_ID, method='lMOMA')
        f2 = WYIELD(BIOMASS_ID, SUCC_ID)
        from mewpy.problems import RKOProblem
        problem = RKOProblem(self.model, [f1, f2], max_candidate_size=6)
        from mewpy.optimization import EA
        ea = EA(problem, max_generations=2)
        ea.run()
        self.assertEqual(ea.get_population_size(), 10)

    def test_OUProblem(self):
        """Tests OU problems
        """
        from mewpy.optimization.evaluation import BPCY_FVA, TargetFlux, ModificationType
        f1 = BPCY_FVA(BIOMASS_ID, SUCC_ID, method='lMOMA')
        f2 = TargetFlux(SUCC_ID)
        f3 = ModificationType()
        from mewpy.problems import ROUProblem
        problem = ROUProblem(self.model, [f1, f2, f3], max_candidate_size=6)
        from mewpy.optimization import EA
        ea = EA(problem, max_generations=1)
        ea.run()
        self.assertEqual(ea.get_population_size(), 10)


class TestOptJMetal(TestOptInspyred):
    """ Unittests for JMetalPy based optimizations.
    """

    def setUp(self):
        """Sets up the the model
        """
        from reframed.io.sbml import load_cbmodel
        self.model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.optimization.settings import set_default_population_size
        set_default_population_size(10)
        from mewpy.optimization import set_default_engine, get_available_engines
        if len(get_available_engines()):
            set_default_engine('jmetal')


class TestGERMOptInspyred(unittest.TestCase):
    """ Unittests for Inspyred based optimizations using germ models.
    """

    def setUp(self):
        """Sets up the the model
        """
        from mewpy.io import read_sbml
        self.model = read_sbml(EC_CORE_MODEL2, regulatory=False, warnings=False)
        from mewpy.optimization.settings import set_default_population_size
        set_default_population_size(10)
        from mewpy.optimization import set_default_engine, get_available_engines
        if len(get_available_engines()):
            set_default_engine('inspyred')

    def test_engine(self):
        """Assert the availability of optimization engines
        """
        from mewpy.optimization import get_available_engines
        eng = get_available_engines()
        self.assertGreater(len(eng), 0)

    def test_KOProblem(self):
        """Tests KO problems
        """
        from mewpy.optimization.evaluation import BPCY, WYIELD
        f1 = BPCY(BIOMASS_ID, SUCC_ID, method='fba')
        f2 = WYIELD(BIOMASS_ID, SUCC_ID, method='fba')
        from mewpy.problems import RKOProblem
        problem = RKOProblem(self.model, [f1, f2], max_candidate_size=6)
        from mewpy.optimization import EA
        ea = EA(problem, max_generations=2)
        ea.run()
        ea.dataframe()
        self.assertEqual(ea.get_population_size(), 10)

    def test_OUProblem(self):
        """Tests OU problems
        """
        from mewpy.optimization.evaluation import BPCY_FVA, TargetFlux, ModificationType
        f1 = BPCY_FVA(BIOMASS_ID, SUCC_ID)
        f2 = TargetFlux(SUCC_ID, method='fba')
        f3 = ModificationType()
        from mewpy.problems import ROUProblem
        problem = ROUProblem(self.model, [f1, f2, f3], max_candidate_size=6)
        from mewpy.optimization import EA
        ea = EA(problem, max_generations=1)
        ea.run()
        self.assertEqual(ea.get_population_size(), 10)


class TestGERMOptJMetal(TestGERMOptInspyred):
    """ Unittests for JMetalPy based optimizations using germ models.
    """

    def setUp(self):
        """Sets up the the model
        """
        from mewpy.io import read_sbml
        self.model = read_sbml(EC_CORE_MODEL2, regulatory=False, warnings=False)
        from mewpy.optimization.settings import set_default_population_size
        set_default_population_size(10)
        from mewpy.optimization import set_default_engine, get_available_engines
        if len(get_available_engines()):
            set_default_engine('jmetal')


if __name__ == '__main__':
    unittest.main()
