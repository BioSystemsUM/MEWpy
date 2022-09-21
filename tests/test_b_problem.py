import unittest
from pathlib import Path

import pytest

MODELS_PATH = 'tests/data/'
EC_CORE_MODEL = MODELS_PATH + 'e_coli_core.xml.gz'

OPTRAM_MODEL = MODELS_PATH + 'yeast_7.6-optram.xml'
OPTRAM_GENES = MODELS_PATH + 'mgene.csv'
OPTRAM_TFS = MODELS_PATH + 'TFnames.csv'
OPTRAM_REGNET = MODELS_PATH + 'regnet.csv'

EC_CORE_MODEL2 = Path(__file__).parent.joinpath('data', 'e_coli_core.xml')
EC_CORE_REG_MODEL = Path(__file__).parent.joinpath('data', 'e_coli_core_trn.csv')

MIN_GROWTH = 0.1

MAX_TRIES = 3


class TestRKOP(unittest.TestCase):

    def setUp(self):
        from reframed.io.sbml import load_cbmodel
        model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.problems import RKOProblem
        self.problem = RKOProblem(model, [])

    def test_targets(self):
        target = self.problem.target_list
        self.assertGreater(len(target), 0)

    def test_generator(self):
        import random
        candidate = self.problem.generator(random)
        self.assertGreater(len(candidate), 0)

    def test_decode(self):
        import random
        candidate = self.problem.generator(random)
        solution = self.problem.decode(candidate)
        n_candidate = self.problem.encode(solution)
        self.assertEqual(candidate, n_candidate)

    def test_to_constraints(self):
        """
        """
        import random
        ispass = False
        tries = 0
        constraints = []
        while not ispass and tries < MAX_TRIES:
            tries += 1
            candidate = self.problem.generator(random)
            solution = self.problem.decode(candidate)
            constraints = self.problem.solution_to_constraints(solution)
            if len(constraints) > 0:
                ispass = True
        self.assertGreater(len(constraints), 0)

    def test_simul_constraints(self):
        import random
        candidate = self.problem.generator(random)
        solution = self.problem.decode(candidate)
        constraints = self.problem.solution_to_constraints(solution)
        self.problem.simulator.simulate(constraints=constraints)


class TestROUP(TestRKOP):

    def setUp(self):
        from reframed.io.sbml import load_cbmodel
        model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.problems import ROUProblem
        self.problem = ROUProblem(model, [])


class TestGKOP(TestRKOP):

    def setUp(self):
        from reframed.io.sbml import load_cbmodel
        model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.problems import GKOProblem
        self.problem = GKOProblem(model, [])


class TestGOUP(TestRKOP):

    def setUp(self):
        from reframed.io.sbml import load_cbmodel
        model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.problems import GOUProblem
        self.problem = GOUProblem(model, [])


class TestOptRAM(TestRKOP):

    def setUp(self):
        from mewpy.problems import OptRamProblem, load_optram
        regnet = load_optram(OPTRAM_GENES, OPTRAM_TFS, OPTRAM_REGNET, gene_prefix='G_')
        from reframed.io.sbml import load_cbmodel
        model = load_cbmodel(OPTRAM_MODEL)
        self.problem = OptRamProblem(model, [], regnet)

    @pytest.mark.xfail
    def test_decode(self):
        import random
        candidate = self.problem.generator(random)
        solution = self.problem.decode(candidate)
        n_candidate = self.problem.encode(solution)
        self.assertEqual(candidate, n_candidate)

    def test_simul_constraints(self):
        """ Can not be run with a community cplex.
        """
        pass

    def test_to_constraints(self):
        """ Can not be run with a community cplex.
        """
        pass


class TestOptORF(unittest.TestCase):

    def setUp(self):

        _BIOMASS_ID = 'Biomass_Ecoli_core'
        _O2 = 'EX_o2_e'
        _GLC = 'EX_glc__D_e'
        _FUM = 'EX_fum_e'
        _AC = 'EX_ac_e'
        _GLU = 'EX_glu__L_e'
        _LAC = 'EX_lac__D_e'
        _SUC = 'EX_succ_e'

        from mewpy.io import read_model, Engines, Reader

        metabolic_reader = Reader(Engines.MetabolicSBML, EC_CORE_MODEL2)
        regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                                   EC_CORE_REG_MODEL,
                                   sep=',',
                                   id_col=0,
                                   rule_col=2,
                                   aliases_cols=[1],
                                   header=0)

        model = read_model(metabolic_reader, regulatory_reader)

        envcond = {'EX_glc__D_e': (-10, 100000.0)}
        from mewpy.simulation import get_simulator
        sim = get_simulator(model, envcond=envcond)
        sim.objective = _BIOMASS_ID

        _PRODUCT_ID = "EX_succ_e"
        from mewpy.optimization.evaluation import BPCY, WYIELD
        evaluator_1 = BPCY(_BIOMASS_ID, _PRODUCT_ID, method='pFBA')
        evaluator_2 = WYIELD(_BIOMASS_ID, _PRODUCT_ID)
        from mewpy.problems import OptORFProblem
        self.problem = OptORFProblem(model, [evaluator_1, evaluator_2], candidate_max_size=6)

    def test_targets(self):
        target = self.problem.target_list
        self.assertGreater(len(target), 0)

    def test_generator(self):
        import random
        candidate = self.problem.generator(random)
        self.assertGreater(len(candidate), 0)


if __name__ == '__main__':
    unittest.main()
