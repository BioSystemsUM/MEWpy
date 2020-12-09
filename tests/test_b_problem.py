import unittest

MODELS_PATH = 'tests/data/'
EC_CORE_MODEL = MODELS_PATH + 'e_coli_core.xml.gz'

OPTRAM_MODEL = MODELS_PATH + 'yeast_7.6-optram.xml'
OPTRAM_GENES = MODELS_PATH + 'mgene.csv'
OPTRAM_TFS = MODELS_PATH + 'TFnames.csv'
OPTRAM_REGNET = MODELS_PATH + 'regnet.csv'

OPTORF_REG = MODELS_PATH + 'core_TRN_v2.csv'
OPTORF_ALIASES = MODELS_PATH + 'core_TRN_rfba_aliases.csv'


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
        candidate = self.problem.generator(random, None)
        self.assertGreater(len(candidate), 0)

    def test_decode(self):
        import random
        candidate = self.problem.generator(random, None)
        solution = self.problem.decode(candidate)
        self.assertGreater(len(solution), 0)

    def test_to_constraints(self):
        """
        """
        import random
        ispass = False
        tries = 0
        while not ispass and tries < MAX_TRIES:
            tries += 1
            candidate = self.problem.generator(random, None)
            solution = self.problem.decode(candidate)
            constraints = self.problem.solution_to_constraints(solution)
            if len(constraints) > 0:
                ispass = True
        self.assertGreater(len(constraints), 0)

    def test_simul_constraints(self):
        import random
        candidate = self.problem.generator(random, None)
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
        from mewpy.regulation.optram import OptRamProblem, load_optram
        regnet = load_optram(OPTRAM_GENES, OPTRAM_TFS, OPTRAM_REGNET, gene_prefix='G_')
        from reframed.io.sbml import load_cbmodel
        model = load_cbmodel(OPTRAM_MODEL)
        self.problem = OptRamProblem(model, [], regnet)

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
        import cobra.test
        model = cobra.test.create_test_model("textbook")
        _BIOMASS_ID = 'Biomass_Ecoli_core'

        envcond = {'EX_glc__D_e': (-10, 100000.0)}
        from mewpy.simulation import get_simulator
        sim = get_simulator(model, envcond=envcond)
        sim.objective = _BIOMASS_ID
        from mewpy.regulation import RFBAModel
        rfba = RFBAModel.from_tabular_format(OPTORF_REG, model, sim,
                                             sep=',', id_col=1, rule_col=2, aliases_cols=[0], header=0)
        rfba.update_aliases_from_tabular_format_file(OPTORF_ALIASES, id_col=1, aliases_cols=[0])

        initial_state = {var: 1 for var in rfba.targets}
        initial_state.update({_BIOMASS_ID: 0.1})
        rfba.initial_state = initial_state

        _PRODUCT_ID = "EX_succ_e"
        from mewpy.optimization.evaluation import BPCY, WYIELD
        evaluator_1 = BPCY(_BIOMASS_ID, _PRODUCT_ID, method='pFBA')
        evaluator_2 = WYIELD(_BIOMASS_ID, _PRODUCT_ID)
        from mewpy.regulation.optorf import OptOrfProblem
        self.problem = OptOrfProblem(model, [evaluator_1, evaluator_2], rfba, candidate_max_size=6)

    def test_targets(self):
        target = self.problem.target_list
        self.assertGreater(len(target), 0)

    def test_generator(self):
        import random
        candidate = self.problem.generator(random, None)
        self.assertGreater(len(candidate), 0)


if __name__ == '__main__':
    unittest.main()
