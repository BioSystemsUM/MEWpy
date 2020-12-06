import unittest

MODELS_PATH = 'tests/data/'
EC_CORE_MODEL = MODELS_PATH + 'e_coli_core.xml.gz'

OPTRAM_MODEL = MODELS_PATH + 'yeast_7.6-optram.xml'
OPTRAM_GENES = MODELS_PATH + 'mgene.csv'
OPTRAM_TFS = MODELS_PATH + 'TFnames.csv'
OPTRAM_REGNET = MODELS_PATH + 'regnet.csv'

MIN_GROWTH = 0.1


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
        import random
        candidate = self.problem.generator(random, None)
        solution = self.problem.decode(candidate)
        constraints = self.problem.solution_to_constraints(solution)
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

if __name__ == '__main__':
    unittest.main()
