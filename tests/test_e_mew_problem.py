import unittest
from pathlib import Path

EC_CORE_MODEL = Path(__file__).parent.joinpath('data', 'e_coli_core.xml')

MIN_GROWTH = 0.1

MAX_TRIES = 3


class TestRKOP(unittest.TestCase):

    def setUp(self):
        from mewpy.io import read_sbml
        model = read_sbml(EC_CORE_MODEL, regulatory=False, warnings=False)
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
        from mewpy.io import read_sbml
        model = read_sbml(EC_CORE_MODEL, regulatory=False, warnings=False)
        from mewpy.problems import ROUProblem
        self.problem = ROUProblem(model, [])


class TestGKOP(TestRKOP):

    def setUp(self):
        from mewpy.io import read_sbml
        model = read_sbml(EC_CORE_MODEL, regulatory=False, warnings=False)
        from mewpy.problems import GKOProblem
        self.problem = GKOProblem(model, [])


class TestGOUP(TestRKOP):

    def setUp(self):
        from mewpy.io import read_sbml
        model = read_sbml(EC_CORE_MODEL, regulatory=False, warnings=False)
        from mewpy.problems import GOUProblem
        self.problem = GOUProblem(model, [])


if __name__ == '__main__':
    unittest.main()
