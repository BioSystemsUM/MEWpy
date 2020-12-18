from inspyred.ec.emo import Pareto
from ...util.process import Evaluable


class IntTuppleBounder(object):
    """
    A bounder for (int,int,...) representations

    :param lower_bound: The integers lower bound.
    :param upper_bound: The integers upper bound.

    """

    def __init__(self, lower_bound, upper_bound):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.range = [self.upper_bound[i] - self.lower_bound[i] +
                      1 for i in range(len(self.lower_bound))]

    def __call__(self, candidate, args):
        bounded_candidate = set()
        for c in candidate:
            al = []
            for i in range(len(c)):
                v = c[i] % self.range[i] + self.lower_bound[i]
                al.append(v)
            bounded_candidate.add(tuple(al))
        return bounded_candidate


class InspyredProblem(Evaluable):
    """Inspyred EA builder helper.

        :param problem: the optimization problem.
    """

    def __init__(self, problem):
        self.problem = problem

    def evaluate(self, solution):
        """Evaluates a single solution

            :param solution: The individual to be evaluated.
            :returns: A list with a fitness value or a Pareto object.

        """
        p = self.problem.evaluate_solution(solution)
        # single objective
        if self.problem.number_of_objectives == 1:
            return p[0]
        # multi objective
        else:
            return Pareto(p)

    def evaluator(self, candidates, *args):
        """
        Evaluator
        Note: shoudn't be dependent on args to ease multiprocessing

        :param candidates: A list of candidate solutions.
        :returns: A list of Pareto fitness values or a list of fitness values.

        """
        fitness = []
        for candidate in candidates:
            p = self.evaluate(candidate)
            fitness.append(p)
        return fitness
