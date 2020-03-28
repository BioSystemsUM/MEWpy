from inspyred.ec.emo import Pareto


class IntTuppleBounder(object):
    """
    A bounder for (int,int,...) representations
    """

    def __init__(self, lower_bound, upper_bound):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.range = [self.upper_bound[i]-self.lower_bound[i] +
                      1 for i in range(len(self.lower_bound))]

    def __call__(self, candidate, args):
        bounded_candidate = set()
        for c in candidate:
            l = []
            for i in range(len(c)):
                v = c[i] % self.range[i] + self.lower_bound[i]
                l.append(v)
            bounded_candidate.add(tuple(l))
        return bounded_candidate


class InspyredProblem:

    def __init__(self, problem):
        self.problem = problem

    def evaluate(self, solution):
        p = self.problem.evaluate_solution(solution)
        # single objective
        if self.problem.number_of_objectives == 1:
            return p[0]
        # multi objective
        else:
            return Pareto(p)

    def evaluator(self, candidates, args):
        """
        Evaluator 
        Note: shoudn't be dependent on args to ease multiprocessing
        returns a list of Pareto fitness values of a candidate list
        """
        fitness = []
        for candidate in candidates:
            p = self.evaluate(candidate)
            fitness.append(p)
        return fitness
