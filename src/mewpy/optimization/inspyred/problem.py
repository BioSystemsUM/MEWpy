# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
##############################################################################
Inspyred Problems

Authors: Vitor Pereira
##############################################################################
"""
from inspyred.ec.emo import Pareto
from mewpy.util.process import Evaluable


class IntTuppleBounder(object):
    """
    A bounder for (int,int,...) representations

    :param lower_bound: The integers lower bound.
    :param upper_bound: The integers upper bound.

    """

    def __init__(self, lower_bound:int, upper_bound:int):
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

    def __init__(self, problem, directions):
        self.problem = problem
        self.direction = directions

    def evaluate(self, solution):
        """Evaluates a single solution

            :param solution: The individual to be evaluated.
            :returns: A list with a fitness value or a Pareto object.

        """
        p = self.problem.evaluate_solution(solution)
        # single objective
        if self.problem.number_of_objectives == 1:
            return p[0]*self.direction[0]
        # multi objective
        else:
            v = [a*b for a, b in zip(p, self.direction)]
            return Pareto(v)

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
