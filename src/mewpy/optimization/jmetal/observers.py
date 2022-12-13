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
Obverser module for EA optimization based on jmetalpy

Authors: Vitor Pereira
##############################################################################
"""
import copy
import logging
from typing import List, TypeVar

import numpy
from mewpy.visualization.plot import StreamingPlot
from ..ea import non_dominated_population, Solution

S = TypeVar('S')
LOGGER = logging.getLogger('mewpy')


class VisualizerObserver():

    def __init__(self,
                 reference_front: List[S] = None,
                 reference_point: list = None,
                 display_frequency: float = 1.0,
                 non_dominated=False,
                 axis_labels=None,
                 nevaluations=None) -> None:
        self.figure = None
        self.display_frequency = display_frequency
        self.reference_point = reference_point
        self.reference_front = reference_front
        self.non_dominated = non_dominated
        self.axis_labels = axis_labels
        self.nevaluations = nevaluations

    def update(self, *args, **kwargs):
        evaluations = kwargs['EVALUATIONS']
        solutions = kwargs['SOLUTIONS']
        obj_directions = kwargs['PROBLEM'].obj_directions

        if solutions:
            population = [Solution(s.variables, s.objectives) for s in solutions]

            # check if reference point has changed
            # reference_point = kwargs.get('REFERENCE_POINT', None)
            # negative fitness values are converted to positive
            for i in range(len(population)):
                obj = population[i].fitness
                population[i].fitness = [(-1*obj[k]*obj_directions[k]) for k in range(len(obj))]

            nds = non_dominated_population(population)
            ds = None

            if not self.non_dominated:
                ds = list(set(population)-set(nds))

            if self.figure is None:
                self.figure = StreamingPlot(axis_labels=self.axis_labels)
                self.figure.plot(nds, dominated=ds)
            else:
                text = str(evaluations)
                self.figure.update(nds, dominated=ds, text=text)

            self.figure.ax.set_title(
                'Eval: {}'.format(evaluations), fontsize=13)


class PrintObjectivesStatObserver():

    def __init__(self, frequency: float = 1.0) -> None:
        """ Show the number of evaluations, best fitness and computing time.

        :param frequency: Display frequency. """
        self.display_frequency = frequency
        self.first = True

    def fitness_statistics(self, solutions, problem):
        """Return the basic statistics of the population's fitness values.
        :param list solutions: List of solutions.
        :param problem: The jMetalPy problem.
        :returns: A statistics dictionary.
        """
        def minuszero(value):
            return round(value, 6)

        stats = {}
        first = solutions[0].objectives
        # number of objectives
        n = len(first)
        for i in range(n):
            direction = problem.obj_directions[i]
            factor = 1 if direction == problem.MINIMIZE else -1
            f = [(factor * p.objectives[i]) for p in solutions]
            if direction == problem.MAXIMIZE:
                worst_fit = min(f)
                best_fit = max(f)
            else:
                worst_fit = max(f)
                best_fit = min(f)
            med_fit = numpy.median(f)
            avg_fit = numpy.mean(f)
            std_fit = numpy.std(f)
            stats['obj_{}'.format(i)] = {'best': minuszero(best_fit), 'worst': minuszero(worst_fit),
                                         'mean': minuszero(avg_fit), 'median': minuszero(med_fit), 'std': minuszero(std_fit)}
        return stats

    def stats_to_str(self, stats, evaluations, title=False):
        if title:
            title = "Eval(s)|"
        values = " {0:>6}|".format(evaluations)

        for key in stats:
            s = stats[key]
            if title:
                title = title + "     Worst      Best    Median   Average   Std Dev|"
            values = values + "  {0:.6f}  {1:.6f}  {2:.6f}  {3:.6f}  {4:.6f}|".format(s['worst'],
                                                                                      s['best'],
                                                                                      s['median'],
                                                                                      s['mean'],
                                                                                      s['std'])
        if title:
            return title + "\n" + values
        else:
            return values

    def update(self, *args, **kwargs):
        evaluations = kwargs['EVALUATIONS']
        solutions = kwargs['SOLUTIONS']
        problem = kwargs['PROBLEM']
        if (evaluations % self.display_frequency) == 0 and solutions:
            if type(solutions) == list:
                stats = self.fitness_statistics(solutions, problem)
                message = self.stats_to_str(stats, evaluations, self.first)
                self.first = False
            else:
                fitness = solutions.objectives
                res = abs(fitness[0])
                message = 'Evaluations: {}\tFitness: {}'.format(
                    evaluations, res)
            print(message)
