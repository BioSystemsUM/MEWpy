from jmetal.core.observer import Observer
from jmetal.lab.visualization import StreamingPlot
from mewpy.optimization.ea import non_dominated_population
from typing import List, TypeVar
import logging
import numpy
import copy

S = TypeVar('S')
LOGGER = logging.getLogger('mewpy')


class VisualizerObserver(Observer):

    def __init__(self,
                 reference_front: List[S] = None,
                 reference_point: list = None,
                 display_frequency: float = 1.0,
                 non_dominated=True) -> None:
        self.figure = None
        self.display_frequency = display_frequency
        self.reference_point = reference_point
        self.reference_front = reference_front
        self.non_dominated = non_dominated

    def update(self, *args, **kwargs):
        evaluations = kwargs['EVALUATIONS']
        solutions = kwargs['SOLUTIONS']

        if solutions:
            if self.figure is None:

                axis_labels = None
                problem = kwargs['PROBLEM']
                if problem and problem.obj_labels:
                    axis_labels = problem.obj_labels

                self.figure = StreamingPlot(reference_point=self.reference_point,
                                            reference_front=self.reference_front,
                                            axis_labels=axis_labels)
                self.figure.plot(solutions)

            if (evaluations % self.display_frequency) == 0:
                # check if reference point has changed
                reference_point = kwargs.get('REFERENCE_POINT', None)
                # negative fitness values are converted to positive
                population = copy.copy(solutions)
                if self.non_dominated:
                    population = non_dominated_population(population)
                for i in range(len(population)):
                    obj = [abs(x) for x in population[i].objectives]
                    population[i].objectives = obj

                if reference_point:
                    self.reference_point = reference_point
                    self.figure.update(population, reference_point)
                else:
                    self.figure.update(population)

                self.figure.ax.set_title(
                    'Eval: {}'.format(evaluations), fontsize=13)


class PrintObjectivesStatObserver(Observer):

    def __init__(self, frequency: float = 1.0) -> None:
        """ Show the number of evaluations, best fitness and computing time.

        :param frequency: Display frequency. """
        self.display_frequency = frequency
        self.first = True

    def fitness_statistics(self, solutions,problem):
        """Return the basic statistics of the population's fitness values.
        :param list solutions: List of solutions.
        :param problem: The jMetalPy problem.
        :returns: A statistics dictionary.
        """

        stats = {}
        first = solutions[0].objectives
        # number of objectives
        n = len(first)
        for i in range(n):
            direction = problem.obj_directions[i]
            factor = 1 if direction== problem.MINIMIZE else -1
            f = [(factor * p.objectives[i]) for p in solutions]
            if direction==problem.MAXIMIZE:
                worst_fit = min(f)
                best_fit = max(f)
            else:
                worst_fit = max(f)
                best_fit = min(f)
            med_fit = numpy.median(f)
            avg_fit = numpy.mean(f)
            std_fit = numpy.std(f)
            stats['obj_{}'.format(i)] = {'best': best_fit, 'worst': worst_fit,
                                         'mean': avg_fit, 'median': med_fit, 'std': std_fit}
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
            return title+"\n"+values
        else:
            return values

    def update(self, *args, **kwargs):
        evaluations = kwargs['EVALUATIONS']
        solutions = kwargs['SOLUTIONS']
        problem = kwargs['PROBLEM']
        if (evaluations % self.display_frequency) == 0 and solutions:
            if type(solutions) == list:
                stats = self.fitness_statistics(solutions,problem)
                message = self.stats_to_str(stats, evaluations, self.first)
                self.first = False
            else:
                fitness = solutions.objectives
                res = abs(fitness[0])
                message = 'Evaluations: {}\tFitness: {}'.format(
                    evaluations, res)
            print(message)
