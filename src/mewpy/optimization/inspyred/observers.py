from inspyred.ec.emo import Pareto
from mewpy.optimization.ea import Solution
from mewpy.visualization.plot import StreamingPlot
from mewpy.utils.utilities import non_dominated_population
import math
import numpy


def fitness_statistics(population):
    """Return the basic statistics of the population's fitness values.

    Arguments:

    - *population* -- the population of individuals 

    """

    stats = {}

    population.sort(reverse=True)
    first = population[0].fitness

    if isinstance(first, Pareto):
        n = len(first.values)
        for i in range(n):
            f = [p.fitness.values[i] for p in population]
            worst_fit = min(f)
            best_fit = max(f)
            med_fit = numpy.median(f)
            avg_fit = numpy.mean(f)
            std_fit = numpy.std(f)
            stats['obj_{}'.format(i)] = {'best': best_fit, 'worst': worst_fit,
                                         'mean': avg_fit, 'median': med_fit, 'std': std_fit}
    else:
        worst_fit = population[-1].fitness
        best_fit = population[0].fitness
        f = [p.fitness for p in population]
        med_fit = numpy.median(f)
        avg_fit = numpy.mean(f)
        std_fit = numpy.std(f)
        stats['obj'] = {'best': best_fit, 'worst': worst_fit,
                        'mean': avg_fit, 'median': med_fit, 'std': std_fit}

    return stats


def results_observer(population, num_generations, num_evaluations, args):
    """
    Print the output of the evolutionary computation to a file with the follow fields:
    - number of generation
    - fitness of candidate
    - the solution candidates
    - the solution encoded candidates

    Args:
        population (list): the population of Individuals
        num_generations (int): the number of elapsed generations
        num_evaluations (int): the number of evaluations already performed
        args (dict): a dictionary of keyword arguments
    """

    stats = fitness_statistics(population)
    title = "Gen    Eval|"
    values = "{0:>4} {1:>6}|".format(num_generations, num_evaluations)

    for key in stats:
        s = stats[key]
        title = title + "     Worst      Best    Median   Average   Std Dev|"
        values = values + "  {0:.6f}  {1:.6f}  {2:.6f}  {3:.6f}  {4:.6f}|".format(s['worst'],
                                                                                  s['best'],
                                                                                  s['median'],
                                                                                  s['mean'],
                                                                                  s['std'])
    if num_generations == 0:
        print(title)
    print(values)


class VisualizerObserver():
    """

    """

    def __init__(self, reference_front=None, reference_point=None, display_frequency=1, axis_labels=None, non_dominated=True, print_stats=True):
        self.figure = None
        self.display_frequency = display_frequency
        self.reference_point = reference_point
        self.reference_front = reference_front
        self.print_stats = print_stats
        self.axis_labels = axis_labels
        self.non_dominated = non_dominated

    def _convertPopulation(self, population):
        p = []
        for i in range(len(population)):
            if isinstance(population[i].fitness, list):
                obj = population[i].fitness
            else:
                obj = [population[i].fitness]
            val = population[i].candidate
            solution = Solution(val, obj)
            p.append(solution)
        return p

    def update(self, population, num_generations, num_evaluations, args):
        generations = num_generations
        evaluations = num_evaluations

        if population:
            population = self._convertPopulation(population)
            if self.non_dominated:
                pop = non_dominated_population(population)
            else:
                pop = population

            if self.figure is None:
                self.figure = StreamingPlot(axis_labels=self.axis_labels)
                solutions = []
                for i in range(len(pop)):
                    obj = pop[i].fitness
                    solutions.append(obj)
                self.figure.plot(solutions)

            if (generations % self.display_frequency) == 0:
                solutions = []
                for i in range(len(pop)):
                    obj = pop[i].fitness
                    solutions.append(obj)
                self.figure.update(solutions)
                self.figure.ax.set_title(
                    'Eval: {}'.format(evaluations), fontsize=13)

            if self.print_stats:
                results_observer(population, num_generations,
                                 num_evaluations, args)
