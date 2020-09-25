from inspyred.ec.emo import Pareto
from mewpy.visualization.plot import StreamingPlot
import numpy


def fitness_statistics(population):
    """Return the basic statistics of the population's fitness values.

    :param population: A population of individuals.

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

    :param population: (list) the population of Individuals.
    :param num_generations: (int) the number of elapsed generations.
    :param num_evaluations: (int) the number of evaluations already performed.
    :param args: (dict) a dictionary of keyword arguments.

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

    def __init__(self, reference_front=None, reference_point=None, display_frequency=1, axis_labels=None,
                 non_dominated=True, print_stats=True):
        self.figure = None
        self.display_frequency = display_frequency
        self.reference_point = reference_point
        self.reference_front = reference_front
        self.print_stats = print_stats
        self.axis_labels = axis_labels
        self.non_dominated = non_dominated

    def update(self, population, num_generations, num_evaluations, args):
        generations = num_generations
        evaluations = num_evaluations

        if population:
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


def non_dominated_population(population, maximize=True):
    """ The non dominated solutions from the population.

    :param population: A list of individuals.
    :param maximize: (bool) Optimization direction.
    :returns: a list of non-dominated solutions.

    """
    population.sort(reverse=True)
    non_dominated = []
    for i in range(len(population)-1):
        individual = population[i]
        j = 0
        dominates = True
        while j < len(population) and dominates:
            if dominance_test(individual, population[j], maximize=maximize) == -1:
                dominates = False
            else:
                j += 1
        if dominates:
            non_dominated.append(individual)

    result = non_dominated
    return result


def dominance_test(solution1, solution2, maximize=True):
    """
    Testes Pareto dominance.

    :param solution1: The first solution.
    :param solution2: The second solution
    :param maximize: (bool) maximization (True) or minimization (False)

    :returns: 1, if the first solution dominates the second;
            -1, if the second solution dominates the first;
            0, if non of the solutions dominates the other.

    """
    best_is_one = 0
    best_is_two = 0

    if isinstance(solution1.fitness, Pareto):
        values1 = solution1.fitness.values
        values2 = solution2.fitness.values
    else:
        values1 = [solution1.fitness]
        values2 = [solution2.fitness]

    for i in range(len(values1)):
        value1 = values1[i]
        value2 = values2[i]
        if value1 != value2:
            if value1 < value2:
                best_is_two = 1
            if value1 > value2:
                best_is_one = 1

    if best_is_one > best_is_two:
        if maximize:
            result = 1
        else:
            result = -1
    elif best_is_two > best_is_one:
        if maximize:
            result = -1
        else:
            result = 1
    else:
        result = 0

    return result
