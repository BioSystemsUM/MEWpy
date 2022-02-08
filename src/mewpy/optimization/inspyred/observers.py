import numpy
from inspyred.ec.emo import Pareto

from mewpy.visualization.plot import StreamingPlot
from ..ea import Solution, non_dominated_population


def fitness_statistics(population):
    """Return the basic statistics of the population's fitness values.

    :param population: A population of individuals.

    """

    def minuszero(value):
        return round(value, 6)
            
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
        stats['obj'] = {'best': minuszero(best_fit), 'worst': minuszero(worst_fit),
                        'mean': minuszero(avg_fit), 'median': minuszero(med_fit), 'std': minuszero(std_fit)}

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
                 non_dominated=False, print_stats=True):
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
        p = []
        for s in population:
            if isinstance(s.fitness, Pareto):
                a = Solution(s.candidate, s.fitness.values)
            else:
                a = Solution(s.candidate, [s.fitness])
            p.append(a)

        nds = non_dominated_population(p)
        ds = None

        if not self.non_dominated:
            ds = list(set(p)-set(nds))

        if self.figure is None:
            self.figure = StreamingPlot(axis_labels=self.axis_labels)
            self.figure.plot(nds, dominated=ds)

        if (generations % self.display_frequency) == 0:
            self.figure.update(nds, dominated=ds)
            self.figure.ax.set_title(
                'Eval: {}'.format(evaluations), fontsize=13)

        if self.print_stats:
            results_observer(population, num_generations,
                             num_evaluations, args)
