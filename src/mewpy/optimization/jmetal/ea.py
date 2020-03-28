from jmetal.algorithm.singleobjective.genetic_algorithm import GeneticAlgorithm
from jmetal.algorithm.multiobjective import SPEA2, NSGAII
from jmetal.util.termination_criterion import StoppingByEvaluations
from jmetal.operator import BinaryTournamentSelection
from jmetal.util.evaluator import MultiprocessEvaluator
from mewpy.optimization.ea import AbstractEA, Solution as MSolution
from mewpy.problems import Strategy
from mewpy.optimization.jmetal.problem import JMetalKOProblem, JMetalOUProblem
from mewpy.optimization.jmetal.observers import PrintObjectivesStatObserver, VisualizerObserver
from mewpy.utils.constants import EAConstants
from mewpy.utils.process import cpu_count
from mewpy.optimization.jmetal.operators import *
from random import Random
from time import time


# MOEA alternatives
moea_map = {
    'NSGAII': NSGAII,
    'SPEA2': SPEA2
}


class EA(AbstractEA):
    """
    EA running helper

    arguments
        *problem*: the optimization problem
        *initial_population* (list): the EA initial population
        *max_generations* (int): the number of iterations of the EA (stopping criteria) 
    """

    def __init__(self, problem, initial_population=[], max_generations=EAConstants.MAX_GENERATIONS, mp=True, visualizer=False):

        super(EA, self).__init__(problem, initial_population=initial_population,
                                 max_generations=max_generations, visualizer=visualizer)

        self.crossover = UniformCrossoverKO(0.8, self.problem.candidate_max_size) if self.problem.strategy == Strategy.KO else UniformCrossoverOU(
            0.5, self.problem.candidate_max_size)
        mutators = []
        if self.problem.strategy == Strategy.KO:
            self.ea_problem = JMetalKOProblem(self.problem)
            mutators.append(GrowMutationKO(
                1.0, max_size=self.problem.candidate_max_size))
            mutators.append(ShrinkMutation(
                1.0, min_size=self.problem.candidate_min_size))
            mutators.append(SingleMutationKO(1.0))
        else:
            self.ea_problem = JMetalOUProblem(self.problem)
            mutators.append(GrowMutationOU(
                1.0, max_size=self.problem.candidate_max_size))
            mutators.append(ShrinkMutation(
                1.0, min_size=self.problem.candidate_min_size))
            mutators.append(SingleMutationOU(1.0))
        self.mutation = MutationContainer(0.3, mutators=mutators)

    def _run_so(self):
        """ Runs a single objective EA optimization ()
        """
        max_evaluations = self.max_generations * 100

        algorithm = GeneticAlgorithm(
            problem=self.ea_problem,
            population_size=100,
            offspring_population_size=100,
            mutation=self.mutation,
            crossover=self.crossover,
            selection=BinaryTournamentSelection(),
            termination_criterion=StoppingByEvaluations(
                max_evaluations=max_evaluations)
        )

        algorithm.observable.register(observer=PrintObjectivesStatObserver())
        algorithm.run()

        result = algorithm.solutions
        return result

    def _run_mo(self):
        """ Runs a multi objective EA (SPEA or NSGAII) optimization
        """
        max_evaluations = self.max_generations * 100
        ncpu = cpu_count()
        if EAConstants.JMETAL_MOEA and EAConstants.JMETAL_MOEA in moea_map.keys():
            f = moea_map[EAConstants.JMETAL_MOEA]
        else:
            f = moea_map['SPEA2']
        algorithm = f(
            population_evaluator=MultiprocessEvaluator(ncpu),
            problem=self.ea_problem,
            population_size=100,
            offspring_population_size=100,
            mutation=self.mutation,
            crossover=self.crossover,
            termination_criterion=StoppingByEvaluations(
                max_evaluations=max_evaluations)
        )

        if self.visualizer:
            algorithm.observable.register(observer=VisualizerObserver())
        algorithm.observable.register(observer=PrintObjectivesStatObserver())

        algorithm.run()
        result = algorithm.solutions
        return result

    def _convertPopulation(self, population):
        p = []
        for i in range(len(population)):
            # Corrects fitness values for maximization problems
            # TODO: verify each objective individualy
            if self.problem.is_maximization:
                obj = [abs(x) for x in population[i].objectives]
            else:
                obj = [x for x in population[i].objectives]
            val = set(population[i].variables[:])
            const = self.problem.decode(val)
            solution = MSolution(val, obj, const)
            p.append(solution)
        return p
