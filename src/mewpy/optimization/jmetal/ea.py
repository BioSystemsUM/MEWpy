from jmetal.algorithm.singleobjective.genetic_algorithm import GeneticAlgorithm
from jmetal.algorithm.multiobjective.nsgaiii import NSGAIII, UniformReferenceDirectionFactory
from jmetal.algorithm.multiobjective.nsgaii import NSGAII
from jmetal.algorithm.multiobjective.spea2 import SPEA2
from jmetal.util.termination_criterion import StoppingByEvaluations
from jmetal.operator import BinaryTournamentSelection
from mewpy.utils.process import MultiProcessorEvaluator
from mewpy.optimization.ea import AbstractEA, Solution 
from mewpy.optimization.jmetal.problem import JMetalKOProblem, JMetalOUProblem
from mewpy.optimization.jmetal.observers import PrintObjectivesStatObserver, VisualizerObserver
from mewpy.utils.constants import EAConstants
from mewpy.utils.process import cpu_count
from mewpy.optimization.jmetal.operators import (ShrinkMutation,GrowMutationKO,GrowMutationOU,UniformCrossoverKO,
        UniformCrossoverOU,SingleMutationKO,SingleMutationOU,SingleMutationOULevel,MutationContainer)
from random import Random
from time import time


# MOEA alternatives
moea_map = {
    'NSGAII': NSGAII,
    'SPEA2': SPEA2,
    'NSGAIII': NSGAIII
}


class EA(AbstractEA):
    """
    EA running helper for JMetal.

    
    :param problem: The optimization problem.
    :param initial_population: (list) The EA initial population.
    :param max_generations: (int) The number of iterations of the EA (stopping criteria). 
    """

    def __init__(self, problem, initial_population=[], max_generations=EAConstants.MAX_GENERATIONS, mp=True, visualizer=False,algorithm = None):

        super(EA, self).__init__(problem, initial_population=initial_population,
                                 max_generations=max_generations, mp=mp, visualizer=visualizer)
        self.algorithm_name = algorithm
        from mewpy.problems import Strategy
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
        if self.algorithm_name in moea_map.keys():
            f = moea_map[self.algorithm_name]
        else:
            if self.ea_problem.number_of_objectives > 2:
                self.algorithm_name== 'NSGAIII'
            else:
                f = moea_map['SPEA2']

        print(f"Running {self.algorithm_name}")
        if self.algorithm_name  and self.algorithm_name== 'NSGAIII':
            algorithm = NSGAIII(
                population_evaluator=MultiProcessorEvaluator(self.ea_problem.evaluate,ncpu),
                problem=self.ea_problem,
                population_size=100,
                mutation=self.mutation,
                crossover=self.crossover,
                termination_criterion=StoppingByEvaluations(max_evaluations=max_evaluations),
                reference_directions=UniformReferenceDirectionFactory(self.ea_problem.number_of_objectives, n_points=99)
            )

        else:
            algorithm = f(
                population_evaluator=MultiProcessorEvaluator(self.ea_problem.evaluate,ncpu),
                problem=self.ea_problem,
                population_size=100,
                offspring_population_size=100,
                mutation=self.mutation,
                crossover=self.crossover,
                termination_criterion=StoppingByEvaluations(
                    max_evaluations=max_evaluations),
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
            values = self.problem.translate(val)
            const = self.problem.decode(val)
            solution = Solution(values, obj, const)
            p.append(solution)
        return p
