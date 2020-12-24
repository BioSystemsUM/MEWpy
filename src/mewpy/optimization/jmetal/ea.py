from jmetal.algorithm.multiobjective import NSGAII, SPEA2
from jmetal.algorithm.multiobjective.nsgaiii import NSGAIII
from jmetal.algorithm.multiobjective.nsgaiii import UniformReferenceDirectionFactory
from jmetal.algorithm.singleobjective import GeneticAlgorithm, SimulatedAnnealing
from jmetal.operator import BinaryTournamentSelection
from jmetal.util.termination_criterion import StoppingByEvaluations

from .operators import (ShrinkMutation, GrowMutationKO, GrowMutationOU, UniformCrossoverKO,
                        UniformCrossoverOU, SingleMutationKO, SingleMutationOU, SingleMutationOULevel,
                        MutationContainer)
from .observers import PrintObjectivesStatObserver, VisualizerObserver
from .problem import JMetalKOProblem, JMetalOUProblem
from ..ea import AbstractEA, Solution
from ...util.constants import EAConstants
from ...util.process import get_evaluator
from .settings import get_population_size

# SOEA alternatives
soea_map = {
    'GA': GeneticAlgorithm,
    'SA': SimulatedAnnealing
}
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

    def __init__(self, problem, initial_population=[], max_generations=EAConstants.MAX_GENERATIONS, mp=True,
                 visualizer=False, algorithm=None, **kwargs):

        super(EA, self).__init__(problem, initial_population=initial_population,
                                 max_generations=max_generations, mp=mp, visualizer=visualizer, **kwargs)

        self.algorithm_name = algorithm

        from mewpy.problems import Strategy

        mutators = []

        if self.problem.strategy == Strategy.KO:
            self.crossover = UniformCrossoverKO(0.8, self.problem.candidate_max_size)
            self.ea_problem = JMetalKOProblem(self.problem, self.initial_population)
            mutators.append(GrowMutationKO(
                1.0, max_size=self.problem.candidate_max_size))
            mutators.append(ShrinkMutation(
                1.0, min_size=self.problem.candidate_min_size))
            mutators.append(SingleMutationKO(1.0))
        else:
            self.crossover = UniformCrossoverOU(0.5, self.problem.candidate_max_size)
            self.ea_problem = JMetalOUProblem(self.problem, self.initial_population)
            mutators.append(GrowMutationOU(
                1.0, max_size=self.problem.candidate_max_size))
            mutators.append(ShrinkMutation(
                1.0, min_size=self.problem.candidate_min_size))
            mutators.append(SingleMutationOU(1.0))
            mutators.append(SingleMutationOULevel(1.0))

        self.mutation = MutationContainer(0.3, mutators=mutators)
        self.population_size = kwargs.get('population_size', get_population_size())
        self.max_evaluations = self.max_generations * self.population_size

    def get_population_size(self):
        return self.population_size

    def _run_so(self):
        """ Runs a single objective EA optimization ()
        """
        self.ea_problem.reset_initial_population_counter()
        if self.algorithm_name == 'SA':
            print("Running SA")
            self.mutation.probability = 1.0
            algorithm = SimulatedAnnealing(
                problem=self.ea_problem,
                mutation=self.mutation.probability,
                termination_criterion=StoppingByEvaluations(max_evaluations=self.max_evaluations)
            )

        else:
            print("Running GA")
            algorithm = GeneticAlgorithm(
                problem=self.ea_problem,
                population_size=self.population_size,
                offspring_population_size=self.population_size,
                mutation=self.mutation,
                crossover=self.crossover,
                selection=BinaryTournamentSelection(),
                termination_criterion=StoppingByEvaluations(
                    max_evaluations=self.max_evaluations)
            )

        algorithm.observable.register(observer=PrintObjectivesStatObserver())
        algorithm.run()

        result = algorithm.solutions
        return result

    def _run_mo(self):
        """ Runs a multi objective EA optimization
        """
        self.ea_problem.reset_initial_population_counter()
        if self.algorithm_name in moea_map.keys():
            f = moea_map[self.algorithm_name]
        else:
            if self.ea_problem.number_of_objectives > 2:
                self.algorithm_name == 'NSGAIII'
            else:
                f = moea_map['SPEA2']

        args = {
            'problem': self.ea_problem,
            'population_size': self.population_size,
            'mutation': self.mutation,
            'crossover': self.crossover,
            'termination_criterion': StoppingByEvaluations(max_evaluations=self.max_evaluations)
        }

        if self.mp:
            args['population_evaluator'] = get_evaluator(self.ea_problem)

        print(f"Running {self.algorithm_name}")
        if self.algorithm_name == 'NSGAIII':
            args['reference_directions'] = UniformReferenceDirectionFactory(self.ea_problem.number_of_objectives,
                                                                            n_points=self.population_size-1)
            algorithm = NSGAIII(**args)
        else:
            args['offspring_population_size'] = self.population_size
            algorithm = f(**args)

        if self.visualizer:
            algorithm.observable.register(observer=VisualizerObserver())
        algorithm.observable.register(observer=PrintObjectivesStatObserver())

        algorithm.run()
        result = algorithm.solutions
        return result

    def _convertPopulation(self, population):
        """Converts a population represented in Inpyred format to
        MEWpy solution format.

        :param list population: A list of solutions.

        :returns: A MEWpy list of solutions.
        """
        p = []
        for i in range(len(population)):
            # Corrects fitness values for maximization problems
            # TODO: verify each objective individualy
            if self.problem.is_maximization:
                obj = [abs(x) for x in population[i].objectives]
            else:
                obj = [x for x in population[i].objectives]
            val = set(population[i].variables[:])
            values = self.problem.decode(val)
            const = self.problem.solution_to_constraints(values)
            solution = Solution(values, obj, const)
            p.append(solution)
        return p
