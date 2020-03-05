from jmetal.algorithm.singleobjective.genetic_algorithm import GeneticAlgorithm
from jmetal.algorithm.multiobjective.spea2 import SPEA2
from jmetal.algorithm.multiobjective.nsgaii import NSGAII
from jmetal.util.termination_criterion import StoppingByEvaluations
from jmetal.operator import BinaryTournamentSelection
from jmetal.util.evaluator import MultiprocessEvaluator
from mewpy.optimization.jmetal.observers import PrintObjectivesStatObserver, VisualizerObserver
from mewpy.utils.constants import EAConstants
from mewpy.utils.process import cpu_count
from mewpy.optimization.problem import Strategy
from mewpy.optimization.jmetal.operators import *
from random import Random
from time import time
import jmetal


# MOEA alternatives
moea_map = {
        'NSGAII' : NSGAII,
        'SPEA2' : SPEA2
    }


class EA:
    """
    EA running helper

    arguments
        *problem*: the optimization problem
        *initial_population* (list): the EA initial population
        *max_generations* (int): the number of iterations of the EA (stopping criteria) 
    """

    def __init__(self, problem, initial_population=[], max_generations = EAConstants.MAX_GENERATIONS, mp = True, visualizer = False ):

        self.problem = problem
        self.initial_population = initial_population
        self.max_generations = max_generations
        self.final_population = None
        self.visualizer = visualizer
        
        self.crossover = UniformCrossoverKO(0.8,self.problem.candidate_max_size) if self.problem.strategy == Strategy.KO else UniformCrossoverOU(0.5,self.problem.candidate_max_size)
        mutators = []
        if self.problem.strategy == Strategy.KO:
            mutators.append(GrowMutationKO(1.0,max_size = self.problem.candidate_max_size))
            mutators.append(ShrinkMutation(1.0, min_size = self.problem.candidate_min_size))
            mutators.append(SingleMutationKO(1.0))
        else:
            mutators.append(GrowMutationOU(1.0,max_size = self.problem.candidate_max_size))
            mutators.append(ShrinkMutation(1.0, min_size = self.problem.candidate_min_size))
            mutators.append(SingleMutationOU(1.0))
        
        self.mutation = MutationContainer(0.3, mutators= mutators)
    
    
    



    def run(self):
        """ Runs the optimization for the defined problem.
            The number of objectives is defined to be the number of evaluation
            functions in fevalution. If there are more than one objective, 
            NSGAII is used as optimization engine. 
        """
        if self.problem.fevaluation is None or len(self.problem.fevaluation) == 0:
            raise ValueError("At leat one objective should be provided.")
        
         # builds the target list
        self.problem.pre_process()
        
        if self.problem.number_of_objectives == 1:
            return self._run_so()
        else:
            return self._run_mo()



    def _run_so(self):
        """ Runs a single objective EA optimization ()
        """
        max_evaluations = self.max_generations * 100

        algorithm = GeneticAlgorithm(
            problem=self.problem,
            population_size=100,
            offspring_population_size=100,
            mutation = self.mutation,
            crossover= self.crossover,
            selection=BinaryTournamentSelection(),
            termination_criterion=StoppingByEvaluations(max=max_evaluations)
        )

        algorithm.observable.register(observer=PrintObjectivesStatObserver())
        algorithm.run()
        
        result = algorithm.solutions
        if self.problem.is_maximization:
            result = self.fitness_min_to_max(result)
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
            problem=self.problem,
            population_size=100,
            offspring_population_size=100,
            mutation=self.mutation,
            crossover=self.crossover,
            termination_criterion=StoppingByEvaluations(max=max_evaluations)
        )
        
        if self.visualizer:
            algorithm.observable.register(observer=VisualizerObserver())
        algorithm.observable.register(observer=PrintObjectivesStatObserver())
        
        algorithm.run()
        result = algorithm.solutions
        if self.problem.is_maximization:
            result = self.fitness_min_to_max(result)
        return result



    def fitness_min_to_max(self, population, inline = True):
        """ Corrects fitness values for maximization problems
            TODO: verify each objective individualy
        """
        if not inline:
            p = copy.copy(population)
        else:
            p = population

        for i in range(len(p)):
                obj = [ abs(x) for x in p[i].objectives]
                p[i].objectives = obj
        return p