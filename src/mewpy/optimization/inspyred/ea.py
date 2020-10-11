from mewpy.utils.process import MultiProcessorEvaluator, cpu_count
from mewpy.utils.constants import EAConstants, ModelConstants
from mewpy.optimization.ea import AbstractEA, Solution
from mewpy.optimization.inspyred.problem import InspyredProblem
from mewpy.optimization.inspyred import operators as op
from mewpy.optimization.inspyred import observers
from random import Random
from time import time
import inspyred



SOEA={
    'GA':inspyred.ec.EvolutionaryComputation,
    'SA':inspyred.ec.SA
}

class EA(AbstractEA):
    """
    EA running helper

    :param problem: the optimization problem.
    :param initial_population: (list) the EA initial population.
    :param max_generations: (int) the number of iterations of the EA (stopping criteria).
    """

    def __init__(self, problem, initial_population=[], max_generations=EAConstants.MAX_GENERATIONS, mp=True,
                 visualizer=False, algorithm=None):

        super(EA, self).__init__(problem, initial_population=initial_population,
                                 max_generations=max_generations, mp=mp, visualizer=visualizer)

        self.algorithm_name = algorithm
        self.ea_problem = InspyredProblem(self.problem)
        from mewpy.problems import Strategy
        if self.problem.strategy == Strategy.OU:
            self.variators = [op.uniform_crossover_OU,
                              op.grow_mutation_OU,
                              op.shrink_mutation,
                              op.single_mutation_OU
                              ]
        elif self.problem.strategy == Strategy.KO:
            self.variators = [op.uniform_crossover_KO,
                              op.grow_mutation_KO,
                              op.shrink_mutation,
                              op.single_mutation_KO
                              ]
        else:
            raise ValueError("Unknow strategy")

        # needs to be defined elsewhere
        self.args = {
            'num_selected': 100,
            'max_generations': self.max_generations,
            # operators probabilities
            'gs_mutation_rate': 0.1,
            'mutation_rate': 0.1,
            'crossover_rate': 0.9,
            # candidate size
            'candidate_min_size': self.problem.candidate_min_size,
            'candidate_max_size': self.problem.candidate_max_size
        }
        if self.problem.number_of_objectives == 1:
            self.args['tournament_size'] = 7

    def _run_so(self):
        """ Runs a single objective EA optimization
        """
        prng = Random()
        prng.seed(time())

        if self.mp:
            nmp = cpu_count()
            if ModelConstants.RESET_SOLVER:
                mp_evaluator = MultiProcessorEvaluator(
                    self.ea_problem.evaluate, nmp)
            else:
                try:
                    from mewpy.utils.process import RayEvaluator
                    mp_evaluator = RayEvaluator(self.ea_problem, nmp)
                except ImportError:
                    Warning(
                        "Multiprocessing with persistente solver requires ray (pip install ray). Linux only")
                    mp_evaluator = MultiProcessorEvaluator(
                        self.ea_problem.evaluate, nmp)
            self.evaluator = mp_evaluator.evaluate
        else:
            self.evaluator = self.ea_problem.evaluator

        if self.algorithm_name == 'SA':
            ea = inspyred.ec.SA(prng)
            print("Running SA")
        else:
            ea = inspyred.ec.EvolutionaryComputation(prng)
            print("Running GA")
        ea.selector = inspyred.ec.selectors.tournament_selection

        ea.variator = self.variators
        ea.observer = observers.results_observer
        ea.replacer = inspyred.ec.replacers.truncation_replacement
        ea.terminator = inspyred.ec.terminators.generation_termination

        final_pop = ea.evolve(generator=self.problem.generator,
                              evaluator=self.evaluator,
                              pop_size=100,
                              seeds=self.initial_population,
                              maximize=self.problem.is_maximization,
                              bounder=self.problem.bounder,
                              **self.args
                              )
        self.final_population = final_pop
        return final_pop

    def _run_mo(self):
        """ Runs a multi objective EA (NSGAII) optimization
        """
        prng = Random()
        prng.seed(time())

        if self.mp:
            nmp = cpu_count()
            if ModelConstants.RESET_SOLVER:
                mp_evaluator = MultiProcessorEvaluator(
                    self.ea_problem.evaluate, nmp)
            else:
                try:
                    from mewpy.utils.process import RayEvaluator
                    mp_evaluator = RayEvaluator(self.ea_problem, nmp)
                except ImportError:
                    Warning(
                        "Multiprocessing with persistente solver requires ray (pip install ray).")
                    mp_evaluator = MultiProcessorEvaluator(
                        self.ea_problem.evaluate, nmp)
            self.evaluator = mp_evaluator.evaluate
        else:
            self.evaluator = self.ea_problem.evaluator

        ea = inspyred.ec.emo.NSGA2(prng)
        print("Running NSGAII")
        ea.variator = self.variators
        ea.terminator = inspyred.ec.terminators.generation_termination
        if self.visualizer:
            axis_labels = [f.short_str() for f in self.problem.fevaluation]
            observer = observers.VisualizerObserver(axis_labels=axis_labels)
            ea.observer = observer.update
        else:
            ea.observer = observers.results_observer

        final_pop = ea.evolve(generator=self.problem.generator,
                              evaluator=self.evaluator,
                              pop_size=100,
                              seeds=self.initial_population,
                              maximize=self.problem.is_maximization,
                              bounder=self.problem.bounder,
                              **self.args
                              )

        self.final_population = final_pop
        return final_pop

    def _convertPopulation(self, population):
        p = []
        for i in range(len(population)):
            if self.problem.number_of_objectives == 1:
                obj = [population[i].fitness]
            else:
                obj = population[i].fitness
            val = population[i].candidate
            values = self.problem.translate(val)
            const = self.problem.decode(val)
            solution = Solution(values, obj, const)
            p.append(solution)
        return p
