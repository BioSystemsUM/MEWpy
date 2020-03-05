from mewpy.utils.process import MultiProcessorEvaluator, cpu_count , DaskEvaluator 
from mewpy.utils.utilities import non_dominated_population
from mewpy.utils.constants import EAConstants, ModelConstants
from mewpy.optimization.problem import Strategy
import mewpy.optimization.operators as op
import mewpy.optimization.observers as observers
from random import Random
from time import time
import inspyred



class EA:
    """
    EA running helper

    arguments
        *problem*: the optimization problem
        *initial_population* (list): the EA initial population
        *max_generations* (int): the number of iterations of the EA (stopping criteria) 
    """

    def __init__(self, problem, initial_population=[], max_generations=EAConstants.MAX_GENERATIONS, mp = True ,visualizer=False):

        self.problem = problem
        self.initial_population = initial_population
        self.max_generations = max_generations
        self.visualizer = visualizer
        self.mp = mp
        self.final_population = None

        self.ou_variators = [op.uniform_crossover_OU,
                             op.grow_mutation_OU,
                             op.shrink_mutation,
                             op.single_mutation_OU
                             ]

        self.ko_variators = [op.uniform_crossover_KO,
                             op.grow_mutation_KO,
                             op.shrink_mutation,
                             op.single_mutation_KO
                             ]

        


        ## needs to be defined elsewhere
        self.args = { 
                'num_selected' : 100,
                'max_generations' :self.max_generations,
                # operators probabilities
                'gs_mutation_rate' : 0.1,
                'mutation_rate' : 0.1,
                'crossover_rate': 0.9,
                # candidate size
                'candidate_min_size':self.problem.candidate_min_size,
                'candidate_max_size':self.problem.candidate_max_size
                }
        if self.problem.number_of_objectives == 1:
            self.args['tournament_size'] = 7




    
    def non_dominated_final_population(self):
        """
        returns the non dominated solutions from the final population.
        """
        return non_dominated_population(self.final_population,self.problem.is_maximization)


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
        """ Runs a single objective EA optimization
        """
        prng = Random()
        prng.seed(time())


        if self.mp:
            nmp = cpu_count()
            if ModelConstants.RESET_SOLVER:
                mp_evaluator = MultiProcessorEvaluator(self.problem.evaluate , nmp)
            else:
                try:
                    from mewpy.utils.process import RayEvaluator
                    mp_evaluator = RayEvaluator(self.problem, nmp)
                except:
                    Warning("Multiprocessing with persistente solver requires ray (pip install ray).")
                    mp_evaluator = MultiProcessorEvaluator(self.problem.evaluate , nmp)
            self.evaluator = mp_evaluator.evaluate
        else:
            self.evaluator = self.problem.evaluator


        ea = inspyred.ec.EvolutionaryComputation(prng)
        ea.selector = inspyred.ec.selectors.tournament_selection

        if self.problem.strategy == Strategy.KO:
            ea.variator = self.ko_variators
        else:
            ea.variator = self.ou_variators
        ea.observer = observers.results_observer
        ea.replacer = inspyred.ec.replacers.truncation_replacement
        ea.terminator = inspyred.ec.terminators.generation_termination
        
        final_pop = ea.evolve(generator=self.problem.generator,
                              evaluator= self.evaluator,
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
                mp_evaluator = MultiProcessorEvaluator(self.problem.evaluate , nmp)
            else:
                try:
                    from mewpy.utils.process import RayEvaluator
                    mp_evaluator = RayEvaluator(self.problem, nmp)
                except:
                    Warning("Multiprocessing with persistente solver requires ray (pip install ray).")
                    mp_evaluator = MultiProcessorEvaluator(self.problem.evaluate , nmp)
            self.evaluator = mp_evaluator.evaluate
        else:
            self.evaluator = self.problem.evaluator


        ea = inspyred.ec.emo.NSGA2(prng)

        if self.problem.strategy == Strategy.KO:
            ea.variator = self.ko_variators
        else:
            ea.variator = self.ou_variators

        ea.terminator = inspyred.ec.terminators.generation_termination
        if self.visualizer:
            axis_labels = [f.short_str() for f in self.problem.fevaluation]
            observer = observers.VisualizerObserver(axis_labels=axis_labels)
            ea.observer = observer.update
        else:
            ea.observer = observers.results_observer
    
        final_pop = ea.evolve(generator=self.problem.generator,
                              evaluator= self.evaluator,
                              pop_size=100,
                              seeds=self.initial_population,
                              maximize=self.problem.is_maximization,
                              bounder=self.problem.bounder,
                              **self.args  
                              )

        self.final_population = final_pop

        return final_pop
