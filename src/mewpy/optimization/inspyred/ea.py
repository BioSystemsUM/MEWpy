# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
""" 
##############################################################################
EA Module for inspyred

Authors: Vitor Pereira
##############################################################################
"""
from random import Random
from time import time

import inspyred
from .settings import get_population_size, KO, PARAMETERS, OU
from .problem import InspyredProblem
from .observers import results_observer, VisualizerObserver
from .terminator import generation_termination
from .operators import OPERATORS
from ..ea import AbstractEA, Solution
from mewpy.util.constants import EAConstants
from mewpy.util.process import get_evaluator, cpu_count


SOEA = {
    'GA': inspyred.ec.EvolutionaryComputation,
    'SA': inspyred.ec.SA
}


class EA(AbstractEA):
    """
    EA running helper

    :param problem: the optimization problem.
    :param initial_population: (list) the EA initial population.
    :param max_generations: (int) the number of iterations of the EA (stopping criteria).
    """

    def __init__(self, problem, initial_population=[], max_generations=EAConstants.MAX_GENERATIONS, mp=True,
                 visualizer=False, algorithm=None, **kwargs):

        super(EA, self).__init__(problem, initial_population=initial_population,
                                 max_generations=max_generations, mp=mp, visualizer=visualizer, **kwargs)

        self.algorithm_name = algorithm
        self.directions = [1 if f.maximize else -1 for f in self.problem.fevaluation]
        self.ea_problem = InspyredProblem(self.problem, self.directions)
        # operators
        if self.problem.operators:
            self.variators = [OPERATORS[x] for x in self.problem.operators.keys()]
        elif self.problem.strategy == 'OU':
            self.variators = OU['variators']
        elif self.problem.strategy == 'KO':
            self.variators = KO['variators']
        else:
            raise ValueError("Unknow strategy")

        self.population_size = kwargs.get('population_size', get_population_size())

        # parameters
        self.args = PARAMETERS.copy()
        if self.problem.operators_param:
            self.args.update(self.problem.operators_param)
        self.args['num_selected'] = self.population_size
        self.args['max_generations'] = max_generations,
        self.args['candidate_min_size'] = self.problem.candidate_min_size
        self.args['candidate_max_size'] = self.problem.candidate_max_size
        if self.problem.number_of_objectives != 1:
            self.args.pop('tournament_size')
        self.seeds = [self.problem.encode(s) for s in initial_population]

    def get_population_size(self):
        return self.population_size

    def _run_so(self):
        """ Runs a single objective EA optimization
        """
        prng = Random()
        prng.seed(time())

        if self.mp:
            self.evaluator = get_evaluator(self.ea_problem, n_mp=cpu_count())
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
        ea.observer = results_observer
        ea.replacer = inspyred.ec.replacers.truncation_replacement
        ea.terminator = generation_termination
        self.algorithm = ea

        setattr(ea, 'directions', self.directions)

        final_pop = ea.evolve(generator=self.problem.generator,
                              evaluator=self.evaluator,
                              pop_size=self.population_size,
                              seeds=self.seeds,
                              maximize=True,
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
            self.evaluator = get_evaluator(self.ea_problem, n_mp=cpu_count())
        else:
            self.evaluator = self.ea_problem.evaluator

        ea = inspyred.ec.emo.NSGA2(prng)
        print("Running NSGAII")
        ea.variator = self.variators
        ea.terminator = generation_termination
        if self.visualizer:
            axis_labels = [f.short_str() for f in self.problem.fevaluation]
            observer = VisualizerObserver(axis_labels=axis_labels)
            ea.observer = observer.update
        else:
            ea.observer = results_observer

        setattr(ea, 'directions', self.directions)

        self.algorithm = ea
        final_pop = ea.evolve(generator=self.problem.generator,
                              evaluator=self.evaluator,
                              pop_size=self.population_size,
                              seeds=self.seeds,
                              maximize=True,
                              bounder=self.problem.bounder,
                              **self.args
                              )

        self.final_population = final_pop
        return final_pop

    def _convertPopulation(self, population):
        """Converts a population represented in Inpyred format to
        MEWpy solution format.

        :param list population: A list of solutions.

        :returns: A MEWpy list of solutions.
        """
        p = []
        for i in range(len(population)):
            if self.problem.number_of_objectives == 1:
                obj = [population[i].fitness * self.directions[0]]
            else:
                obj = [ a*b for a,b in zip(population[i].fitness.values,self.directions)]
            val = population[i].candidate
            values = self.problem.decode(val)
            const = self.problem.solution_to_constraints(values)
            solution = Solution(values, obj, const)
            p.append(solution)
        return p

    def _get_current_population(self):
        """Dumps the population for gracefull exit."""
        pop = self.algorithm.population
        cv = self._convertPopulation(pop)
        return cv
