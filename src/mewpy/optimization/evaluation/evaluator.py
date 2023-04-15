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
Abstract evaluators 

Author: Vitor Pereira
##############################################################################
"""
from abc import ABCMeta, abstractmethod
import math

class EvaluationFunction:
    __metaclass__ = ABCMeta

    def __init__(self, maximize:bool=True,
                 worst_fitness:float=0.0):
        """This abstract class should be extended by all evaluation functions.

        :param maximize: Wether to maximize (True) or minimize (False), defaults to True
        :type maximize: bool, optional
        :param worst_fitness: The worst fitness value, defaults to 0.0
        :type worst_fitness: float, optional
        """
        self.worst_fitness = worst_fitness
        self.maximize = maximize

    @abstractmethod
    def get_fitness(self, simul_results, candidate, **kwargs):
        """Evaluates a candidate

        :param simul_results: (dic) A dictionary of phenotype SimulationResult objects
        :param candidate:  Candidate beeing evaluated
        :returns: A fitness value.

        """
        raise NotImplementedError

    @abstractmethod
    def method_str(self):
        raise NotImplementedError

    def short_str(self):
        return self.method_str

    def __str__(self):
        return self.method_str()

    @abstractmethod
    def required_simulations(self):
        return None

    @property
    def no_solution(self):
        """
        Value to be retuned for wost case evaluation
        """
        if self.worst_fitness is not None:
            res = self.worst_fitness
        elif self.maximize:
            res = -math.inf
        else:
            res = math.inf
        return res

    def __call__(self, simulationResult, candidate, **kwargs):
        return self.get_fitness(simulationResult, candidate, **kwargs)


class PhenotypeEvaluationFunction(EvaluationFunction):

    def __init__(self, maximize=True, worst_fitness=0.0):
        super(PhenotypeEvaluationFunction, self).__init__(maximize=maximize,
                                                          worst_fitness=worst_fitness)


class KineticEvaluationFunction(EvaluationFunction):

    def __init__(self, maximize=True, worst_fitness=0.0):
        super(KineticEvaluationFunction, self).__init__(maximize=maximize,
                                                        worst_fitness=worst_fitness)
