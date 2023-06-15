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
General evaluators 

Author: Vitor Pereira
##############################################################################
"""
from .evaluator import (PhenotypeEvaluationFunction,
                        KineticEvaluationFunction, 
                        EvaluationFunction)
from functools import reduce
import warnings
import numpy as np

from typing import List, Dict


class AggregatedSum(PhenotypeEvaluationFunction, KineticEvaluationFunction):
    
    def __init__(self,
                 fevaluation: List[EvaluationFunction], 
                 tradeoffs: List[float]=None, 
                 maximize: bool=True):
        """
        Aggredated sum evaluation function. Used to converte MOEAs into Single Objective EAs.

        :param fevaluation: (list) List of evaluation functions.
        :param tradeoffs: (list) Tradeoff values for each evaluation function. If None, all functions have \
            the same associated weight.
    """
        super(AggregatedSum, self).__init__(
            maximize=maximize, worst_fitness=0.0)
        self.fevaluation = fevaluation
        if tradeoffs and len(tradeoffs) == len(fevaluation):
            self.tradeoffs = tradeoffs
        else:
            self.tradeoffs = [1 / len(self.fevaluation)] * \
                (len(self.fevaluation))

    def required_simulations(self):
        methods = []
        for f in self.fevaluation:
            methods.extend(f.required_simulations())
        return list(set(methods))

    def get_fitness(self, simul_results, candidate, **kwargs):
        """Evaluates a candidate

        :param simul_results: (dic) A dictionary of phenotype SimulationResult objects
        :param candidate:  Candidate beeing evaluated
        :returns: A fitness value.

        """
        res = []
        for f in self.fevaluation:
            res.append(f.get_fitness(simul_results, candidate, **kwargs))
        # return sum(map(lambda x, y: x * y, f, self.tradeoffs))
        return np.dot(res, self.tradeoffs)

    def short_str(self):
        return "Agg"

    def method_str(self):
        return "Aggregated Sum = " + reduce(lambda a, b: a + " " + b, [f.method_str() for f in self.fevaluation], "")


class CandidateSize(PhenotypeEvaluationFunction, KineticEvaluationFunction):
   
    def __init__(self, maximize:bool=False):
        """
        Maximize/minimize the number of modifications.

        :param (bool) maximize: Optimization sense. Default False (minimize)
        """
        super(CandidateSize, self).__init__(maximize=maximize, worst_fitness=0.0)

    def get_fitness(self, simulResult, candidate, **kwargs):
        return len(candidate)

    def required_simulations(self):
        """
        If the evaluation function requires a pre-simulation to compute fitness values
        """
        return []

    def short_str(self):
        return "Size"

    def method_str(self):
        return "Minimize/maximize the number of alterations"


class MinCandSize(CandidateSize):

    def __init__(self, maximize=False):
        super(MinCandSize, self).__init__(maximize=maximize, worst_fitness=0.0)
        warnings.warn("This class will soon be depricated. Use CandidateSize instead.")


class ModificationType(PhenotypeEvaluationFunction, KineticEvaluationFunction):
    

    def __init__(self, 
                 penalizations:Dict[str,float]={'KO': 5, 'UE': 2, 'OE': 0}, 
                 maximize:bool=True):
        """
        This Objective function favors solutions with deletions, under expression
        and over expression, in this same order.

        :param penalizations: Penalizations by modification type, defaults to {'KO': 5, 'UE': 2, 'OE': 0}
        :type penalizations: _type_, optional
        :param maximize: optimization direction: True (maximize) False (minimize), defaults to True
        :type maximize: bool, optional
        """
        super(ModificationType, self).__init__(maximize=maximize, worst_fitness=0.0)
        self.penalizations = penalizations
        

    def get_fitness(self, simulResult, candidate, **kwargs):
        sum = 0
        for v in candidate.values():
            if v == 0:
                sum += self.penalizations['KO']
            elif v < 1:
                sum += self.penalizations['UE']
            else:
                sum += self.penalizations['OE']
        return sum / len(candidate)

    def required_simulations(self):
        """
        If the evaluation function requires a pre-simulation to compute fitness values
        """
        return []

    def short_str(self):
        return "ModificationType"

    def method_str(self):
        return "ModificationType"
