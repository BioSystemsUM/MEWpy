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
Phenotype evaluators 

Author: Vitor Pereira
##############################################################################
"""
from .evaluator import (PhenotypeEvaluationFunction,
                        KineticEvaluationFunction, 
                        EvaluationFunction)
from mewpy.simulation.simulation import SimulationMethod, SStatus
from mewpy.solvers.ode import ODEStatus
from mewpy.simulation import get_simulator
from mewpy.util.constants import ModelConstants, EAConstants
import numpy as np
import math
import warnings

from typing import Dict, Union, List, Tuple

class TargetFlux(PhenotypeEvaluationFunction,KineticEvaluationFunction):

    def __init__(self, 
                 reaction:str, 
                 biomass:str=None, 
                 maximize:bool=True,
                 min_biomass_value:float=None,
                 min_biomass_per:float=0.0,
                 method:Union[str,SimulationMethod]=SimulationMethod.pFBA):
        """ Target Flux evaluation function.
        
        The fitness value is the flux value of the identified reaction.
        If the reaction parameter is None, the fitness value is the optimization objective value.
        Additional parameters include a minimum of allowed biomass value computed from the min_biomass_per
        and reference flux values

        :param reaction: (str) The reaction identifier whose flux value is to be used as fitness. Default None \
            in which case the model objective is considered.
        :param biomass: (str) The biomass reaction identifier.
        :param maximize: (boolean) The optimization direction. Default True for maximization.
        :param min_biomass_value: (float) The minimum biomass value.
        :param min_biomass_per: (float) Minimum biomass percentage. Only used if no min_biomass_value is provided.
        """
        super(TargetFlux, self).__init__(maximize=maximize, worst_fitness=0.0)
        self.reaction = reaction
        self.biomass = biomass
        self.min_biomass_value = min_biomass_value
        self.min_biomass_per = min_biomass_per
        self.method = method
        self.kinetic = False

    def get_fitness(self, simul_results, candidate, **kwargs):
        """Evaluates a candidate

        :param simul_results: (dic) A dictionary of phenotype SimulationResult objects
        :param candidate:  Candidate beeing evaluated
        :returns: A fitness value.

        """
        if self.kinetic:
            sim = simul_results
        else:
            sim = simul_results[self.method] if self.method in simul_results.keys(
            ) else None

        if not sim or sim.status not in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL, ODEStatus.OPTIMAL):
            return self.no_solution

        # only executed once if required
        if self.biomass and self.min_biomass_value is None:
            self.min_biomass_value = 0.0
            if self.min_biomass_per > 0.0:
                simulation = get_simulator(
                    sim.model, envcond=sim.envcond, constraints=sim.model_constraints)
                result = simulation.simulate(objective={self.biomass: 1})
                self.min_biomass_value = self.min_biomass_per * \
                    result.fluxes[self.biomass]

        if self.biomass and sim.fluxes[self.biomass] < self.min_biomass_value:
            res = self.no_solution
        elif self.reaction and self.reaction in sim.fluxes.keys():
            res = sim.fluxes[self.reaction]
        else:
            res = sim.objective_value

        return res

    def _repr_latex_(self):
        sense = '\\max' if self.maximize else '\\min'
        return "$$ %s %s $$" % (sense, self.reaction.replace('_', '\\_'))

    def required_simulations(self):
        """
        If the evaluation function requires a pre-simulation to compute fitness values
        """
        return [self.method]

    def short_str(self):
        return "TargetFlux"

    def method_str(self):
        return "TargetFlux {} with at least {} of biomass ({})".format(self.reaction, self.min_biomass_per,
                                                                       self.biomass)


class WYIELD(PhenotypeEvaluationFunction):
    def __init__(self,
                 biomassId: str,
                 productId: str, 
                 maximize: bool =True,
                 **kwargs):
        """ Weighted Yield (WYIELD) objective function, a linear combination of the target 
            product minimum and maximum FVA under the introduced metabolic modifications.

        :param biomassId: (str) Biomass reaction identifier.
        :param productId: (str) Target product reaction identifier.

        kwargs options:

        :param min_biomass_value: (float) Minimum biomass value (default None, in which case the min_biomass_per is used).
        :param min_biomass_per: (float) Instead of defining explicitly a minimum biomass value, a percentage of the wild \
            type biomass is used. Only used when no min_biomass_value is defined (default min_biomass_per: 0.10).
        :param alpha: (float) Tradeoff between the Max and min FVA of the target product (alpha Max + (1-alpha) min). \
            Must be in range [0,1]  (default alpha: 0.3).
        :param scale: (boolean) Defines if the WYIELD is devided by the biomass of the simulated result (default false).
    """

        super(WYIELD, self).__init__(maximize=maximize, worst_fitness=0.0)
        self.biomassId = biomassId
        self.productId = productId
        # parameters
        self.min_biomass_value = kwargs.get('min_biomass_value', None)
        self.min_biomass_per = kwargs.get('min_biomass_per', 0.1)
        self.scale = kwargs.get('scale', False)
        self.method = SimulationMethod.FBA
        self.alpha = kwargs.get('alpha', 0.3)
        self.obj_frac = kwargs.get('obj_frac', 0.99)
        if self.alpha > 1 or self.alpha < 0:
            warnings.warn(
                "The value of the tradeoff parameter alpha should be in range 0 to 1. Setting default value.")
            self.alpha = 0.3

    def get_fitness(self, simul_results, candidate, **kwargs):
        """Evaluates a candidate

        :param simul_results: (dic) A dictionary of phenotype SimulationResult objects
        :param candidate:  Candidate beeing evaluated
        :returns: A fitness value.

        """

        sim = simul_results[self.method] if self.method in simul_results.keys(
        ) else None
        if not sim or sim.status not in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL):
            return self.no_solution

        fvaMaxProd = 0.0
        fvaMinProd = 0.0
        scalefactor = kwargs.get('scalefactor', None)

        model = sim.model
        ssFluxes = sim.fluxes

        simulation = get_simulator(
            model, envcond=sim.envcond, constraints=sim.model_constraints)

        if not ssFluxes:
            return self.no_solution
        ids = list(ssFluxes.keys())
        if self.biomassId not in ids or self.productId not in ids:
            raise ValueError(
                "Reaction ids are not present in the fluxes distribution.")

        biomassFluxValue = ssFluxes[self.biomassId] * self.obj_frac

        try:
            # computed only once
            if self.min_biomass_value is None or self.min_biomass_value < 0.0:
                solution = simulation.simulate(
                    objective={self.biomassId: 1}, scalefactor=scalefactor)
                wtBiomassValue = solution.fluxes[self.biomassId]
                minBiomass = wtBiomassValue * self.min_biomass_per
                self.min_biomass_value = minBiomass
            else:
                minBiomass = self.min_biomass_value

            constraints = {}
            constraints.update(sim.simulation_constraints)
            # add biomass constraint
            constraints[self.biomassId] = (biomassFluxValue, ModelConstants.REACTION_UPPER_BOUND)

            # only need to simulate FVA max if alpha is larger than 0, otherwise it will always be zero
            if (self.alpha > 0):
                fvaMaxResult = simulation.simulate(
                    objective={self.productId: 1}, constraints=constraints, scalefactor=scalefactor)
                if fvaMaxResult.status == SStatus.OPTIMAL:
                    fvaMaxProd = fvaMaxResult.fluxes[self.productId]
                else:
                    return self.no_solution

            # only need to simulate FVA min if alpha is lesser than 1, otherwise it will always be zero
            if (self.alpha < 1):
                fvaMinResult = simulation.simulate(objective={
                    self.productId: 1}, constraints=constraints, maximize=False, scalefactor=scalefactor)
                if fvaMinResult.status == SStatus.OPTIMAL:
                    fvaMinProd = fvaMinResult.fluxes[self.productId]
                else:
                    return self.no_solution

            res = self.no_solution
            if EAConstants.DEBUG:
                print(f"WYIELD FVA max: {fvaMaxProd} min:{fvaMinProd}")
            if biomassFluxValue > minBiomass:
                res = (self.alpha * fvaMaxProd + (1 - self.alpha) * fvaMinProd)
                if self.scale:
                    res = res / biomassFluxValue
            return res
        except Exception:
            return self.no_solution

    def _repr_latex_(self):
        sense = '\\max' if self.maximize else '\\min'
        return "$$ %s \\left( %f \\times FVA_{max}(%s) + (1-%f) \\times FVA_{min}(%s) \\right) $$" % (sense,
                                                                                                      self.alpha,
                                                                                                      self.productId.replace(
                                                                                                          '_', '\\_'),
                                                                                                      self.alpha,
                                                                                                      self.productId.replace(
                                                                                                          '_', '\\_'),
                                                                                                      )

    def required_simulations(self):
        return [self.method]

    def short_str(self):
        return "WYIELD"

    def method_str(self):
        return "WYIELD biomass: {} product: {} ".format(self.biomassId, self.productId)


class BPCY(PhenotypeEvaluationFunction):

    def __init__(self,
                 biomass: str,
                 product: str,
                 uptake: str=None,
                 maximize: bool=True,
                 **kwargs):
        """
        This class implements the "Biomass-Product Coupled Yield" objective function. The fitness is given by the equation:
        (biomass_flux * product_flux)/ uptake_flux

        :param biomass: (str) Biomass reaction identifier
        :param product: (str) Target product reaction identifier
        :param uptake: (str) (optional) Reaction of uptake. If no substract is defined, ie uptake is None, a substract \
            flux value of 1.0 is considered.

        kargs options:

        :param method: (SimulationMethod) The simulation method. Default Node in which case received simulation results \
            are used to compute the biomass product coupled yield.
        :param reference: (dic) Wild type reference values when MOMA, lMOMA or ROOM are defined as method. \
            If not provided, wild type reference values will be computed.
        """
        super(BPCY, self).__init__(maximize=maximize, worst_fitness=0.0)
        self.biomassId = biomass
        self.productId = product
        self.uptakeId = uptake
        self.method = kwargs.get('method', SimulationMethod.pFBA)
        self.reference = kwargs.get('reference', None)
        self.worst_fitness = 0.0

    def get_fitness(self, simul_results, candidate, **kwargs):
        """Evaluates a candidate

        :param simul_results: (dic) A dictionary of phenotype SimulationResult objects
        :param candidate:  Candidate beeing evaluated
        :returns: A fitness value.

        """

        sim = simul_results[self.method] if self.method in simul_results.keys(
        ) else None
        if not sim or sim.status not in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL):
            return self.no_solution

        ssFluxes = sim.fluxes

        ids = list(ssFluxes.keys())
        if self.biomassId not in ids or self.productId not in ids:
            raise ValueError(
                "Biomass or product reactions ids are not present in the fluxes distribution.")

        if self.uptakeId and self.uptakeId in ids:
            uptake = abs(ssFluxes[self.uptakeId])
        else:
            uptake = 1.0

        if uptake == 0.0:
            return self.no_solution
        if EAConstants.DEBUG:
            try:
                print("BPCY Bionamss: {} product: {}".format(ssFluxes[self.biomassId], ssFluxes[self.productId]))
            except Exception:
                print("BPCY No Fluxes")
        return (ssFluxes[self.biomassId] * ssFluxes[self.productId]) / uptake

    def _repr_latex_(self):
        sense = '\\max' if self.maximize else '\\min'
        if self.uptakeId:
            return "$$ %s \\frac{%s \\times %s}{%s} $$" % (sense,
                                                           self.biomassId.replace('_', '\\_'),
                                                           self.productId.replace('_', '\\_'),
                                                           self.uptakeId.replace('_', '\\_'))
        else:
            return "$$ %s \\left( %s \\times %s \\right) $$" % (sense,
                                                                self.biomassId.replace('_', '\\_'),
                                                                self.productId.replace('_', '\\_'))

    def required_simulations(self):
        return [self.method]

    def short_str(self):
        return "BPCY"

    def method_str(self):
        if self.uptakeId:
            return "BPCY (" + self.biomassId + " * " + self.productId + ") / " + self.uptakeId
        else:
            return "BPCY (" + self.biomassId + " * " + self.productId + ")"


class BPCY_FVA(PhenotypeEvaluationFunction):
    
    def __init__(self,
                 biomass: str, 
                 product: str, 
                 uptake: str=None, 
                 maximize: bool=True, 
                 **kwargs):
        """
        This class implements the "Biomass-Product Coupled Yield" objective function with FVA as defined in
        "OptRAM: In-silico strain design via integrative regulatory-metabolic network modeling".
        It combines BPCY with WYIELD objective functions.

        The fitness is given by the equation:
        ((biomass_flux * product_flux) * / uptake_flux ) * (1-log((range)/(target))) where range=(FVA_max-FVA_min)/2
        and target= (FVA_max+FVA_min)/2

        :param biomass: (str) Biomass reaction identifier.
        :param product: (str) Target product reaction identifier.
        :param uptake: (str) (optional) Reaction of uptake. If no substract is defined, ie uptakeId is None, a substract \
            flux value of 1.0 is considered.

        kwargs options:

        :param method: (SimulationMethod) The simulation method. Default Node in which case received simulation results \
            are used to compute the biomass product coupled yield.
        :param reference: (dic) Wild type reference values when MOMA, lMOMA or ROOM are defined as method. \
            If not provided, wild type reference values will be computed.
    """
        super(BPCY_FVA, self).__init__(maximize=maximize, worst_fitness=0.0)
        self.biomassId = biomass
        self.productId = product
        self.uptakeId = uptake
        self.method = kwargs.get('method', SimulationMethod.pFBA)
        self.reference = kwargs.get('reference', None)
        self.worst_fitness = 0.0
        self.obj_frac = kwargs.get('obj_frac', 0.99)

    def get_fitness(self, simul_results, candidate, **kwargs):
        """Evaluates a candidate.

        :param simul_results: (dic) A dictionary of phenotype SimulationResult objects
        :param candidate:  Candidate beeing evaluated
        :returns: A fitness value.

        """
        sim = simul_results[self.method] if self.method in simul_results.keys(
        ) else None
        if not sim or sim.status not in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL):
            return self.no_solution

        ssFluxes = sim.fluxes

        ids = list(ssFluxes.keys())
        if self.biomassId not in ids or self.productId not in ids:
            raise ValueError(
                "Biomass or product reactions ids are not present in the fluxes distribution.")

        if self.uptakeId and self.uptakeId in ids:
            uptake = abs(ssFluxes[self.uptakeId])
        else:
            uptake = 1.0

        if uptake == 0.0:
            return self.no_solution

        # computes target FVA min and max
        simulation = get_simulator(
            sim.model, envcond=sim.envcond, constraints=sim.model_constraints)
        v_min = 0
        v_max = 1
        constraints = sim.simulation_constraints
        # add biomass constraint
        biomassFluxValue = ssFluxes[self.biomassId] * self.obj_frac
        constraints[self.biomassId] = (biomassFluxValue, ModelConstants.REACTION_UPPER_BOUND)

        fvaMaxResult = simulation.simulate(
            objective={self.productId: 1}, constraints=constraints)
        v_max = fvaMaxResult.fluxes[self.productId]

        if not v_max:
            return self.worst_fitness

        fvaMinResult = simulation.simulate(
            objective={self.productId: 1}, constraints=constraints, maximize=False)
        v_min = fvaMinResult.fluxes[self.productId]

        if abs(v_max) == abs(v_min):
            return (ssFluxes[self.biomassId] * ssFluxes[self.productId]) / uptake
        else:
            return ((ssFluxes[self.biomassId] * ssFluxes[self.productId]) / uptake) * (
                1 - math.log(abs((v_max - v_min) / (v_max + v_min))))

    def required_simulations(self):
        return [self.method]

    def short_str(self):
        return "BPCY_FVA"

    def method_str(self):
        if self.uptakeId:
            return "BPCY_FVA (" + self.biomassId + " * " + self.productId + ") / " + self.uptakeId
        else:
            return "BPCY_FVA (" + self.biomassId + " * " + self.productId + ")"


class CNRFA(PhenotypeEvaluationFunction):

    def __init__(self, 
                 reactions:List[str],
                 threshold:float=0.1,
                 maximize:bool=True,
                 **kwargs):
        """Counts the Number of Reaction Fluxes Above a specified value.

        :param reactions: List of reactions
        :type reactions: List[str]
        :param threshold: the threshold, defaults to 0.1
        :type threshold: float, optional
        :param maximize: optimization direction: True (maximize) False (minimize), defaults to True
        :type maximize: bool, optional
        """
        super(CNRFA, self).__init__(maximize=maximize, worst_fitness=0)
        self.reactions = reactions
        self.method = kwargs.get('method', SimulationMethod.pFBA)
        self.theshold = threshold

    def get_fitness(self, simul_results, candidate, **kwargs):
        """Evaluates a candidate

        :param simul_results: (dic) A dictionary of phenotype SimulationResult objects
        :param candidate:  Candidate beeing evaluated
        :returns: A fitness value.

        """
        sim = simul_results[self.method] if self.method in simul_results.keys(
        ) else None
        if not sim or sim.status not in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL):
            return self.no_solution

        count = 0
        for rxn in self.reactions:
            if sim.fluxes[rxn]> self.theshold:
                count +=1
        return count

    def required_simulations(self):
        return [self.method]

    def short_str(self):
        return "CNRFA"

    def method_str(self):
        return "Count N Reaction Fluxes Above"



class MolecularWeight(PhenotypeEvaluationFunction):
    
    def __init__(self,
                 reactions:List[str],
                 maximize:bool=False, **kwargs):
        """ 
        Minimizes the sum of molecular weights of the products of a set of 
        reactions (g/gDW/h).

        :param reactions: List of reactions
        :type reactions: List[str]
        :param maximize: optimization direction: True (maximize) False (minimize), defaults to True
        :type maximize: bool, optional
        """
        super(MolecularWeight, self).__init__(maximize=maximize, worst_fitness=np.inf)
        self.reactions = reactions
        self.method = kwargs.get('method', SimulationMethod.pFBA)
        # sum of molar masses of product compounds for the reactions
        self.__mw = None

    def compute_rxnmw(self, model):
        from mewpy.util.constants import atomic_weights
        self.__mw = {}
        simulator = get_simulator(model)
        for rx in self.reactions:
            p = simulator.get_products(rx)
            if not p:
                p = simulator.get_substrates(rx)
            rmw = 0
            for m, v in p.items():
                elem = simulator.metabolite_elements(m)
                w = 0
                for e, n in elem.items():
                    try:
                        w += atomic_weights[e] * n
                    except:
                        pass

                rmw += abs(v) * w

            self.__mw[rx] = rmw

    def get_fitness(self, simul_results, candidate, **kwargs):
        try:
            sim = simul_results[self.method]
        except Exception:
            sim = None

        if sim.status not in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL):
            return self.no_solution

        if self.__mw is None:
            self.compute_rxnmw(sim.model)

        fitness = 0
        for rx in self.reactions:
            fitness += self.__mw[rx] * abs(sim.fluxes[rx])

        if fitness > 0:
            return fitness * 0.001
        else:
            return self.no_solution

    def required_simulations(self):
        """
        If the evaluation function requires a pre-simulation to compute fitness values
        """
        return [self.method]

    def short_str(self):
        return "Molecular Weight"

    def method_str(self):
        return "MW"


class FluxDistance(EvaluationFunction):

    def __init__(self, fluxes: Dict[str, float],
                 maximize: bool = False,
                 worst_fitness: float = 1000):
        """Minimizes the distance to a flux distribution

        :param fluxes: A dictionaty of flux distribution
        :type fluxes: Dict[str,float]
        :param maximize: optimization direction, defaults to False
        :type maximize: bool, optional
        :param worst_fitness: fitness to assign to bad solutions, defaults to inf
        :type worst_fitness: float, optional
        """
        super().__init__(maximize, worst_fitness)
        self.fluxes = fluxes
        self.method = 'pFBA'

    def get_fitness(self, simul_results, candidate, **kwargs):
        """Evaluates a candidate

        :param simul_results: (dic) A dictionary of phenotype SimulationResult objects
        :param candidate:  Candidate beeing evaluated
        :returns: A fitness value.
        """
        sim = simul_results[self.method] if self.method in simul_results.keys() else None

        if not sim or sim.status not in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL) or not sim.fluxes:
            return self.no_solution

        sum = 0
        for rxn in self.fluxes:
            sum += (sim.fluxes[rxn]-self.fluxes[rxn])**2
        return math.sqrt(sum)

    def required_simulations(self):
        return [self.method]

    def short_str(self):
        return "DISTFLUX"

    def method_str(self):
        return "Distance to a flux distribution"


class TargetFluxWithConstraints(EvaluationFunction):
    
    def __init__(self,
                 reaction:str,
                 constraints:Dict[str,Union[float,Tuple[float,float]]],
                 maximize:bool=False,
                 ):
        """_summary_

        :param reaction: _description_
        :type reaction: str
        :param constraints: _description_
        :type constraints: Dict[str,Union[float,Tuple[float,float]]]
        :param maximize: _description_, defaults to False
        :type maximize: bool, optional
        """
        super().__init__(maximize=maximize, worst_fitness=1000)
        self.reaction = reaction
        self.constraints = constraints

    def get_fitness(self, simul_results, candidate, **kwargs):
        """Evaluates a candidate

        :param simul_results: (dic) A dictionary of phenotype SimulationResult objects
        :param candidate:  Candidate beeing evaluated
        :returns: A fitness value.

        """
        simulator = kwargs.get('simulator', None)
        if simulator is None:
            return self.no_solution

        res = simulator.simulate(method='FBA', constraints=self.constraints)
        if res.status in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL):
            return res.fluxes[self.reaction]
        else:
            return self.no_solution

    def required_simulations(self):
        """
        If the evaluation function requires a pre-simulation to compute fitness values
        """
        return []

    def short_str(self):
        return "Target Flux with constraints"

    def method_str(self):
        return self.short_str()
