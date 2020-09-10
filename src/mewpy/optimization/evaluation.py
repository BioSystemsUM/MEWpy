from abc import ABCMeta, abstractmethod
from itertools import chain
from collections import OrderedDict
from mewpy.utils.constants import EAConstants
from mewpy.simulation import get_simulator, SimulationMethod, SStatus
from mewpy.utils.ode import ODEStatus
import numpy as np
from functools import reduce
import math
import warnings


class EvaluationFunction:
   
    __metaclass__ = ABCMeta

    def __init__(self, maximize=True, worst_fitness=0.0):
        """This abstract class should be extended by all evaluation functions.
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
        
        return

    @abstractmethod
    def method_str(self):
        return

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
            res = EAConstants.MIN_THRESHOLD
        else:
            res = EAConstants.MAX_THRESHOLD
        return res

    def __call__(self, simulationResult, candidate, **kwargs):
        return self.get_fitness(simulationResult, candidate, **kwargs)


class PhenotypeEvaluationFunction(EvaluationFunction):
    
    def __init__(self, maximize=True, worst_fitness=0.0):
        super(PhenotypeEvaluationFunction, self).__init__(maximize=maximize, worst_fitness=0.0)


class KineticEvaluationFunction(EvaluationFunction):
    
    def __init__(self, maximize=True, worst_fitness=0.0):
        super(KineticEvaluationFunction, self).__init__(maximize=maximize, worst_fitness=0.0)


class TargetFlux(PhenotypeEvaluationFunction):
    """ Target Flux evaluation function. 
    The fitness value is the flux value of the identified reaction.
    If the reaction parameter is None, the fitness value is the optimization objective value.
    Additional parameters include a minimum of allowed biomass value computed from the min_biomass_per and reference flux values

    
    :param reaction: (str) The reaction identifier whose flux value is to be used as fitness. Default None in which case the model objective is considered.
    :param biomass: (str) The biomass reaction identifier.
    :param maximize: (boolean) The optimization direction. Default True for maximization.
    :param min_biomass_value: (float) The minimum biomass value.
    :param min_biomass_per: (float) Minimum biomass percentage. Only used if no min_biomass_value is provided.
    
    """

    def __init__(self, reaction, biomass=None, maximize=True, min_biomass_value=None, min_biomass_per=0.0, method=SimulationMethod.pFBA):
        super(TargetFlux, self).__init__(maximize=maximize, worst_fitness=0.0)
        self.reaction = reaction
        self.biomass = biomass
        self.min_biomass_value = min_biomass_value
        self.min_biomass_per = min_biomass_per
        self.method = method

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

    def required_simulations(self):
        """
        If the evaluation function requires a pre-simulation to compute fitness values
        """
        return [self.method]

    def short_str(self):
        return "TargetFlux"

    def method_str(self):
        return "TargetFlux {} with at least {} of biomass ({})".format(self.reaction, self.min_biomass_per, self.biomass)


class WYIELD (PhenotypeEvaluationFunction):
    """ Weighted Yield (WYIELD) objective function, a linear combination of the target product minimum and maximum FVA under the introduced metabolic modifications.

    :param biomassId: (str) Biomass reaction identifier.
    :param productId: (str) Target product reaction identifier.
    
    kwargs options:
    
    :param min_biomass_value: (float) Minimum biomass value (default None, in which case the min_biomass_per is used).
    :param min_biomass_per: (float) Instead of defining explicitly a minimum biomass value, a percentage of the wild type biomass is used. Only used when no min_biomass_value is defined (default min_biomass_per: 0.10).                          
    :param alpha: (float) Tradeoff between the Max and min FVA of the target product (alpha Max + (1-alpha) min). Must be in range [0,1]  (default alpha: 0.3).
    :param scale: (boolean) Defines if the WYIELD is devided by the biomass of the simulated result (default false).  
    
    """

    def __init__(self, biomassId, productId, maximize=True, **kwargs):
        super(WYIELD, self).__init__(maximize=maximize, worst_fitness=0.0)
        self.biomassId = biomassId
        self.productId = productId
        # parameters
        self.min_biomass_value = kwargs.get('min_biomass_value', None)
        self.min_biomass_per = kwargs.get('min_biomass_per', 0.1)
        self.scale = kwargs.get('scale', False)
        self.method = SimulationMethod.FBA
        self.alpha = kwargs.get('alpha', 0.3)
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
                "Reaction ids are not present in the fluxes distribution. Please check if the objective function ids are correct.")

        biomassFluxValue = ssFluxes[self.biomassId] * 0.999

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

            constraints = sim.simulation_constraints
            # add biomass constraint
            constraints[self.biomassId] = (biomassFluxValue, 100000.0)

            # only need to simulate FVA max if alpha is larger than 0, otherwise it will always be zero
            if(self.alpha > 0):
                fvaMaxResult = simulation.simulate(
                    objective={self.productId: 1}, constraints=constraints, scalefactor=scalefactor)
                if fvaMaxResult.status == SStatus.OPTIMAL:
                    fvaMaxProd = fvaMaxResult.fluxes[self.productId]
                else:
                    return self.no_solution

            # only need to simulate FVA min if alpha is lesser than 1, otherwise it will always be zero
            if(self.alpha < 1):
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
        except:
            return self.no_solution

    def required_simulations(self):
        return [self.method]

    def short_str(self):
        return "WYIELD"

    def method_str(self):
        return "WYIELD biomass: {} product: {} ".format(self.biomassId, self.productId)


class BPCY(PhenotypeEvaluationFunction):
    """
    This class implements the "Biomass-Product Coupled Yield" objective function. The fitness is given by the equation:
    (biomass_flux * product_flux)/ uptake_flux

    
    :param biomass: (str) Biomass reaction identifier
    :param product: (str) Target product reaction identifier
    :param uptake: (str) (optional) Reaction of uptake. If no substract is defined, ie uptake is None, a substract flux value of 1.0 is considered.  

    kargs options:
    
    :param method: (SimulationMethod) The simulation method. Default Node in which case received simulation results are used to compute the biomass product coupled yield.
    :param reference: (dic) Wild type reference values when MOMA, lMOMA or ROOM are defined as method. If not provided, wild type reference values will be computed.
    
    """

    def __init__(self, biomass, product, uptake=None,  maximize=True, **kwargs):
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
            except:
                print("BPCY No Fluxes")
        return (ssFluxes[self.biomassId] * ssFluxes[self.productId]) / uptake

    def required_simulations(self):
        return [self.method]

    def short_str(self):
        return "BPCY"

    def method_str(self):
        if self.uptakeId:
            return "BPCY (" + self.biomassId + " * " + self.productId + ") / " + self.uptakeId
        else:
            return "BPCY (" + self.biomassId + " * " + self.productId + ")"


class BPCY_FVA (PhenotypeEvaluationFunction):
    """
    This class implements the "Biomass-Product Coupled Yield" objective function with FVA as defined in 
    "OptRAM: In-silico strain design via integrative regulatory-metabolic network modeling". It combines BPCY with WYIELD objective functions. 

    The fitness is given by the equation:
    ((biomass_flux * product_flux) * / uptake_flux ) * (1-log((range)/(target))) where range=(FVA_max-FVA_min)/2 and target= (FVA_max+FVA_min)/2
    
    :param biomass: (str) Biomass reaction identifier.
    :param product: (str) Target product reaction identifier.
    :param uptake: (str) (optional) Reaction of uptake. If no substract is defined, ie uptakeId is None, a substract flux value of 1.0 is considered.  

    kwargs options:
    
    :param method: (SimulationMethod) The simulation method. Default Node in which case received simulation results are used to compute the biomass product coupled yield.
    :param reference: (dic) Wild type reference values when MOMA, lMOMA or ROOM are defined as method. If not provided, wild type reference values will be computed.
    
    """

    def __init__(self, biomass, product, uptake=None,  maximize=True, **kwargs):
        super(BPCY_FVA, self).__init__(maximize=maximize, worst_fitness=0.0)
        self.biomassId = biomass
        self.productId = product
        self.uptakeId = uptake
        self.method = kwargs.get('method', SimulationMethod.pFBA)
        self.reference = kwargs.get('reference', None)
        self.worst_fitness = 0.0

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
        biomassFluxValue = ssFluxes[self.biomassId] * 0.999
        constraints[self.biomassId] = (biomassFluxValue, 100000.0)

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
            return ((ssFluxes[self.biomassId] * ssFluxes[self.productId]) / uptake) * (1-math.log(abs((v_max-v_min)/(v_max+v_min))))

    def required_simulations(self):
        return [self.method]

    def short_str(self):
        return "BPCY_FVA"

    def method_str(self):
        if self.uptakeId:
            return "BPCY_FVA (" + self.biomassId + " * " + self.productId + ") / " + self.uptakeId
        else:
            return "BPCY_FVA (" + self.biomassId + " * " + self.productId + ")"


class AggregatedSum(PhenotypeEvaluationFunction,KineticEvaluationFunction):
    """
    Aggredated sum evaluation function. Used to converte MOEAs into Single Objective EAs. 

        
    :param fevaluation: (list) List of evaluation functions.
    :param tradeoffs: (list) Tradeoff values for each evaluation function. If None, all functions have the same associated weight.

    """

    def __init__(self, fevaluation, tradeoffs=None, maximize=True):
        super(AggregatedSum, self).__init__(
            maximize=maximize, worst_fitness=0.0)
        self.fevaluation = fevaluation
        if tradeoffs and len(tradeoffs) == len(fevaluation):
            self.tradeoffs = tradeoffs
        else:
            self.tradeoffs = [1/len(self.fevaluation)] * \
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
        return "Aggregated Sum = " + reduce(lambda a, b: a+" "+b, [f.method_str() for f in self.fevaluation], "")


class MinCandSize(PhenotypeEvaluationFunction, KineticEvaluationFunction):
    """
    This class implements the "minimization of number of reactions" objective function. The fitness is given by
    1 - size(candidate)/ candidate_max_size, where the candidate_max_size is the maximum size that a candidate can have
    during optimization.

    :param maxCandidateSize: (int) Maximum size allowed for candidate.
    
    """

    def __init__(self, candidate_max_size=EAConstants.MAX_SOLUTION_SIZE, maximize=True):
        super(MinCandSize, self).__init__(maximize=maximize, worst_fitness=0.0)
        self.candidate_max_size = candidate_max_size

    def get_fitness(self, simulResult, candidate, **kwargs):
        return 1 - len(candidate)/self.candidate_max_size

    def required_simulations(self):
        """
        If the evaluation function requires a pre-simulation to compute fitness values
        """
        return []

    def short_str(self):
        return "MinCandSize"

    def method_str(self):
        return "Minimizes the number of alterations"



"""
Kinetic specific evaluation functions
"""


class KineticTargetFlux(KineticEvaluationFunction):
    """ Target Flux evaluation function. 
    The fitness value is the flux value of the identified reaction.

    
    :param reaction: (str) The reaction identifier whose flux value is to be used as fitness. Default None in which case the model objective is considered.
    
    """

    def __init__(self, reaction, maximize = True):
        super(KineticTargetFlux, self).__init__(maximize=maximize, worst_fitness=0.0)
        self.reaction = reaction
       
    def get_fitness(self, simul_results, candidate, **kwargs):
        if simul_results.status == ODEStatus.ERROR:
            return self.worst_fitness
        else:
            try :
                return simul_results.fluxes[self.reaction]
            except Exception as e:
                print(e)
                return self.worst_fitness

    def required_simulations(self):
        """
        If the evaluation function requires a pre-simulation to compute fitness values
        """
        return []

    def short_str(self):
        return "TargetFlux"

    def method_str(self):
        return "TargetFlux"

