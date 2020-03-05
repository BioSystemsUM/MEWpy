"""
##################################################################

GECKO model optimization test using JMetalPy

##################################################################
"""
from mewpy.model.gecko import GeckoModel
from mewpy.simulation.reframed import GeckoSimulation
from mewpy.simulation import SimulationMethod
from mewpy.simulation.simulation import SimulationResult
from mewpy.optimization.evaluation import BPCY,WYIELD,TargetFlux
from mewpy.optimization.jmetal.problem import GeckoRKOProblem, GeckoROUProblem
from mewpy.optimization.jmetal.ea import EA
import mewpy.utils.utilities as utl
from collections import OrderedDict
from time import time






def geckomo(display=False, filename = None):
    """ Gecko OU MO example
    """
    model = GeckoModel('single-pool', biomass_reaction_id='r_2111')
    model.set_objective({'r_2111': 1.0 , 'r_4041':0.0})
    envcond = OrderedDict()
    
    # wild type reference values
    simulation = GeckoSimulation(model,envcond= envcond)
    res = simulation.simulate(method=SimulationMethod.pFBA)
    reference = res.fluxes
    # min_biomass = res.fluxes['r_2111'] * 0.1

    additional_constraints = OrderedDict()
    # imposes a minimal biomass 
    # additional_constraints["r_2111"]=(min_biomass,10000)


    # the evaluation (objective) functions
    evaluator_1 = WYIELD("r_2111", "r_1913", parsimonious = True)
    evaluator_2 = BPCY("r_2111", "r_1913", uptake = "r_1714_REV", method=SimulationMethod.lMOMA ,reference=reference)
    evaluator_3 = TargetFlux("r_1913", biomass = "r_2111")
    
    # The optimization problem
    # Notes:
    #  - A scale factor for the LP can be defined by setting the 'scalefactor' acordingly. 
    #  - The scale factor is only used in the solver context and all results are scale free.   
    problem = GeckoRKOProblem(model, 
                              fevaluation=[evaluator_1,evaluator_2,evaluator_3], 
                              envcond = envcond, 
                              constraints = additional_constraints,
                              reference = reference, 
                              scalefactor = None, 
                              candidate_max_size = 30)

    # A new instance of the EA optimizer
    ea = EA(problem, max_generations = 500)
    # runs the optimization
    final_pop = ea.run()
    # optimization results
    if display:
        individual = max(final_pop)
        best = list(problem.decode(individual.candidate).keys())
        print('Best Solution: \n{0}'.format(str(best)))
    # save final population to file
    if filename:
       print("Simplifying and saving solutions to file") 
       utl.population_to_csv(problem,final_pop,filename, simplify= True)




if __name__ == '__main__':
 
    millis = int(round(time() * 1000))
    geckomo(display=True, filename="gecko_SPEA_KO_{}.csv".format(millis))
    