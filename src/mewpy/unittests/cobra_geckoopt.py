"""
##################################################################

GECKO model over COBRApy optimization test

##################################################################
"""
from geckopy.gecko import GeckoModel
from mewpy.simulation.cobra import GeckoSimulation
from mewpy.simulation import SimulationMethod
from mewpy.simulation.simulation import SimulationResult
from mewpy.problems.gecko import GeckoRKOProblem, GeckoROUProblem
from mewpy.optimization.evaluation import BPCY, WYIELD, TargetFlux
from mewpy.optimization import EA
import mewpy.utils.utilities as utl
from collections import OrderedDict
from time import time

import warnings


ITERATIONS = 100


def gecko_ko(compound, display=False, filename=None):
    """ Gecko OU MO example
    """

    warnings.filterwarnings("ignore")

    model = GeckoModel('single-pool', biomass_reaction_id='r_2111')
    model.solver = 'cplex'
    model.objective = {model.reactions.get_by_id('r_2111'): 1.0}
    envcond = OrderedDict()

    # the evaluation (objective) functions
    evaluator_1 = WYIELD("r_2111", compound, parsimonious=True)
    evaluator_2 = BPCY("r_2111", compound, uptake="r_1714_REV",
                       method=SimulationMethod.lMOMA)

    # The optimization problem
    # Notes:
    #  - A scale factor for the LP can be defined by setting the 'scalefactor' acordingly.
    #  - The scale factor is only used in the solver context and all results are scale free.
    problem = GeckoRKOProblem(model,
                              fevaluation=[evaluator_1, evaluator_2],
                              envcond=envcond,
                              scalefactor=None,
                              candidate_max_size=30)

    # A new instance of the EA optimizer
    ea = EA(problem, max_generations=ITERATIONS, mp=True)
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
        utl.population_to_csv(problem, final_pop, filename, simplify=False)


def gecko_ou(compound, display=False, filename=None):
    """ Gecko OU MO example
    """

    model = GeckoModel('single-pool', biomass_reaction_id='r_2111')
    model.objective = 'r_2111'
    envcond = OrderedDict()

    # wild type reference values
    simulation = GeckoSimulation(model, envcond=envcond)
    res = simulation.simulate(method=SimulationMethod.pFBA)
    reference = res.fluxes

    # the evaluation (objective) functions
    evaluator_1 = WYIELD("r_2111", compound, parsimonious=False)
    evaluator_2 = BPCY("r_2111", compound, uptake="r_1714_REV",
                       method=SimulationMethod.lMOMA, reference=reference)

    # The optimization problem
    # Notes:
    #  - A scale factor for the LP can be defined by setting the 'scalefactor' acordingly.
    #  - The scale factor is only used in the solver context and all results are scale free.
    problem = GeckoROUProblem(model,
                              fevaluation=[evaluator_1, evaluator_2],
                              envcond=envcond,
                              reference=reference,
                              scalefactor=None,
                              candidate_max_size=30)

    # A new instance of the EA optimizer
    ea = EA(problem, max_generations=ITERATIONS, mp=False)
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
        utl.population_to_csv(problem, final_pop, filename, simplify=False)


if __name__ == '__main__':

    N_EXP = 1

    compounds = {'SUC': 'r_2056',
                 'TYR': 'r_1913',
                 'PHE': 'r_1903',
                 'TRY': 'r_1912'}

    for k, v in compounds.items():
        for _ in range(N_EXP):
            millis = int(round(time() * 1000))
            gecko_ou(v, display=False,
                     filename="gecko_{}_OU_{}.csv".format(k, millis))

    for k, v in compounds.items():
        for _ in range(N_EXP):
            millis = int(round(time() * 1000))
            gecko_ko(v, display=False,
                     filename="gecko_{}_KO_{}.csv".format(k, millis))
