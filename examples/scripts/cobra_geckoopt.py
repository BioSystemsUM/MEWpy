"""
##################################################################

GECKO model over COBRApy optimization test

##################################################################
"""
import warnings
from collections import OrderedDict
from time import time
from geckopy.gecko import GeckoModel
from mewpy.optimization import EA
from mewpy.optimization.evaluation import BPCY, WYIELD
from mewpy.problems.gecko import GeckoKOProblem, GeckoOUProblem
from mewpy.simulation import SimulationMethod


# This should be inscreased. Only for illustrative purposes
ITERATIONS = 10


def gecko_ko(compound, display=False, filename=None):
    """ GECKO enzyme deletion example.
    It runs a multi objective optimization for the increased production of a certain compound on yeast.
    The GECKO model is the yeast model companion from the GECKO paper "Improving the phenotype predictions
    of a yeast genome‐scale metabolic model by incorporating enzymatic constraints"
    https://doi.org/10.15252/msb.20167411.
    Runs over the GECKO original implementation.

    :param compound: A target reaction identifier.
    :param display: Prints the best solution.
    :param filename: If given, saves the results as csv to filename.

    """

    warnings.filterwarnings("ignore")

    model = GeckoModel('single-pool', biomass_reaction_id='r_2111')
    model.objective = {model.reactions.get_by_id('r_2111'): 1.0}
    envcond = OrderedDict()

    # the evaluation (objective) functions
    evaluator_1 = WYIELD("r_2111", compound)
    evaluator_2 = BPCY("r_2111", compound, uptake="r_1714_REV",
                       method=SimulationMethod.lMOMA)

    # The optimization problem
    problem = GeckoKOProblem(model,
                             fevaluation=[evaluator_1, evaluator_2],
                             envcond=envcond,
                             scalefactor=None,
                             candidate_max_size=30)

    # A new instance of the EA optimizer
    ea = EA(problem, max_generations=ITERATIONS)
    # runs the optimization
    ea.run()
    # optimization results
    if display:
        print(ea.dataframe())
    # save final population to file
    if filename:
        print("Saving solutions to file")
        ea.dataframe().to_csv(filename)


def gecko_ou(compound, display=False, filename=None):
    """ GECKO enzyme over/under expressoion example.
    It runs a multi objective optimization for the increased production of a certain compound on yeast.
    The GECKO model is the yeast model companion from the GECKO paper "Improving the phenotype predictions
    of a yeast genome‐scale metabolic model by incorporating enzymatic
    constraints" https://doi.org/10.15252/msb.20167411.
    Runs over the GECKO original implementation.

    :param compound: A target reaction identifier.
    :param display: Prints the best solution.
    :param filename: If given, saves the results as csv to filename.

    """

    model = GeckoModel('single-pool', biomass_reaction_id='r_2111')
    model.objective = 'r_2111'
    envcond = OrderedDict()

    # the evaluation (objective) functions
    evaluator_1 = WYIELD("r_2111", compound)
    evaluator_2 = BPCY("r_2111", compound, uptake="r_1714_REV",
                       method='lMOMA')

    # The optimization problem
    problem = GeckoOUProblem(model,
                             fevaluation=[evaluator_1, evaluator_2],
                             envcond=envcond,
                             candidate_max_size=30)

    # A new instance of the EA optimizer
    ea = EA(problem, max_generations=ITERATIONS)
    # runs the optimization
    ea.run()
    # optimization results
    if display:
        print(ea.dataframe())
    # save final population to file
    if filename:
        print("Saving solutions to file")
        ea.dataframe().to_csv(filename)


if __name__ == '__main__':

    N_EXP = 1

    compounds = {'SUC': 'r_2056'}

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
