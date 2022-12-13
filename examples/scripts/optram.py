import os
from time import time

from mewpy.optimization import EA
from mewpy.optimization.evaluation import BPCY, WYIELD
from mewpy.problems.optram import OptRamProblem, load_optram
from mewpy.simulation import SimulationMethod


def test():
    """ An example on to use OptRAM optimization problems to find
    regulatory modification for the increased production of tryptophan.
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(dir_path, '../../../examples/models/optram/')
    gene_file = os.path.join(PATH, 'mgene.csv')
    ft_file = os.path.join(PATH, 'TFnames.csv')
    matrix_file = os.path.join(PATH, 'regnet.csv')
    model_file = os.path.join(PATH, 'yeast_7.6-optram.xml')

    BIOMASS_ID = 'r_2111'
    PRODUCT_ID = 'r_1912'
    GLC = 'r_1714'

    envcond = {GLC: (-10, 0)}

    # adds the prefix 'G_' to genes. Only for REFRAMED models
    regnet = load_optram(gene_file, ft_file, matrix_file, gene_prefix='G_')
    # the general objective is to maximize the target
    from reframed.io.sbml import load_cbmodel
    model = load_cbmodel(model_file)
    model.set_objective({BIOMASS_ID: 1})

    evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method=SimulationMethod.lMOMA)
    evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID, parsimonious=True)

    # OptRAM problem
    problem = OptRamProblem(model, [evaluator_1, evaluator_2],
                            regnet, envcond=envcond, candidate_min_size=10, candidate_max_size=30)

    print('Target List:', problem.target_list)
    print("\n\n")
    print('Metabolic Genes', problem.simulator.genes)
    print("\n\n")

    ea = EA(problem, max_generations=3, mp=True)
    ea.run()
    millis = int(round(time() * 1000))
    filename = "OPTRAM{}_OU_{}.csv".format('TRP', millis)
    df = ea.dataframe()
    df.to_csv(filename)


if __name__ == "__main__":
    for _ in range(1):
        test()
