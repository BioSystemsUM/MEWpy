from mewpy.regulation.optram import OptRAMRegModel, OptRamProblem, load_optram
from mewpy.simulation import SimulationMethod, get_simulator
from mewpy.optimization.evaluation import BPCY_FVA, BPCY, WYIELD
from mewpy.optimization import EA
from time import time
import os


def test():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(dir_path, '../../../examples/models/optram/')
    gene_file = os.path.join(PATH, 'mgene.csv')
    ft_file = os.path.join(PATH, 'TFnames.csv')
    matrix_file = os.path.join(PATH, 'regnet.csv')
    model_file = os.path.join(PATH, 'yeast_7.6-optram.xml')

    BIOMASS_ID = 'r_2111'
    PRODUCT_ID = 'r_1912'
    GLC = 'r_1714'

    envcond = {GLC:(-10,0)}

    # adds the prefix 'G_' to genes. Only for REFRAMED models
    regnet = load_optram(gene_file, ft_file, matrix_file, gene_prefix='G_')
    #from cobra.io import read_sbml_model
    #model = read_sbml_model(model_file)
    
    # the general objective is to maximize the target
    from reframed.io.sbml import load_cbmodel
    model = load_cbmodel(model_file)
    model.set_objective({BIOMASS_ID:1})
    
    evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method=SimulationMethod.lMOMA)
    evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID, parsimonious=True)
    
    # OptRAM problem
    problem = OptRamProblem(model, [evaluator_1, evaluator_2],
                            regnet, envcond = envcond, candidate_min_size=10, candidate_max_size=30)

    #print(problem.target_list)
    #print("\n\n")
    #print(problem.simulator.genes)    
    #print("\n\n")
    #print(set(problem.target_list).intersection(set(problem.simulator.genes)))


    ea = EA(problem, max_generations=3, mp=True)
    final_pop = ea.run()
    import mewpy.utils.utilities as utl
    millis = int(round(time() * 1000))
    filename = "OPTRAM{}_OU_{}.csv".format('TRP', millis)
    utl.population_to_csv(problem, final_pop, filename, simplify=False)


if __name__ == "__main__":
    for _ in range(1):
        test()
