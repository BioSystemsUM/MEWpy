from mewpy.regulation.optram import *
from reframed.io.sbml import load_cbmodel
from mewpy.simulation import SimulationMethod
from mewpy.optimization.evaluation import BPCY_FVA, BPCY,WYIELD
from mewpy.optimization.ea import EA
import os






def test():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(dir_path,'../../../examples/models/')
    gene_file = os.path.join(PATH, 'mgene.csv')
    ft_file = os.path.join(PATH, 'TFnames.csv')
    matrix_file = os.path.join(PATH, 'regnet.csv')
    model_file = os.path.join(PATH, 'iMM904SL_v6.xml')
    
    BIOMASS_ID = 'R_biomass_SC5_notrace'
    PRODUCT_ID = 'R_EX_succ_e_'
    GLC = 'R_EX_glc_e_'

    
    # adds the prefix 'G_' to genes. Only for REFRAMED models
    regnet = load_optram(gene_file,ft_file,matrix_file, gene_prefix = 'G_')
    
    model = load_cbmodel(model_file, flavor='cobra')
    # the general objective is to maximize the target
    old_obj = model.get_objective().keys()
    obj = { key: 0 for key in old_obj if key != BIOMASS_ID }
    obj[BIOMASS_ID] = 1
    model.set_objective(obj)
        
    evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method = SimulationMethod.lMOMA )
    evaluator_2 = WYIELD(BIOMASS_ID,PRODUCT_ID, parsimonious = True)
    # evaluation function used in the OptRAM paper. 
    evaluator_3 = BPCY_FVA(BIOMASS_ID, PRODUCT_ID, uptake= GLC , method = SimulationMethod.pFBA)
    
    problem = OptRamProblem(model,[evaluator_1, evaluator_2, evaluator_3],regnet, candidate_min_size = 4 , candidate_max_size = 6)
    #print(problem.target_list)
    
    #A = set(model.genes.keys())
    #B = set(regnet.genes.keys())
    #C = A.intersection(B)
    #print(A-B)
    #print(B-A)
    #print(len(C))
   
    
    
    ea = EA(problem, max_generations= 100, mp = True)
    final_pop = ea.run()
    
    
if __name__ == "__main__":
    
    test()