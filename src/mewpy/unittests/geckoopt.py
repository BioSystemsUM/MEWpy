"""
##################################################################

GECKO E. coli model optimization

##################################################################
"""
from mewpy.model.gecko import GeckoModel
from mewpy.simulation.reframed import GeckoSimulation
from mewpy.simulation import SimulationMethod
from mewpy.simulation.simulation import SimulationResult
from mewpy.problems.gecko import GeckoKOProblem, GeckoOUProblem
from mewpy.optimization.evaluation import BPCY, WYIELD, TargetFlux
from mewpy.optimization import EA, set_default_engine
import mewpy.utils.utilities as utl
from collections import OrderedDict
from time import time


ITERATIONS = 100




def ec_gecko_ko(compound, display=False, filename=None):
    """ GECKO enzyme deletion example. 
    It runs a multi objective optimization for the increased production of a certain compound on E. Coli.
    The GECKO model is the yeast model companion from the GECKO paper "Improving the phenotype predictions
    of a yeast genome‐scale metabolic model by incorporating enzymatic constraints" https://doi.org/10.15252/msb.20167411. 
    Runs over the MEWpy implementation.
    
    :param compound: A target reaction identifier.
    :param display: Prints the best solution.
    :param filename: If given, saves the results as csv to filename.
         
    """

    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(dir_path, '../../../examples/models/gecko')
    DATA_FILE = os.path.join(PATH, 'eciML1515_batch.xml')
    from reframed.io.sbml import load_cbmodel
    m = load_cbmodel(DATA_FILE)
    
    model = GeckoModel(m, biomass_reaction_id='R_BIOMASS_Ec_iML1515_core_75p37M',protein_pool_exchange_id='R_prot_pool_exchange',reaction_prefix='R_')
    model.set_objective({'R_BIOMASS_Ec_iML1515_core_75p37M': 1.0})
    
    
    envcond = OrderedDict()

   
    # the evaluation (objective) functions
    
    evaluator_1 = BPCY("R_BIOMASS_Ec_iML1515_core_75p37M", compound,
                       method=SimulationMethod.lMOMA)
    evaluator_2 = WYIELD("R_BIOMASS_Ec_iML1515_core_75p37M", compound)
    # The optimization problem
    problem = GeckoKOProblem(model,
                              fevaluation=[evaluator_1, evaluator_2],
                              envcond=envcond,
                              prot_prefix='R_draw_prot_',
                              candidate_max_size=6)
    
    
    # A new instance of the EA optimizer
    ea = EA(problem, max_generations=ITERATIONS)
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
    

def ec_gecko_ou(compound, display=False, filename=None):
    """ GECKO enzyme over/under expression example. 
    It runs a multi objective optimization for the increased production of a certain compound on E. Coli.
    The GECKO model is the yeast model companion from the GECKO paper "Improving the phenotype predictions
    of a yeast genome‐scale metabolic model by incorporating enzymatic constraints" https://doi.org/10.15252/msb.20167411. 
    Runs over the MEWpy implementation.
    
    :param compound: A target reaction identifier.
    :param display: Prints the best solution.
    :param filename: If given, saves the results as csv to filename.
         
    """
    # define the default solver
    #from reframed.solvers import set_default_solver
    #set_default_solver('gurobi')

    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(dir_path, '../../../examples/models/gecko/')
    DATA_FILE = os.path.join(PATH, 'eciML1515_batch.xml')
    from reframed.io.sbml import load_cbmodel
    m = load_cbmodel(DATA_FILE)
    model = GeckoModel(m, biomass_reaction_id='R_BIOMASS_Ec_iML1515_core_75p37M',protein_pool_exchange_id='R_prot_pool_exchange',reaction_prefix='R_')
    model.set_objective({'R_BIOMASS_Ec_iML1515_core_75p37M': 1.0})
    
    # change protein pool bound (suggested by Leslie)
    model.reactions['R_prot_pool_exchange'].ub = 0.26 
    # define environmental consitions
    envcond = OrderedDict()

    # the evaluation (objective) functions
    evaluator_1 = BPCY("R_BIOMASS_Ec_iML1515_core_75p37M", compound,
                       method=SimulationMethod.MOMA)
    # FVA MAX is strangely very high... changing the default alpha (0.3) to compensate..
    evaluator_2 = WYIELD("R_BIOMASS_Ec_iML1515_core_75p37M", compound, alpha= 0.01 )

    evaluator_3 = TargetFlux(compound)
    
    # The optimization problem
    problem = GeckoOUProblem(model,
                              fevaluation=[evaluator_1, evaluator_2,evaluator_3],
                              envcond=envcond,
                              prot_prefix='R_draw_prot_',
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





if __name__ == '__main__':
    """
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
    """
    #from mewpy.utils.constants import EAConstants
    #EAConstants.DEBUG = True
    from reframed.solvers import set_default_solver
    set_default_solver('gurobi')
    ec_gecko_ou('R_EX_tyr__L_e',filename='gecko_ec_tyr.csv')