from mewpy.model.smoment import SMomentModel
from mewpy.simulation import get_simulator
from mewpy.problems import GeckoOUProblem
from mewpy.utils.utilities import population_to_csv
import os


def test1(compound = 'R_EX_tyr__L_e'):
    """
    AutoPACMEN, "Automatic construction of metabolic models with enzyme constraints" (https://doi.org/10.1186/s12859-019-3329-9), is able to construct GECKO like models.
    This example optimizes the production of a compound using the E. coli autoPACMEN model where enzymes are defined as pseudo reactants.
    
    """
    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/autopacmen/')
    DATA_FILE = os.path.join(PATH, "iJO1366_2019_06_25_GECKO.xml")

    model = SMomentModel(DATA_FILE, enzyme_reaction_prefix='R__TG_ER_')

    print(model.proteins)

    sim = get_simulator(model)
    sim.set_objective("R_BIOMASS_Ec_iJO1366_core_53p95M")
    solution = sim.simulate()
    print(solution)
    print('Wild type tyrosine production :',solution.fluxes['R_EX_tyr__L_e'])
    
    from mewpy.optimization.evaluation import BPCY,WYIELD
    from mewpy.simulation import SimulationMethod  
    from collections import OrderedDict

    envcond = OrderedDict()
    envcond.update({'R_EX_glc__D_e': (-10.0, 100000.0)})


   
    # the evaluation (objective) functions
    evaluator_1 = BPCY("R_BIOMASS_Ec_iJO1366_core_53p95M", compound, method=SimulationMethod.lMOMA)
    evaluator_2 = WYIELD("R_BIOMASS_Ec_iJO1366_core_53p95M", compound)
    
    
    problem = GeckoOUProblem(model,
                              fevaluation=[evaluator_1, evaluator_2],
                              envcond=envcond,
                              candidate_max_size=6,
                              prot_prefix='R__TG_ER_')

    # A new instance of the EA optimizer
    from mewpy.optimization import EA
    ea = EA(problem, max_generations=500)
    # runs the optimization
    final_pop = ea.run()

    import mewpy.utils.utilities as utl
    from time import time
    
    millis = int(round(time() * 1000))
    filename = "sMOMEMT{}_OU_{}.csv".format(compound, millis)
    utl.population_to_csv(problem, final_pop, filename, simplify=False)
    
    


def test2(compoud = 'R_EX_tyr__L_e',filename = None):
    """
    AutoPACMEN, "Automatic construction of metabolic models with enzyme constraints" (https://doi.org/10.1186/s12859-019-3329-9), is able to construct GECKO like models.
    This example optimizes the production of a compound using the E. coli autoPACMEN model where enzymes are defined as pseudo reactants.
    The model defines a linear constraint over the protein pool as reactant, by adding the protein pool as a metabolite in the stochiometric matrix.
    Therefore, the model may be treated as a regular GSM. 

    """
    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/autopacmen/')
    DATA_FILE = os.path.join(PATH, "iJO1366_sMOMENT_2019_06_25.xml")
    from reframed.io.sbml import load_cbmodel
    model = load_cbmodel(DATA_FILE)
    
    #from cobra.io import read_sbml_model
    #model = read_sbml_model(DATA_FILE)
    # increases the enzyme pool upper bound
    #model.reactions.get_by_id('R_ER_pool_TG_').upper_bound = 1

    sim = get_simulator(model)
    s = sim.simulate()
    print(s)
    #print(s.fluxes)
    print('Wildtype tyrosine production :',s.fluxes['R_EX_tyr__L_e'])
    print('Pool:',s.fluxes['R_ER_pool_TG_'])
    
    fva = sim.FVA(reactions=['R_ER_pool_TG_'])
    print(fva)
    
    #for r_id, reaction in model.reactions.items():
    #    if reaction.gpr:
    #        print(str(reaction.gpr))

    #for rx in model.reactions:
    #        reaction = model.reactions.get_by_id(rx.id)
    #        if reaction.gene_reaction_rule:
    #           print(str(reaction.gene_reaction_rule))
        

    # implements a knockout over genes that encode enzymes
    BIOMASS_ID = 'R_BIOMASS_Ec_iJO1366_WT_53p95M'
    PRODUCT_ID = compoud
    from mewpy.optimization.evaluation import BPCY,WYIELD,TargetFlux
    from mewpy.simulation import SimulationMethod  
    evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method=SimulationMethod.lMOMA)
    evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)
    # evaluator_3 = TargetFlux(PRODUCT_ID)
    from mewpy.problems.genes import GOUProblem
    problem = GOUProblem(model, fevaluation=[
                         evaluator_1, evaluator_2], max_candidate_size =30)
    print(problem.target_list)
    from mewpy.optimization import EA
    ea = EA(problem, max_generations=20, mp=True)
    final_pop = ea.run()
    if filename:
        print("Simplifying and saving solutions to file")
        population_to_csv(problem, final_pop, filename, simplify=False)



if  __name__ == "__main__":
    for i in range(1):
        test1()
