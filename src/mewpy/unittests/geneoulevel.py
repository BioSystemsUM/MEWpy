"""
Finds expression levels for modifications suggested in the literature
"""
from mewpy.simulation.reframed import GeckoSimulation, Simulation
from mewpy.simulation import SimulationMethod
from mewpy.optimization.evaluation import WYIELD, BPCY, TargetFlux
from mewpy.optimization import EA, set_default_engine
from mewpy.model.gecko import GeckoModel
import mewpy.utils.utilities as utl
from collections import OrderedDict
from reframed.io.sbml import load_cbmodel
from reframed.core.cbmodel import CBModel
from random import Random
from time import time
import os


# Number of iterations 
ITERATIONS = 10 



def load_ec():
    # E. Coli
    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/ec/')
    DATA_FILE = os.path.join(PATH, "iJO1366SL.xml")
    BIOMASS_ID = 'R_Ec_biomass_iJO1366_core_53p95M'
    O2 = 'R_EX_o2_LPAREN_e_RPAREN_'
    GLC = 'R_EX_glc_LPAREN_e_RPAREN_'

    model = load_cbmodel(DATA_FILE, flavor='cobra')

    old_obj = model.get_objective().keys()
    obj = {key: 0 for key in old_obj if key != BIOMASS_ID}
    obj[BIOMASS_ID] = 1
    model.set_objective(obj)

    envcond = OrderedDict()
    envcond.update({GLC: (-10.0, 100000.0), O2: (-9.66, 100000.0)})
    return {'model': model, 'biomass': BIOMASS_ID, 'envcond': envcond}


def load_yeast():
    # Yeast
    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/yeast/')
    DATA_FILE = os.path.join(PATH, "iMM904SL_v6.xml")
    BIOMASS_ID = 'R_biomass_SC5_notrace'
    O2 = 'R_EX_o2_e_'
    GLC = 'R_EX_glc_e_'

    model = load_cbmodel(DATA_FILE, flavor='cobra')

    old_obj = model.get_objective().keys()
    obj = {key: 0 for key in old_obj if key != BIOMASS_ID}
    obj[BIOMASS_ID] = 1
    model.set_objective(obj)

    envcond = OrderedDict()
    envcond.update({GLC: (-10.0, 999999.0), O2: (-12.25, 100000.0)})

    simulation = Simulation(model, envcond=envcond)
    res = simulation.simulate(method=SimulationMethod.pFBA)
    reference = res.fluxes

    return {'model': model, 'biomass': BIOMASS_ID, 'envcond': envcond, 'reference': reference, 'non_target': []}


def cb_ou(product, modification_targets, chassis='ec', display=False, filename=None):
    
    if chassis == 'ec':
        conf = load_ec()
    elif chassis == 'ys':
        conf = load_yeast()
    else:
        raise ValueError

    BIOMASS_ID = conf['biomass']
    PRODUCT_ID = product
    model = conf['model']
    envcond = conf['envcond']
    evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method=SimulationMethod.lMOMA)
    evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)
    from mewpy.problems import GOUProblem
    problem = GOUProblem(model, fevaluation=[
                         evaluator_1, evaluator_2], envcond=envcond, candidate_min_size=len(modification_targets),candidate_max_size=len(modification_targets), target = modification_targets)

    ea = EA(problem, max_generations=ITERATIONS, visualizer=False)
    final_pop = ea.run()

    if display:
        individual = max(final_pop)
        best = list(problem.decode(individual.candidate).keys())
        print('Best Solution: \n{0}'.format(str(best)))

    if filename:
        print("Simplifying and saving solutions to file")
        utl.population_to_csv(problem, final_pop, filename, simplify=False)




if __name__ == '__main__':

    from mewpy.utils.constants import EAConstants
    # Define the number of CPUs to be used in parallel evaluations.
    EAConstants.NUM_CPUS = 10
    
    # change the EA engine and MOEA
    from mewpy.optimization import set_preferred_EA, set_default_engine
    set_preferred_EA('SPEA2')
    set_default_engine('jmetal')

    compounds_EC = {"PHE": "R_EX_phe_DASH_L_LPAREN_e_RPAREN_",
                    "TYR": "R_EX_tyr_DASH_L_LPAREN_e_RPAREN_",
                    "TRP": "R_EX_trp_DASH_L_LPAREN_e_RPAREN_"}

    compounds_YS = {"PHE": "R_EX_phe_L_e_",
                    "TYR": "R_EX_tyr_L_e_",
                    "TRY": "R_EX_trp_L_e_"
                   }
    
    # Can also be used to find expression values for
    # modifications found in the literature.
    millis = int(round(time() * 1000))
    modification_targets = ["G_YBR166C","G_YNL241C","G_YBR249C","G_YDR380W"]
    cb_ou("R_EX_tyr_L_e_", modification_targets ,chassis='ys', filename="CBMODEL_{}_OU_{}_.csv".format("Tyr", millis))
