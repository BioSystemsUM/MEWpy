"""
CB model optimization test Set of using inspyred EA module
"""
from mewpy.simulation import SimulationMethod, get_simulator
from mewpy.optimization.evaluation import WYIELD, BPCY
from mewpy.optimization import EA, set_default_engine
import mewpy.utils.utilities as utl
from collections import OrderedDict
from reframed.io.sbml import load_cbmodel
from time import time
import os


ITERATIONS = 300
set_default_engine('jmetal')


def load_ec():
    # E. Coli
    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/ec/')
    DATA_FILE = os.path.join(PATH, "iJO1366SL.xml")
    NON_TARGET_FILE = os.path.join(
        PATH, "nontargets#RK#iJO1366SL#[lim-aerobic#glucose].txt")
    BIOMASS_ID = 'R_Ec_biomass_iJO1366_core_53p95M'
    O2 = 'R_EX_o2_LPAREN_e_RPAREN_'
    GLC = 'R_EX_glc_LPAREN_e_RPAREN_'

    model = load_cbmodel(DATA_FILE, flavor='cobra')

    old_obj = model.get_objective().keys()
    obj = {key: 0 for key in old_obj if key != BIOMASS_ID}
    obj[BIOMASS_ID] = 1
    model.set_objective(obj)

    non_target = [O2, GLC, 'R_ATPM']
    with open(NON_TARGET_FILE) as f:
        line = f.readline()
        while line:
            non_target.append(line.strip())
            line = f.readline()

    envcond = OrderedDict()
    envcond.update({GLC: (-10.0, 100000.0), O2: (-9.66, 100000.0)})

    simulation = get_simulator(model, envcond=envcond)
    res = simulation.simulate(method=SimulationMethod.pFBA)
    reference = res.fluxes

    return {'model': model, 'biomass': BIOMASS_ID, 'envcond': envcond, 'reference': reference, 'non_target': non_target}


def load_ec2():
    # E. Coli
    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/ec/')
    DATA_FILE = os.path.join(PATH, "iML1515.xml")
    NON_TARGET_FILE = os.path.join(
        PATH, "nontargets#RK#iJO1366SL#[lim-aerobic#glucose].txt")
    BIOMASS_ID = 'R_BIOMASS_Ec_iML1515_core_75p37M'
    O2 = 'R_EX_o2_e'
    GLC = 'R_EX_glc__D_e'

    model = load_cbmodel(DATA_FILE, flavor='cobra')

    old_obj = model.get_objective().keys()
    obj = {key: 0 for key in old_obj if key != BIOMASS_ID}
    obj[BIOMASS_ID] = 1
    model.set_objective(obj)

    non_target = [O2, GLC, 'R_ATPM']
    with open(NON_TARGET_FILE) as f:
        line = f.readline()
        while line:
            non_target.append(line.strip())
            line = f.readline()

    envcond = OrderedDict()
    envcond.update({GLC: (-10.0, 100000.0), O2: (-9.66, 100000.0)})

    simulation = get_simulator(model, envcond=envcond)
    res = simulation.simulate(method=SimulationMethod.pFBA)
    reference = res.fluxes

    res = simulation.simulate()
    print(res)

    return {'model': model, 'biomass': BIOMASS_ID, 'envcond': envcond, 'reference': reference, 'non_target': non_target}


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

    simulation = get_simulator(model, envcond=envcond)
    res = simulation.simulate(method=SimulationMethod.pFBA)
    reference = res.fluxes

    return {'model': model, 'biomass': BIOMASS_ID, 'envcond': envcond, 'reference': reference, 'non_target': []}


def cb_ou(product, chassis='ec', display=False, filename=None):
    "CBModel Reaction KO SO example"
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
    reference = conf['reference']

    evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method=SimulationMethod.lMOMA)
    evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)


    from mewpy.problems import GOUProblem
    problem = GOUProblem(model, fevaluation=[
                         evaluator_1, evaluator_2], envcond=envcond, reference=reference,
                         candidate_min_size=4, candidate_max_size=6,
                         operators=("lambda x,y: min(x,y)", "lambda x,y: max(x,y)"),
                         product = PRODUCT_ID)

    ea = EA(problem, max_generations=ITERATIONS, visualizer=False,algorithm='NSGAIII')
    final_pop = ea.run()

    if display:
        individual = max(final_pop)
        best = list(problem.decode(individual.candidate).keys())
        print('Best Solution: \n{0}'.format(str(best)))

    if filename:
        print("Simplifying and saving solutions to file")
        utl.population_to_csv(problem, final_pop, filename, simplify=False)


def cb_ko(product, chassis='ec', display=False, filename=None):
    "CBModel Reaction KO SO example"
    if chassis == 'ec':
        conf = load_ec()
    elif chassis == 'ys':
        conf = load_yeast()
    else:
        raise ValueError

    BIOMASS_ID = conf['biomass']
    PRODUCT_ID = product
    model = conf['model']
    non_target = conf['non_target']
    envcond = conf['envcond']
    reference = conf['reference']

    evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method=SimulationMethod.lMOMA)
    evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)
    from mewpy.problems.genes import GKOProblem
    problem = GKOProblem(model, fevaluation=[
                         evaluator_1, evaluator_2], non_target=non_target, envcond=envcond, reference=reference)

    ea = EA(problem, max_generations=ITERATIONS, mp=True)
    final_pop = ea.run()

    if display:
        individual = max(final_pop)
        best = list(problem.decode(individual.candidate).keys())
        print('Best Solution: \n{0}'.format(str(best)))

    if filename:
        print("Simplifying and saving solutions to file")
        utl.population_to_csv(problem, final_pop, filename, simplify=False)


if __name__ == '__main__':
    RUNS = 10
    compounds_EC = {"PHE": "R_EX_phe_DASH_L_LPAREN_e_RPAREN_",
                    "TYR": "R_EX_tyr_DASH_L_LPAREN_e_RPAREN_",
                    "TRP": "R_EX_trp_DASH_L_LPAREN_e_RPAREN_"}

    #compounds_EC = {"TYR": "R_EX_tyr__L_e"}
    compounds_YS = {# "PHE": "R_EX_phe_L_e_",
                    "TYR": "R_EX_tyr_L_e_"
                    # "TRY": "R_EX_trp_L_e_"
                   }

    """
    for k, v in compounds_EC.items():
        for i in range(RUNS):
            millis = int(round(time() * 1000))
            cb_ko(v, filename="CBMODEL_{}_KO_{}.csv".format(k, millis))
    """
    for k, v in compounds_EC.items():
        for i in range(RUNS):
            millis = int(round(time() * 1000))
            cb_ou(v, filename="CBMODEL_{}_OU_{}.csv".format(k, millis))
    """
    for k, v in compounds_YS.items():
        for i in range(RUNS):
            millis = int(round(time() * 1000))
            cb_ko(v, chassis='ys', filename="CBMODEL_{}_KO_{}.csv".format(k, millis))
    """

    # from reframed.solvers import set_default_solver
    # set_default_solver('gurobi')
    """
    EAConstants.DEBUG = True
    for k, v in compounds_EC.items():
        for i in range(RUNS):
            millis = int(round(time() * 1000))
            cb_ou(v, chassis='ec', filename="CBMODEL_{}_OU_{}_.csv".format(k, millis))
    """
