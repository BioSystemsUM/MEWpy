import os

from mewpy.optimization import EA
from mewpy.optimization.evaluation import BPCY, WYIELD
from mewpy.problems import OptORFProblem


def optorf_imc():
    """

    From Kim, J., Reed, J.L. OptORF: Optimal metabolic and regulatory perturbations for metabolic engineering of
    microbial strains. BMC Syst Biol 4, 53 (2010). https://doi.org/10.1186/1752-0509-4-53

    In the simulations, we excluded transport reactions for acetate, carbon dioxide, formate, phosphate,
    and water from consideration as eliminating transport may be challenging. In addition, ATP synthase deletion was
    excluded from consideration since the deletion resulted in a high variability in ethanol production at the
    predicted optimal growth condition. Equivantly, the deletion of focA, focB, and atp operon were excluded from the
    OptORF simulations. The OptORF approach was applied to identify metabolic engineering strategies for
    overproduction of ethanol or higher alcohols (i.e. c j = 1 for j = desired alcohol secretion) in glucose minimal
    media. Maximum glucose uptake rate (GUR) and oxygen uptake rate (OUR) are specified in order to simulate
    anaerobic growth conditions (GUR = 18.5 mmol/gDW/hr, OUR = 0 mmol/gDW/hr) [21]. A minimal growth rate was set to
    0.1 hr-1 for all simulations.

    :return:
    """
    DIR = os.path.dirname(os.path.realpath(__file__))
    cbm_model_f = os.path.join(DIR, "../models/regulation/iJR904_srfba.xml")
    reg_model_f = os.path.join(DIR, '../models/regulation/imc1010.csv')

    _BIOMASS_ID = 'BiomassEcoli'
    _O2 = 'EX_o2_e'
    _GLC = 'EX_glc_DASH_D_e'
    _CO2 = 'EX_co2_e'
    _FE2 = "EX_fe2_e"
    _H = "EX_h_e"
    _H2O = "EX_h2o_e"
    _K = "EX_k_e"
    _NA1 = "EX_na1_e"
    _NH4 = "EX_nh4_e"
    _PI = "EX_pi_e"
    _SO4 = "EX_so4_e"
    _SUCC = "EX_succ_e"
    _ETOH = "EX_etoh_e"

    _GLY = "EX_gly_e"
    _DALA = "EX_ala_DASH_D_e"

    from mewpy.io import read_model, Engines, Reader

    metabolic_reader = Reader(Engines.MetabolicSBML, cbm_model_f)
    regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                               reg_model_f,
                               sep=';',
                               id_col=1,
                               rule_col=4,
                               aliases_cols=[0, 2, 3],
                               header=0)

    model = read_model(metabolic_reader, regulatory_reader)
    model.objective = {_BIOMASS_ID: 1}

    envcond = {_GLC: (-18.5, 100000.0),
               # _SUCC: (0.0, 100000.0),
               # _NH4: (-10.0, 100000.0),
               _O2: (0, 100000.0),
               # _CO2: (-15.0, 100000.0),
               # _PI: (-15.0, 100000.0),
               # _SO4: (-10.0, 100000.0),
               # _H: (-10.0, 100000.0),
               # _H2O: (-55.0, 100000.0)
               }

    for rxn, bds in envcond.items():
        model.get(rxn).bounds = bds

    model.get(_BIOMASS_ID).lower_bound = 0.1

    _PRODUCT_ID = _ETOH

    evaluator_1 = BPCY(_BIOMASS_ID, _PRODUCT_ID)
    evaluator_2 = WYIELD(_BIOMASS_ID, _PRODUCT_ID)

    problem = OptORFProblem(model, [evaluator_1, evaluator_2], candidate_max_size=6)

    ea = EA(problem, max_generations=10, mp=True)
    final_pop = ea.run()

    from mewpy.util.io import population_to_csv

    filename = "OPTORF_{}_KO_{}.csv".format(_PRODUCT_ID, "iJR904_srfba")
    population_to_csv(problem, final_pop, filename, simplify=False)


def optorf_ec():

    DIR = os.path.dirname(os.path.realpath(__file__))
    cbm_model_f = os.path.join(DIR, "../models/regulation/e_coli_core.xml")
    reg_model_f = os.path.join(DIR, '../models/regulation/e_coli_core_trn.csv')

    _BIOMASS_ID = 'Biomass_Ecoli_core'
    _O2 = 'EX_o2_e'
    _GLC = 'EX_glc__D_e'
    _FUM = 'EX_fum_e'
    _AC = 'EX_ac_e'
    _GLU = 'EX_glu__L_e'
    _LAC = 'EX_lac__D_e'
    _SUC = 'EX_succ_e'

    from mewpy.io import read_model, Engines, Reader

    metabolic_reader = Reader(Engines.MetabolicSBML, cbm_model_f)
    regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                               reg_model_f,
                               sep=',',
                               id_col=0,
                               rule_col=2,
                               aliases_cols=[1],
                               header=0)

    model = read_model(metabolic_reader, regulatory_reader)
    model.objective = {_BIOMASS_ID: 1}

    model.get(_GLC).bounds = (-10.0, 100000.0)
    model.get(_BIOMASS_ID).lower_bound = 0.1

    _PRODUCT_ID = _SUC

    evaluator_1 = BPCY(_BIOMASS_ID, _PRODUCT_ID)
    evaluator_2 = WYIELD(_BIOMASS_ID, _PRODUCT_ID)

    problem = OptORFProblem(model, [evaluator_1, evaluator_2], candidate_max_size=10)

    ea = EA(problem, max_generations=100, mp=True)
    final_pop = ea.run()

    from mewpy.util.io import population_to_csv

    filename = "OPTORF_{}_KO_{}.csv".format(_PRODUCT_ID, "ec")
    population_to_csv(problem, final_pop, filename, simplify=False)


if __name__ == '__main__':
    optorf_ec()
    # optorf_imc()
