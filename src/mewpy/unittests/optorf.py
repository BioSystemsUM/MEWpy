from mewpy.optimization.evaluation import BPCY, WYIELD
from mewpy.optimization import EA
from mewpy.simulation import SimulationMethod
from mewpy.regulation import RFBAModel
from mewpy.regulation.optorf import OptOrfProblem
import time
import os


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

    from reframed.io.sbml import load_cbmodel
    from mewpy.simulation.reframed import Simulation

    DIR = os.path.dirname(os.path.realpath(__file__))
    cbm_model_f = os.path.join(DIR,"../../../examples/models/regulation/iJR904_srfba.xml")
    reg_model_f = os.path.join(DIR,'../../../examples/models/regulation/imc1010_v6.csv')
    aliases_f = os.path.join(DIR,'../../../examples/models/regulation/imc1010_rfba_aliases.csv')
    # env_cond_f = "../../../examples/models/regulation/imc1010_env_cond.xlsx"

    _BIOMASS_ID = 'R_BiomassEcoli'
    _O2 = 'R_EX_o2_e'
    _GLC = 'R_EX_glc_DASH_D_e'
    _CO2 = 'R_EX_co2_e'
    _FE2 = "R_EX_fe2_e"
    _H = "R_EX_h_e"
    _H2O = "R_EX_h2o_e"
    _K = "R_EX_k_e"
    _NA1 = "R_EX_na1_e"
    _NH4 = "R_EX_nh4_e"
    _PI = "R_EX_pi_e"
    _SO4 = "R_EX_so4_e"
    _SUCC = "R_EX_succ_e"
    _ETOH = "R_EX_etoh_e"

    _GLY = "R_EX_gly_e"
    _DALA = "R_EX_ala_DASH_D_e"

    model = load_cbmodel(cbm_model_f)
    model.set_objective({_BIOMASS_ID: 1})

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

    simulation = Simulation(model, envcond=envcond)

    rfba = RFBAModel.from_tabular_format(reg_model_f, model, simulation,
                                         sep=';', id_col=1, rule_col=4, aliases_cols=[0, 2, 3], header=0)
    rfba.update_aliases_from_tabular_format_file(aliases_f, sep=';', id_col=0, aliases_cols=[1])

    initial_state = {var: 1 for var in rfba.targets}
    initial_state.update({_BIOMASS_ID: 0.1})
    rfba.initial_state = initial_state

    _PRODUCT_ID = _ETOH

    evaluator_1 = BPCY(_BIOMASS_ID, _PRODUCT_ID, method=SimulationMethod.pFBA)
    evaluator_2 = WYIELD(_BIOMASS_ID, _PRODUCT_ID)

    problem = OptOrfProblem(model, [evaluator_1, evaluator_2], rfba, candidate_max_size=6)
    
    ea = EA(problem, max_generations=100, mp=True)
    final_pop = ea.run()

    import mewpy.utils.utilities as utl

    filename = "OPTORF{}_KO_{}.csv".format(_PRODUCT_ID, "iJR904_srfba" )
    utl.population_to_csv(problem, final_pop, filename, simplify=False)


def optorf_ec():
    
    import cobra.test
    from mewpy.simulation.cobra import Simulation

    DIR = os.path.dirname(os.path.realpath(__file__))
    reg_model_f = os.path.join(DIR,'../../../examples/models/regulation/core_TRN_v2.csv')
    aliases_f = os.path.join(DIR,'../../../examples/models/regulation/core_TRN_rfba_aliases.csv')

    _BIOMASS_ID = 'Biomass_Ecoli_core'
    _O2 = 'EX_o2_e'
    _GLC = 'EX_glc__D_e'
    _FUM = 'EX_fum_e'
    _AC = 'EX_ac_e'
    _GLU = 'EX_glu__L_e'
    _LAC = 'EX_lac__D_e'
    _SUC = 'EX_succ_e'
    
    model = cobra.test.create_test_model("textbook")
    model.objective = _BIOMASS_ID

    envcond = {_GLC: (-10.0, 100000.0)}
    
    simulation = Simulation(model, envcond=envcond)
    
    

    rfba = RFBAModel.from_tabular_format(reg_model_f, model, simulation,
                                         sep=',', id_col=1, rule_col=2, aliases_cols=[0], header=0)
    rfba.update_aliases_from_tabular_format_file(aliases_f, id_col=1, aliases_cols=[0])
    

    initial_state = {var: 1 for var in rfba.targets}
    initial_state.update({_BIOMASS_ID: 0.1})
    rfba.initial_state = initial_state

    _PRODUCT_ID = _SUC

    evaluator_1 = BPCY(_BIOMASS_ID, _PRODUCT_ID, method=SimulationMethod.pFBA)
    evaluator_2 = WYIELD(_BIOMASS_ID, _PRODUCT_ID)


    problem = OptOrfProblem(model, [evaluator_1, evaluator_2], rfba, candidate_max_size=6)
    
    print(len(problem.target_list))


    ea = EA(problem, max_generations=100, mp=True)
    final_pop = ea.run()

    import mewpy.utils.utilities as utl

    filename = "OPTORF{}_KO_{}.csv".format(_PRODUCT_ID, "ec" )
    utl.population_to_csv(problem, final_pop, filename, simplify=False)




if __name__ == '__main__':
    optorf_ec()
