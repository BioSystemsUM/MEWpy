
from mewpy.model.kinetic import load_ODEModel
from mewpy.simulation.kinetic import KineticSimulation
from mewpy.problems.kinetic import KineticOUProblem
from mewpy.optimization.evaluation import KineticTargetFlux
from mewpy.optimization import EA
from mewpy.utils.ode import SolverMethod
import os




def test1():

    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/kinetic/')
    """
    DATA_FILE = os.path.join(PATH, "chassagnole2002_mMs.xml")
    odemodel = load_ODEModel(DATA_FILE)
    """

    DATA_FILE = os.path.join(PATH, "Jahan2016_chemostat_fixed_july_mMs.xml")
    mapParamReacs = {"vE_6PGDH": ["v6PGDH_max"], "vE_Ack": ["vAck_max"], "vE_Ack_medium": ["vAck_max"],
                        "vE_Cya": ["vCya_max"], "vE_Eda": ["vEda_max"], "vE_Edd": ["vEdd_max"], "vE_Fum": ["Fum"],
                        "vE_G6PDH": ["vG6PDH_max"], "vE_MDH": ["MDH"], "vE_Pgi": ["vPgi_max"],
                        "vE_Pgl": ["vPgl_max"], "vE_Pta": ["vPta_max"], "vE_R5PI": ["vR5PI_max"], "vE_Ru5P": ["vRu5P_max"],
                        "vE_Tal": ["vTal_max"], "vE_TktA": ["vTktA_max"], "vE_TktB": ["vTktB_max"],
                        "vE_cAMPdegr": ["vcAMPdegr_max"], "vNonPTS": ["vNonPTS_max"], "vNonPTS_medium": ["vNonPTS_max"],
                        "vPTS4": ["vPTS4_max"], "vPTS4_medium": ["vPTS4_max"], "vE_AceKki": ["AceK"],
                        "vE_AceKph": ["AceK"], "vE_Acs": ["Acs"], "vE_Acs_medium": ["Acs"], "vE_CS": ["CS"],
                        "vE_Fba": ["Fba"], "vE_Fbp": ["Fbp"], "vE_GAPDH": ["GAPDH"], "vE_Glk": ["Glk"],
                        "vE_ICDH": ["ICDH"], "vE_Icl": ["Icl"], "vE_MS": ["MS"], "vE_Mez": ["Mez"], "vE_PDH": ["PDH"],
                        "vE_Pck": ["Pck"], "vE_Pfk": ["Pfk"], "vE_Ppc": ["Ppc"], "vE_Pps": ["Pps"], "vE_Pyk": ["Pyk"],
                        "vE_SDH": ["SDH"], "vE_aKGDH": ["aKGDH"]}

    model = load_ODEModel(DATA_FILE, map= mapParamReacs)
    
    print("reactions:",model.reactions)
    #factors = None
    #factors = {'CS': 32, 'v6PGDH_max': 0.0625, 'vEdd_max': 8, 'ICDH': 0.125, 'MS': 8, 'vPta_max': 0, 'SDH': 0.03125, 'vPTS4_max': 32}
    factors = {'ICDH': 0.03125, 'vAck_max': 0, 'SDH': 0.0625, 'vPTS4_max': 32}
    s = model.build_ode(factors)
    #print(s)
    res = KineticSimulation(model,method= SolverMethod.LSODE)
    res.simulate(factors = factors)
    #print("-------------------------------")



def test2():

    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/kinetic/')
    DATA_FILE = os.path.join(PATH, "Jahan2016_chemostat_fixed_july_mMs.xml")
    mapParamReacs = {"vE_6PGDH": ["v6PGDH_max"], "vE_Ack": ["vAck_max"], "vE_Ack_medium": ["vAck_max"],
                        "vE_Cya": ["vCya_max"], "vE_Eda": ["vEda_max"], "vE_Edd": ["vEdd_max"], "vE_Fum": ["Fum"],
                        "vE_G6PDH": ["vG6PDH_max"], "vE_MDH": ["MDH"], "vE_Pgi": ["vPgi_max"],
                        "vE_Pgl": ["vPgl_max"], "vE_Pta": ["vPta_max"], "vE_R5PI": ["vR5PI_max"], "vE_Ru5P": ["vRu5P_max"],
                        "vE_Tal": ["vTal_max"], "vE_TktA": ["vTktA_max"], "vE_TktB": ["vTktB_max"],
                        "vE_cAMPdegr": ["vcAMPdegr_max"], "vNonPTS": ["vNonPTS_max"], "vNonPTS_medium": ["vNonPTS_max"],
                        "vPTS4": ["vPTS4_max"], "vPTS4_medium": ["vPTS4_max"], "vE_AceKki": ["AceK"],
                        "vE_AceKph": ["AceK"], "vE_Acs": ["Acs"], "vE_Acs_medium": ["Acs"], "vE_CS": ["CS"],
                        "vE_Fba": ["Fba"], "vE_Fbp": ["Fbp"], "vE_GAPDH": ["GAPDH"], "vE_Glk": ["Glk"],
                        "vE_ICDH": ["ICDH"], "vE_Icl": ["Icl"], "vE_MS": ["MS"], "vE_Mez": ["Mez"], "vE_PDH": ["PDH"],
                        "vE_Pck": ["Pck"], "vE_Pfk": ["Pfk"], "vE_Ppc": ["Ppc"], "vE_Pps": ["Pps"], "vE_Pyk": ["Pyk"],
                        "vE_SDH": ["SDH"], "vE_aKGDH": ["aKGDH"]}

    model = load_ODEModel(DATA_FILE, map= mapParamReacs)
    f1 = KineticTargetFlux("vD_SUC")
    problem = KineticOUProblem(model,[f1])
    ea = EA(problem, max_generations=10, mp=True)
    final_pop = ea.run()

if __name__ == "__main__":
    test2()