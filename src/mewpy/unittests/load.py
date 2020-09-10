from mewpy.simulation import get_simulator
import os


def load_reframed():

    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/optram/')
    DATA_FILE = os.path.join(PATH, "yeast_7.6-optram.xml")

    from reframed.io.sbml import load_cbmodel
    model = load_cbmodel(DATA_FILE, flavor='cobra')
    model.summary()
    simul = get_simulator(model)
    # int(simul.reactions)


def load_ec_gecko():

    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/autopacmen/')
    DATA_FILE = os.path.join(PATH, "iJO1366_2019_06_25_GECKO.xml")

    from mewpy.model.gecko import GeckoModel
    from reframed.io.sbml import load_cbmodel
    cbmodel = load_cbmodel(DATA_FILE)
    # ='R_ER_pool_TG_'
    model = GeckoModel(cbmodel, biomass_reaction_id='R_BIOMASS_Ec_iJO1366_WT_53p95M',
                       protein_reaction_id='R_PROTRS_TG_1', common_protein_pool_id='M_prot_pool')

    simul = get_simulator(model)
    for rxn in simul.reactions:
        if "R_PROT" in rxn:
            print(rxn)


def load_cobra():

    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/optram/')
    DATA_FILE = os.path.join(PATH, "yeast_7.6-optram.xml")

    from cobra.io import read_sbml_model
    model = read_sbml_model(DATA_FILE)
    simul = get_simulator(model)
    simul.summary()


def load_gecko():
    from mewpy.model.gecko import GeckoModel
    model = GeckoModel('single-pool')
    simul = get_simulator(model)
    simul.summary()


if __name__ == "__main__":
    load_cobra()
    load_gecko()
    load_reframed()
    load_ec_gecko()
