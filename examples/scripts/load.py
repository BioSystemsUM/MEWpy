'''
Simple example/tests for loading models
'''

import os
from mewpy.simulation import get_simulator


def load_reframed():
    """
    Loads a model with REFRAMED
    """
    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../models/optram/')
    DATA_FILE = os.path.join(PATH, "yeast_7.6-optram.xml")

    from reframed.io.sbml import load_cbmodel
    model = load_cbmodel(DATA_FILE, flavor='cobra')
    model.summary()
    simul = get_simulator(model)
    simul.summary()


def load_ec_gecko():
    """ Loads a GECKO like model from AUTOPACMEN
    """
    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../models/autopacmen/')
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
    """Load a model using COBRApy
    """
    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../models/optram/')
    DATA_FILE = os.path.join(PATH, "yeast_7.6-optram.xml")

    from cobra.io import read_sbml_model
    model = read_sbml_model(DATA_FILE)
    simul = get_simulator(model)
    simul.summary()


def load_gecko():
    """Loads yeast GECKO model using REFRAMED
    """
    from mewpy.model.gecko import GeckoModel
    model = GeckoModel('single-pool')
    simul = get_simulator(model)
    simul.summary()


def load_germ_model():
    """
    Loads a GERM model
    """
    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../models/optram/')
    DATA_FILE = os.path.join(PATH, "yeast_7.6-optram.xml")

    from mewpy.io import read_sbml
    model = read_sbml(DATA_FILE, metabolic=True, regulatory=False)
    simul = get_simulator(model)
    simul.summary()


if __name__ == "__main__":
    load_cobra()
    load_gecko()
    load_reframed()
    load_ec_gecko()
