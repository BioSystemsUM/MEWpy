from mewpy.simulation import get_simulator
import os


def load_reframed():

    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/')
    DATA_FILE = os.path.join(PATH, "iJO1366SL.xml")

    from reframed.io.sbml import load_cbmodel
    model = load_cbmodel(DATA_FILE, flavor='cobra')

    simul = get_simulator(model)
    print(simul.reactions)


def load_cobra():

    DIR = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(DIR, '../../../examples/models/')
    DATA_FILE = os.path.join(PATH, "iJO1366SL.xml")

    from cobra.io import read_sbml_model
    model = read_sbml_model(DATA_FILE)

    simul = get_simulator(model)
    print(simul.reactions)


load_cobra()
load_reframed()
