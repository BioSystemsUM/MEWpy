from mewpy.io.sbml import load_ODEModel
from mewpy.simulation.kinetic import KineticSimulation
from scipy.integrate import odeint, solve_ivp
from collections import OrderedDict
from numpy import linspace, array, dot, isnan
import os


DIR = os.path.dirname(os.path.realpath(__file__))
PATH = os.path.join(DIR, '../models/kinetic/')
DATA_FILE = os.path.join(PATH, "chassagnole2002.xml")


model = load_ODEModel(DATA_FILE)
sim = KineticSimulation(model, timeout=True)
res = sim.simulate()
print(res.fluxes)