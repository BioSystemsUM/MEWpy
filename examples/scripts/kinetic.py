from mewpy.io.sbml import load_ODEModel
from mewpy.simulation.kinetic import KineticSimulation
from scipy.integrate import odeint, solve_ivp
from collections import OrderedDict
from numpy import linspace, array, dot, isnan
import os


#from mewpy.solvers import set_default_ode_solver
#set_default_ode_solver('scikits')

DIR = os.path.dirname(os.path.realpath(__file__))
PATH = os.path.join(DIR, '../models/kinetic/')
DATA_FILE = os.path.join(PATH, "chassagnole2002.xml")
# DATA_FILE = os.path.join(PATH, "Jahan2016_chemostat_fixed.xml")


model = load_ODEModel(DATA_FILE)
factors = {'cpep':2, 'cglcex':3}
print(model.build_ode(factors=factors))
#print()
#print(model.concentrations)
#print(len(model.concentrations))
#sim = KineticSimulation(model, timeout=True)
#res = sim.simulate()
#print(res.fluxes)
