# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
Author: Vitor Pereira
"""
from mewpy.io.sbml import load_ODEModel
from mewpy.simulation.kinetic import KineticSimulation
from mewpy.solvers import set_default_ode_solver
from mewpy.problems import KineticKOProblem, KineticOUProblem
from mewpy.optimization.evaluation import TargetFlux, ModificationType
from mewpy.optimization import EA
import os

# to verify the solvers available on the environment call 
# >>mewpy.info()
set_default_ode_solver('scipy')

# Load the kinetic model
DIR = os.path.dirname(os.path.realpath(__file__))
PATH = os.path.join(DIR, '../models/kinetic/')
DATA_FILE = os.path.join(PATH, "chassagnole2002.xml")
model = load_ODEModel(DATA_FILE)

def simulation():
    "Simulation example"
    # model parameters
    print("Parameters:", model.get_parameters(exclude_compartments=True))
    # initial concentrstions
    print("Metabolite concentrations:",model.concentrations)
    print()
    # Builds a simulator for kinetic models
    # - timeout is set in seconds
    # - time points may be defined as a span, an array with the starting and end point 
    #   (the number of time steps is defined in SolverConfigurations, see bellow)
    #   or as time points [1,2,3,4,5,6,7,8,10]
    #   from mewpy.solvers.ode import SolverConfigurations
    #   SolverConfigurations.N_STEPS = 10
    sim = KineticSimulation(model, timeout=60, t_points=[0,1])
    
    # simulation with defaults from the model
    res = sim.simulate()
    print('Kinetic simulation:')
    print('final rates:', res.fluxes)
    print('t=',res.t)
    print('y=',res.y)
    print()
    # Changing the time points
    print('Kinetic simulation with time points:')
    res = sim.simulate(t_points=[1,2,3,4,5,6,7,8,9,10])
    print('t=',res.t)
    print('y=',res.y)
    print()
    # simulate with user defined initial concentrations
    consc ={k:5*v for k,v in model.concentrations.items()}
    print("Modified metabolite concentrations:",consc)
    r1 = res.fluxes
    # simulate with user defined folds of vMAXs
    # 0 folds are KO 
    factors = {'rmaxPGI': 0, 
               'rmaxPGM': 0,
               'rmaxTKb': 0,
               'rmaxTrpSynth': 0,
               'rmaxPK': 0,
               'rmaxPGDH': 0,
               'rmaxG1PAT': 0
               }
    print("Modified kinetic paramenter by a factor", factors)
    res = sim.simulate(factors=factors)
    r2 = res.fluxes
   
    # pretty printing of results
    import pandas as pd
    df = pd.DataFrame([r1, r2])
    print(df)


def optimization():
    print('\n\nOptimization')

    # Optimization objective functions
    f1 = TargetFlux('vsersynth')
    f2 = ModificationType()
    
    # optimization problems:

    # Single objective
    problem = KineticKOProblem(model,[f1], t_points=[0, 1e9])
    ea = EA(problem,max_generations=10)
    ea.run()    

    # Multi objective
    problem = KineticOUProblem(model,[f1,f2],t_points=[0, 1e9])
    ea = EA(problem,max_generations=10)
    ea.run()
    

if __name__ == '__main__':
    simulation()
    #optimization()
