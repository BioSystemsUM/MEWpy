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
    # timeout is set in seconds
    sim = KineticSimulation(model, timeout=60,tSteps=[0,1])
    
    # simulation with defaults from the model
    res = sim.simulate()
    print('Kinetic simulation:')
    print('final rates:', res.fluxes)
    print('t=',res.t)
    print('y=',res.y)

    # simulate with user defined initial concentrations
    consc ={k:5*v for k,v in model.concentrations.items()}
    print("Modified metabolite concentrations:",consc)
    res = sim.simulate(initcon=consc)
    r1 = res.fluxes
    print(res.fluxes)
    
    # simulate with user defined folds of vMAXs
    # 0 folds are KO 
    factors = {'vPGI_rmaxPGI': 0, 
               'vPGM_rmaxPGM': 0,
               'vTKB_rmaxTKb': 0,
               'vTRPSYNTH_rmaxTrpSynth': 0,
               'vPK_rmaxPK': 0,
               'vPGDH_rmaxPGDH': 0,
               'vG1PAT_rmaxG1PAT': 0
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
    problem = KineticKOProblem(model,[f1], tSteps=[0, 1])
    ea = EA(problem,max_generations=10)
    ea.run()    

    # Multi objective
    problem = KineticOUProblem(model,[f1,f2],tSteps=[0, 1])
    ea = EA(problem,max_generations=10)
    ea.run()
    

if __name__ == '__main__':
    simulation()
    optimization()
