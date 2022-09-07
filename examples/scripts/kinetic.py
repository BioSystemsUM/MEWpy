from mewpy.io.sbml import load_ODEModel
from mewpy.simulation.kinetic import KineticSimulation
from mewpy.solvers import set_default_ode_solver
from mewpy.problems import KineticKOProblem, KineticOUProblem
from mewpy.optimization.evaluation import TargetFlux, ModificationType
from mewpy.optimization import EA
import os

set_default_ode_solver('scipy')

DIR = os.path.dirname(os.path.realpath(__file__))
PATH = os.path.join(DIR, '../models/kinetic/')
DATA_FILE = os.path.join(PATH, "chassagnole2002.xml")
model = load_ODEModel(DATA_FILE)

def simulation():
    print("Parameters:", model.get_parameters(exclude_compartments=True))
    print("Metabolite concentrations:",model.concentrations)
    print()
    sim = KineticSimulation(model, timeout=True,tSteps=[0,1])
    res = sim.simulate()
    print('Kinetic simulation:')
    print('final rates:', res.fluxes)
    print('t=',res.t)
    print('y=',res.y)
    consc ={k:5*v for k,v in model.concentrations.items()}
    print("Modified metabolite concentrations:",consc)
    res = sim.simulate(initcon=consc)
    r1 = res.fluxes
    print(res.fluxes)
    
    factors = {'vPGI_rmaxPGI': 0, 'vPGM_rmaxPGM': 0, 'vTKB_rmaxTKb': 0, 'vTRPSYNTH_rmaxTrpSynth': 0, 'vPK_rmaxPK': 0, 'vPGDH_rmaxPGDH': 0, 'vG1PAT_rmaxG1PAT': 0}
    print("Modidief kinetic paramenter by a factor", factors)
    res = sim.simulate(factors=factors)
    r2 = res.fluxes
    import pandas as pd
    df = pd.DataFrame([r1,r2])
    print(df)


def optimization():
    print('\n\nOptimization')
    f1 = TargetFlux('vsersynth')
    f2 = ModificationType()
    problem = KineticKOProblem(model,[f1],tSteps=[0,1])
    ea = EA(problem,max_generations=10)
    ea.run()    

    problem = KineticOUProblem(model,[f1,f2],tSteps=[0,1])
    ea = EA(problem,max_generations=10)
    ea.run()
    

if __name__ == '__main__':
    simulation()
    optimization()
