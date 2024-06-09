"""
##############################################################################
Regularized Flux Balance Analysis for communities
Author: Vitor Pereira
##############################################################################   
"""
from mewpy.solvers import solver_instance
from mewpy.simulation import get_simulator, SStatus, Simulator
from mewpy.solvers.solution import to_simulation_result
from warnings import warn

from . import CommunityModel

def regComFBA(cmodel, objective=None, maximize=True, constraints=None, obj_frac=0.99):
    """ Run a Regularized Flux Balance Analysis simulation:

    Arguments:
        model (CommunityModel): a constraint-based model
        objective (dict: objective coefficients (optional)
        minimize (bool): minimize objective function (False by default)
        constraints (dict): environmental or additional constraints (optional)
        
    Returns:
        Solution: solution
    """
    if isinstance(cmodel, CommunityModel):
        sim = cmodel.get_community_model()  
    elif isinstance(cmodel,Simulator): 
        sim = cmodel
    else:
        sim = get_simulator(cmodel)

    if not objective:
        objective = sim.objective
        if len(objective) == 0:
            warn('Model objective undefined.')

    solver = solver_instance(sim)

    if not constraints:
        constraints = {}

    if not objective:
        objective = sim.get_objective()

    pre_solution = sim.simulate(objective,maximize=maximize,constraints=constraints)
    if pre_solution.status != SStatus.OPTIMAL:
        return pre_solution

    solver.add_constraint('obj', objective, '>',
                              obj_frac * pre_solution.objective_value)

    solver.update()
   
    org_bio=list(sim.organisms_biomass.values())
    qobjective = {(rid,rid):1 for rid in org_bio}

    solution = solver.solve(quadratic=qobjective, minimize=True, constraints=constraints)

    result = to_simulation_result(sim, solution.fobj, constraints, sim, solution, regComFBA )
    
    return result
