# Regularized Flux Balance Analysis for communities
from mewpy.solvers import solver_instance
from mewpy.simulation import get_simulator, SStatus, SimulationResult
from mewpy.com.com import CommunityModel
from mewpy.solvers.solution import to_simulation_result
from warnings import warn

def regComFBA(model:CommunityModel, objective=None, minimize=False, constraints=None, alpha=0.9):
    """ Run a Regularized Flux Balance Analysis simulation:

    Arguments:
        model (CommunityModel): a constraint-based model
        objective (dict: objective coefficients (optional)
        minimize (bool): minimize objective function (False by default)
        constraints (dict): environmental or additional constraints (optional)
        
    Returns:
        Solution: solution
    """
    sim = get_simulator(model)

    if not objective:
        objective = sim.objective
        if len(objective) == 0:
            warn('Model objective undefined.')

    solver = solver_instance(sim)

    if not constraints:
        constraints = {}

    if not objective:
        objective = model.get_objective()

    pre_solution = sim.simulate(objective,minimize=minimize,constraints=constraints)
    if pre_solution.status != SStatus.OPTIMAL:
        return pre_solution

    solver.add_constraint('obj', objective, '>',
                              alpha * pre_solution.objective_value)

    solver.update()
   
    org_bio=list(model.get_organisms_biomass().values())
    qobjective = {(rid,rid):1 for rid in org_bio}

    solution = solver.solve(quadratic=qobjective, minimize=True, constraints=constraints)

    

    result = to_simulation_result(model, solution.fobj, constraints, sim, solution)
    
    return result
