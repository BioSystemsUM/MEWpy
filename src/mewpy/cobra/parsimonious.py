from math import inf
from mewpy.solvers.solver import VarType
from mewpy.solvers import solver_instance
from mewpy.simulation import get_simulator, SStatus
from mewpy.simulation.simulation import Simulator


def epFBA(model, objective=None, reactions=None, constraints=None, obj_frac=None):
    """
    Modified versions of the Parsimonious Flux Balance Analysis allowing to minimize
    the sum of enzyme usage instead of the sum of reaction flux rates when a model includes 
    enzymatic constraints, such as GECKO and sMOMENT formulations.

    If the model defines protein constraints, and no set of reactions are defined,
    the objective will be to minimize enzyme usage, otherwise the objective is to
    minimize the sum of metabolic fluxes.

    :param model: a COBRAPY or REFRAMED model, or an instance of Simulator
    :param objective: The linear objective function as a dict of reaction and coefficient, defaults to None
    :type objective: dict, optional
    :param reactions: list of reactions whose sum of fluxes is to be minimized, defaults to None
    :type reactions: list, optional
    :param constraints: constraints to be imposed, defaults to None
    :type constraints: dict, optional
    :param obj_frac: fraction of the objective, defaults to None
    :type obj_frac: float, optional
    :return: a solver solution
    :rtype: mewpy.solver.Solution
    """

    if isinstance(model, Simulator):
        sim = model
    else:
        sim = get_simulator(model)

    if not objective:
        objective = sim.objective

    solver = solver_instance(sim)

    if not constraints:
        constraints = {}

    # update with simulation constraints if any
    constraints.update(sim.environmental_conditions)
    # constraints.update(sim._constraints)

    # make irreversible
    for r_id in sim.reactions:
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
            pos, neg = r_id + '+', r_id + '-'
            solver.add_variable(pos, 0, inf, update=False)
            solver.add_variable(neg, 0, inf, update=False)
    solver.update()

    for r_id in sim.reactions:
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
            pos, neg = r_id + '+', r_id + '-'
            solver.add_constraint(
                'c' + pos, {r_id: -1, pos: 1}, '>', 0, update=False)
            solver.add_constraint(
                'c' + neg, {r_id: 1, neg: 1}, '>', 0, update=False)

    solver.update()

    # add biomass constraint
    pre_solution = sim.simulate(constraints=constraints)

    if pre_solution.status != SStatus.OPTIMAL:
        return pre_solution

    if obj_frac is None:
        solver.add_constraint('obj', objective, '=', pre_solution.objective_value)
    else:
        solver.add_constraint('obj', objective, '>',
                              obj_frac * pre_solution.objective_value)

    solver.update()

    if not reactions:
        try:
            proteins = sim.proteins
            if proteins:
                reactions = [f"{sim.protein_prefix}{protein}" for protein in proteins]
        except Exception:
            reactions = sim.reactions

    sobjective = dict()

    if isinstance(reactions, dict):
        for k, v in reactions.items():
            sobjective[f"{sim.protein_prefix}{k}"] = v
    else:
        for r_id in reactions:
            lb, _ = sim.get_reaction_bounds(r_id)
            if lb < 0:
                pos, neg = r_id + '+', r_id + '-'
                sobjective[pos] = 1
                sobjective[neg] = 1
            else:
                sobjective[r_id] = 1

    solution = solver.solve(sobjective, minimize=True, constraints=constraints)

    return solution
