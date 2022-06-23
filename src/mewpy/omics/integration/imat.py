from math import inf

from mewpy.solvers.solver import VarType
from mewpy.solvers import solver_instance
from mewpy.simulation import get_simulator
from mewpy.simulation.simulation import Simulator
from mewpy.solvers.solution import to_simulation_result
from .. import Preprocessing, ExpressionSet


def iMAT(model, expr, constraints=None, cutoff=(25, 75),
        condition=0, epsilon=1):

    if isinstance(model, Simulator):
        sim = model
    else:
        sim = get_simulator(model)

    if isinstance(expr, ExpressionSet):
        pp = Preprocessing(sim, expr)
        coeffs, _ = pp.percentile(condition, cutoff=cutoff)
        low_coeffs, high_coeffs = coeffs
    else:
        low_coeffs, high_coeffs = expr

    solver = solver_instance(sim)

    if not constraints:
        constraints = {}

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

    objective = list()

    for r_id, val in high_coeffs.items():
        lb, ub = sim.get_reaction_bounds(r_id)
        pos_cons = (lb-epsilon)
        neg_cons = (ub+epsilon)
        pos, neg = 'y_' + r_id + '+', 'y_' + r_id + '-'
        objective.append(pos)
        solver.add_variable(pos, 0, 1, vartype=VarType.BINARY, update=True)
        solver.add_constraint(
            'c' + pos, {r_id: 1, pos: pos_cons}, '>', lb, update=False)
        objective.append(neg)
        solver.add_variable(neg, 0, 1, vartype=VarType.BINARY, update=True)
        solver.add_constraint(
            'c' + neg, {r_id: 1, neg: neg_cons}, '<', ub, update=False)

    solver.update()

    for r_id, val in low_coeffs.items():
        lb, ub = sim.get_reaction_bounds(r_id)
        x_var = 'x_' + r_id
        objective.append(x_var)
        solver.add_variable(x_var, 0, 1, vartype=VarType.BINARY, update=True)
        solver.add_constraint(
            'c' + x_var + '_pos', {r_id: 1, x_var: lb}, '>', lb, update=False)
        solver.add_constraint(
            'c' + x_var + '_neg', {r_id: 1, x_var: ub}, '<', ub, update=False)

    solver.update()

    object = {x: 1 for x in objective}

    solution = solver.solve(object, minimize=False, constraints=constraints)
    
    res = to_simulation_result(model, None, constraints, sim, solution)
    return res