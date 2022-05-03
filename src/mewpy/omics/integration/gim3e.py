
from math import inf
import copy
from ...solvers.solver import VarType
from ...solvers import solver_instance
from ...simulation import get_simulator
from ...simulation.simulation import Simulator
from .. import Preprocessing


def GIM3E(model, expr, metabolites=[], biomass=None, constraints=None,
          opt_frac=0.9, penalty=1.01, tolerance=1E-9, inline=False):

    if not inline:
        model = copy.deepcopy(model)

    if isinstance(model, Simulator):
        sim = model
    else:
        sim = get_simulator(model)

    solver = solver_instance(sim)

    turnover_mets = []

    for met in sim.metabolites:
        compar = sim.get_metabolite(met)['compartment']
        turnover_met = 'TM_' + met
        sim.add_metabolite(turnover_met, turnover_met, compar)
        turnover_mets.append(turnover_met)

    for met in turnover_mets:
        sink = 'SK_' + met
        if met[3:] in metabolites:
            sim.add_reaction(sink, {met: 1}, lb=(penalty * tolerance), ub=1)
        else:
            sim.add_reaction(sink, {met: 1}, lb=0, ub=1)

    if biomass is None:
        try:
            biomass = list(sim.objective.keys())[0]
        except Exception:
            raise ValueError(
                "A biomass reaction identifier is required or "
                "needs to be set as the model objective")

    wt_solution = sim.simulate(constraints=constraints)

    constraints[biomass] = (opt_frac * wt_solution.fluxes[biomass], inf)

    pp = Preprocessing(sim, expr)

    express = pp.reactions_expression()

    _, exp_max = expr.minmax()

    for r_id in sim.reactions:
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
            pos, neg = r_id + '+', r_id + '-'
            solver.add_variable(pos, 0, 1, vartype=VarType.BINARY, update=False)
            solver.add_variable(neg, 0, 1, vartype=VarType.BINARY, update=False)
            solver.add_constraint(
                'c' + pos, {pos: 1, neg: 1}, '=', 1, update=False)
    solver.update()

    rxn_coef = {}

    for r_id in sim.reactions:
        penal_coeff = exp_max - express[r_id]
        rxn_coef[r_id] = penal_coeff
        solver.add_constraint('c' + r_id, {r_id: 1}, '<', penal_coeff, update=False)
    solver.update()

    # optimizar penalties

    sol_cons = solver.solve(rxn_coef, minimize=True, contraints=constraints)

    solver.add_constraint('penal_cons', {'penal_const': 1}, '<', (penalty * sol_cons), update=False)
    solver.update()

    # Final Optimization
    solution = solver.solve({biomass: 1}, minimize=False)

    return solution
