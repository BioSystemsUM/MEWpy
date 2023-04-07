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
##############################################################################
Adapted from REFRAMED
##############################################################################
"""
from mewpy.solvers import solver_instance
from mewpy.solvers.solver import VarType
from mewpy.solvers.solution import Status
from mewpy.simulation import get_simulator
from mewpy.util.utilities import molecular_weight
from warnings import warn
from math import inf


def minimal_medium(model, exchange_reactions=None, direction=-1, min_mass_weight=False, min_growth=1,
                   max_uptake=100, max_compounds=None, n_solutions=1, validate=True, abstol=1e-6,
                   warnings=True, milp=True, use_pool=False, pool_gap=None, biomass_reaction=None,solver=None):
    """ Minimal medium calculator. Determines the minimum number of medium components for the organism to grow.
    Options: simply minimize the total number of components;
    minimize nutrients by molecular weight (as implemented by Zarecki et al, 2014)
    
    :param model: model
    :param exchange_reactions: list of exchange reactions (if not provided all model exchange reactions are used)
    :param (int) direction (int): direction of uptake reactions (negative or positive, default: -1)
    :param min_mass_weight (bool): minimize by molecular weight of compounds (default: False)
    :param min_growth (float): minimum growth rate (default: 1)
    :param max_uptake (float): maximum uptake rate (default: 100)
    :param max_compounds (int): limit maximum number of compounds (optional)
    :param n_solutions (int): enumerate multiple solutions (default: 1)
    :param validate (bool): validate solution using FBA (for debugging purposes, default: False)
    :param abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)
    :param warnings (bool): print warnings (default: False)
    :param milp (bool): minimize total number of compounds, otherwise use total flux (default: True)
    :param use_pool (bool): get multiple solutions from solution pool, otherwise enumerate (default: False)
    :param pool_gap (float): pool gap value when using solution pool (default: solver defined)
    
    :return:
        list: minimal set of exchange reactions
        Solution: solution(s) from solver
    """

    sim = get_simulator(model)

    def warn_wrapper(message):
        if warnings:
            warn(message)

    if exchange_reactions is None:
        exchange_reactions = sim.get_exchange_reactions()
    
    if solver is None:
        solver = solver_instance(sim)

    if not milp and max_compounds is not None:
        raise RuntimeError("max_compounds can only be used with MILP formulation")

    if not milp and n_solutions > 1:
        raise RuntimeError("n_solutions can only be used with MILP formulation")

    if milp:
        for r_id in exchange_reactions:
            solver.add_variable('y_' + r_id, 0, 1, vartype=VarType.BINARY, update=False)
    else:
        for r_id in exchange_reactions:
            solver.add_variable('f_' + r_id, 0, max_uptake, update=False)

    solver.update()

    if milp:
        for r_id in exchange_reactions:
            if direction < 0:
                solver.add_constraint('c_' + r_id, {r_id: 1, 'y_' + r_id: max_uptake}, '>', 0, update=False)
            else:
                solver.add_constraint('c_' + r_id, {r_id: 1, 'y_' + r_id: -max_uptake}, '<', 0, update=False)

        if max_compounds:
            lhs = {'y_' + r_id: 1 for r_id in exchange_reactions}
            solver.add_constraint('max_cmpds', lhs, '<', max_compounds, update=False)

    else:
        for r_id in exchange_reactions:
            if direction < 0:
                solver.add_constraint('c_' + r_id, {r_id: 1, 'f_' + r_id: 1}, '>', 0, update=False)
            else:
                solver.add_constraint('c_' + r_id, {r_id: 1, 'f_' + r_id: -1}, '<', 0, update=False)

    solver.update()

    valid_reactions = []

    if min_mass_weight:
        objective = {}

        multiple_compounds =[]
        no_compounds = []
        no_formula = []
        invalid_formulas = []

        for r_id in exchange_reactions:

            if direction < 0:
                compounds = sim.get_substrates(r_id)
            else:
                compounds = sim.get_products(r_id)

            if len(compounds) > 1:
                multiple_compounds.append(r_id)
                continue

            if len(compounds) == 0:
                no_compounds.append(r_id)
                continue
            
            metabolite = sim.get_metabolite(list(compounds.keys())[0])
            
            if not metabolite.formula:
                no_formula.append(metabolite.id)
                continue

            weight = molecular_weight(metabolite.formula)
            
            if weight is None:
                invalid_formulas.append(metabolite.id)
                continue

            if milp:
                objective['y_' + r_id] = weight
            else:
                objective['f_' + r_id] = weight

            valid_reactions.append(r_id)

        if multiple_compounds:
            warn_wrapper(f"Reactions ignored (multiple compounds): {multiple_compounds}")
        if no_compounds:
            warn_wrapper(f"Reactions ignored (no compounds): {no_compounds}")
        if multiple_compounds:
            warn_wrapper(f"Compounds ignored (no formula): {no_formula}")
        if invalid_formulas:
            warn_wrapper(f"Compounds ignored (invalid formula): {invalid_formulas}")

    else:
        if milp:
            objective = {'y_' + r_id: 1 for r_id in exchange_reactions}
        else:
            objective = {'f_' + r_id: 1 for r_id in exchange_reactions}

        valid_reactions = exchange_reactions

    result, ret_sols = None, None

    if direction < 0:
        constraints = {r_id: (-max_uptake if r_id in valid_reactions else 0, sim.get_reaction(r_id).ub)
                       for r_id in exchange_reactions}
    else:
        constraints = {r_id: (sim.get_reaction(r_id).lb, max_uptake if r_id in valid_reactions else 0)
                       for r_id in exchange_reactions}

    if biomass_reaction is None:
        biomass_reaction = list(sim.objective.keys())[0]
    
    constraints[biomass_reaction] = (min_growth, inf)

    if n_solutions == 1:

        solution = solver.solve(objective, minimize=True, constraints=constraints, get_values=exchange_reactions)
        
        if solution.status != Status.OPTIMAL:
            warn_wrapper('No solution found')
            result, ret_sols = None, solution
        else:
            medium = get_medium(solution, exchange_reactions, direction, abstol)

            if validate:
                validate_solution(model, medium, exchange_reactions, direction, min_growth, max_uptake)

            result, ret_sols = medium, solution

    elif use_pool:
        solutions = solver.solve(objective, minimize=True, constraints=constraints, get_values=exchange_reactions,
                                 pool_size=n_solutions, pool_gap=pool_gap)

        if solutions is None:
            result, ret_sols = [], []
        else:
            media = [get_medium(solution, exchange_reactions, direction, abstol) for solution in solutions]
            result, ret_sols = media, solutions
    else:
        media = []
        solutions = []

        for i in range(0, n_solutions):
            if i > 0:
                constr_id = f"iteration_{i}"
                previous_sol = {'y_' + r_id: 1 for r_id in medium}
                solver.add_constraint(constr_id, previous_sol, '<', len(previous_sol) - 1)

            solution = solver.solve(objective, minimize=True, constraints=constraints, get_values=exchange_reactions)

            if solution.status != Status.OPTIMAL:
                break

            medium = get_medium(solution, exchange_reactions, direction, abstol)
            media.append(medium)
            solutions.append(solution)

            result, ret_sols = media, solutions

    return result, ret_sols


def get_medium(solution, exchange, direction, abstol):
    return set(r_id for r_id in exchange
                 if (direction < 0 and solution.values[r_id] < -abstol
                     or direction > 0 and solution.values[r_id] > abstol))


def validate_solution(model, medium, exchange_reactions, direction, min_growth, max_uptake):
    sim = get_simulator(model)
    if direction == -1:
        constraints = {r_id: (-max_uptake, inf) if r_id in medium else (0, inf) for r_id in exchange_reactions}
    else:
        constraints = {r_id: (-inf, max_uptake) if r_id in medium else (-inf, 0) for r_id in exchange_reactions}

    sol = sim.simulate(constraints=constraints)

    if sol.objective_value < min_growth:
        warn('Solution appears to be invalid.')
