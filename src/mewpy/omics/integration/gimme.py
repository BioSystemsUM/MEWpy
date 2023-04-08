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
GIMME algorithm

Author: Vitor Pereira
Contributors: Paulo Carvalhais
##############################################################################
"""
from math import inf
from copy import deepcopy
from mewpy.solvers.solution import to_simulation_result
from mewpy.solvers import solver_instance
from mewpy.simulation import get_simulator
from mewpy.cobra.util import convert_to_irreversible
from .. import Preprocessing, ExpressionSet



def GIMME(model, expr, biomass=None, condition=0, cutoff=25, growth_frac=0.9,
          constraints=None, parsimonious=False, inline = False, **kwargs):
    """ Run a GIMME simulation [1]_.

    Arguments:
        model: a REFRAMED or COBRApy model or a MEWpy Simulator.
        expr (ExpressionSet): transcriptomics data.
        biomass: the biomass reaction identifier
        condition: the condition to use in the simulation\
            (default:0, the first condition is used if more than one.)
        cutoff (int): percentile cuttof (default: 25).
        growth_frac (float): minimum growth requirement (default: 0.9)
        constraints (dict): additional constraints
        parsimonious (bool): compute a parsimonious solution (default: False)
        inline (bool): returns a tissue specific model

    Returns:
        Solution: solution

    References
    ----------
    .. [1] Becker, S. and Palsson, B. O. (2008).
           Context-specific metabolic networks are consistent with experiments.
           PLoS Computational Biology, 4(5), e1000082.
           doi:10.1371/journal.pcbi.1000082
    """
    if not inline:
        sim = get_simulator(model)
    else:
        sim = get_simulator(deepcopy(model))
        
    if isinstance(expr, ExpressionSet):
        pp = Preprocessing(sim, expr)
        coeffs, _ = pp.percentile(condition, cutoff=cutoff)
    else:
        coeffs = expr
    
    solver = solver_instance(sim)

    if biomass is None:
        try:
            biomass = list(sim.objective.keys())[0]
        except Exception:
            raise ValueError(
                "A biomass reaction identifier is required or "
                "needs to be set as the model objective")

    wt_solution = sim.simulate(constraints=constraints)

    if not constraints:
        constraints = {}
    # add growth constraint
    constraints[biomass] = (growth_frac * wt_solution.fluxes[biomass], inf)

    # make model irreversible
    if not inline:
        for r_id in sim.reactions:
            lb, _ = sim.get_reaction_bounds(r_id)
            if lb < 0:
                pos, neg = r_id + '_p', r_id + '_n'
                solver.add_variable(pos, 0, inf, update=False)
                solver.add_variable(neg, 0, inf, update=False)
        solver.update()
    
        for r_id in sim.reactions:
            lb, _ = sim.get_reaction_bounds(r_id)
            if lb < 0:
                pos, neg = r_id + '_p', r_id + '_n'
                solver.add_constraint(
                    'c' + pos, {r_id: -1, pos: 1}, '>', 0, update=False)
                solver.add_constraint(
                    'c' + neg, {r_id: 1, neg: 1}, '>', 0, update=False)
        solver.update()

    else:
        convert_to_irreversible(sim,inline=True)


    objective = dict()
    for r_id, val in coeffs.items():
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
            pos, neg = r_id + '_p', r_id + '_n'
            objective[pos] = val
            objective[neg] = val
        else:
            objective[r_id] = val

    solution = solver.solve(objective, minimize=True, constraints=constraints)

    if parsimonious:
        pre_solution = solution

        solver.add_constraint('obj', objective, '=', pre_solution.fobj)
        objective = dict()

        for r_id in sim.reactions:
            lb, _ = sim.get_reaction_bounds(r_id)
            if lb < 0:
                pos, neg = r_id + '_p', r_id + '_n'
                objective[pos] = 1
                objective[neg] = 1
            else:
                objective[r_id] = 1

        solution = solver.solve(objective, minimize=True,
                                constraints=constraints)
        solver.remove_constraint('obj')
        solution.pre_solution = pre_solution

    for r_id in sim.reactions:
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
            pos, neg = r_id + '_p', r_id + '_n'
            del solution.values[pos]
            del solution.values[neg]

    res = to_simulation_result(model, biomass, constraints, sim, solution)
    if hasattr(solution,'pre_solution'):
        res.pre_solution = solution.pre_solution
    
    return res

