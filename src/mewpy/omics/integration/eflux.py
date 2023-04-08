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
E-Flux algorithm

Author: Vitor Pereira
Contributors: Paulo Carvalhais
##############################################################################
"""
from mewpy.simulation import get_simulator
from .. import Preprocessing, ExpressionSet


def eFlux(model, expr, condition=0, scale_rxn=None, scale_value=1,
          constraints=None, parsimonious=False, max_exp = None, **kwargs):
    """ Run an E-Flux simulation (Colijn et al, 2009).

    
    :param model: a REFRAMED or COBRApy model or a MEWpy Simulator.
    :param expr (ExpressionSet): transcriptomics data.
    :param condition: the condition to use in the simulation\
                (default:0, the first condition is used if more than one.)
    :param scale_rxn (str): reaction to scale flux vector (optional)
    :param scale_value (float): scaling factor (mandatory if scale_rxn\
            is specified)
    :param constraints (dict): additional constraints (optional)
    :param parsimonious (bool): compute a parsimonious solution (default: False)

    :return: Solution: solution
    """

    sim = get_simulator(model)

    if isinstance(expr, ExpressionSet):
        pp = Preprocessing(sim, expr)
        rxn_exp = pp.reactions_expression(condition)
    else:
        rxn_exp = expr

    if max_exp is None:
        max_exp = max(rxn_exp.values())

    bounds = {}

    for r_id in sim.reactions:
        val = rxn_exp[r_id] / max_exp if r_id in rxn_exp else 1
        lb, ub = sim.get_reaction_bounds(r_id)
        lb2 = -val if lb < 0 else 0
        ub2 = val if ub > 0 else 0
        bounds[r_id] = (lb2, ub2)

    if constraints:
        for r_id, x in constraints.items():
            lb, ub = x if isinstance(x, tuple) else (x, x)
            lb2 = -1 if lb < 0 else 0
            ub2 = 1 if ub > 0 else 0
            bounds[r_id] = (lb2, ub2)

    if parsimonious:
        sol = sim.simulate(constraints=bounds, method='pFBA')
    else:
        sol = sim.simulate(constraints=bounds)

    if scale_rxn is not None:

        if sol.fluxes[scale_rxn] != 0:
            k = abs(scale_value / sol.fluxes[scale_rxn])
        else:
            k = 0

        for r_id, val in sol.fluxes.items():
            sol.fluxes[r_id] = val * k

    return sol
