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
COBRA utility module

Authors: Vitor Pereira
##############################################################################
"""
from mewpy.simulation import get_simulator
from mewpy.simulation.simulation import Simulator
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from cobra import Model
    from reframed.core.cbmodel import CBModel


def convert_to_irreversible(model: Union[Simulator, "Model", "CBModel"]):
    """Split reversible reactions into two irreversible reactions
    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.
    
    :param model: A COBRApy or REFRAMED Model or an instance of 
        mewpy.simulation.simulation.Simulator

    TODO: Add anotations
    """
    if isinstance(model, Simulator):
        sim = model
    else:
        sim = get_simulator(model)

    objective = sim.objective.copy()

    for r_id in sim.reactions:
        lb, ub = sim.get_reaction_bounds(r_id)
        if lb < 0 and ub > 0:
            rxn = sim.get_reaction(r_id)
            rev_rxn_id = r_id+"_REV"

            rev_rxn = dict()
            rev_rxn['name'] = rxn.name + " reverse"
            rev_rxn['lb'] = 0
            rev_rxn['ub'] = -rxn.lb
            rev_rxn['gpr'] = rxn.gpr
            sth = {k: v * -1 for k, v in rxn.stoichiometry.items()}
            rev_rxn['stoichiometry'] = sth
            rev_rxn['reversible'] = False

            sim.add_reaction(rev_rxn_id, **rev_rxn)
            sim.set_reaction_bounds(r_id, 0, rxn.ub, False)

            if r_id in objective:
                objective[rev_rxn_id] = -objective[r_id]

    sim.objective = objective
    return sim


def split_isozymes(model):

    if isinstance(model, Simulator):
        sim = model
    else:
        sim = get_simulator(model)
        
    mapping = dict()

    for r_id in sim.reactions:
        rx_d = sim.get_reaction(r_id)
        gpr = rx_d.gpr

        if gpr is not None:
            mapping[r_id] = []
            for i, protein in enumerate(gpr.proteins):
                r_id_new = '{}_pc{}'.format(reaction.id, i+1)
                mapping[r_id].append(r_id_new)
                reaction_new = Reaction(r_id_new, reaction.name, reaction.reversible,
                                        reaction.stoichiometry, reaction.regulators)
                model.add_reaction(reaction_new)
                model.set_flux_bounds(r_id_new, *model.bounds[r_id])
                model.set_reaction_objective(r_id_new, model.objective[r_id])
                gpr_new = GPRAssociation()
                gpr_new.proteins.append(protein)
                model.set_gpr_association(r_id_new, gpr_new)
            model.remove_reaction(r_id)
    return mapping
