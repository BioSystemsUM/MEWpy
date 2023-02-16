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
from mewpy.util.parsing import isozymes
from copy import deepcopy, copy
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from cobra import Model
    from reframed.core.cbmodel import CBModel


def convert_to_irreversible(model: Union[Simulator, "Model", "CBModel"], inline=False):
    """Split reversible reactions into two irreversible reactions
    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.
    
    :param model: A COBRApy or REFRAMED Model or an instance of 
        mewpy.simulation.simulation.Simulator
    """
    if isinstance(model, Simulator):
        if inline:
            sim = model
        else:
            sim = deepcopy(model)
    else:
        if inline:
            sim = get_simulator(model)
        else:    
            sim = get_simulator(deepcopy(model))

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
            rev_rxn['annotations'] = copy(rxn.annotations)

            sim.add_reaction(rev_rxn_id, **rev_rxn)
            sim.set_reaction_bounds(r_id, 0, rxn.ub, False)

            if r_id in objective:
                objective[rev_rxn_id] = -objective[r_id]

    sim.objective = objective
    return sim



def split_isozymes(model: Union[Simulator, "Model", "CBModel"], inline=False):
    """Splits reactions with isozymes into separated reactions

    :param model: A COBRApy or REFRAMED Model or an instance of 
        mewpy.simulation.simulation.Simulator
    :param (boolean) inline: apply the modifications to the same of generate a new model. Default generates a new model.
    :return: 
    :rtype: _type_
    """

    if isinstance(model, Simulator):
        if inline:
            sim = model
        else:
            sim = deepcopy(model)
    else:
        if inline:
            sim = get_simulator(model)
        else:
            sim = get_simulator(deepcopy(model))

    objective = sim.objective
    mapping = dict()
    newobjective = {}

    for r_id in sim.reactions:
        rxn = sim.get_reaction(r_id)
        gpr = rxn.gpr

        if gpr is not None and len(gpr.strip()) > 0:
            proteins = isozymes(gpr)
            mapping[r_id] = []
            for i, protein in enumerate(proteins):
                r_id_new = '{}_No{}'.format(r_id, i+1)
                mapping[r_id].append(r_id_new)

                rxn_new = dict()
                rxn_new['name'] = '{} No{}'.format(rxn.name, i+1)
                rxn_new['lb'] = rxn.lb
                rxn_new['ub'] = rxn.ub
                rxn_new['gpr'] = protein
                rxn_new['stoichiometry'] = rxn.stoichiometry.copy()
                rxn_new['annotations'] = copy(rxn.annotations)
                
                sim.add_reaction(r_id_new, **rxn_new)
            sim.remove_reaction(r_id)

    # set the objective
    for k, v in objective.items():
        if k in mapping.keys():
            for r in mapping[k]:
                newobjective[r] = v
        else:
            newobjective[k] = v

    sim.objective = newobjective

    return sim, mapping


def genes_to_species(model: Union[Simulator, "Model", "CBModel"], inline=False):

    if isinstance(model, Simulator):
        if inline:
            sim = model
        else:
            sim = deepcopy(model)
    else:
        if inline:
            sim = get_simulator(model)
        else:
            sim = get_simulator(deepcopy(model))


    for gene in sim.genes:
        sim.add_metabolite(gene,name=gene,compartment="")