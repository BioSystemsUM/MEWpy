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
from mewpy.simulation import get_simulator, Simulator
from mewpy.util.parsing import isozymes, build_tree, Boolean
from mewpy.util.constants import ModelConstants
from copy import deepcopy, copy
from math import inf
from tqdm import tqdm
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from cobra import Model
    from reframed.core.cbmodel import CBModel


def convert_gpr_to_dnf(model) -> None:
    """
    Convert all existing GPR associations to DNF.
    """
    sim = get_simulator(model)
    for rxn_id in tqdm(sim.reactions):
        rxn = sim.get_reaction(rxn_id)
        if not rxn.gpr:
            continue
        tree = build_tree(rxn.gpr, Boolean)
        gpr = tree.to_infix()
        # TODO: update the gpr
        
        return gpr 
        
        
def convert_to_irreversible(model: Union[Simulator, "Model", "CBModel"], inline: bool = False):
    """Split reversible reactions into two irreversible reactions
    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    :param model: A COBRApy or REFRAMED Model or an instance of 
        mewpy.simulation.simulation.Simulator
    :return: a irreversible model simulator, a reverse mapping.
    :rtype:(Simulator,dict)
    """
    
    sim = get_simulator(deepcopy(model))

    objective = sim.objective.copy()
    irrev_map=dict()
    for r_id in tqdm(sim.reactions, "Converting to irreversible"):
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
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
            
            irrev_map[r_id] = rev_rxn_id
            
            if r_id in objective:
                objective[rev_rxn_id] = -objective[r_id]

    sim.objective = objective
    return sim, irrev_map


def split_isozymes(model: Union[Simulator, "Model", "CBModel"], inline: bool = False):
    """Splits reactions with isozymes into separated reactions

    :param model: A COBRApy or REFRAMED Model or an instance of 
        mewpy.simulation.simulation.Simulator
    :param (boolean) inline: apply the modifications to the same of generate a new model. Default generates a new model.
    :return: a simulator and a mapping from original to splitted reactions
    :rtype: (Simulator, dict)
    """

    sim = get_simulator(deepcopy(model))

    objective = sim.objective
    mapping = dict()
    newobjective = {}

    for r_id in tqdm(sim.reactions, "Splitting isozymes"):
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


def __enzime_constraints(model: Union[Simulator, "Model", "CBModel"],
                         prot_mw=None,
                         enz_kcats=None,
                         c_compartment: str = 'c',
                         inline: bool = False):
    """Auxiliary method to add enzyme constraints to a model

    :param model: A model or simulator
    :type model: A COBRApy or REFRAMED Model or an instance of 
        mewpy.simulation.simulation.Simulator
    :param data: Protein MW and Kcats 
    :type data: None
    :param c_compartment: The compartment where gene/proteins pseudo species are to be added.
        Defaults to 'c'
    :type c_compartment: str, optional
    :param (boolean) inline: apply the modifications to the same of generate a new model.
        Default generates a new model.
    :type inline: bool, optional
    :return: a new enzyme constrained model
    :rtype: Simulator
    """

    if inline:
        sim = get_simulator(model)
    else:
        sim = deepcopy(get_simulator(model))
    
    objective = sim.objective

    if prot_mw is None:
        prot_mw = dict()
        for gene in sim.genes:
            prot_mw[gene] = {'protein': gene[len(sim._g_prefix):], 'mw': 1}

    if enz_kcats is None:
        enz_kcats = dict()
        for gene in sim.genes:
            enz_kcats[gene] = dict()
            rxns = sim.get_gene(gene).reactions
            for rxn in rxns:
                enz_kcats[gene][rxn] = {'protein': gene[len(sim._g_prefix):], 'kcat': 1}


    # Add protein pool and species
    common_protein_pool_id = sim._m_prefix+'prot_pool_c'
    pool_reaction = sim._r_prefix+'prot_pool_exchange'

    sim.add_metabolite(common_protein_pool_id,
                       name='prot_pool [cytoplasm]',
                       compartment=c_compartment)

    sim.add_reaction(pool_reaction,
                     name='protein pool exchange',
                     stoichiometry={common_protein_pool_id: 1},
                     lb=0,
                     ub=inf,
                     reversible=False,
                     reaction_type='EX'
                     )

    # Add gene/protein species and draw protein pseudo-reactions
    # MW in kDa, [kDa = g/mmol]
    gene_meta = dict()
    for gene in tqdm(sim.genes, "Adding gene species"):
        mw = prot_mw[gene]
        m_prot_id = f"prot_{mw['protein']}_{c_compartment}"
        m_name = f"prot_{mw['protein']} {c_compartment}"
        sim.add_metabolite(m_prot_id,
                           name=m_name,
                           compartment=c_compartment)

        gene_meta[gene] = m_prot_id

        r_prot_id = f"draw_prot_{mw['protein']}"
        sim.add_reaction(r_prot_id,
                         name=r_prot_id,
                         stoichiometry={common_protein_pool_id: -1*mw['mw'],
                                        m_prot_id: 1},
                         lb=0,
                         ub=inf,
                         reversible=False,
                         )

    # Add enzymes to reactions stoichiometry.
    # 1/Kcats in per hour. Considering kcats in per second.
    for rxn_id in tqdm(sim.reactions, "Adding proteins usage to reactions"):
        rxn = sim.get_reaction(rxn_id)
        if rxn.gpr:
            s = rxn.stoichiometry
            genes = build_tree(rxn.gpr, Boolean).get_operands()
            for g in genes:
                # TODO: mapping of (gene, reaction ec) to kcat
                try:
                    s[gene_meta[g]] = -1/(enz_kcats[g][rxn_id]['kcat']*3600)
                except Exception:
                    s[gene_meta[g]] = -1/(ModelConstants.DEFAULT_KCAT*3600)
            sim.update_stoichiometry(rxn_id, s)
    sim.objective = objective
    return sim


def add_enzyme_constraints(model: Union[Simulator, "Model", "CBModel"],
                           prot_mw=None,
                           enz_kcats=None,
                           c_compartment: str='c',
                           inline: bool=False):
    """Adds enzyme constraints to a model.

    :param model: A model or simulator
    :type model: A COBRApy or REFRAMED Model or an instance of 
        mewpy.simulation.simulation.Simulator
    :param data: Protein MW and Kcats 
    :type data: None
    :param c_compartment: The compartment where gene/proteins pseudo species are to be added.
        Defaults to 'c'
    :type c_compartment: str, optional
    :param (boolean) inline: apply the modifications to the same of generate a new model.
        Default generates a new model.
    :type inline: bool, optional
    :return: a new enzyme constrained model
    :rtype: Simulator
    """
    sim, _ = convert_to_irreversible(model, inline)
    sim, _ = split_isozymes(sim, True)
    sim = __enzime_constraints(sim,
                               prot_mw=prot_mw,
                               enz_kcats=enz_kcats,
                               c_compartment=c_compartment,
                               inline=True)
    return sim
