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

Author: Vitor Pereira
##############################################################################   
"""
from mewpy.solvers import solver_instance
from mewpy.solvers.solver import VarType
from mewpy.solvers.solution import Status
from mewpy.simulation import Environment
from mewpy.cobra.medium import minimal_medium
from mewpy.util.constants import ModelConstants
from mewpy.util import AttrDict

from warnings import warn
from collections import Counter
from itertools import combinations, chain
from math import isinf, inf


def sc_score(community, environment=None, min_growth=0.1, n_solutions=100, verbose=True, abstol=1e-6,
             use_pool=True):
    """
    Calculate frequency of community species dependency on each other
    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)
    Args:
        community (CommunityModel): microbial community
        environment (Environment): metabolic environment (optional)
        min_growth (float): minimum growth rate (default: 0.1)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)
        n_solutions (int): number of alternative solutions to calculate (default: 100)
    Returns:
        dict: Keys are dependent organisms, values are dictionaries with required organism frequencies
    """

    community = community.copy()
    community.add_compartments = True
    comm_model = community.get_community_model()
    if environment:
        environment.apply(comm_model, inplace=True, warning=False)

    for b in community.organisms_biomass.values():
        comm_model.set_reaction_bounds(b, 0, ModelConstants.REACTION_UPPER_BOUND, False)

    solver = solver_instance(comm_model)

    for org_id in community.organisms:
        org_var = 'y_{}'.format(org_id)
        solver.add_variable(org_var, 0, 1, vartype=VarType.BINARY, update=False)

    solver.update()

    bigM = ModelConstants.REACTION_UPPER_BOUND
    for org_id, sim in community.organisms.items():
        org_var = 'y_{}'.format(org_id)
        rxns = set(sim.reactions)-set(sim.get_exchange_reactions())
        for rxn in rxns:
            r_id = community.reaction_map[(org_id, rxn)]
            if r_id == community.organisms_biomass[org_id]:
                continue
            solver.add_constraint('c_{}_lb'.format(r_id), {r_id: 1, org_var: bigM}, '>', 0, update=False)
            solver.add_constraint('c_{}_ub'.format(r_id), {r_id: 1, org_var: -bigM}, '<', 0, update=False)

    solver.update()

    scores = AttrDict()

    for org_id, biomass_id in community.organisms_biomass.items():
        other = {o for o in community.organisms if o != org_id}
        solver.add_constraint('COM_Biomass', {biomass_id: 1}, '>', min_growth)
        objective = {"y_{}".format(o): 1.0 for o in other}

        if not use_pool:
            previous_constraints = []
            donors_list = []
            failed = False

            for i in range(n_solutions):
                sol = solver.solve(objective, minimize=True, get_values=list(objective.keys()))

                if sol.status != Status.OPTIMAL:
                    failed = i == 0
                    break

                donors = [o for o in other if sol.values["y_{}".format(o)] > abstol]
                donors_list.append(donors)

                previous_con = 'iteration_{}'.format(i)
                previous_constraints.append(previous_con)
                previous_sol = {"y_{}".format(o): 1 for o in donors}
                solver.add_constraint(previous_con, previous_sol, '<', len(previous_sol) - 1)

            solver.remove_constraints(['COM_Biomass'] + previous_constraints)

            if not failed:
                donors_list_n = float(len(donors_list))
                donors_counter = Counter(chain(*donors_list))
                scores[org_id] = {o: donors_counter[o] / donors_list_n for o in other}
            else:
                if verbose:
                    warn('SCS: Failed to find a solution for growth of ' + org_id)
                scores[org_id] = None

        else:
            sols = solver.solve(objective, minimize=True, get_values=list(objective.keys()),
                                pool_size=n_solutions, pool_gap=0.5)
            solver.remove_constraint('COM_Biomass')

            if len(sols) == 0:
                scores[org_id] = None
                if verbose:
                    warn('SCS: Failed to find a solution for growth of ' + org_id)
            else:
                donor_count = [o for sol in sols for o in other if sol.values["y_{}".format(o)] > abstol]
                donor_count = Counter(donor_count)
                scores[org_id] = {o: donor_count[o] / len(sols) for o in other}

    return scores


def mu_score(community, environment=None, min_mol_weight=False, min_growth=0.1, max_uptake=10.0,
             abstol=1e-6, validate=False, n_solutions=100, pool_gap=0.5, verbose=True):
    """
    Calculate frequency of metabolite requirement for species growth
    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)
    Args:
        community (CommunityModel): microbial community
        environment (Environment): metabolic environment
        min_mol_weight (bool): Prefer smaller compounds (default: False)
        min_growth (float): minimum growth rate (default: 0.1)
        max_uptake (float): maximum uptake rate (default: 10)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)
        validate (bool): validate solution using FBA (for debugging purposes, default: False)
        n_solutions (int): number of alternative solutions to calculate (default: 100)
    Returns:
        dict: Keys are organism names, values are dictionaries with metabolite frequencies 
        dict: Extra information
    """
    community.add_compartments = True
    sim = community.get_community_model()
    
    def ex_met(r_id,original=False):
        met = list(sim.get_reaction_metabolites(r_id).keys())[0]
        if original:
            for k,v in community.metabolite_map.items():
                if v==met:
                    return k[1]
        else:    
            return met
    
    if environment:
        environment.apply(sim, inplace=True, warning=False)

    max_uptake = max_uptake * len(community.organisms)
    scores = AttrDict()
    
    solver = solver_instance(sim)

    for org_id in community.organisms:
        org_ex = community.organisms[org_id].get_exchange_reactions()
        exchange_rxns = [community.reaction_map[(org_id,rx_id)] for rx_id in org_ex]
        biomass_reaction = community.organisms_biomass[org_id]

        medium_list, sols = minimal_medium(sim, exchange_reactions=exchange_rxns,
                                           min_mass_weight=min_mol_weight, min_growth=min_growth,
                                           n_solutions=n_solutions, max_uptake=max_uptake, validate=validate,
                                           abstol=abstol, use_pool=True, pool_gap=pool_gap,
                                           warnings=False,solver=solver,biomass_reaction=biomass_reaction)

        if medium_list:
            counter = Counter(chain(*medium_list))

            scores[org_id] = {ex_met(ex,True): counter[ex] / len(medium_list)
                              for ex in exchange_rxns}
        else:
            if verbose:
                warn('MUS: Failed to find a minimal growth medium for ' + org_id)
            scores[org_id] = None

    return scores


def mp_score(community, environment=None, abstol=1e-3):
    """
    Discover metabolites which species can produce in community
    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)
    Args:
        community (CommunityModel): community object
        environment (Environment): Metabolic environment in which the SMETANA score is colulated
        min_growth (float): minimum growth rate (default: 0.1)
        max_uptake (float): maximum uptake rate (default: 10)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)
    Returns:
        dict: Keys are model names, values are list with produced compounds
        dict: Extra information
    """
    community.add_compartments = True
    sim = community.get_community_model()
    
    def ex_met(r_id,original=False):
        met = list(sim.get_reaction_metabolites(r_id).keys())[0]
        if original:
            for k,v in community.metabolite_map.items():
                if v==met:
                    return k[1]
        else:    
            return met
    
    if environment:
        environment.apply(sim, inplace=True, warning=False)
        env_compounds = environment.get_compounds(fmt_func=lambda x: x[5:-5])
    else:
        env_compounds = set()

    for org_id in community.organisms:
        org_ex = community.organisms[org_id].get_exchange_reactions()
        exchange_rxns = [community.reaction_map[(org_id,rx_id)] for rx_id in org_ex]
        
        for r_id in exchange_rxns:
            rxn = sim.get_reaction(r_id)
            if isinf(rxn.ub):
                sim.set_reaction_bounds(r_id,rxn.lb,1000)
                
    solver = solver_instance(sim)

    scores = AttrDict()

    for org_id in community.organisms:
        org_ex = community.organisms[org_id].get_exchange_reactions()
        exchange_rxns = [community.reaction_map[(org_id,rx_id)] for rx_id in org_ex]
    
        scores[org_id] = {}

        remaining = [r_id for r_id in exchange_rxns if ex_met(r_id) not in env_compounds]

        while len(remaining) > 0:
            sol = solver.solve(linear={r_id: 1 for r_id in remaining}, minimize=False, get_values=remaining)

            if sol.status != Status.OPTIMAL:
                break

            blocked = [r_id for r_id in remaining if sol.values[r_id] < abstol]

            if len(blocked) == len(remaining):
                break

            for r_id in remaining:
                if sol.values[r_id] >= abstol:
                    scores[org_id][ex_met(r_id,True)] = 1

            remaining = blocked

        for r_id in remaining:
            sol = solver.solve(linear={r_id: 1}, minimize=False, get_values=False)
            

            if sol.status == Status.OPTIMAL and sol.fobj > abstol:
                scores[org_id][ex_met(r_id,True)] = 1
            else:
                scores[org_id][ex_met(r_id,True)] = 0

    return scores


def mip_score(community, environment=None, min_mol_weight=False, min_growth=0.1, direction=-1, max_uptake=10,
              validate=False, verbose=True, use_lp=False, exclude=None):
    """
    Implements the metabolic interaction potential (MIP) score as defined in (Zelezniak et al, 2015).
    Args:
        community (CommunityModel): microbial community model
        environment (Environment): Metabolic environment in which the SMETANA score is calculated
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        extracellular_id (str): extracellular compartment id
        min_mol_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 0.1)
        max_uptake (float): maximum uptake rate (default: 10)
        validate (bool): validate solution using FBA (for debugging purposes, default: False)
    Returns:
        float: MIP score
    """
    community.add_compartments = True
    sim = community.get_community_model()
    
    def ex_met(r_id,trim=False):
        met = list(sim.get_reaction_metabolites(r_id).keys())[0]
        if trim:
            return met[len(sim._m_prefix):].split('_')[0]
        else:
            return met
    # revisit for noninteracting
    noninteracting = community.copy()
    
    exch_reactions = set(sim.get_exchange_reactions())
    max_uptake = max_uptake * len(community.organisms)

    if environment:
        environment.apply(noninteracting.merged, inplace=True, warning=False)
        exch_reactions &= set(environment)

    noninteracting_medium, sol1 = minimal_medium(noninteracting.get_community_model(), exchange_reactions=exch_reactions,
                                                 direction=direction, min_mass_weight=min_mol_weight,
                                                 min_growth=min_growth, max_uptake=max_uptake, validate=validate,
                                                 warnings=False, milp=(not use_lp))
    if noninteracting_medium is None:
        if verbose:
            warn('MIP: Failed to find a valid solution for non-interacting community')
        return None, None

    # anabiotic environment is limited to non-interacting community minimal media
    noninteracting_env = Environment.from_reactions(noninteracting_medium, max_uptake=max_uptake)
    noninteracting_env.apply(sim, inplace=True)

    interacting_medium, sol2 = minimal_medium(sim, direction=direction, exchange_reactions=noninteracting_medium,
                                              min_mass_weight=min_mol_weight, min_growth=min_growth, milp=(not use_lp),
                                              max_uptake=max_uptake, validate=validate, warnings=False)

    if interacting_medium is None:
        if verbose:
            warn('MIP: Failed to find a valid solution for interacting community')
        return None, None

    if exclude is not None:
        exclude_rxns = {'R_EX_M_{}_e_pool'.format(x) for x in exclude}
        interacting_medium = set(interacting_medium) - exclude_rxns
        noninteracting_medium = set(noninteracting_medium) - exclude_rxns

    score = len(noninteracting_medium) - len(interacting_medium)

    noninteracting_medium = [ex_met(r_id,True) for r_id in noninteracting_medium]
    interacting_medium = [ex_met(r_id,True) for r_id in interacting_medium]

    extras = {
        'noninteracting_medium': noninteracting_medium,
        'interacting_medium': interacting_medium
    }

    return score, extras


def mro_score(community, environment=None, direction=-1, min_mol_weight=False, min_growth=0.1, max_uptake=10,
              validate=False, verbose=True, use_lp=False, exclude=None):
    """
    Implements the metabolic resource overlap (MRO) score as defined in (Zelezniak et al, 2015).
    Args:
        community (CommunityModel): microbial community model
        environment (Environment): Metabolic environment in which the SMETANA score is colulated
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        extracellular_id (str): extracellular compartment id
        min_mol_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 0.1)
        max_uptake (float): maximum uptake rate (default: 10)
    Returns:
        float: MRO score
    """
    community.add_compartments = True
    sim = community.get_community_model()
    
    def ex_met(r_id,trim=False):
        met = list(sim.get_reaction_metabolites(r_id).keys())[0]
        if trim:
            return met[len(sim._m_prefix):].split('_')[0]
        else:
            return met
    
    exch_reactions = set(sim.get_exchange_reactions())
    max_uptake = max_uptake * len(community.organisms)

    if environment:
        environment.apply(sim, inplace=True, warning=False)
        exch_reactions &= set(environment)

    medium, sol = minimal_medium(sim, exchange_reactions=exch_reactions, direction=direction,
                                 min_mass_weight=min_mol_weight, min_growth=min_growth, max_uptake=max_uptake,
                                 validate=validate,  warnings=False, milp=(not use_lp),
                                 biomass_reaction=community.biomass)
    if sol.status != Status.OPTIMAL:
        if verbose:
            warn('MRO: Failed to find a valid solution for community')
        return None, None

    interacting_env = Environment.from_reactions(medium, max_uptake=max_uptake)
    interacting_env.apply(sim, inplace=True)

    if exclude is None:
        exclude = set()

    medium = {ex_met(x,True) for x in medium} - exclude
    
    individual_media = AttrDict()

    for org_id in community.organisms:
        biomass_reaction = community.organisms_biomass[org_id]
        
        org_ex = community.organisms[org_id].get_exchange_reactions()
        org_interacting_exch = [community.reaction_map[(org_id,rx_id)] for rx_id in org_ex]

        medium_i, sol = minimal_medium(sim, exchange_reactions=org_interacting_exch, direction=direction,
                                     min_mass_weight=min_mol_weight, min_growth=min_growth, max_uptake=max_uptake,
                                     validate=validate, warnings=False, milp=(not use_lp),
                                     biomass_reaction=biomass_reaction)
        
        if sol.status != Status.OPTIMAL:
            warn('MRO: Failed to find a valid solution for: ' + org_id)
            return None, None

        individual_media[org_id] = {ex_met(r,True) for r in medium_i} - exclude

    pairwise = {(o1, o2): individual_media[o1] & individual_media[o2] for o1, o2 in combinations(community.organisms, 2)}

    numerator = sum(map(len, pairwise.values())) / len(pairwise) if len(pairwise) != 0 else 0
    denominator = sum(map(len, individual_media.values())) / len(individual_media) if len(individual_media) != 0 else 0
    score = numerator / denominator if denominator != 0 else None

    extras = AttrDict({
        'community_medium': medium,
        'individual_media': individual_media
    })

    return score, extras


def minimal_environment(community, aerobic=None, min_mol_weight=False, min_growth=0.1, max_uptake=10,
                        validate=False, verbose=True, use_lp=False, biomass_reaction=None):

    sim = community.get_community_model()
    def ex_by_comp(compound):
        for rx in sim.get_exchange_reactions():
            mets = list(sim.get_reaction_metabolites(rx).keys())
            if len(mets)!=1:
                continue
            formula = sim.get_metabolite(mets[0]).formula
            if formula==compound:
                    return rx
        return None

    exch_reactions = set(sim.get_exchange_reactions())
    r_h2o= ex_by_comp('H2O')
    if r_h2o is not None:
        exch_reactions -= {r_h2o}
        sim.set_reaction_bounds(r_h2o, -inf, inf)

    
    if aerobic is not None:
        r_o2 = ex_by_comp('O2')
        exch_reactions -= {r_o2}
        if aerobic:
            sim.set_reaction_bounds(r_o2, -max_uptake, inf)
        else:
            sim.set_reaction_bounds(r_o2, 0, inf)
    
    ex_rxns, sol = minimal_medium(sim, exchange_reactions=exch_reactions,
        min_mass_weight=min_mol_weight, min_growth=min_growth, milp=(not use_lp),
        max_uptake=max_uptake, validate=validate, warnings=False,biomass_reaction=biomass_reaction)

    if ex_rxns is None:
        if verbose:
            warn('Failed to find a medium for interacting community.')
        return None
    else:
        if aerobic is not None and aerobic and r_o2 is not None:
            ex_rxns |= {r_o2}
        env = Environment.from_reactions(ex_rxns, max_uptake=max_uptake)
        if r_h2o is not None:
            env[r_h2o] = (-inf, inf)
        return env