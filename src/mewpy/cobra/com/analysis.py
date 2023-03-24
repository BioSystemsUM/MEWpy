from mewpy.solvers import solver_instance
from mewpy.solvers.solver import VarType
from mewpy.solvers.solution import Status
from mewpy.simulation import Environment
from mewpy.cobra.medium import minimal_medium

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

    if environment:
        environment.apply(community.merged, inplace=True, warning=False)

    for b in community.organisms_biomass_reactions.values():
        community.merged.reactions[b].lb = 0

    solver = solver_instance(community.merged)

    for org_id in community.organisms:
        org_var = 'y_{}'.format(org_id)
        solver.add_variable(org_var, 0, 1, vartype=VarType.BINARY, update=False)

    solver.update()

    bigM = 1000
    for org_id, rxns in community.organisms_reactions.items():
        org_var = 'y_{}'.format(org_id)
        for r_id in rxns:
            if r_id == community.organisms_biomass_reactions[org_id]:
                continue
            solver.add_constraint('c_{}_lb'.format(r_id), {r_id: 1, org_var: bigM}, '>', 0, update=False)
            solver.add_constraint('c_{}_ub'.format(r_id), {r_id: 1, org_var: -bigM}, '<', 0, update=False)

    solver.update()

    scores = {}

    for org_id, biomass_id in community.organisms_biomass_reactions.items():
        other = {o for o in community.organisms if o != org_id}
        solver.add_constraint('COM_Biomass', {community.organisms_biomass_reactions[org_id]: 1}, '>', min_growth)
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
        community (Community): microbial community
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

    if environment:
        environment.apply(community.merged, inplace=True, warning=False)

    max_uptake = max_uptake * len(community.organisms)
    scores = {}
    solver = solver_instance(community.merged)

    for org_id in community.organisms:
        exchange_rxns = community.organisms_exchange_reactions[org_id]
        biomass_reaction = community.organisms_biomass_reactions[org_id]
        community.merged.biomass_reaction = biomass_reaction

        medium_list, sols = minimal_medium(community.merged, exchange_reactions=list(exchange_rxns.keys()),
                                           min_mass_weight=min_mol_weight, min_growth=min_growth,
                                           n_solutions=n_solutions, max_uptake=max_uptake, validate=validate,
                                           abstol=abstol, use_pool=True, pool_gap=pool_gap, solver=solver,
                                           warnings=False)

        if medium_list:
            counter = Counter(chain(*medium_list))

            scores[org_id] = {cnm.original_metabolite: counter[ex] / len(medium_list)
                              for ex, cnm in exchange_rxns.items()}
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
        community (Community): community object
        environment (Environment): Metabolic environment in which the SMETANA score is colulated
        min_growth (float): minimum growth rate (default: 0.1)
        max_uptake (float): maximum uptake rate (default: 10)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)
    Returns:
        dict: Keys are model names, values are list with produced compounds
        dict: Extra information
    """

    if environment:
        environment.apply(community.merged, inplace=True, warning=False)
        env_compounds = environment.get_compounds(fmt_func=lambda x: x[5:-5])
    else:
        env_compounds = set()

    for exchange_rxns in community.organisms_exchange_reactions.values():
        for r_id in exchange_rxns.keys():
            rxn = community.merged.reactions[r_id]
            if isinf(rxn.ub):
                rxn.ub = 1000

    solver = solver_instance(community.merged)

    scores = {}

    for org_id, exchange_rxns in community.organisms_exchange_reactions.items():
        scores[org_id] = {}

        remaining = [r_id for r_id, cnm in exchange_rxns.items() if cnm.original_metabolite not in env_compounds]

        while len(remaining) > 0:
            sol = solver.solve(linear={r_id: 1 for r_id in remaining}, minimize=False, get_values=remaining)

            if sol.status != Status.OPTIMAL:
                break

            blocked = [r_id for r_id in remaining if sol.values[r_id] < abstol]

            if len(blocked) == len(remaining):
                break

            for r_id in remaining:
                if sol.values[r_id] >= abstol:
                    cnm = exchange_rxns[r_id]
                    scores[org_id][cnm.original_metabolite] = 1

            remaining = blocked

        for r_id in remaining:
            sol = solver.solve(linear={r_id: 1}, minimize=False, get_values=False)
            cnm = exchange_rxns[r_id]

            if sol.status == Status.OPTIMAL and sol.fobj > abstol:
                scores[org_id][cnm.original_metabolite] = 1
            else:
                scores[org_id][cnm.original_metabolite] = 0

    return scores


def mip_score(community, environment=None, min_mol_weight=False, min_growth=0.1, direction=-1, max_uptake=10,
              validate=False, verbose=True, use_lp=False, exclude=None):
    """
    Implements the metabolic interaction potential (MIP) score as defined in (Zelezniak et al, 2015).
    Args:
        community (Community): microbial community model
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

    noninteracting = community.copy(copy_models=False, interacting=False)
    exch_reactions = set(community.merged.get_exchange_reactions())
    max_uptake = max_uptake * len(community.organisms)

    if environment:
        environment.apply(noninteracting.merged, inplace=True, warning=False)
        exch_reactions &= set(environment)

    noninteracting_medium, sol1 = minimal_medium(noninteracting.merged, exchange_reactions=exch_reactions,
                                                 direction=direction, min_mass_weight=min_mol_weight,
                                                 min_growth=min_growth, max_uptake=max_uptake, validate=validate,
                                                 warnings=False, milp=(not use_lp))
    if noninteracting_medium is None:
        if verbose:
            warn('MIP: Failed to find a valid solution for non-interacting community')
        return None, None

    # anabiotic environment is limited to non-interacting community minimal media
    noninteracting_env = Environment.from_reactions(noninteracting_medium, max_uptake=max_uptake)
    noninteracting_env.apply(community.merged, inplace=True)

    interacting_medium, sol2 = minimal_medium(community.merged, direction=direction, exchange_reactions=noninteracting_medium,
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

    noninteracting_medium = [r_id[7:-7] for r_id in noninteracting_medium]
    interacting_medium = [r_id[7:-7] for r_id in interacting_medium]

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
        community (Community): microbial community model
        environment (Environment): Metabolic environment in which the SMETANA score is colulated
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        extracellular_id (str): extracellular compartment id
        min_mol_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 0.1)
        max_uptake (float): maximum uptake rate (default: 10)
    Returns:
        float: MRO score
    """

    exch_reactions = set(community.merged.get_exchange_reactions())
    max_uptake = max_uptake * len(community.organisms)

    if environment:
        environment.apply(community.merged, inplace=True, warning=False)
        exch_reactions &= set(environment)

    medium, sol = minimal_medium(community.merged, exchange_reactions=exch_reactions, direction=direction,
                                 min_mass_weight=min_mol_weight, min_growth=min_growth, max_uptake=max_uptake,
                                 validate=validate,  warnings=False, milp=(not use_lp))

    if sol.status != Status.OPTIMAL:
        if verbose:
            warn('MRO: Failed to find a valid solution for community')
        return None, None

    interacting_env = Environment.from_reactions(medium, max_uptake=max_uptake)
    interacting_env.apply(community.merged, inplace=True)

    if exclude is None:
        exclude = set()

    medium = {x[7:-7] for x in medium} - exclude
    individual_media = {}
    solver = solver_instance(community.merged)

    for org_id in community.organisms:
        biomass_reaction = community.organisms_biomass_reactions[org_id]
        community.merged.biomass_reaction = biomass_reaction
        org_interacting_exch = community.organisms_exchange_reactions[org_id]

        medium_i, sol = minimal_medium(community.merged, exchange_reactions=org_interacting_exch, direction=direction,
                                     min_mass_weight=min_mol_weight, min_growth=min_growth, max_uptake=max_uptake,
                                     validate=validate, solver=solver, warnings=False, milp=(not use_lp))

        if sol.status != Status.OPTIMAL:
            warn('MRO: Failed to find a valid solution for: ' + org_id)
            return None, None

        individual_media[org_id] = {org_interacting_exch[r].original_metabolite[2:-2] for r in medium_i} - exclude

    pairwise = {(o1, o2): individual_media[o1] & individual_media[o2] for o1, o2 in combinations(community.organisms, 2)}

    numerator = sum(map(len, pairwise.values())) / len(pairwise) if len(pairwise) != 0 else 0
    denominator = sum(map(len, individual_media.values())) / len(individual_media) if len(individual_media) != 0 else 0
    score = numerator / denominator if denominator != 0 else None

    extras = {
        'community_medium': medium,
        'individual_media': individual_media
    }

    return score, extras


def minimal_environment(community, aerobic=None, min_mol_weight=False, min_growth=0.1, max_uptake=10,
                        validate=False, verbose=True, use_lp=False):

    exch_reactions = set(community.merged.get_exchange_reactions())

    exch_reactions -= {"R_EX_M_h2o_e_pool"}
    community.merged.set_flux_bounds("R_EX_M_h2o_e_pool", -inf, inf)

    if aerobic is not None:
        exch_reactions -= {"R_EX_M_o2_e_pool"}
        if aerobic:
            community.merged.set_flux_bounds("R_EX_M_o2_e_pool", -max_uptake, inf)
        else:
            community.merged.set_flux_bounds("R_EX_M_o2_e_pool", 0, inf)

    ex_rxns, sol = minimal_medium(community.merged, exchange_reactions=exch_reactions,
        min_mass_weight=min_mol_weight, min_growth=min_growth, milp=(not use_lp),
        max_uptake=max_uptake, validate=validate, warnings=False)

    if ex_rxns is None:
        if verbose:
            warn('Failed to find a medium for interacting community.')
        return None
    else:
        if aerobic is not None and aerobic:
            ex_rxns |= {"R_EX_M_o2_e_pool"}
        env = Environment.from_reactions(ex_rxns, max_uptake=max_uptake)
        env["R_EX_M_h2o_e_pool"] = (-inf, inf)
        return env