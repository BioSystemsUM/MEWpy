from collections import namedtuple
from dataclasses import dataclass
from typing import TYPE_CHECKING, Tuple, Union, Dict, Callable

from math import log, exp

from mewpy.mew.algebra import And, Or
from mewpy.solvers.solution import Status
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.solvers import Solution
    from mewpy.mew.lp import LinearProblem
    from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel


def decode_solver_solution(solution: 'Solution') -> Tuple[float, str]:
    """
    It decodes the solution of the solver and returns the objective value.

    :param solution: the solution of the solver
    :return: the objective value and the status of the solution
    """
    sol_status = solution.status.value

    if solution.status == Status.OPTIMAL:
        sol_f_obj = solution.fobj

    elif solution.status == Status.UNBOUNDED or solution.status == Status.INF_OR_UNB:

        sol_f_obj = ModelConstants.REACTION_UPPER_BOUND

    else:
        sol_f_obj = 0.0

    return sol_f_obj, sol_status


def run_method_and_decode(method: 'LinearProblem',
                          objective: Union[str, Dict[str, float]] = None,
                          constraints: Dict[str, Tuple[float, float]] = None,
                          **kwargs) -> Tuple[float, str]:
    """
    It runs a method and decodes the objective value and status returned by the solver.
    :param method: the method to be run
    :param objective: an alternative temporary objective function
    :param constraints: alternative temporary constraints
    :param kwargs: additional arguments to be passed to the method
    :return: the objective value and the status of the solution
    """
    solver_kwargs = {'get_values': False}

    if objective:
        if hasattr(objective, 'keys'):
            solver_kwargs['linear'] = objective.copy()
        else:
            solver_kwargs['linear'] = {str(objective): 1.0}

    if constraints:
        solver_kwargs['constraints'] = constraints

    if 'minimize' in kwargs:
        solver_kwargs['minimize'] = kwargs['minimize']

    solution = method.optimize(to_solver=True, solver_kwargs=solver_kwargs, **kwargs)
    objective_value, status = decode_solver_solution(solution=solution)
    return objective_value, status


# ---------------------------------
# CoRegFlux utils
# ---------------------------------
CoRegMetabolite = namedtuple('CoRegMetabolite', ('id', 'concentration', 'exchange'))
CoRegBiomass = namedtuple('CoRegBiomass', ('id', 'biomass_yield'))


@dataclass
class CoRegResult:
    objective_value: float = None
    values: Dict[str, float] = None
    metabolites: Dict[str, CoRegMetabolite] = None
    biomass: CoRegBiomass = None
    state: Dict[str, float] = None
    constraints: Dict[str, Tuple[float, float]] = None


def build_metabolites(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                      metabolites: Dict[str, float]) -> Dict[str, CoRegMetabolite]:
    res = {}
    for metabolite, concentration in metabolites.items():
        exchange = model.get(metabolite).exchange_reaction.id

        res[metabolite] = CoRegMetabolite(id=metabolite,
                                          concentration=concentration,
                                          exchange=exchange)
    return res


def build_biomass(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                  biomass: float) -> CoRegBiomass:
    variable = next(iter(model.objective))
    return CoRegBiomass(id=variable.id, biomass_yield=biomass)


def concentration_to_lb(concentration,
                        biomass,
                        time_step):
    return concentration / (biomass * time_step)


def lb_soft_plus(soft_plus, coefficient):
    return -1 * log(1 + exp(soft_plus + abs(coefficient)))


def ub_soft_plus(soft_plus, coefficient):
    return log(1 + exp(soft_plus + coefficient))


def euler_step_biomass(old_biomass_yield, growth_rate, time_step):
    return old_biomass_yield * exp(growth_rate * time_step)


def euler_step_metabolites(metabolite_concentration,
                           metabolite_rate,
                           old_biomass_yield,
                           growth_rate,
                           time_step):
    next_concentration = metabolite_rate / growth_rate * old_biomass_yield * (1 - exp(growth_rate * time_step))
    return metabolite_concentration - next_concentration


def biomass_yield_to_rate(biomass):
    return 1 * biomass


def metabolites_constraints(constraints: Dict[str, Tuple[float, float]],
                            metabolites: Dict[str, CoRegMetabolite],
                            biomass: CoRegBiomass,
                            time_step: float):
    constraints = constraints.copy()
    for metabolite in metabolites.values():

        next_lb = concentration_to_lb(concentration=metabolite.concentration,
                                      biomass=biomass.biomass_yield,
                                      time_step=time_step)

        rxn = metabolite.exchange

        if next_lb < ModelConstants.TOLERANCE:
            next_lb = 0

        prev_lb, prev_ub = constraints[rxn]

        if prev_lb == 0:
            next_lb = -next_lb

        else:

            if abs(prev_lb) < abs(next_lb):
                next_lb = prev_lb

            else:
                next_lb = -abs(next_lb)

        constraints[rxn] = (next_lb, prev_ub)

    return constraints


def continuous_gpr(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                   state: Dict[str, float],
                   scale: bool = False):
    operators = {And: min, Or: max}

    states = {}
    for reaction in model.yield_reactions():

        if reaction.gpr.is_none:
            continue

        if not set(reaction.genes).issubset(state):
            continue

        states[reaction.id] = reaction.gpr.evaluate(values=state, operators=operators, missing_value=0)

    if scale:
        _max_state = max(states.values())

        return {key: round((val / _max_state), 6) * 1000 for key, val in states.items()}

    return states


def gene_state_constraints(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                           constraints: Dict[str, Tuple[float, float]],
                           state: Dict[str, float],
                           soft_plus: float = 0,
                           tolerance: float = ModelConstants.TOLERANCE,
                           scale: bool = False):
    constraints = constraints.copy()

    reactions_state = continuous_gpr(model=model, state=state, scale=scale)
    # find the reactions which lower bound was changed by the rules
    # find the reactions which upper bound was changed by the rules
    # evaluate the soft plus over the bounds computed using the continuous version of the gpr rules
    # gene_state_bounds does not contain reactions without gpr rules
    for reaction, coefficient in reactions_state.items():

        old_lb, old_ub = constraints[reaction]

        if coefficient > tolerance:

            if abs(old_lb) > tolerance:
                new_lb = lb_soft_plus(soft_plus=soft_plus, coefficient=coefficient)
            else:
                new_lb = old_lb

            if old_ub > tolerance:
                new_ub = ub_soft_plus(soft_plus=soft_plus, coefficient=coefficient)
            else:
                new_ub = old_ub

        else:
            new_lb = old_lb
            new_ub = old_ub

        constraints[reaction] = (new_lb, new_ub)

    return constraints


def system_state_update(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                        flux_state: Dict[str, float],
                        metabolites: Dict[str, CoRegMetabolite],
                        biomass: CoRegBiomass,
                        time_step: float,
                        biomass_fn: Callable = None) -> Tuple[CoRegBiomass, Dict[str, CoRegMetabolite]]:

    growth_rate = flux_state[biomass.id]
    old_biomass_yield = biomass.biomass_yield

    if biomass_fn is not None:
        growth_rate = biomass_fn(growth_rate)

    if growth_rate == 0:
        growth_rate = ModelConstants.TOLERANCE

    biomass_yield = euler_step_biomass(old_biomass_yield=old_biomass_yield,
                                       growth_rate=growth_rate,
                                       time_step=time_step)

    next_biomass = build_biomass(model, biomass_yield)

    next_metabolites = {}

    for met_id, met in metabolites.items():

        concentration = euler_step_metabolites(metabolite_concentration=met.concentration,
                                               metabolite_rate=flux_state[met.exchange],
                                               old_biomass_yield=old_biomass_yield,
                                               growth_rate=growth_rate,
                                               time_step=time_step)

        if concentration < 0:
            concentration = 0

        next_metabolites[met_id] = CoRegMetabolite(id=met_id,
                                                   concentration=concentration,
                                                   exchange=met.exchange)

    return next_biomass, next_metabolites
