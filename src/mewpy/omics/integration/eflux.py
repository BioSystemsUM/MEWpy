from mewpy.simulation import get_simulator
from mewpy.simulation.simulation import Simulator
from .. import Preprocessing, ExpressionSet


def eFlux(model, expr, condition=0, scale_rxn=None, scale_value=1,
          constraints=None, parsimonious=False, max_exp = None, **kwargs):
    """ Run an E-Flux simulation (Colijn et al, 2009).

    Arguments:
        model: a REFRAMED or COBRApy model or a MEWpy Simulator.
        expr (ExpressionSet): transcriptomics data.
        condition: the condition to use in the simulation\
                (default:0, the first condition is used if more than one.)
            scale_rxn (str): reaction to scale flux vector (optional)
        scale_value (float): scaling factor (mandatory if scale_rxn\
            is specified)
        constraints (dict): additional constraints (optional)
        parsimonious (bool): compute a parsimonious solution (default: False)

    Returns:
        Solution: solution
    """
    if isinstance(model, Simulator):
        sim = model
    else:
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
