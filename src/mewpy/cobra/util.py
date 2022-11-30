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