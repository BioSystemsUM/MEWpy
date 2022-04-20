from math import inf
import copy

import numpy as np

# from ...solvers import solver_instance
from ...simulation import get_simulator
from ...simulation.simulation import Simulator
from .. import Preprocessing


def MADE(model, expr, p_vals, fold_change, biomass=None, constraints=None,
         opt_frac=0.9, weight='log', inline=False, condition=0):

    if isinstance(model, Simulator):
        sim = model
    else:
        sim = get_simulator(model)

    # solver = solver_instance(sim)

    if biomass is None:
        try:
            biomass = list(sim.objective.keys())[0]
        except Exception:
            raise ValueError(
                "A biomass reaction identifier is required or "
                "needs to be set as the model objective")

    wt_solution = sim.simulate(constraints=constraints)

    constraints[biomass] = (opt_frac * wt_solution.fluxes[biomass], inf)

    # Define weighting type
    weights = ['log', 'linear', 'unit', 'none']
    # weighting = ''
    if weight in weights:
        # weighting = weight
        pass
    else:
        pass
        # weighting = 'log'

    assert len(expr) > 2, 'MADE requires at least three inputs'

    if p_vals is None:
        p_vals = np.ones(np.shape(fold_change))
        # weighting = 'none'
    else:
        assert np.shape(p_vals) == np.shape(fold_change), 'Fold_change and P_values must have the same dimentions'

    # Transitions and conditions
    transi = np.shape(fold_change)[1]
    condi = transi + 1

    models = {}
    for i in range(len(condi)):
        models[i] = copy.deepcopy(model)

    n = 1
    expression_sets = {}
    for exp in expr:
        n += 1
        exp_id = 'expression_' + str(n)
        pp = Preprocessing(sim, exp)
        coeffs, thresh = pp.percentile(condition)
        expression_sets[exp_id] = coeffs

    dif = {}
    binary_states = {}

    for i in range(len(expression_sets.keys())-1):
        dif[i] = {}
        binary_states[i] = {}
        expr_1 = expression_sets['expression_' + str(i+1)]
        expr_2 = expression_sets['expression_' + str(i+2)]
        for rxn_id in expr_1.keys():
            if expr_1[rxn_id] < expr_2[rxn_id]:
                dif[i][rxn_id] = -1
            if expr_1[rxn_id] > expr_2[rxn_id]:
                dif[i][rxn_id] = +1
            if expr_1[rxn_id] == expr_2[rxn_id]:
                dif[i][rxn_id] = 0

    for key in binary_states.keys():
        for rxn, val in dif[key].items():
            # if val
            pass
