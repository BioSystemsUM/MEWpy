from typing import Union, TYPE_CHECKING, List, Dict, Callable
from warnings import warn

import pandas as pd

from mewpy.mew.algebra import Symbolic
from mewpy.mew.variables import Regulator
from .analysis_utils import decode_solver_solution
from .milp_bool import milpBool
from .sim_bool import SimBool

if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


def slim_milp_bool(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                   milp_bool: milpBool = None,
                   objective: Union[str, Dict[str, Union[float, int]]] = None,
                   minimize: bool = False,
                   initial_state: Dict[str, Union[float, int]] = None,
                   get_values: bool = True) -> Union[float, Dict[str, float]]:
    """
    A slim Mixed-Integer Boolean simulation (milpBool) of a regulatory model.
    A slim analysis produces a single and simple solution for the model. This method returns a dictionary with the
    solution of the model.

    Fundamentals of the milpBool procedure:
        - A linear problem based on the model regulatory interactions
        - Linear constraints are generated from the regulatory events using mixed-integer linear constraints
        - Regulators are binary variables of the linear problem
        - The linear problem is solved using a mixed-integer linear solver

    :param model: A regulatory model to be simulated
    :param milp_bool: A milpBool object to be used for the simulation. If None, a new milpBool object is created
    :param objective: The objective function of the linear problem. If None, no objective function is used
    :param minimize: If True, the objective function is minimized. If False, the objective function is maximized
    :param initial_state: A dictionary with the initial state of the model. If None, the default initial state is used
    :param get_values: If True, the values of the regulators are returned.
    If False, only the objective function value is returned
    :return: A dictionary with the solution of the model or the objective function value
    """

    if not milp_bool:
        milp_bool = milpBool(model, build=True, attach=False)

    # noinspection PyProtectedMember
    old_linear_obj = milp_bool._linear_objective.copy()
    # noinspection PyProtectedMember
    old_quadratic_obj = milp_bool._quadratic_objective.copy()
    old_sense = milp_bool.minimize

    sol = milp_bool.optimize(initial_state=initial_state,
                             objective=objective,
                             minimize=minimize,
                             to_solver=True,
                             get_values=get_values)

    milp_bool.set_objective(linear=old_linear_obj, quadratic=old_quadratic_obj, minimize=old_sense)

    if get_values:
        return sol.values

    sol, _ = decode_solver_solution(solution=sol, minimize=minimize)
    return sol


def slim_sim_bool(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                  sim_bool: SimBool = None,
                  initial_state: Dict[str, Union[float, int]] = None) -> Dict[str, float]:
    """
    A slim Boolean simulation (simBool) of a regulatory model.
    A slim analysis produces a single and simple solution for the model. This method returns a dictionary with the
    solution of the model.

    Fundamentals of the simBool procedure:
        - A Boolean network based on the model regulatory interactions
        - The regulatory events are evaluated using Boolean algebra and regulator coefficients are used as weights

    :param model: A regulatory or metabolic-regulatory model to be simulated
    :param sim_bool: A simBool object to be used for the simulation. If None, a new simBool object is created
    :param initial_state: A dictionary with the initial state of the model. If None, the default initial state is used
    :return: A dictionary with the solution of the model
    """
    if not sim_bool:
        sim_bool = SimBool(model, build=True, attach=False)

    return sim_bool.optimize(initial_state=initial_state, to_dict=True)


def single_regulator_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                              regulators: List[Union[str, Regulator]] = None,
                              initial_state: Dict[str, Union[float, int]] = None,
                              to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], pd.DataFrame]:
    """
    Single regulator deletion for a regulatory model.
    The regulatory table of each regulator is evaluated using a slim analysis. The results are stored in a dictionary
    or a pandas DataFrame.

    :param model: A regulatory or metabolic-regulatory model to be simulated
    :param regulators: A list of regulators to be deleted. If None, all regulators are deleted
    :param initial_state: A dictionary with the initial state of the model. If None, the default initial state is used
    :param to_dict: If True, the results are returned as a dictionary.
    If False, the results are returned as a pandas DataFrame
    :return: A dictionary or a pandas DataFrame with the results of the analysis
    """
    if not initial_state:
        initial_state = {}

    reg_state = {reg.id: reg.coefficient.default_coefficient
                 for reg in model.yield_regulators()}

    reg_state.update(initial_state)

    res = {}

    if not regulators:
        regulators = model.yield_regulators()

    for regulator in regulators:

        if isinstance(regulator, str):

            regulator_obj = model.get(regulator)

        else:
            regulator_obj = regulator

        reg_state[regulator_obj.id] = 0.0

        dfs = []

        for interaction in regulator.yield_interactions():

            for coefficient, regulatory_event in interaction.regulatory_events.items():

                if not regulatory_event.is_none:

                    values = {key: val for key, val in reg_state.items()
                              if key in regulatory_event.variables}

                    df = regulatory_event.truth_table(values=values,
                                                      default_coefficients=True,
                                                      coefficient=coefficient,
                                                      default=0.0)

                    df.index = [interaction.target.id] * df.shape[0]

                    dfs.append(df)

        df = pd.concat(dfs)

        df = df[['result'] + [col for col in df.columns if col != 'result']]

        del reg_state[regulator_obj.id]

        res[regulator_obj.id] = df

    if to_dict:
        return res

    return pd.concat(res.values(), keys=res.keys())


def regulatory_events(model,
                      initial_state: Dict[str, Union[float, int]] = None,
                      default_coefficients: bool = True,
                      operators: Dict[Symbolic, Callable] = None,
                      decoder: dict = None) -> pd.DataFrame:
    """
    Regulatory events for a regulatory model.
    The regulatory table of each regulator is evaluated using a slim analysis. The results are stored in a dictionary
    or a pandas DataFrame.

    :param model: A regulatory or metabolic-regulatory model to be simulated
    :param initial_state: A dictionary with the initial state of the model. If None, the default initial state is used
    :param default_coefficients: If True, it only calculates results for the default coefficients. If False, it calculates
    results for all coefficients.
    :param operators: A dictionary with custom operators to be used in the evaluation of the regulatory events
    :param decoder: A dictionary with the decoder to be used in the evaluation of the regulatory events
    :return: A pandas DataFrame with the results of the analysis
    """
    if not default_coefficients:
        warn('Attention! Active states was set to False without an initial state. '
             'Note that computing all possible states may take some time for large models')

    dfs = []

    for interaction in model.yield_interactions():

        for coefficient, regulatory_event in interaction.regulatory_events.items():

            if not regulatory_event.is_none:
                df = regulatory_event.truth_table(values=initial_state,
                                                  default_coefficients=default_coefficients,
                                                  coefficient=coefficient,
                                                  operators=operators,
                                                  default=0.0,
                                                  decoder=decoder)

                df.index = [interaction.target.id] * df.shape[0]

                dfs.append(df)

    df = pd.concat(dfs)

    df = df[['result'] + [col for col in df.columns if col != 'result']]

    return df
