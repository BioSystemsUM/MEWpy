from typing import Union, TYPE_CHECKING, List, Dict, Callable
from warnings import warn

# TODO: this module depends on pandas dataframes. Should it be set as package requirement?
from pandas import DataFrame, concat

from .milp_bool import milpBool
from .sim_bool import SimBool
from .analysis_utils import decode_solver_solution
from mewpy.algebra import Symbolic
from mewpy.variables import Regulator

if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


# TODO: type hinting and documentation
def slim_milp_bool(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                   milp_bool: milpBool = None,
                   objective: Union[str, Dict[str, Union[float, int]]] = None,
                   minimize: bool = False,
                   initial_state: Dict[str, Union[float, int]] = None,
                   get_values: bool = False) -> Union[int, float, None]:

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

    return decode_solver_solution(solution=sol, minimize=minimize, status=False), sol.values


def slim_sim_bool(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                  sim_bool: SimBool = None,
                  initial_state: Dict[str, Union[float, int]] = None) -> Union[int,
                                                                               float,
                                                                               None]:
    if not sim_bool:
        sim_bool = SimBool(model, build=True, attach=False)

    return sim_bool.optimize(initial_state=initial_state, to_dict=True)


def single_regulator_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                              regulators: List[Union[str, Regulator]] = None,
                              initial_state: Dict[str, Union[float, int]] = None,
                              to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], DataFrame]:
    """

    Single regulator deletions

    :param model:
    :param regulators:
    :param initial_state:
    :param to_dict:

    :return:
    """

    if not initial_state:
        initial_state = {}

    reg_state = {reg.id: reg.coefficient.active_coefficient
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
                                                      active_states=True,
                                                      coefficient=coefficient,
                                                      default=0.0)

                    df.index = [interaction.target.id] * df.shape[0]

                    dfs.append(df)

        df = concat(dfs)

        df = df[['result'] + [col for col in df.columns if col != 'result']]

        del reg_state[regulator_obj.id]

        res[regulator_obj.id] = df

    if to_dict:
        return res

    return concat(res.values(), keys=res.keys())


def regulatory_events(model,
                      initial_state: Dict[str, Union[float, int]] = None,
                      active_states: bool = True,
                      operators: Dict[Symbolic, Callable] = None,
                      decoder: dict = None):

    """

    Regulatory events of a model

    :param model:
    :param initial_state:
    :param active_states:
    :param operators:
    :param decoder:
    :return:
    """

    if not active_states:
        warn('Attention! Active states was set to False without an initial state. '
             'Note that computing all possible states may take some time for large models')

    dfs = []

    for interaction in model.yield_interactions():

        for coefficient, regulatory_event in interaction.regulatory_events.items():

            if not regulatory_event.is_none:
                df = regulatory_event.truth_table(values=initial_state,
                                                  active_states=active_states,
                                                  coefficient=coefficient,
                                                  operators=operators,
                                                  default=0.0,
                                                  decoder=decoder)

                df.index = [interaction.target.id] * df.shape[0]

                dfs.append(df)

    df = concat(dfs)

    df = df[['result'] + [col for col in df.columns if col != 'result']]

    return df
