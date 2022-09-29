from typing import Union, TYPE_CHECKING, Dict, Callable, Any, Type, Sequence
from warnings import warn

import pandas as pd

from mewpy.mew.algebra import Symbolic

if TYPE_CHECKING:
    from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel


def regulatory_truth_table(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                           interactions: Sequence[str] = None,
                           initial_state: Dict[str, float] = None,
                           strategy: str = 'max',
                           operators: Union[Dict[Type[Symbolic], Callable], Dict[Type[Symbolic], Any]] = None,
                           decoder: dict = None) -> pd.DataFrame:
    """
    The regulatory truth table of a regulatory model contains the evaluation of all regulatory events.
    RegulatoryModel's interactions are evaluated using an initial state or regulators coefficients.
    The results are stored in a dictionary or a pandas DataFrame.

    :param model: A regulatory or metabolic-regulatory model to be simulated
    :param interactions: A list of interactions to be evaluated. If None, all interactions are evaluated (default).
    :param initial_state: A dictionary with the initial state of the model. If None, the default initial state is used
    :param strategy: The truth table can be calculated using the maximum or minimum value
    in the variables' coefficients. Otherwise, the truth table is calculated using all variables' coefficients.
    :param operators: A dictionary with custom operators to be used in the evaluation of the regulatory events
    :param decoder: A dictionary with the decoder to be used in the evaluation of the regulatory events
    :return: A pandas DataFrame with the results of the analysis
    """
    if not interactions:
        interactions = model.yield_interactions()
    else:
        interactions = [model.interactions[interaction] for interaction in interactions]

    if initial_state is None and strategy == 'all':
        warn('Attention! Missing initial state and calculating "all" coefficients may take some time for large models!')

    dfs = []
    for interaction in interactions:

        for coefficient, regulatory_event in interaction.regulatory_events.items():

            if not regulatory_event.is_none:
                df = regulatory_event.truth_table(values=initial_state,
                                                  strategy=strategy,
                                                  coefficient=coefficient,
                                                  operators=operators,
                                                  decoder=decoder)

                df.index = [interaction.target.id] * df.shape[0]

                dfs.append(df)

    df = pd.concat(dfs)
    result_col = df.pop('result')
    df = pd.concat([result_col, df], axis=1)
    return df
