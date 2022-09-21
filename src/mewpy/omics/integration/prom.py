from typing import Union, TYPE_CHECKING

from mewpy.omics import ExpressionSet
from mewpy.omics.preprocessing import quantile_preprocessing_pipeline, target_regulator_interaction_probability

if TYPE_CHECKING:
    from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel


def PROM(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
         expr: ExpressionSet,
         regulator: str):
    """ Run a PROM simulation. Consult mewpy.mew.analysis.PROM for more information.

    Arguments:
        model: an integrated Metabolic-Regulatory model aka MEW model.
        expr (ExpressionSet): transcriptomics data.
        regulator (str): the regulator to be knocked out in the simulation.

    Returns:
        Solution: solution
    """
    expression_df = expr.dataframe
    quantile_expression, binary_expression = quantile_preprocessing_pipeline(expression_df)
    initial_state, _ = target_regulator_interaction_probability(model,
                                                                expression=quantile_expression,
                                                                binary_expression=binary_expression)
    from mewpy.mew.analysis import PROM as PROMLP
    prom = PROMLP(model).build()
    return prom.optimize(initial_state=initial_state, regulators=regulator, to_solver=True)
