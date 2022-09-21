from typing import Union, TYPE_CHECKING, Dict, Sequence

from mewpy.util.constants import ModelConstants

from mewpy.omics import ExpressionSet

if TYPE_CHECKING:
    from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel


def CoRegFlux(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
              expr: ExpressionSet,
              condition: str,
              dynamic: bool = False,
              metabolites: Dict[str, float] = None,
              growth_rate: float = None,
              time_steps: Sequence[float] = None,
              soft_plus: float = 0,
              tolerance: float = ModelConstants.TOLERANCE,
              scale: bool = False):
    """ Run a CoRegFlux simulation. Consult mewpy.mew.analysis.CoRegFlux for more information.

    Arguments:
        model: an integrated Metabolic-Regulatory model aka MEW model.
        expr (ExpressionSet): transcriptomics data.
        condition (str): the condition to be simulated.
        dynamic (bool): whether to run a dynamic simulation or not.
        metabolites (dict): initial metabolite concentrations.
        growth_rate (float): initial growth rate.
        time_steps (list): time steps for dynamic simulations.
        soft_plus (float): soft plus parameter.
        tolerance (float): solver tolerance.
        scale (bool): whether to scale the flux vector or not.

    Returns:
        Solution: solution
    """
    expression_df = expr.dataframe
    initial_state = expression_df[condition].to_dict()

    from mewpy.mew.analysis import CoRegFlux as CoRegFluxLP
    co_reg_flux = CoRegFluxLP(model).build()
    return co_reg_flux.optimize(to_solver=True,
                                initial_state=initial_state,
                                dynamic=dynamic,
                                metabolites=metabolites,
                                growth_rate=growth_rate,
                                time_steps=time_steps,
                                soft_plus=soft_plus,
                                tolerance=tolerance,
                                scale=scale)
