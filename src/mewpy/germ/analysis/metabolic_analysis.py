from collections import defaultdict
from typing import Union, TYPE_CHECKING, Dict, Tuple, Sequence, Optional

import pandas as pd

from mewpy.util.constants import ModelConstants
from .analysis_utils import run_method_and_decode
from .fba import FBA
from .pfba import pFBA

if TYPE_CHECKING:
    from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel


def slim_fba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
             objective: Union[str, Dict[str, float]] = None,
             constraints: Dict[str, Tuple[float, float]] = None) -> Optional[float]:
    """
    A Flux Balance Analysis simulation of a metabolic model.
    A slim analysis produces a single and simple solution for the model. This method returns the objective value for the
    FBA simulation.

    Fundamentals of the FBA procedure:
        - A linear problem based on the mass balance constraints
        - Reactions are linear variables constrained by their bounds
        - The objective function is a linear combination of the reactions
        - The objective function is solved using a linear solver

    :param model: a metabolic model to be simulated
    :param objective: the objective function to be used for the simulation.
    If not provided, the default objective is used.
    :param constraints: additional constraints to be used for the simulation.
    :return: the objective value for the simulation
    """
    fba = FBA(model).build()

    objective_value, _ = run_method_and_decode(method=fba, objective=objective, constraints=constraints)
    return objective_value


def slim_pfba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
              objective: Union[str, Dict[str, float]] = None,
              constraints: Dict[str, Tuple[float, float]] = None) -> Optional[float]:
    """
    A parsimonious Flux Balance Analysis simulation of a metabolic model.
    A slim analysis produces a single and simple solution for the model. This method returns the objective value for the
    pFBA simulation.

    Fundamentals of the pFBA procedure:
        - A linear problem based on the mass balance constraints
        - Reactions are linear variables constrained by their bounds
        - The objective function is a linear combination of the reactions plus the sum of the absolute values of the
        reactions
        - The objective function is solved using a linear solver by minimizing the sum of the absolute values of the
        reactions

    :param model: a metabolic model to be simulated
    If not provided, a new instance will be created.
    :param objective: the objective function to be used for the simulation.
    If not provided, the default objective is used.
    :param constraints: additional constraints to be used for the simulation.
    :return: the objective value for the simulation
    """
    pfba = pFBA(model).build()

    objective_value, _ = run_method_and_decode(method=pfba, objective=objective, constraints=constraints)
    return objective_value


def fva(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
        fraction: float = 1.0,
        reactions: Sequence[str] = None,
        objective: Union[str, Dict[str, float]] = None,
        constraints: Dict[str, Tuple[float, float]] = None) -> pd.DataFrame:
    """
    Flux Variability Analysis (FVA) of a metabolic model.
    FVA is a method to determine the minimum and maximum fluxes for each reaction in a metabolic model.
    It can be used to identify the reactions that are limiting the growth of a cell.
    In MEWpy, FVA is performed by solving a linear problem for each reaction in the model.

    :param model: a metabolic model to be simulated
    :param fraction: the fraction of the optimal solution to be used as the upper bound for the objective function
    (default: 1.0)
    :param reactions: the reactions to be simulated (default: all reactions in the model)
    :param objective: the objective function to be used for the simulation (default: the default objective)
    :param constraints: additional constraints to be used for the simulation (default: None)
    :return: a pandas DataFrame with the minimum and maximum fluxes for each reaction
    """
    if not reactions:
        reactions = model.reactions.keys()

    if not constraints:
        constraints = {}

    if objective:
        if hasattr(objective, 'keys'):
            obj = next(iter(objective.keys()))
        else:
            obj = str(objective)

    else:
        obj = next(iter(model.objective)).id

    _fba = FBA(model).build()
    objective_value, _ = run_method_and_decode(method=_fba, objective=objective, constraints=constraints)
    constraints[obj] = (fraction * objective_value, ModelConstants.REACTION_UPPER_BOUND)

    fba = FBA(model).build()

    result = defaultdict(list)
    for rxn in reactions:
        min_val, _ = run_method_and_decode(method=fba, objective={rxn: 1.0}, constraints=constraints, minimize=True)
        result[rxn].append(min_val)

        max_val, _ = run_method_and_decode(method=fba, objective={rxn: 1.0}, constraints=constraints, minimize=False)
        result[rxn].append(max_val)

    return pd.DataFrame.from_dict(data=result, orient='index', columns=['minimum', 'maximum'])


def single_gene_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                         genes: Sequence[str] = None,
                         constraints: Dict[str, Tuple[float, float]] = None) -> pd.DataFrame:
    """
    Single gene deletion analysis of a metabolic model.
    Single gene deletion analysis is a method to determine the effect of deleting each gene in a metabolic model.
    It can be used to identify the genes that are essential for the growth of a cell.
    In MEWpy, single gene deletion analysis is performed by solving a linear problem for each gene in the model.
    A gene knockout can switch off reactions associated with the gene, only if the gene is essential for the reaction.

    :param model: a metabolic model to be simulated
    :param genes: the genes to be simulated (default: all genes in the model)
    :param constraints: additional constraints to be used for the simulation (default: None)
    :return: a pandas DataFrame with the fluxes for each gene
    """
    if not constraints:
        constraints = {}

    if not genes:
        genes = model.yield_genes()
    else:
        genes = [model.genes[gene] for gene in genes if gene in model.genes]

    fba = FBA(model).build()
    wt_objective_value, wt_status = run_method_and_decode(method=fba, constraints=constraints)

    state = {gene.id: max(gene.coefficients) for gene in model.yield_genes()}

    result = {}
    for gene in genes:

        gene_coefficient = state.pop(gene.id, 0.0)
        state[gene.id] = 0.0

        gene_constraints = {}
        for reaction in gene.yield_reactions():

            if reaction.gpr.is_none:
                continue

            gpr_eval = reaction.gpr.evaluate(values=state)

            if gpr_eval:
                continue

            gene_constraints[reaction.id] = (0.0, 0.0)

        if gene_constraints:
            solution, status = run_method_and_decode(method=fba, constraints={**constraints, **gene_constraints})
            result[gene.id] = [solution, status]

        else:
            result[gene.id] = [float(wt_objective_value), str(wt_status)]

        state[gene.id] = gene_coefficient

    return pd.DataFrame.from_dict(data=result, orient='index', columns=['growth', 'status'])


def single_reaction_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                             reactions: Sequence[str] = None,
                             constraints: Dict[str, Tuple[float, float]] = None) -> pd.DataFrame:
    """
    Single reaction deletion analysis of a metabolic model.
    Single reaction deletion analysis is a method to determine the effect of deleting each reaction
    in a metabolic model.
    It can be used to identify the reactions that are essential for the growth of a cell.
    In MEWpy, single reaction deletion analysis is performed by solving a linear problem for each reaction in the model.

    :param model: a metabolic model to be simulated
    :param reactions: the reactions to be simulated (default: all reactions in the model)
    :param constraints: additional constraints to be used for the simulation (default: None)
    :return: a pandas DataFrame with the fluxes for each reaction
    """
    if not reactions:
        reactions = model.reactions.keys()

    if not constraints:
        constraints = {}

    fba = FBA(model).build()

    result = {}
    for reaction in reactions:
        reaction_constraints = {reaction: (0.0, 0.0)}
        solution, status = run_method_and_decode(method=fba, constraints={**constraints, **reaction_constraints})
        result[reaction] = [solution, status]

    return pd.DataFrame.from_dict(data=result, orient='index', columns=['growth', 'status'])
