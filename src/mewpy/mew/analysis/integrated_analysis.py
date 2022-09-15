from collections import defaultdict
from typing import Union, TYPE_CHECKING, List, Dict, Tuple, Optional, Sequence

import pandas as pd

from mewpy.util.constants import ModelConstants
from .analysis_utils import run_method_and_decode
from .rfba import RFBA
from .srfba import SRFBA

if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


INTEGRATED_ANALYSIS_METHODS = {'rfba': RFBA,
                               'srfba': SRFBA}


def slim_rfba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
              initial_state: Dict[str, float] = None,
              objective: Union[str, Dict[str, float]] = None,
              constraints: Dict[str, Tuple[float, float]] = None) -> Optional[float]:
    """
    A Regulatory Flux Balance Analysis (RFBA) simulation of an integrated Metabolic-Regulatory model.
    A slim analysis produces a single and simple solution for the model. This method returns the objective value of
    the RFBA simulation.

    Fundamentals of the RFBA procedure:
        - A linear problem based on reactions
        - A synchronous simulation of the metabolic (GPRs) and regulatory networks is performed to obtain the
        metabolic state (boundaries of the reactions) and the regulatory state (boundaries of the regulators)
        - The metabolic state is used to constrain the reactions of the model
        - A FBA simulation is performed to obtain the fluxes of the model

    :param model: an integrated metabolic-regulatory model to be simulated
    :param initial_state: the initial state of the model. If None, the default initial state is used
    :param objective: the objective of the simulation. If None, the default objective is used
    :param constraints: the constraints of the model. If None, the default constraints are used
    :return: the objective value of the simulation
    """
    rfba = RFBA(model).build()

    objective_value, _ = run_method_and_decode(method=rfba, objective=objective, constraints=constraints,
                                               initial_state=initial_state)
    return objective_value


def slim_srfba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
               initial_state: Dict[str, float] = None,
               objective: Union[str, Dict[str, float]] = None,
               constraints: Dict[str, Tuple[float, float]] = None) -> Optional[float]:
    """
    A Synchronous Regulatory Flux Balance Analysis (SRFBA) simulation of an integrated Metabolic-Regulatory model.
    A slim analysis produces a single and simple solution for the model. This method returns the objective value of
    the SRFBA simulation.

    Fundamentals of the SRFBA procedure:
        - A linear problem based on reactions, GPRs, and regulatory interactions using mixed-integer constraints
        - Reactions, genes, and regulators are constrained by the bounds or coefficients
        - SRFBA is solved using a mixed-integer solver

    :param model: an integrated metabolic-regulatory model to be simulated
    :param initial_state: the initial state of the model. If None, the default initial state is used
    :param objective: the objective of the simulation. If None, the default objective is used
    :param constraints: the constraints of the model. If None, the default constraints are used
    :return: the objective value of the simulation
    """
    srfba = SRFBA(model).build()

    objective_value, _ = run_method_and_decode(method=srfba, objective=objective, constraints=constraints,
                                               initial_state=initial_state)
    return objective_value


def ifva(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
         fraction: float = 1.0,
         reactions: Sequence[str] = None,
         objective: Union[str, Dict[str, float]] = None,
         constraints: Dict[str, Tuple[float, float]] = None,
         initial_state: Dict[str, float] = None,
         method: str = 'srfba',
         to_dict: bool = False) -> Union[Dict[str, List[float]], pd.DataFrame]:
    """
    Integrated Flux Variability Analysis (iFVA) of an integrated Metabolic-Regulatory model.
    iFVA is a flux variability analysis method that considers:
        - Metabolites
        - Reactions
        - Genes
        - Regulatory interactions
        - Target genes
        - Regulators
    It can be used to identify the reactions, genes, and regulators that are limiting the growth of a cell.
    In MEWpy, FVA is performed by solving a linear problem for each reaction in the model.
    The method can be either RFBA or SRFBA.

    :param model: an integrated metabolic-regulatory model to be simulated
    :param fraction: the fraction of the optimal solution to be used as the upper bound for the objective function
    (default: 1.0)
    :param reactions: the reactions to be simulated. If None, all reactions are simulated (default: None)
    :param objective: the objective of the simulation. If None, the default objective is used (default: None)
    :param constraints: additional constraints to be added to the model. If None, no additional constraints are added
    :param method: the method to be used for the simulation. Available methods: 'rfba', 'srfba'. Default: 'srfba'
    :param initial_state: the initial state of the model. If None, the default initial state is used (default: None)
    :param to_dict: if True, the results are returned as a dictionary. If False, the results are returned as a DataFrame
    :return: a dictionary or a DataFrame with the results of the simulation
    """
    if not reactions:
        reactions = model.reactions.keys()

    if objective:
        if hasattr(objective, 'keys'):
            obj = next(iter(objective.keys()))
        else:
            obj = str(objective)

    else:
        obj = next(iter(model.objective)).id

    if not constraints:
        constraints = {}

    LP = INTEGRATED_ANALYSIS_METHODS[method]

    _lp = LP(model).build()
    objective_value, _ = run_method_and_decode(method=_lp, objective=objective, constraints=constraints,
                                               initial_state=initial_state)
    constraints[obj] = (fraction * objective_value, ModelConstants.REACTION_UPPER_BOUND)

    lp = LP(model).build()

    result = defaultdict(list)
    for rxn in reactions:
        min_val, _ = run_method_and_decode(method=lp, objective={rxn: 1.0}, constraints=constraints, minimize=True)
        result[rxn].append(min_val)

        max_val, _ = run_method_and_decode(method=lp, objective={rxn: 1.0}, constraints=constraints, minimize=False)
        result[rxn].append(max_val)

    if to_dict:
        return result

    return pd.DataFrame.from_dict(data=result, orient='index', columns=['minimum', 'maximum'])


def isingle_gene_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                          genes: Sequence[str] = None,
                          constraints: Dict[str, Tuple[float, float]] = None,
                          initial_state: Dict[str, float] = None,
                          method: str = 'srfba',
                          to_dict: bool = False) -> Union[Dict[str, List[float]], pd.DataFrame]:
    """
    Integrated single gene deletion analysis of an integrated Metabolic-Regulatory model.
    Integrated single gene deletion analysis is a method to determine the effect of deleting each gene
    in a model.
    It can be used to identify the genes that are essential for the growth of a cell.
    In MEWpy, single gene deletion analysis is performed by solving a linear problem for each gene in the model.
    A gene knockout can switch off reactions associated with the gene, only if the gene is essential for the reaction.
    The methods RFBA and SRFBA can determine if the gene is essential for the growth of a cell.

    :param model: an integrated metabolic-regulatory model to be simulated
    :param genes: the genes to be simulated. If None, all genes are simulated (default: None)
    :param constraints: additional constraints to be added to the model. If None, no additional constraints are added
    :param initial_state: the initial state of the model. If None, the default initial state is used (default: None)
    :param method: the method to be used for the simulation. Available methods: 'rfba', 'srfba'. Default: 'srfba'
    :param to_dict: if True, the results are returned as a dictionary. If False, the results are returned as a DataFrame
    :return: a dictionary or a DataFrame with the results of the simulation
    """
    if not initial_state:
        initial_state = {}

    if not genes:
        genes = model.genes.keys()

    LP = INTEGRATED_ANALYSIS_METHODS[method]
    lp = LP(model).build()

    result = {}
    for gene in genes:

        gene_coefficient = initial_state.pop(gene, None)
        initial_state[gene] = 0.0

        solution, status = run_method_and_decode(method=lp, constraints=constraints, initial_state=initial_state)

        result[gene] = [solution, status]

        if gene_coefficient is not None:
            initial_state[gene] = gene_coefficient
        else:
            initial_state.pop(gene)

    if to_dict:
        return result

    return pd.DataFrame.from_dict(data=result, orient='index', columns=['growth', 'status'])


def isingle_reaction_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                              reactions: Sequence[str] = None,
                              constraints: Dict[str, Tuple[float, float]] = None,
                              initial_state: Dict[str, float] = None,
                              method: str = 'srfba',
                              to_dict: bool = False) -> Union[Dict[str, List[float]], pd.DataFrame]:
    """
    Integrated single reaction deletion analysis of an integrated Metabolic-Regulatory model.
    Integrated single reaction deletion analysis is a method to determine the effect of deleting each reaction
    in a model.
    It can be used to identify the reactions that are essential for the growth of a cell.
    In MEWpy, single reaction deletion analysis is performed by solving a linear problem for each reaction in the model.
    The methods RFBA and SRFBA can determine if the reaction is essential for the growth of a cell.

    :param model: an integrated metabolic-regulatory model to be simulated
    :param reactions: the reactions to be simulated. If None, all reactions are simulated (default: None)
    :param constraints: additional constraints to be added to the model. If None, no additional constraints are added
    :param initial_state: the initial state of the model. If None, the default initial state is used (default: None)
    :param method: the method to be used for the simulation. Available methods: 'rfba', 'srfba'. Default: 'srfba'
    :param to_dict: if True, the results are returned as a dictionary. If False, the results are returned as a DataFrame
    :return: a dictionary or a DataFrame with the results of the simulation
    """
    if not constraints:
        constraints = {}

    if not reactions:
        reactions = model.reactions.keys()

    LP = INTEGRATED_ANALYSIS_METHODS[method]
    lp = LP(model).build()

    result = {}
    for reaction in reactions:

        reaction_constraint = constraints.pop(reaction, None)
        constraints[reaction] = (0.0, 0.0)

        solution, status = run_method_and_decode(method=lp, constraints=constraints, initial_state=initial_state)

        result[reaction] = [solution, status]

        if reaction_constraint is not None:
            constraints[reaction] = reaction_constraint
        else:
            constraints.pop(reaction)

    if to_dict:
        return result

    return pd.DataFrame.from_dict(data=result, orient='index', columns=['growth', 'status'])


def isingle_regulator_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                               regulators: Sequence[str] = None,
                               constraints: Dict[str, Tuple[float, float]] = None,
                               initial_state: Dict[str, float] = None,
                               method: str = 'srfba',
                               to_dict: bool = False) -> Union[Dict[str, List[float]], pd.DataFrame]:
    """
    Integrated single regulator deletion analysis of a regulatory model.
    Integrated single regulator deletion analysis is a method to determine the effect of deleting each regulator
    in a regulatory model.
    It can be used to identify the regulators that are essential for the growth of a cell.
    In MEWpy, single regulator deletion analysis is performed by solving a linear problem
    for each regulator in the model.
    A regulator knockout can switch off reactions associated with the regulator,
    only if the regulator is essential for the interactions and gprs.
    The methods RFBA and SRFBA can determine if the regulator is essential for the growth of a cell.

    :param model: an integrated metabolic-regulatory model to be simulated
    :param regulators: the regulators to be simulated. If None, all regulators are simulated (default: None)
    :param constraints: additional constraints to be added to the model. If None, no additional constraints are added
    :param initial_state: the initial state of the model. If None, the default initial state is used (default: None)
    :param method: the method to be used for the simulation. Available methods: 'rfba', 'srfba'. Default: 'srfba'
    :param to_dict: if True, the results are returned as a dictionary. If False, the results are returned as a DataFrame
    :return: a dictionary or a DataFrame with the results of the simulation
    """
    if not initial_state:
        initial_state = {}

    if not regulators:
        regulators = model.regulators.keys()

    LP = INTEGRATED_ANALYSIS_METHODS[method]
    lp = LP(model).build()

    result = {}
    for regulator in regulators:

        regulator_coefficient = initial_state.pop(regulator, None)
        initial_state[regulator] = 0.0

        solution, status = run_method_and_decode(method=lp, constraints=constraints, initial_state=initial_state)

        result[regulator] = [solution, status]

        if regulator_coefficient is not None:
            initial_state[regulator] = regulator_coefficient
        else:
            initial_state.pop(regulator)

    if to_dict:
        return result

    return pd.DataFrame.from_dict(data=result, orient='index', columns=['growth', 'status'])
