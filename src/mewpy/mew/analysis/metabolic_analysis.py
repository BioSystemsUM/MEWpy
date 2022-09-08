from collections import defaultdict
from typing import Union, TYPE_CHECKING, List, Dict, Tuple

import pandas as pd

from mewpy.util.constants import ModelConstants
from mewpy.mew.variables import Reaction, Gene
from .fba import FBA, pFBA, milpFBA
from .analysis_utils import decode_solver_solution

if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


def _fba(method, objective, minimize, constraints):
    """
    Performs a FBA simulation.
    Internal use only.
    :param method: a FBA-like method
    :param objective: objective
    :param minimize: minimize
    :param constraints: constraints
    :return: decoded solution
    """

    # noinspection PyProtectedMember
    old_linear_obj = method._linear_objective.copy()
    # noinspection PyProtectedMember
    old_quadratic_obj = method._quadratic_objective.copy()
    old_sense = method.minimize

    sol = method.optimize(objective=objective,
                          minimize=minimize,
                          constraints=constraints,
                          to_solver=True,
                          get_values=False)

    method.set_objective(linear=old_linear_obj, quadratic=old_quadratic_obj, minimize=old_sense)

    return decode_solver_solution(solution=sol, minimize=minimize)


def slim_fba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
             fba: FBA = None,
             objective: Union[str, Dict[str, Union[float, int]]] = None,
             minimize: bool = False,
             constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None) -> Union[int, float, None]:
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
    :param fba: a FBA instance to be used for the simulation (optional).
    If not provided, a new instance will be created.
    :param objective: the objective function to be used for the simulation.
    If not provided, the default objective is used.
    :param minimize: whether to minimize the objective function (default: False)
    :param constraints: additional constraints to be used for the simulation.
    :return: the objective value for the simulation
    """
    if not fba:
        fba = FBA(model, build=True, attach=False)

    sol, _ = _fba(method=fba, objective=objective, minimize=minimize, constraints=constraints)
    return sol


def slim_pfba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
              pfba: pFBA = None,
              objective: Union[str, Dict[str, Union[float, int]]] = None,
              minimize: bool = False,
              constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None) -> Union[int, float, None]:
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
    :param pfba: a pFBA instance to be used for the simulation (optional).
    If not provided, a new instance will be created.
    :param objective: the objective function to be used for the simulation.
    If not provided, the default objective is used.
    :param minimize: Whether to minimize the objective function (default: False)
    :param constraints: additional constraints to be used for the simulation.
    :return: the objective value for the simulation
    """
    if not pfba:
        pfba = pFBA(model, build=True, attach=False)

    sol, _ = _fba(method=pfba, objective=objective, minimize=minimize, constraints=constraints)
    return sol


def slim_milp_fba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                  milp_fba: milpFBA = None,
                  objective: Union[str, Dict[str, Union[float, int]]] = None,
                  minimize: bool = False,
                  constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None) -> Union[int,
                                                                                                       float,
                                                                                                       None]:
    """
    A Mixed-Integer Flux Balance Analysis (milpFBA) simulation of a metabolic model.
    A slim analysis produces a single and simple solution for the model. This method returns the objective value for the
    milpFBA simulation.

    Fundamentals of the FBA procedure:
        - A linear problem based on the mass balance constraints plus the GPR constraints (linearized)
        using mixed-integer variables
        - Reactions are linear variables constrained by their bounds and the gene variables are constrained by their
        coefficients
        - The objective function is a linear combination of the reactions
        - The objective function is solved using a mixed-integer solver

    :param model: a metabolic model to be simulated
    :param milp_fba: a milpFBA instance to be used for the simulation (optional).
    If not provided, a new instance will be created.
    :param objective: the objective function to be used for the simulation.
    If not provided, the default objective is used.
    :param minimize: whether to minimize the objective function (default: False)
    :param constraints: additional constraints to be used for the simulation.
    :return: the objective value for the simulation
    """
    if not milp_fba:
        milp_fba = milpFBA(model, build=True, attach=False)

    sol, _ = _fba(method=milp_fba, objective=objective, minimize=minimize, constraints=constraints)
    return sol


def _inputs_processing(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                       method: str = None,
                       fraction: float = None,
                       reactions: List[Union[str, Reaction]] = None,
                       genes: List[Union[str, Gene]] = None,
                       objective: Union[str, Dict[str, Union[float, int]]] = None,
                       constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None):
    """
    A method to process the inputs for the FVA and FVA slim methods.
    :param model: a metabolic model to be simulated
    :param method: the method to be used for the simulation
    :param fraction: the fraction of the optimal solution to be used as the objective function
    :param reactions: a list of reactions to be used for the simulation
    :param genes: a list of genes to be used for the simulation
    :param objective: the objective function to be used for the simulation
    :param constraints: additional constraints to be used for the simulation
    :return: the processed inputs
    """
    methods_map = {'fba': FBA,
                   'milp': milpFBA,
                   'pfba': pFBA}

    if not constraints:
        constraints = {}

    if objective:

        if isinstance(objective, str):

            rxn = model.get(objective)

            if not rxn:
                ValueError(f'The {objective} objective is not available in the {model} model')

            objective = {objective: 1.0}

        elif isinstance(objective, dict):

            initial_obj = objective

            objective = {}

            for key in initial_obj:

                if isinstance(key, tuple):
                    ValueError('Quadratic objectives are not allowed')

                rxn = model.get(key)

                if not rxn:
                    ValueError(f'The {key} objective is not available in the {model} model')

        else:
            raise TypeError(f'Objective wrong type {type(objective)}')

    else:
        objective = {rxn.id: coef for rxn, coef in model.objective.items()}

        if len(objective) > 1:
            raise ValueError('linear objectives having more than one variable or quadratic objectives are not allowed')

    obj_rxn = next(iter(objective.keys()))

    _FBA = methods_map.get(method, None)

    if _FBA is None:
        raise ValueError(f'Simulation method {method} is not valid. Available methods: {methods_map}')

    fba = _FBA(model, build=True, attach=False)
    fba.set_objective(linear=objective, minimize=False)

    if fraction is not None:

        if fraction > 0.0:
            sol, _ = _fba(method=fba, objective=objective, minimize=False, constraints=constraints)

            constraints[obj_rxn] = (fraction * sol, ModelConstants.REACTION_UPPER_BOUND)

    if not reactions:
        reactions = model.reactions.keys()

    if not genes:
        genes = model.genes.keys()

    return fba, reactions, genes, constraints


def fva(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
        method: str = 'fba',
        fraction: float = 1.0,
        reactions: List[Union[str, Reaction]] = None,
        objective: Union[str, Dict[str, Union[float, int]]] = None,
        constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None,
        to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], pd.DataFrame]:
    """
    Flux Variability Analysis (FVA) of a metabolic model.
    FVA is a method to determine the minimum and maximum fluxes for each reaction in a metabolic model.
    It can be used to identify the reactions that are limiting the growth of a cell.
    In MEWpy, FVA is performed by solving a linear problem for each reaction in the model.
    The method can be either FBA, MILP-FBA or pFBA.

    :param model: a metabolic model to be simulated
    :param method: the method to be used for the simulation (default: 'fba').
    Available methods: 'fba', 'milp', 'pfba'
    :param fraction: the fraction of the optimal solution to be used as the upper bound for the objective function
    (default: 1.0)
    :param reactions: the reactions to be simulated (default: all reactions in the model)
    :param objective: the objective function to be used for the simulation (default: the default objective)
    :param constraints: additional constraints to be used for the simulation (default: None)
    :param to_dict: whether to return the results as a dictionary (default: False)
    :return: a dictionary or a pandas DataFrame with the minimum and maximum fluxes for each reaction
    """
    fba, reactions, _, constraints = _inputs_processing(model=model,
                                                        method=method,
                                                        fraction=fraction,
                                                        reactions=reactions,
                                                        genes=None,
                                                        objective=objective,
                                                        constraints=constraints)

    res = defaultdict(list)

    for rxn in reactions:

        if isinstance(rxn, str):

            rxn_objective = rxn

        else:
            rxn_objective = rxn.id

        sol, _ = _fba(method=fba, objective=rxn_objective, minimize=True, constraints=constraints)

        res[rxn_objective].append(sol)

        sol, _ = _fba(method=fba, objective=rxn_objective, minimize=False, constraints=constraints)

        res[rxn_objective].append(sol)

    if to_dict:
        return res

    return pd.DataFrame.from_dict(data=res, orient='index', columns=['minimum', 'maximum'])


def _milp_gene_deletion(genes, constraints, fba):
    """
    A method to perform gene deletions using MILP-FBA.
    :param genes:
    :param constraints:
    :param fba:
    :return:
    """
    res = {}

    for gene in genes:

        if isinstance(gene, str):

            gene_id = gene

        else:
            gene_id = gene.id

        constraints[gene_id] = (0.0, 0.0)

        sol, status = _fba(method=fba, objective=None, minimize=False, constraints=constraints)

        res[gene_id] = [sol, status]

        del constraints[gene_id]

    return res


def _fba_gene_deletion(model, genes, constraints, fba):
    """
    A method to perform gene deletions using FBA.
    :param model:
    :param genes:
    :param constraints:
    :param fba:
    :return:
    """
    values = {gene.id: 1.0 for gene in model.yield_genes()}

    res = {}

    for gene in genes:

        if isinstance(gene, str):

            gene_obj = model.get(gene)

        else:
            gene_obj = gene

        values[gene_obj.id] = 0.0

        to_remove = []
        for rxn in gene_obj.yield_reactions():

            if not rxn.gpr.is_none:

                gpr_eval = rxn.gpr.evaluate(values=values)

                if not gpr_eval:
                    constraints[rxn.id] = (0.0, 0.0)
                    to_remove.append(rxn.id)

        sol, status = _fba(method=fba, objective=None, minimize=False, constraints=constraints)

        res[gene_obj.id] = [sol, status]

        for rxn in to_remove:
            del constraints[rxn]

    return res


def single_gene_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                         method: str = 'fba',
                         genes: List[Union[str, Gene]] = None,
                         objective: Union[str, Dict[str, Union[float, int]]] = None,
                         constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None,
                         to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], pd.DataFrame]:
    """
    Single gene deletion analysis of a metabolic model.
    Single gene deletion analysis is a method to determine the effect of deleting each gene in a metabolic model.
    It can be used to identify the genes that are essential for the growth of a cell.
    In MEWpy, single gene deletion analysis is performed by solving a linear problem for each gene in the model.
    A gene knockout can switch off reactions associated with the gene, only if the gene is essential for the reaction.
    The methods FBA, MILP-FBA or pFBA can determine if the gene is essential for the growth of a cell.

    :param model: a metabolic model to be simulated
    :param method: the method to be used for the simulation (default: 'fba').
    Available methods: 'fba', 'milp', 'pfba'
    :param genes: the genes to be simulated (default: all genes in the model)
    :param objective: the objective function to be used for the simulation (default: the default objective)
    :param constraints: additional constraints to be used for the simulation (default: None)
    :param to_dict: whether to return the results as a dictionary (default: False)
    :return: a dictionary or a pandas DataFrame with the fluxes for each gene
    """
    fba, _, genes, constraints = _inputs_processing(model=model,
                                                    method=method,
                                                    fraction=None,
                                                    reactions=None,
                                                    genes=genes,
                                                    objective=objective,
                                                    constraints=constraints)

    if method == 'milp':

        res = _milp_gene_deletion(genes=genes, constraints=constraints, fba=fba)

    else:

        res = _fba_gene_deletion(model=model, genes=genes, constraints=constraints, fba=fba)

    if to_dict:
        return res

    return pd.DataFrame.from_dict(data=res, orient='index', columns=['growth', 'status'])


def single_reaction_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                             method: str = 'fba',
                             reactions: List[Union[str, Reaction]] = None,
                             objective: Union[str, Dict[str, Union[float, int]]] = None,
                             constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None,
                             to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], pd.DataFrame]:
    """
    Single reaction deletion analysis of a metabolic model.
    Single reaction deletion analysis is a method to determine the effect of deleting each reaction in a metabolic model.
    It can be used to identify the reactions that are essential for the growth of a cell.
    In MEWpy, single reaction deletion analysis is performed by solving a linear problem for each reaction in the model.
    The methods FBA, MILP-FBA or pFBA can determine if the reaction is essential for the growth of a cell.

    :param model: a metabolic model to be simulated
    :param method: the method to be used for the simulation (default: 'fba').
    Available methods: 'fba', 'milp', 'pfba'
    :param reactions: the reactions to be simulated (default: all reactions in the model)
    :param objective: the objective function to be used for the simulation (default: the default objective)
    :param constraints: additional constraints to be used for the simulation (default: None)
    :param to_dict: whether to return the results as a dictionary (default: False)
    :return: a dictionary or a pandas DataFrame with the fluxes for each reaction
    """
    fba, reactions, _, constraints = _inputs_processing(model=model,
                                                        method=method,
                                                        fraction=None,
                                                        reactions=reactions,
                                                        genes=None,
                                                        objective=objective,
                                                        constraints=constraints)

    res = {}

    for rxn in reactions:

        if isinstance(rxn, str):

            rxn_id = rxn

        else:
            rxn_id = rxn.id

        constraints[rxn_id] = (0.0, 0.0)

        sol, status = _fba(method=fba, objective=None, minimize=False, constraints=constraints)

        res[rxn_id] = [sol, status]

        del constraints[rxn_id]

    if to_dict:
        return res

    return pd.DataFrame.from_dict(data=res, orient='index', columns=['growth', 'status'])
