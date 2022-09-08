from collections import defaultdict
from typing import Union, TYPE_CHECKING, List, Dict, Tuple

import pandas as pd

from mewpy.util.constants import ModelConstants
from mewpy.mew.variables import Regulator, Reaction, Gene
from .rfba import RFBA
from .srfba import SRFBA
from .analysis_utils import decode_solver_solution

if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


def _rfba(method, objective, minimize, initial_state, constraints):
    # noinspection PyProtectedMember
    old_linear_obj = method._linear_objective.copy()
    # noinspection PyProtectedMember
    old_quadratic_obj = method._quadratic_objective.copy()
    old_sense = method.minimize

    sol = method.optimize(objective=objective,
                          minimize=minimize,
                          initial_state=initial_state,
                          constraints=constraints,
                          to_solver=True,
                          get_values=False)

    method.set_objective(linear=old_linear_obj, quadratic=old_quadratic_obj, minimize=old_sense)

    return decode_solver_solution(solution=sol, minimize=minimize)


def slim_rfba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
              rfba: RFBA = None,
              objective: Union[str, Dict[str, Union[float, int]]] = None,
              minimize: bool = False,
              initial_state: Dict[str, Union[float, int]] = None,
              constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None) -> Union[int, float, None]:
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
    :param rfba: a RFBA object to be used for the simulation. If None, a new RFBA object is created and used
    :param objective: the objective of the simulation. If None, the default objective is used
    :param minimize: if True, the objective is minimized. If False, the objective is maximized
    :param initial_state: the initial state of the model. If None, the default initial state is used
    :param constraints: the constraints of the model. If None, the default constraints are used
    :return: the objective value of the simulation
    """
    if not rfba:
        rfba = RFBA(model, build=True, attach=False)

    sol, _ = _rfba(method=rfba, objective=objective, minimize=minimize, initial_state=initial_state,
                   constraints=constraints)
    return sol


def slim_srfba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
               srfba: SRFBA = None,
               objective: Union[str, Dict[str, Union[float, int]]] = None,
               minimize: bool = False,
               initial_state: Dict[str, Union[float, int]] = None,
               constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None) -> Union[int, float, None]:
    """
    A Synchronous Regulatory Flux Balance Analysis (SRFBA) simulation of an integrated Metabolic-Regulatory model.
    A slim analysis produces a single and simple solution for the model. This method returns the objective value of
    the SRFBA simulation.

    Fundamentals of the SRFBA procedure:
        - A linear problem based on reactions, GPRs, and regulatory interactions using mixed-integer constraints
        - Reactions, genes, and regulators are constrained by the bounds or coefficients
        - SRFBA is solved using a mixed-integer solver

    :param model: an integrated metabolic-regulatory model to be simulated
    :param srfba: a SRFBA object to be used for the simulation. If None, a new SRFBA object is created and used
    :param objective: the objective of the simulation. If None, the default objective is used
    :param minimize: if True, the objective is minimized. If False, the objective is maximized
    :param initial_state: the initial state of the model. If None, the default initial state is used
    :param constraints: the constraints of the model. If None, the default constraints are used
    :return: the objective value of the simulation
    """
    if not srfba:
        srfba = SRFBA(model, build=True, attach=False)

    sol, _ = _rfba(method=srfba, objective=objective, minimize=minimize, initial_state=initial_state,
                   constraints=constraints)
    return sol


def _inputs_processing(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                       method: str = None,
                       fraction: float = None,
                       reactions: List[Union[str, Reaction]] = None,
                       genes: List[Union[str, Gene]] = None,
                       regulators: List[Union[str, Regulator]] = None,
                       objective: Union[str, Dict[str, Union[float, int]]] = None,
                       initial_state: Dict[str, Union[float, int]] = None,
                       constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None):
    """
    Processes the inputs of the analysis methods.
    :param model: an integrated metabolic-regulatory model to be simulated
    :param method: the method to be used for the simulation
    :param fraction: the fraction of the objective to be used for the simulation
    :param reactions: the reactions to be simulated
    :param genes: the genes to be simulated
    :param regulators: the regulators to be simulated
    :param objective: the objective of the simulation. If None, the default objective is used
    :param initial_state: the initial state of the model. If None, the default initial state is used
    :param constraints: the constraints of the model. If None, the default constraints are used
    :return: processed inputs
    """
    methods_map = {'rfba': RFBA,
                   'srfba': SRFBA}

    if not constraints:
        constraints = {}

    if not initial_state:
        initial_state = {}

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

    _RFBA = methods_map.get(method, None)

    if _RFBA is None:
        raise ValueError(f'Simulation method {method} is not valid. Available methods: {methods_map}')

    rfba = _RFBA(model, build=True, attach=False)
    rfba.set_objective(linear=objective, minimize=False)

    if fraction is not None:

        if fraction > 0.0:
            sol, _ = _rfba(method=rfba, objective=objective, minimize=False, initial_state=initial_state,
                           constraints=constraints)

            if sol is None:
                sol = 0.0

            constraints[obj_rxn] = (fraction * sol, ModelConstants.REACTION_UPPER_BOUND)

    if not reactions:
        reactions = model.reactions.keys()

    if not genes:
        genes = model.genes.keys()

    if not regulators:
        regulators = model.regulators.keys()

    return rfba, reactions, genes, regulators, initial_state, constraints


def ifva(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
         method: str = 'srfba',
         fraction: float = 1.0,
         reactions: List[Union[str, Reaction]] = None,
         objective: Union[str, Dict[str, Union[float, int]]] = None,
         initial_state: Dict[str, Union[float, int]] = None,
         constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None,
         to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], pd.DataFrame]:
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
    :param method: the method to be used for the simulation. Available methods: 'rfba', 'srfba'. Default: 'srfba'
    :param fraction: the fraction of the optimal solution to be used as the upper bound for the objective function
    (default: 1.0)
    :param reactions: the reactions to be simulated. If None, all reactions are simulated (default: None)
    :param objective: the objective of the simulation. If None, the default objective is used (default: None)
    :param initial_state: the initial state of the model. If None, the default initial state is used (default: None)
    :param constraints: additional constraints to be added to the model. If None, no additional constraints are added
    :param to_dict: if True, the results are returned as a dictionary. If False, the results are returned as a DataFrame
    :return: a dictionary or a DataFrame with the results of the simulation
    """
    rfba, reactions, _, _, initial_state, constraints = _inputs_processing(model=model,
                                                                           method=method,
                                                                           fraction=fraction,
                                                                           reactions=reactions,
                                                                           genes=None,
                                                                           regulators=None,
                                                                           objective=objective,
                                                                           initial_state=initial_state,
                                                                           constraints=constraints)

    res = defaultdict(list)

    for rxn in reactions:

        if isinstance(rxn, str):

            rxn_objective = rxn

        else:
            rxn_objective = rxn.id

        sol, _ = _rfba(method=rfba, objective=rxn_objective, minimize=True, initial_state=initial_state,
                       constraints=constraints)

        res[rxn_objective].append(sol)

        sol, _ = _rfba(method=rfba, objective=rxn_objective, minimize=False, initial_state=initial_state,
                       constraints=constraints)

        res[rxn_objective].append(sol)

    if to_dict:
        return res

    return pd.DataFrame.from_dict(data=res, orient='index', columns=['minimum', 'maximum'])


def isingle_gene_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                          method: str = 'srfba',
                          genes: List[Union[str, Gene]] = None,
                          objective: Union[str, Dict[str, Union[float, int]]] = None,
                          initial_state: Dict[str, Union[float, int]] = None,
                          constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None,
                          to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], pd.DataFrame]:
    """
    Integrated single gene deletion analysis of an integrated Metabolic-Regulatory model.
    Integrated single gene deletion analysis is a method to determine the effect of deleting each gene
    in a model.
    It can be used to identify the genes that are essential for the growth of a cell.
    In MEWpy, single gene deletion analysis is performed by solving a linear problem for each gene in the model.
    A gene knockout can switch off reactions associated with the gene, only if the gene is essential for the reaction.
    The methods RFBA and SRFBA can determine if the gene is essential for the growth of a cell.

    :param model: an integrated metabolic-regulatory model to be simulated
    :param method: the method to be used for the simulation. Available methods: 'rfba', 'srfba'. Default: 'srfba'
    :param genes: the genes to be simulated. If None, all genes are simulated (default: None)
    :param objective: the objective of the simulation. If None, the default objective is used (default: None)
    :param initial_state: the initial state of the model. If None, the default initial state is used (default: None)
    :param constraints: additional constraints to be added to the model. If None, no additional constraints are added
    :param to_dict: if True, the results are returned as a dictionary. If False, the results are returned as a DataFrame
    :return: a dictionary or a DataFrame with the results of the simulation
    """
    rfba, _, genes, _, initial_state, constraints = _inputs_processing(model=model,
                                                                       method=method,
                                                                       fraction=None,
                                                                       reactions=None,
                                                                       genes=genes,
                                                                       regulators=None,
                                                                       objective=objective,
                                                                       initial_state=initial_state,
                                                                       constraints=constraints)

    res = {}

    for gene in genes:

        if isinstance(gene, str):

            gene_id = gene

        else:
            gene_id = gene.id

        initial_state[gene_id] = 0.0

        sol, status = _rfba(method=rfba, objective=None, minimize=False, initial_state=initial_state,
                            constraints=constraints)

        res[gene_id] = [sol, status]

        del initial_state[gene_id]

    if to_dict:
        return res

    return pd.DataFrame.from_dict(data=res, orient='index', columns=['growth', 'status'])


def isingle_reaction_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                              method: str = 'srfba',
                              reactions: List[Union[str, Reaction]] = None,
                              objective: Union[str, Dict[str, Union[float, int]]] = None,
                              initial_state: Dict[str, Union[float, int]] = None,
                              constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None,
                              to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], pd.DataFrame]:
    """
    Integrated single reaction deletion analysis of an integrated Metabolic-Regulatory model.
    Integrated single reaction deletion analysis is a method to determine the effect of deleting each reaction
    in a model.
    It can be used to identify the reactions that are essential for the growth of a cell.
    In MEWpy, single reaction deletion analysis is performed by solving a linear problem for each reaction in the model.
    The methods RFBA and SRFBA can determine if the reaction is essential for the growth of a cell.

    :param model: an integrated metabolic-regulatory model to be simulated
    :param method: the method to be used for the simulation. Available methods: 'rfba', 'srfba'. Default: 'srfba'
    :param reactions: the reactions to be simulated. If None, all reactions are simulated (default: None)
    :param objective: the objective of the simulation. If None, the default objective is used (default: None)
    :param initial_state: the initial state of the model. If None, the default initial state is used (default: None)
    :param constraints: additional constraints to be added to the model. If None, no additional constraints are added
    :param to_dict: if True, the results are returned as a dictionary. If False, the results are returned as a DataFrame
    :return: a dictionary or a DataFrame with the results of the simulation
    """

    rfba, reactions, _, _, initial_state, constraints = _inputs_processing(model=model,
                                                                           method=method,
                                                                           fraction=None,
                                                                           reactions=reactions,
                                                                           genes=None,
                                                                           regulators=None,
                                                                           objective=objective,
                                                                           initial_state=initial_state,
                                                                           constraints=constraints)

    res = {}

    for rxn in reactions:

        if isinstance(rxn, str):

            rxn_id = rxn

        else:
            rxn_id = rxn.id

        constraints[rxn_id] = (0.0, 0.0)

        sol, status = _rfba(method=rfba, objective=None, minimize=False, initial_state=initial_state,
                            constraints=constraints)

        res[rxn_id] = [sol, status]

        del constraints[rxn_id]

    if to_dict:
        return res

    return pd.DataFrame.from_dict(data=res, orient='index', columns=['growth', 'status'])


def isingle_regulator_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                               method: str = 'srfba',
                               regulators: List[Union[str, Regulator]] = None,
                               objective: Union[str, Dict[str, Union[float, int]]] = None,
                               initial_state: Dict[str, Union[float, int]] = None,
                               constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None,
                               to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], pd.DataFrame]:
    """
    Integrated single regulator deletion analysis of a regulatory model.
    Integrated single regulator deletion analysis is a method to determine the effect of deleting each regulator
    in a regulatory model.
    It can be used to identify the regulators that are essential for the growth of a cell.
    In MEWpy, single regulator deletion analysis is performed by solving a linear problem for each regulator in the model.
    A regulator knockout can switch off reactions associated with the regulator,
    only if the regulator is essential for the interactions and gprs.
    The methods RFBA and SRFBA can determine if the regulator is essential for the growth of a cell.

    :param model: an integrated metabolic-regulatory model to be simulated
    :param method: the method to be used for the simulation. Available methods: 'rfba', 'srfba'. Default: 'srfba'
    :param regulators: the regulators to be simulated. If None, all regulators are simulated (default: None)
    :param objective: the objective of the simulation. If None, the default objective is used (default: None)
    :param initial_state: the initial state of the model. If None, the default initial state is used (default: None)
    :param constraints: additional constraints to be added to the model. If None, no additional constraints are added
    :param to_dict: if True, the results are returned as a dictionary. If False, the results are returned as a DataFrame
    :return: a dictionary or a DataFrame with the results of the simulation
    """
    rfba, _, _, regulators, initial_state, constraints = _inputs_processing(model=model,
                                                                            method=method,
                                                                            fraction=None,
                                                                            reactions=None,
                                                                            genes=None,
                                                                            regulators=regulators,
                                                                            objective=objective,
                                                                            initial_state=initial_state,
                                                                            constraints=constraints)

    res = {}

    for regulator in regulators:

        if isinstance(regulator, str):

            regulator_id = regulator

        else:
            regulator_id = regulator.id

        initial_state[regulator_id] = 0.0

        sol, status = _rfba(method=rfba, objective=None, minimize=False, initial_state=initial_state,
                            constraints=constraints)

        res[regulator_id] = [sol, status]

        del initial_state[regulator_id]

    if to_dict:
        return res

    return pd.DataFrame.from_dict(data=res, orient='index', columns=['growth', 'status'])
