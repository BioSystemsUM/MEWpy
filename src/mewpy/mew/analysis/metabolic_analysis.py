from collections import defaultdict
from typing import Union, TYPE_CHECKING, List, Dict, Tuple

# TODO: this module depends on pandas dataframes. Should it be set as package requirement?
# noinspection PyPackageRequirements
from pandas import DataFrame

from mewpy.util.constants import ModelConstants
from mewpy.mew.variables import Reaction, Gene
from .fba import FBA, pFBA, milpFBA
from .analysis_utils import decode_solver_solution

if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


# TODO: missing documentation and typing
def _fba(method, objective, minimize, constraints, status=False):

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

    return decode_solver_solution(solution=sol, minimize=minimize, status=status)


def slim_fba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
             fba: FBA = None,
             objective: Union[str, Dict[str, Union[float, int]]] = None,
             minimize: bool = False,
             constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None) -> Union[int, float, None]:

    if not fba:
        fba = FBA(model, build=True, attach=False)

    return _fba(method=fba,
                objective=objective,
                minimize=minimize,
                constraints=constraints)


def slim_pfba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
              pfba: pFBA = None,
              objective: Union[str, Dict[str, Union[float, int]]] = None,
              minimize: bool = False,
              constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None) -> Union[int, float, None]:
    if not pfba:
        pfba = pFBA(model, build=True, attach=False)

    return _fba(method=pfba,
                objective=objective,
                minimize=minimize,
                constraints=constraints)


def slim_milp_fba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                  milp_fba: milpFBA = None,
                  objective: Union[str, Dict[str, Union[float, int]]] = None,
                  minimize: bool = False,
                  constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None) -> Union[int,
                                                                                                       float,
                                                                                                       None]:
    if not milp_fba:
        milp_fba = pFBA(model, build=True, attach=False)

    return _fba(method=milp_fba,
                objective=objective,
                minimize=minimize,
                constraints=constraints)


def _inputs_processing(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                       method: str = None,
                       fraction: float = None,
                       reactions: List[Union[str, Reaction]] = None,
                       genes: List[Union[str, Gene]] = None,
                       objective: Union[str, Dict[str, Union[float, int]]] = None,
                       constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None):
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
            sol = _fba(method=fba,
                       objective=objective,
                       minimize=False,
                       constraints=constraints)

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
        to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], DataFrame]:
    """

    Flux Variability Analysis (FVA)

    :param model:
    :param method:
    :param fraction:
    :param reactions:
    :param objective:
    :param constraints:
    :param to_dict:

    :return:
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

        sol = _fba(method=fba,
                   objective=rxn_objective,
                   minimize=True,
                   constraints=constraints)

        res[rxn_objective].append(sol)

        sol = _fba(method=fba,
                   objective=rxn_objective,
                   minimize=False,
                   constraints=constraints)

        res[rxn_objective].append(sol)

    if to_dict:
        return res

    return DataFrame.from_dict(data=res, orient='index', columns=['minimum', 'maximum'])


def _milp_gene_deletion(genes, constraints, fba):
    res = {}

    for gene in genes:

        if isinstance(gene, str):

            gene_id = gene

        else:
            gene_id = gene.id

        constraints[gene_id] = (0.0, 0.0)

        sol, status = _fba(method=fba,
                           objective=None,
                           minimize=False,
                           constraints=constraints,
                           status=True)

        res[gene_id] = [sol, status]

        del constraints[gene_id]

    return res


def _fba_gene_deletion(model, genes, constraints, fba):
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

        sol, status = _fba(method=fba,
                           objective=None,
                           minimize=False,
                           constraints=constraints,
                           status=True)

        res[gene_obj.id] = [sol, status]

        for rxn in to_remove:
            del constraints[rxn]

    return res


def single_gene_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                         method: str = 'fba',
                         genes: List[Union[str, Gene]] = None,
                         objective: Union[str, Dict[str, Union[float, int]]] = None,
                         constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None,
                         to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], DataFrame]:
    """

    Single gene deletions

    :param model:
    :param method:
    :param genes:
    :param objective:
    :param constraints:
    :param to_dict:

    :return:
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

    return DataFrame.from_dict(data=res, orient='index', columns=['growth', 'status'])


def single_reaction_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                             method: str = 'fba',
                             reactions: List[Union[str, Reaction]] = None,
                             objective: Union[str, Dict[str, Union[float, int]]] = None,
                             constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None,
                             to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], DataFrame]:
    """

    Single reaction deletions

    :param model:
    :param method:
    :param reactions:
    :param objective:
    :param constraints:
    :param to_dict:

    :return:
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

        sol, status = _fba(method=fba,
                           objective=None,
                           minimize=False,
                           constraints=constraints,
                           status=True)

        res[rxn_id] = [sol, status]

        del constraints[rxn_id]

    if to_dict:
        return res

    return DataFrame.from_dict(data=res, orient='index', columns=['growth', 'status'])
