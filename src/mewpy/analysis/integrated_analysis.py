from collections import defaultdict
from typing import Union, TYPE_CHECKING, List, Dict, Tuple

# TODO: this module depends on pandas dataframes. Should it be set as package requirement?
from pandas import DataFrame

from mewpy.util.constants import ModelConstants
from mewpy.variables import Regulator, Reaction, Gene
from .rfba import RFBA
from .srfba import SRFBA
from .analysis_utils import decode_solver_solution

if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


# TODO: missing documentation and typing
def _rfba(method, objective, minimize, initial_state, constraints, status=False):

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

    return decode_solver_solution(solution=sol, minimize=minimize, status=status)


def slim_rfba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
              rfba: RFBA = None,
              objective: Union[str, Dict[str, Union[float, int]]] = None,
              minimize: bool = False,
              initial_state: Dict[str, Union[float, int]] = None,
              constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None) -> Union[int, float, None]:

    if not rfba:
        rfba = RFBA(model, build=True, attach=False)

    return _rfba(method=rfba,
                 objective=objective,
                 minimize=minimize,
                 initial_state=initial_state,
                 constraints=constraints)


def slim_srfba(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
               srfba: SRFBA = None,
               objective: Union[str, Dict[str, Union[float, int]]] = None,
               minimize: bool = False,
               initial_state: Dict[str, Union[float, int]] = None,
               constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None) -> Union[int, float, None]:
    if not srfba:
        srfba = SRFBA(model, build=True, attach=False)

    return _rfba(method=srfba,
                 objective=objective,
                 minimize=minimize,
                 initial_state=initial_state,
                 constraints=constraints)


def _inputs_processing(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                       method: str = None,
                       fraction: float = None,
                       reactions: List[Union[str, Reaction]] = None,
                       genes: List[Union[str, Gene]] = None,
                       regulators: List[Union[str, Regulator]] = None,
                       objective: Union[str, Dict[str, Union[float, int]]] = None,
                       initial_state: Dict[str, Union[float, int]] = None,
                       constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None):
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
            sol = _rfba(method=rfba,
                        objective=objective,
                        minimize=False,
                        initial_state=initial_state,
                        constraints=constraints,
                        status=False)

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
         to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], DataFrame]:
    """

    Integrated Flux Variability Analysis (iFVA)

    :param model:
    :param method:
    :param fraction:
    :param reactions:
    :param objective:
    :param initial_state:
    :param constraints:
    :param to_dict:

    :return:
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

        sol = _rfba(method=rfba,
                    objective=rxn_objective,
                    minimize=True,
                    initial_state=initial_state,
                    constraints=constraints)

        res[rxn_objective].append(sol)

        sol = _rfba(method=rfba,
                    objective=rxn_objective,
                    minimize=False,
                    initial_state=initial_state,
                    constraints=constraints)

        res[rxn_objective].append(sol)

    if to_dict:
        return res

    return DataFrame.from_dict(data=res, orient='index', columns=['minimum', 'maximum'])


def isingle_gene_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                          method: str = 'srfba',
                          genes: List[Union[str, Gene]] = None,
                          objective: Union[str, Dict[str, Union[float, int]]] = None,
                          initial_state: Dict[str, Union[float, int]] = None,
                          constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None,
                          to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], DataFrame]:
    """

    Integrated single gene deletions

    :param model:
    :param method:
    :param genes:
    :param objective:
    :param initial_state:
    :param constraints:
    :param to_dict:

    :return:
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

        sol, status = _rfba(method=rfba,
                            objective=None,
                            minimize=False,
                            initial_state=initial_state,
                            constraints=constraints,
                            status=True)

        res[gene_id] = [sol, status]

        del initial_state[gene_id]

    if to_dict:
        return res

    return DataFrame.from_dict(data=res, orient='index', columns=['growth', 'status'])


def isingle_reaction_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                              method: str = 'srfba',
                              reactions: List[Union[str, Reaction]] = None,
                              objective: Union[str, Dict[str, Union[float, int]]] = None,
                              initial_state: Dict[str, Union[float, int]] = None,
                              constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None,
                              to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], DataFrame]:
    """

    Integrated single reaction deletions

    :param model:
    :param method:
    :param reactions:
    :param objective:
    :param initial_state:
    :param constraints:
    :param to_dict:

    :return:
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

        sol, status = _rfba(method=rfba,
                            objective=None,
                            minimize=False,
                            initial_state=initial_state,
                            constraints=constraints,
                            status=True)

        res[rxn_id] = [sol, status]

        del constraints[rxn_id]

    if to_dict:
        return res

    return DataFrame.from_dict(data=res, orient='index', columns=['growth', 'status'])


def isingle_regulator_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                               method: str = 'srfba',
                               regulators: List[Union[str, Regulator]] = None,
                               objective: Union[str, Dict[str, Union[float, int]]] = None,
                               initial_state: Dict[str, Union[float, int]] = None,
                               constraints: Dict[str, Tuple[Union[float, int], Union[float, int]]] = None,
                               to_dict: bool = False) -> Union[Dict[str, List[Union[float, int]]], DataFrame]:
    """

    Integrated single regulator deletions

    :param model:
    :param method:
    :param regulators:
    :param objective:
    :param initial_state:
    :param constraints:
    :param to_dict:

    :return:
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

        sol, status = _rfba(method=rfba,
                            objective=None,
                            minimize=False,
                            initial_state=initial_state,
                            constraints=constraints,
                            status=True)

        res[regulator_id] = [sol, status]

        del initial_state[regulator_id]

    if to_dict:
        return res

    return DataFrame.from_dict(data=res, orient='index', columns=['growth', 'status'])
