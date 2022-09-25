from collections import defaultdict
from typing import Union, TYPE_CHECKING, Dict, Tuple, Optional, Sequence

import pandas as pd

from mewpy.util.constants import ModelConstants
from .analysis_utils import run_method_and_decode
from .coregflux import CoRegFlux
from .prom import PROM
from .rfba import RFBA
from .srfba import SRFBA

if TYPE_CHECKING:
    from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel

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
        - A linear problem based on mass balance constraints (reactions/metabolites)
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


def slim_prom(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
              initial_state: Dict[Tuple[str, str], float] = None,
              regulator: str = None,
              objective: Union[str, Dict[str, float]] = None,
              constraints: Dict[str, Tuple[float, float]] = None) -> Optional[float]:
    """
    A Probabilistic Regulation of Metabolism (PROM) simulation of an integrated Metabolic-Regulatory model.
    A slim analysis produces a single and simple solution for the model. This method returns the objective value of
    the PROM simulation.

    Fundamentals of the PROM procedure:
        - A linear problem based on mass balance constraints (reactions/metabolites)
        - Reactions are further constrained using the regulator-target probabilities. The probabilities are
        obtained from the regulatory network. Namely, the probability of a gene being active when the its regulator
        is knocked out is used to constrain the associated reaction.
        - A FBA simulation is performed to obtain the flux distribution.

    :param model: an integrated metabolic-regulatory model to be simulated
    :param initial_state: the initial state of the model. If None, the default initial state is used
    :param regulator: the regulator to be knocked out. If None, the first regulator is knocked out
    :param objective: the objective of the simulation. If None, the default objective is used
    :param constraints: the constraints of the model. If None, the default constraints are used
    :return: the objective value of the simulation
    """
    if not regulator:
        regulator = next(iter(model.regulators.keys()))

    prom = PROM(model).build()

    objective_value, _ = run_method_and_decode(method=prom, objective=objective, constraints=constraints,
                                               initial_state=initial_state, regulators=regulator)
    return objective_value


def slim_coregflux(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                   initial_state: Dict[str, float] = None,
                   objective: Union[str, Dict[str, float]] = None,
                   constraints: Dict[str, Tuple[float, float]] = None) -> Optional[float]:
    """
    A CoRegFlux simulation of an integrated Metabolic-Regulatory model.
    A slim analysis produces a single and simple solution for the model. This method returns the objective value of
    the CoRegFlux simulation.

    Fundamentals of the CoRegFlux procedure:
         - A linear problem based on mass balance constraints (reactions/metabolites)
         - A linear regression estimator predicts the expression of target genes
         as function of the co-expression of regulators
         - the predicted expression of the target genes to constrain the bounds of the associated reactions
         - A FBA simulation is performed to obtain the flux distribution.

    :param model: an integrated metabolic-regulatory model to be simulated
    :param initial_state: the initial state of the model. If None, the default initial state is used
    :param objective: the objective of the simulation. If None, the default objective is used
    :param constraints: the constraints of the model. If None, the default constraints are used
    :return: the objective value of the simulation
    """
    co_reg_flux = CoRegFlux(model).build()
    objective_value, _ = run_method_and_decode(method=co_reg_flux, objective=objective, constraints=constraints,
                                               initial_state=initial_state)
    return objective_value


def ifva(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
         fraction: float = 1.0,
         reactions: Sequence[str] = None,
         objective: Union[str, Dict[str, float]] = None,
         constraints: Dict[str, Tuple[float, float]] = None,
         initial_state: Dict[str, float] = None,
         method: str = 'srfba') -> pd.DataFrame:
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
    :return: a DataFrame with the results of the simulation
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

    return pd.DataFrame.from_dict(data=result, orient='index', columns=['minimum', 'maximum'])


def isingle_gene_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                          genes: Sequence[str] = None,
                          constraints: Dict[str, Tuple[float, float]] = None,
                          initial_state: Dict[str, float] = None,
                          method: str = 'srfba') -> pd.DataFrame:
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
    :return: a DataFrame with the results of the simulation
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

    return pd.DataFrame.from_dict(data=result, orient='index', columns=['growth', 'status'])


def isingle_reaction_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                              reactions: Sequence[str] = None,
                              constraints: Dict[str, Tuple[float, float]] = None,
                              initial_state: Dict[str, float] = None,
                              method: str = 'srfba') -> pd.DataFrame:
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
    :return: a DataFrame with the results of the simulation
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

    return pd.DataFrame.from_dict(data=result, orient='index', columns=['growth', 'status'])


def isingle_regulator_deletion(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                               regulators: Sequence[str] = None,
                               constraints: Dict[str, Tuple[float, float]] = None,
                               initial_state: Dict[str, float] = None,
                               method: str = 'srfba') -> pd.DataFrame:
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
    :return: a DataFrame with the results of the simulation
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

    return pd.DataFrame.from_dict(data=result, orient='index', columns=['growth', 'status'])


def _decode_initial_state(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                          state: Dict[str, float]) -> Dict[str, float]:
    """
    Method responsible for retrieving the initial state of the model.
    The initial state is the state of all regulators found in the Metabolic-Regulatory model.
    :param model: the model to be simulated
    :param state: the initial state of the model
    :return: dict of regulatory/metabolic variable keys (regulators) and a value of 0 or 1
    """
    if not state:
        state = {}

    initial_state = {}
    for regulator in model.yield_regulators():
        if regulator.id in state:
            initial_state[regulator.id] = state[regulator.id]

        elif regulator.is_metabolite() and regulator.exchange_reaction:
            if regulator.exchange_reaction.id in state:
                initial_state[regulator.id] = state[regulator.exchange_reaction.id]

            else:
                initial_state[regulator.id] = abs(regulator.exchange_reaction.lower_bound)

        else:
            initial_state[regulator.id] = max(regulator.coefficients)

    return initial_state


def _decode_interactions(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                         state: Dict[str, float]) -> Dict[str, float]:
    """
    Decodes the state of the model to a dictionary with the state of each gene.
    :param model: the model to be decoded
    :param state: the state of the model
    :return: a dictionary with the state of each gene
    """
    target_state = {}
    for interaction in model.yield_interactions():

        for coefficient, event in interaction.regulatory_events.items():
            if event.is_none:
                continue

            result = event.evaluate(state)
            if result:
                target_state[interaction.target.id] = coefficient

            else:
                target_state[interaction.target.id] = 0.0

    return target_state


def _decode_gprs(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                 state: Dict[str, float]) -> Dict[str, float]:
    """
    Decodes the state of the model to a dictionary with the state of each reaction.
    :param model: the model to be decoded
    :param state: the state of the model
    :return: a dictionary with the state of each gene
    """
    reaction_state = {}
    for reaction in model.yield_reactions():
        if reaction.gpr.is_none:
            continue

        result = reaction.gpr.evaluate(state)
        if result:
            reaction_state[reaction.id] = 1.0

        else:
            reaction_state[reaction.id] = 0.0

    return reaction_state


def find_conflicts(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                   strategy: str = 'two-step',
                   constraints: Dict[str, Tuple[float, float]] = None,
                   initial_state: Dict[str, float] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    It finds conflicts between the regulatory and metabolic states of an integrated MEW model.
    Setting up the regulators' initial state in integrated models is a difficult task. Most of the time,
    the initial state is not known and hinders feasible solutions during simulation.
    This method can be used to find regulatory states that affect the growth of the cell.

    One can use two strategies to find conflicts:
        - One-step strategy: it simulates the model with the initial state of the regulators and
        checks if the metabolic state violates essential genes and reactions.
        - Two-step strategy: it simulates the model with the initial state of the regulators to retrieve the regulatory
        state of all regulators. Then, it simulates the model with the regulatory state of the regulators to retrieve
        the metabolic state. Finally, it checks if the metabolic state violates essential genes and reactions.

    NOTE: Although this method finds direct conflicts between regulators and genes affecting cellular growth,
    it might not detect indirect conflicts.
    For example, an active regulator can repress an essential metabolic gene and thus being reported as a conflict.
    However, this regulator can only be active due to an erroneous initial state of another regulator or environmental
    stimuli, metabolite or reaction. That is, this active regulator is also a target being affected by other regulators.
    In these cases, these regulators should be reported as conflicts as well.
    Hence, this method should be run several times with different initial states
    of the regulators to find all conflicts!!!

    :param model: a metabolic-regulatory model to be simulated
    :param strategy: the strategy to be used to find conflicts. Available strategies: 'one-step', 'two-step'.
    one-step performs a single regulatory simulation. two-step performs two simulations. Default: 'two-step'
    :param constraints: additional constraints to be added to the model. If None, no additional constraints are added
    :param initial_state: the initial state of the model. If None, the default initial state is used (default: None)
    :return: a DataFrame with the conflicts between essential genes and regulatory states, and a DataFrame with the
    conflicts between essential reactions and metabolic states
    """
    if not constraints:
        constraints = {}
    else:
        constraints = constraints.copy()

    if not initial_state:
        initial_state = {}
    else:
        initial_state = initial_state.copy()

    # 1. it performs a FBA simulation to find the optimal growth rate
    from mewpy.mew.analysis import FBA
    solution = FBA(model).build().optimize(solver_kwargs={'constraints': constraints})

    if not solution.objective_value:
        raise RuntimeError('FBA solution is not feasible (objective value is 0). To find inconsistencies, '
                           'the metabolic model must be feasible.')

    # 2. it performs an essential genes analysis using FBA
    from mewpy.mew.analysis.metabolic_analysis import single_gene_deletion
    gene_deletion = single_gene_deletion(model, constraints=constraints)
    essential_genes = gene_deletion[gene_deletion['growth'] < ModelConstants.TOLERANCE]

    from mewpy.mew.analysis.metabolic_analysis import single_reaction_deletion
    reaction_deletion = single_reaction_deletion(model, constraints=constraints)
    essential_reactions = reaction_deletion[reaction_deletion['growth'] < ModelConstants.TOLERANCE]

    state = _decode_initial_state(model, initial_state)

    # noinspection PyTypeChecker
    regulatory_state = _decode_interactions(model, state)

    if strategy == 'two-step':
        regulatory_state = {gene: value for gene, value in regulatory_state.items() if model.get(gene).is_regulator()}
        regulatory_state = {**state, **regulatory_state}
        # noinspection PyTypeChecker
        metabolic_state = _decode_interactions(model, regulatory_state)

    else:
        metabolic_state = regulatory_state.copy()

    deleted_genes = set(gene for gene, val in metabolic_state.items() if val < ModelConstants.TOLERANCE)
    gene_conflicts = set.intersection(deleted_genes, set(essential_genes.index))

    reaction_state = _decode_gprs(model, metabolic_state)
    deleted_reactions = set(reaction for reaction, val in reaction_state.items() if val < ModelConstants.TOLERANCE)
    reaction_conflicts = set.intersection(deleted_reactions, set(essential_reactions.index))

    repressed_genes = []
    for gene in gene_conflicts:
        gene = model.genes[gene]

        if gene.is_target():
            # noinspection PyUnresolvedReferences
            repressed_gene = {regulator: regulatory_state[regulator] for regulator in gene.regulators}
            # noinspection PyUnresolvedReferences
            repressed_gene['interaction'] = str(gene.interaction)
        else:
            repressed_gene = {}

        df = pd.DataFrame(repressed_gene, index=[gene.id])
        repressed_genes.append(df)

    if repressed_genes:
        repressed_genes = pd.concat(repressed_genes)
        cols = ['interaction'] + [col for col in repressed_genes.columns if col != 'interaction']
        repressed_genes = repressed_genes[cols]
    else:
        repressed_genes = pd.DataFrame()

    repressed_reactions = []
    for reaction in reaction_conflicts:
        reaction = model.reactions[reaction]

        repressed_reaction = {gene: metabolic_state[gene] for gene in reaction.genes}
        repressed_reaction['gpr'] = reaction.gene_protein_reaction_rule

        df = pd.DataFrame(repressed_reaction, index=[reaction.id])
        repressed_reactions.append(df)

    if repressed_reactions:
        repressed_reactions = pd.concat(repressed_reactions)
        cols = ['gpr'] + [col for col in repressed_reactions.columns if col != 'gpr']
        repressed_reactions = repressed_reactions[cols]
    else:
        repressed_reactions = pd.DataFrame()

    return repressed_genes, repressed_reactions
