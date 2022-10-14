from typing import TYPE_CHECKING, Union, Dict, Sequence, List, Tuple

import numpy as np
import pandas as pd

from mewpy.germ.analysis import FBA
from mewpy.germ.analysis.analysis_utils import biomass_yield_to_rate, \
    CoRegMetabolite, CoRegBiomass, metabolites_constraints, gene_state_constraints, system_state_update, \
    build_metabolites, build_biomass, CoRegResult
from mewpy.germ.solution import ModelSolution, DynamicSolution
from mewpy.germ.variables import Gene, Target
from mewpy.solvers.solution import Solution, Status
from mewpy.solvers.solver import Solver
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel


def _run_and_decode(lp, additional_constraints=None, solver_kwargs=None):
    if not solver_kwargs:
        solver_kwargs = {}

    if additional_constraints:
        solver_kwargs['constraints'] = {**solver_kwargs.get('constraints', {}), **additional_constraints}

    solution = lp.solver.solve(**solver_kwargs)

    if not solution.values:
        return {rxn: 0 for rxn in lp.model.reactions}, 0

    return solution.values, solution.fobj


def result_to_solution(result: CoRegResult, model: 'Model', to_solver: bool = False) -> Union[ModelSolution, Solution]:
    """
    It converts a CoRegResult object to a ModelSolution object.

    :param result: the CoRegResult object
    :param model: the model
    :param to_solver: if True, it returns a Solution object
    :return: the ModelSolution object
    """
    if to_solver:
        return Solution(status=Status.OPTIMAL, fobj=result.objective_value, values=result.values)

    solution = ModelSolution(method='CoRegFlux',
                             x=result.values,
                             objective_value=result.objective_value,
                             status='optimal',
                             model=model)

    solution.metabolites = {key: met.concentration for key, met in result.metabolites.items()}
    solution.biomass = result.biomass.biomass_yield
    solution.constraints = result.constraints
    return solution


class CoRegFlux(FBA):

    def __init__(self,
                 model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                 solver: Union[str, Solver] = None,
                 build: bool = False,
                 attach: bool = False):
        """
        CoRegFlux is aimed at integrating reverse engineered transcriptional regulatory networks and gene-expression
        into metabolic models to improve prediction of phenotypes.
        It builds a linear regression estimator to predict the expression of target genes
        as function of the co-expression of regulators. The influence score of the regulators is used as input
        for the linear regression model.
        Then, it uses the predicted expression of the target genes to constrain the bounds of the associated reactions.
        It infers continuous constraints on the fluxes of the metabolic model.

        Author: Pauline Trébulle, Daniel Trejo-Banos, Mohamed Elati
        For more detail consult: https://dx.doi.org/10.1186%2Fs12918-017-0507-0
        :param model: a GERM model aka an integrated RegulatoryMetabolicModel
        :param solver: A Solver, CplexSolver, GurobiSolver or OptLangSolver instance.
        Alternatively, the name of the solver is also accepted.
        The solver interface will be used to load and solve a linear problem in a given solver.
        If none, a new solver is instantiated. An instantiated solver may be used,
        but it will be overwritten if build is true.
        :param build: Whether to build the linear problem upon instantiation. Default: False
        :param attach: Whether to attach the linear problem to the model upon instantiation. Default: False
        """
        super().__init__(model=model, solver=solver, build=build, attach=attach)

    # ---------------------------------
    # Dynamic simulation
    # ---------------------------------
    def next_state(self,
                   solver_kwargs: Dict = None,
                   state: Dict[str, float] = None,
                   metabolites: Dict[str, CoRegMetabolite] = None,
                   biomass: CoRegBiomass = None,
                   time_step: float = None,
                   soft_plus: float = 0,
                   tolerance: float = ModelConstants.TOLERANCE,
                   scale: bool = False) -> CoRegResult:
        """
        It computes the next state of the system given the current state and the time step.
        :param solver_kwargs: solver arguments
        :param state: current state of the system
        :param metabolites: metabolites constraints
        :param biomass: biomass constraints
        :param time_step: time step
        :param soft_plus: soft plus parameter
        :param tolerance: tolerance
        :param scale: whether to scale the metabolites
        :return: next state of the system
        """
        # Similar to the Simulation_step in the R implementation
        result = CoRegResult()

        constraints = {reaction.id: reaction.bounds for reaction in self.model.yield_reactions()}

        if metabolites:
            # updating coregflux constraints using metabolites concentrations
            constraints = metabolites_constraints(constraints=constraints,
                                                  metabolites=metabolites,
                                                  biomass=biomass,
                                                  time_step=time_step)

        if state:
            # updating coregflux bounds using gene state (predicted with predict_gene_state from gene expression and
            # regulator coregnet influence score)
            constraints = gene_state_constraints(model=self.model,
                                                 constraints=constraints,
                                                 state=state,
                                                 soft_plus=soft_plus,
                                                 tolerance=tolerance,
                                                 scale=scale)

        # retrieve the fba simulation from the inferred constraints
        values, objective_value = _run_and_decode(self, additional_constraints=constraints, solver_kwargs=solver_kwargs)
        result.values = values
        result.objective_value = objective_value

        # Updating the system state by solving an euler step for metabolites and
        # biomass given the previously fba solution, namely the flux state
        next_biomass, next_metabolites = system_state_update(model=self.model,
                                                             flux_state=values,
                                                             metabolites=metabolites,
                                                             biomass=biomass,
                                                             time_step=time_step,
                                                             biomass_fn=biomass_yield_to_rate)

        result.metabolites = next_metabolites
        result.biomass = next_biomass
        result.constraints = constraints
        return result

    def _dynamic_optimize(self,
                          to_solver: bool = False,
                          solver_kwargs: Dict = None,
                          initial_state: Sequence[Dict[str, float]] = None,
                          metabolites: Dict[str, CoRegMetabolite] = None,
                          biomass: CoRegBiomass = None,
                          time_steps: Sequence[float] = None,
                          soft_plus: float = 0,
                          tolerance: float = ModelConstants.TOLERANCE,
                          scale: bool = False) -> Union[DynamicSolution, Dict[float, Solution]]:
        solutions = []

        previous_time_step = 0
        for i_initial_state, time_step in zip(initial_state, time_steps):
            time_step_diff = time_step - previous_time_step

            next_state = self.next_state(solver_kwargs=solver_kwargs,
                                         state=i_initial_state,
                                         metabolites=metabolites,
                                         biomass=biomass,
                                         time_step=time_step_diff,
                                         soft_plus=soft_plus,
                                         tolerance=tolerance,
                                         scale=scale)

            previous_time_step = time_step

            metabolites = next_state.metabolites
            biomass = next_state.biomass

            solution = result_to_solution(result=next_state, model=self.model, to_solver=to_solver)
            solutions.append(solution)

        if to_solver:
            return dict(zip(time_steps, solutions))

        return DynamicSolution(*solutions, time=time_steps)

    def _steady_state_optimize(self,
                               to_solver: bool = False,
                               solver_kwargs: Dict = None,
                               initial_state: Dict[str, float] = None,
                               metabolites: Dict[str, CoRegMetabolite] = None,
                               biomass: CoRegBiomass = None,
                               time_step: float = 1,
                               soft_plus: float = 0,
                               tolerance: float = ModelConstants.TOLERANCE,
                               scale: bool = False) -> Union[ModelSolution, Solution]:

        result = self.next_state(solver_kwargs=solver_kwargs,
                                 state=initial_state,
                                 metabolites=metabolites,
                                 biomass=biomass,
                                 time_step=time_step,
                                 soft_plus=soft_plus,
                                 tolerance=tolerance,
                                 scale=scale)
        return result_to_solution(result=result, model=self.model, to_solver=to_solver)

    def _optimize(self,
                  to_solver: bool = False,
                  solver_kwargs: Dict = None,
                  initial_state: Union[Dict[str, float], Sequence[Dict[str, float]]] = None,
                  metabolites: Dict[str, CoRegMetabolite] = None,
                  biomass: CoRegBiomass = None,
                  time_steps: Sequence[float] = None,
                  soft_plus: float = 0,
                  tolerance: float = ModelConstants.TOLERANCE,
                  scale: bool = False) -> Union[DynamicSolution, ModelSolution, Solution, List[Solution]]:
        """
        CoRegFlux optimization method.
        It supports steady state and dynamic optimization.
        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False
        :param solver_kwargs: Keyword arguments to pass to the solver. See LinearProblem.optimize for details.
        :param initial_state: a dictionary of targets ids and expression predictions to set as initial state.
        :param metabolites: a dictionary of metabolites ids and concentrations to set as initial state
        :param biomass: a float value to set as initial biomass
        :param time_steps: a list of time points to simulate the model over
        :param soft_plus: the soft plus parameter to use for the gene state update
        :param tolerance: the tolerance to use for the gene state update
        :param scale: whether to scale the gene state update
        :return: a ModelSolution instance if dynamic is False,
        a DynamicSolution instance otherwise (if to_solver is False)
        """
        if len(initial_state) == 1:
            return self._steady_state_optimize(to_solver=to_solver,
                                               solver_kwargs=solver_kwargs,
                                               initial_state=initial_state[0],
                                               metabolites=metabolites,
                                               biomass=biomass,
                                               time_step=time_steps[0],
                                               soft_plus=soft_plus,
                                               tolerance=tolerance,
                                               scale=scale)

        return self._dynamic_optimize(to_solver=to_solver,
                                      solver_kwargs=solver_kwargs,
                                      initial_state=initial_state,
                                      metabolites=metabolites,
                                      biomass=biomass,
                                      time_steps=time_steps,
                                      soft_plus=soft_plus,
                                      tolerance=tolerance,
                                      scale=scale)

    def optimize(self,
                 to_solver: bool = False,
                 solver_kwargs: Dict = None,
                 initial_state: Union[Dict[str, float], Sequence[Dict[str, float]]] = None,
                 metabolites: Dict[str, float] = None,
                 growth_rate: float = None,
                 time_steps: Sequence[float] = None,
                 soft_plus: float = 0,
                 tolerance: float = ModelConstants.TOLERANCE,
                 scale: bool = False) -> Union[DynamicSolution, ModelSolution, Solution, List[Solution]]:
        """
        CoRegFlux optimization method.
        It supports steady state and dynamic optimization.
        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False
        :param solver_kwargs: Keyword arguments to pass to the solver. See LinearProblem.optimize for details.
        :param initial_state: a dictionary of targets ids and expression predictions to set as initial state. For dynamic
        optimization, a list of such dictionaries can be provided to set the initial state at each time point.
        :param metabolites: a dictionary of metabolites ids and concentrations to set as initial state
        :param growth_rate: the initial growth rate to set as initial state
        :param time_steps: a list of time points to simulate the model over
        :param soft_plus: the soft plus parameter to use for the gene state update
        :param tolerance: the tolerance to use for the gene state update
        :param scale: whether to scale the gene state update
        :return: a ModelSolution instance if dynamic is False,
        a DynamicSolution instance otherwise (if to_solver is False)
        """
        if not solver_kwargs:
            solver_kwargs = {}
        solver_kwargs['get_values'] = True

        if time_steps is None:
            time_steps = [1]

        if not initial_state:
            initial_state = [{}]

        else:
            if isinstance(initial_state, dict):
                initial_state = [initial_state]

        if len(initial_state) != len(time_steps):
            raise ValueError("The number of time steps must match the number of initial states.")

        if not metabolites:
            metabolites = {}

        if growth_rate is None:
            _, growth_rate = _run_and_decode(self, solver_kwargs=solver_kwargs)

        metabolites = build_metabolites(self.model, metabolites)
        biomass = build_biomass(self.model, growth_rate)
        return self._optimize(to_solver=to_solver,
                              solver_kwargs=solver_kwargs,
                              initial_state=initial_state,
                              metabolites=metabolites,
                              biomass=biomass,
                              time_steps=time_steps,
                              soft_plus=soft_plus,
                              tolerance=tolerance,
                              scale=scale)


# ----------------------------------------------------------------------------------------------------------------------
# Preprocessing using LinearRegression to predict genes expression from regulators co-expression
# Useful for CoRegFlux method
# ----------------------------------------------------------------------------------------------------------------------
def _get_target_regulators(gene: Union['Gene', 'Target'] = None) -> List[str]:
    """
    It returns the list of regulators of a target gene
    :param gene: Target gene
    :return: List of regulators of the target gene
    """
    if gene is None:
        return []

    if gene.is_target():
        return [regulator.id for regulator in gene.yield_regulators()]

    return []


def _filter_influence_and_expression(interactions: Dict[str, List[str]],
                                     influence: pd.DataFrame,
                                     expression: pd.DataFrame,
                                     experiments: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    It filters influence, expression and experiments matrices to keep only the targets and their regulators.
    :param interactions: Dictionary with the interactions between targets and regulators
    :param influence: Influence matrix
    :param expression: Expression matrix
    :param experiments: Experiments matrix
    :return: Filtered influence matrix, filtered expression matrix, filtered experiments matrix
    """
    targets = pd.Index(set(interactions.keys()))
    regulators = pd.Index(set([regulator for regulators in interactions.values() for regulator in regulators]))

    # filter the expression matrix for the target genes only
    expression = expression.loc[expression.index.intersection(targets)].copy()
    influence = influence.loc[influence.index.intersection(regulators)].copy()
    experiments = experiments.loc[experiments.index.intersection(regulators)].copy()
    return influence, expression, experiments


def _predict_experiment(interactions: Dict[str, List[str]],
                        influence: pd.DataFrame,
                        expression: pd.DataFrame,
                        experiment: pd.Series) -> pd.Series:
    try:
        # noinspection PyPackageRequirements
        from sklearn.linear_model import LinearRegression
    except ImportError:
        raise ImportError('The package sklearn is not installed. '
                          'To compute the probability of target-regulator interactions, please install sklearn '
                          '(pip install sklearn).')

    predictions = {}

    for target, regulators in interactions.items():

        if not regulators:
            predictions[target] = np.nan
            continue

        if target not in expression.index:
            predictions[target] = np.nan
            continue

        if not set(regulators).issubset(influence.index):
            predictions[target] = np.nan
            continue

        if not set(regulators).issubset(experiment.index):
            predictions[target] = np.nan
            continue

        # a linear regression model is trained for
        # y = expression of the target gene for all samples
        # x1 = influence score of the regulator 1 in the train data set
        # x2 = influence score of the regulator 2 in the train data set
        # x3 ...

        x = influence.loc[regulators].transpose().to_numpy()
        y = expression.loc[target].to_numpy()

        regressor = LinearRegression()
        regressor.fit(x, y)

        # the expression of the target gene is predicted for the experiment
        x_pred = experiment.loc[regulators].to_frame().transpose().to_numpy()
        predictions[target] = regressor.predict(x_pred)[0]

    return pd.Series(predictions)


def predict_gene_expression(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                            influence: pd.DataFrame,
                            expression: pd.DataFrame,
                            experiments: pd.DataFrame) -> pd.DataFrame:
    """
    It predicts the expression of genes in the experiments set using the co-expression of regulators
    in the expression and influence datasets.
    Adapted from CoRegFlux docs:
        - A GEM model containing GPRs
        - A TRN network containing target genes and the co-activators and co-repressors of each target gene
        - GEM genes and TRN targets must match
        - An influence dataset containing the influence scores of the regulators. Influence score is similar to a
        correlation score between the expression of the regulator and the expression of the target gene.
        Influence dataset format: (rows: regulators, columns: samples). Influence scores are calculated using the
        CoRegNet algorithm in the gene expression dataset.
        - A gene expression dataset containing the expression of the genes. Also called the training dataset.
        The gene expression dataset format: (rows: genes, columns: samples)
        - An experiments dataset containing the influence scores of the regulators in the experiments. These experiments
        are not used for training the linear regression model (the gene state predictor).
        These experiments are condition-specific environmental or genetic conditions that can be used
        to perform phenotype simulations.
        Experiments dataset format: (rows: regulators, columns: experiments/conditions)

    The result is a matrix of predicted gene expression values for each experiment.
    :param model: an integrated Metabolic-Regulatory model aka a GERM model
    :param influence: Influence scores of the regulators in the train data set
    :param expression: Expression of the genes in the train data set
    :param experiments: Influence scores of the regulators in the test data set
    :return: Predicted expression of the genes in the test data set
    """
    # Filtering only the gene expression and influences data of metabolic genes available in the model
    interactions = {target.id: _get_target_regulators(target) for target in model.yield_targets()}
    influence, expression, experiments = _filter_influence_and_expression(interactions=interactions,
                                                                          influence=influence,
                                                                          expression=expression,
                                                                          experiments=experiments)

    predictions = []
    for column in experiments.columns:
        experiment_prediction = _predict_experiment(interactions=interactions,
                                                    influence=influence,
                                                    expression=expression,
                                                    experiment=experiments[column])
        predictions.append(experiment_prediction)

    predictions = pd.concat(predictions, axis=1)
    predictions.columns = experiments.columns
    return predictions.dropna()
