from typing import TYPE_CHECKING, Union, Dict, Sequence, List

from mewpy.mew.analysis import FBA
from mewpy.mew.analysis.analysis_utils import biomass_yield_to_rate, \
    CoRegMetabolite, CoRegBiomass, metabolites_constraints, gene_state_constraints, system_state_update, \
    build_metabolites, build_biomass, CoRegResult
from mewpy.mew.solution import ModelSolution, DynamicSolution
from mewpy.solvers.solution import Solution
from mewpy.solvers.solver import Solver
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel


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
        return Solution(status='Optimal', fobj=result.objective_value, values=result.values)

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
        :param model: a MEW model aka an integrated RegulatoryMetabolicModel
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
                          initial_state: Dict[str, float] = None,
                          metabolites: Dict[str, CoRegMetabolite] = None,
                          biomass: CoRegBiomass = None,
                          time_steps: Sequence[float] = None,
                          soft_plus: float = 0,
                          tolerance: float = ModelConstants.TOLERANCE,
                          scale: bool = False) -> Union[DynamicSolution, Dict[float, Solution]]:
        solutions = []

        previous_time_step = 0
        for time_step in time_steps:
            time_step_diff = time_step - previous_time_step

            next_state = self.next_state(solver_kwargs=solver_kwargs,
                                         state=initial_state,
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
                               soft_plus: float = 0,
                               tolerance: float = ModelConstants.TOLERANCE,
                               scale: bool = False) -> Union[ModelSolution, Solution]:

        result = self.next_state(solver_kwargs=solver_kwargs,
                                 state=initial_state,
                                 metabolites=metabolites,
                                 biomass=biomass,
                                 time_step=1,
                                 soft_plus=soft_plus,
                                 tolerance=tolerance,
                                 scale=scale)
        return result_to_solution(result=result, model=self.model, to_solver=to_solver)

    def _optimize(self,
                  to_solver: bool = False,
                  solver_kwargs: Dict = None,
                  initial_state: Dict[str, float] = None,
                  dynamic: bool = False,
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
        :param initial_state: a dictionary of targets ids and expression predictions to set as initial state
        :param dynamic: If True, the model is simulated over a time course. Default: False
        :param metabolites: a dictionary of metabolites ids and concentrations to set as initial state
        :param biomass: a float value to set as initial biomass
        :param time_steps: a list of time points to simulate the model over
        :param soft_plus: the soft plus parameter to use for the gene state update
        :param tolerance: the tolerance to use for the gene state update
        :param scale: whether to scale the gene state update
        :return: a ModelSolution instance if dynamic is False,
        a DynamicSolution instance otherwise (if to_solver is False)
        """
        if dynamic:
            return self._dynamic_optimize(to_solver=to_solver,
                                          solver_kwargs=solver_kwargs,
                                          initial_state=initial_state,
                                          metabolites=metabolites,
                                          biomass=biomass,
                                          time_steps=time_steps,
                                          soft_plus=soft_plus,
                                          tolerance=tolerance,
                                          scale=scale)

        return self._steady_state_optimize(to_solver=to_solver,
                                           solver_kwargs=solver_kwargs,
                                           initial_state=initial_state,
                                           metabolites=metabolites,
                                           biomass=biomass,
                                           soft_plus=soft_plus,
                                           tolerance=tolerance,
                                           scale=scale)

    def optimize(self,
                 to_solver: bool = False,
                 solver_kwargs: Dict = None,
                 initial_state: Dict[str, float] = None,
                 dynamic: bool = False,
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
        :param initial_state: a dictionary of targets ids and expression predictions to set as initial state
        :param dynamic: If True, the model is simulated over a time course. Default: False
        :param metabolites: a dictionary of metabolites ids and concentrations to set as initial state
        :param growth_rate: the initial growth rate to set as initial state
        :param time_steps: a list of time points to simulate the model over
        :param soft_plus: the soft plus parameter to use for the gene state update
        :param tolerance: the tolerance to use for the gene state update
        :param scale: whether to scale the gene state update
        :return: a ModelSolution instance if dynamic is False,
        a DynamicSolution instance otherwise (if to_solver is False)
        """
        if dynamic and len(time_steps) == 0:
            raise ValueError('Time steps must be provided for dynamic optimization')

        if not solver_kwargs:
            solver_kwargs = {}

        if not initial_state:
            initial_state = {}

        if not metabolites:
            metabolites = {}

        if growth_rate is None:
            _, growth_rate = _run_and_decode(self, solver_kwargs=solver_kwargs)

        metabolites = build_metabolites(self.model, metabolites)
        biomass = build_biomass(self.model, growth_rate)
        return self._optimize(to_solver=to_solver,
                              solver_kwargs=solver_kwargs,
                              initial_state=initial_state,
                              dynamic=dynamic,
                              metabolites=metabolites,
                              biomass=biomass,
                              time_steps=time_steps,
                              soft_plus=soft_plus,
                              tolerance=tolerance,
                              scale=scale)
