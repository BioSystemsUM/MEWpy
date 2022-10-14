from typing import Union, Dict, Tuple, List
from warnings import warn

from mewpy.germ.analysis import FBA
from mewpy.solvers import Solution
from mewpy.solvers.solver import Solver
from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel
from mewpy.util.constants import ModelConstants
from mewpy.germ.solution import ModelSolution, DynamicSolution


def _get_boolean_state_from_reaction_flux(flux_rate: float) -> bool:
    # A reaction can have different flux values between two solutions,
    # though the state remains the same. Thus, non-zero flux stands for ON state
    # while zero flux stands for OFF state
    if abs(flux_rate) > ModelConstants.TOLERANCE:
        return True

    return False


def _find_duplicated_state(state, regulatory_solution, regulatory_reactions, regulatory_metabolites):
    mask = []
    for regulator, value in regulatory_solution.items():

        state_value = state[regulator]

        if regulator in regulatory_reactions or regulator in regulatory_metabolites:
            state_value = _get_boolean_state_from_reaction_flux(state_value)
            value = _get_boolean_state_from_reaction_flux(value)

        mask.append(value == state_value)

    return all(mask)


class RFBA(FBA):

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = False,
                 attach: bool = False):
        """
        Regulatory Flux Balance Analysis (RFBA) of a metabolic-regulatory model.
        Implementation of a steady-state and dynamic versions of RFBA for a integrated metabolic-regulatory model.

        For more details consult Covert et al. 2004 at https://doi.org/10.1038/nature02456

        :param model: a MetabolicModel, RegulatoryModel or GERM model. The model is used to retrieve
        variables and constraints to the linear problem
        :param solver: A Solver, CplexSolver, GurobiSolver or OptLangSolver instance.
        Alternatively, the name of the solver is also accepted.
        The solver interface will be used to load and solve a linear problem in a given solver.
        If none, a new solver is instantiated. An instantiated solver may be used, but it will be overwritten
        if build is true.
        :param build: Whether to build the linear problem upon instantiation. Default: False
        :param attach: Whether to attach the linear problem to the model upon instantiation. Default: False
        """
        self._regulatory_reactions = []
        self._regulatory_metabolites = []
        super().__init__(model=model, solver=solver, build=build, attach=attach)

    def _build(self):
        """
        It builds the linear problem for RFBA. It is a linear problem with the following constraints:
            - Metabolic constraints

        The regulatory layer is not considered in the linear problem. It is only considered in the simulation step.
        :return:
        """
        if self.model.is_metabolic():
            self._build_mass_constraints()

            self._regulatory_reactions = [rxn.id
                                          for rxn in self.model.yield_reactions()
                                          if rxn.is_regulator()]

            self._regulatory_metabolites = [met.id
                                            for met in self.model.yield_metabolites()
                                            if met.is_regulator()]

            self._linear_objective = {var.id: value for var, value in self.model.objective.items()}
            self._minimize = False

    def initial_state(self, state: Dict[str, float] = None) -> Dict[str, float]:
        """
        Method responsible for retrieving the initial state of the model.
        The initial state is the state of all regulators found in the Metabolic-Regulatory model.
        :param state: the initial state of the model
        :return: dict of regulatory/metabolic variable keys (regulators) and a value of 0 or 1
        """
        if not state:
            state = {}

        if not self.model.is_regulatory():
            return state

        initial_state = {}
        for regulator in self.model.yield_regulators():
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

    def decode_regulatory_state(self, state: Dict[str, float]) -> Dict[str, float]:
        """
        It solves the boolean regulatory network for a given specific state.
        It also updates all targets having a valid regulatory interaction associated with it for the resulting state

        :param state: dict of regulatory variable keys (regulators) and a value of 0, 1 or float
        (reactions and metabolites predicates)
        :return: dict of target keys and a value of the resulting state
        """
        if not self.model.is_regulatory():
            return {}

        # Solving regulatory model synchronously for the regulators according to the initial state
        # Targets are associated with a single regulatory interaction
        result = {}
        for interaction in self.model.yield_interactions():

            # solving regulators state only
            if not interaction.target.is_regulator():
                continue

            # an interaction can have multiple regulatory events, namely one for 0 and another for 1
            for coefficient, event in interaction.regulatory_events.items():
                if event.is_none:
                    continue

                eval_result = event.evaluate(values=state)
                if eval_result:
                    result[interaction.target.id] = coefficient
                else:
                    result[interaction.target.id] = 0.0
        return result

    def decode_metabolic_state(self, state: Dict[str, float]) -> Dict[str, float]:
        """
        It solves the boolean regulatory network for a given specific state.
        It also updates all targets having a valid regulatory interaction associated with it for the resulting state

        :param state: dict of regulatory variable keys (regulators) and a value of 0, 1 or float
        (reactions and metabolites predicates)
        :return: dict of target keys and a value of the resulting state
        """
        if not self.model.is_regulatory():
            return {}

        # Solving the whole regulatory model synchronously, as asynchronously would take too much time
        # Targets are associated with a single regulatory interaction
        result = {}
        for interaction in self.model.yield_interactions():

            # an interaction can have multiple regulatory events, namely one for 0 and another for 1
            for coefficient, event in interaction.regulatory_events.items():
                if event.is_none:
                    continue

                eval_result = event.evaluate(values=state)
                if eval_result:
                    result[interaction.target.id] = coefficient
                else:
                    result[interaction.target.id] = 0.0
        return result

    def decode_constraints(self, state: Dict[str, float]) -> Dict[str, Tuple[float, float]]:
        """
        Method responsible for decoding the RFBA metabolic state, namely the state of all metabolic genes associated
        at least with one reaction in the GPRs rule.

        :param state: dict of regulatory/metabolic variable keys (metabolic target) and a value of 0 or 1
        :return: dict of constraints of the resulting metabolic state
        """
        # This method retrieves the constraints associated with a given metabolic/regulatory state

        if not self.model.is_metabolic():
            return {}

        constraints = {}
        for rxn in self.model.yield_reactions():

            if rxn.gpr.is_none:
                continue

            res = rxn.gpr.evaluate(values=state)

            if not res:
                constraints[rxn.id] = (0.0, 0.0)

        return constraints

    def next_state(self, state: Dict[str, float], solver_kwargs: Dict = None) -> Tuple[Dict[str, float],
                                                                                       Solution]:
        """
        Retrieves the next state for the provided state

        Solves the boolean regulatory model using method regulatory_simulation(state) and decodes the metabolic state
        for that state or initial state (according to the flag regulatory_state)

        :param state: dict of regulatory/metabolic variable keys (regulatory and metabolic target) and a value of 0, 1
        or float (reactions and metabolites predicates)
        :param solver_kwargs: solver kwargs
        :return: dict of all regulatory/metabolic variables keys and a value of the resulting state
        """
        state = state.copy()

        if not solver_kwargs:
            solver_kwargs = {}

        if 'constraints' in solver_kwargs:
            constraints = solver_kwargs['constraints'].copy()
        else:
            constraints = {}

        # Regulatory state from a synchronous boolean simulation
        # noinspection PyTypeChecker
        regulatory_state = self.decode_regulatory_state(state=state)

        # Next state is the previous state plus the regulatory state
        next_state = {**state, **regulatory_state}

        # After a simulation of the regulators outputs, the state of the targets are retrieved now
        metabolic_state = self.decode_metabolic_state(state=next_state)

        regulatory_constraints = self.decode_constraints(metabolic_state)
        regulatory_constraints = {**constraints, **regulatory_constraints}

        solver_kwargs['constraints'] = regulatory_constraints
        solver_solution = self.solver.solve(**solver_kwargs)

        if solver_solution.values:
            solver_solution.values = {**metabolic_state, **solver_solution.values}
        else:
            solver_solution.values = {**metabolic_state, **{rxn: 0.0 for rxn in self.model.reactions}}

        # update the next state with the regulatory/metabolic state
        for reaction in self._regulatory_reactions:
            next_state[reaction] = solver_solution.values.get(reaction, 0.0)

        for metabolite in self._regulatory_metabolites:

            reaction = self.model.metabolites[metabolite].exchange_reaction

            if reaction:
                next_state[metabolite] = solver_solution.values.get(reaction.id, 0.0)

        return next_state, solver_solution

    def _dynamic_optimize(self,
                          initial_state: Dict[str, float] = None,
                          iterations: int = 10,
                          to_solver: bool = False,
                          solver_kwargs: Dict = None):
        """
        RFBA model dynamic simulation (until the metabolic-regulatory state is reached).

        First, the boolean regulatory model is solved using the initial state. If there is no initial state ,
        the state of all regulatory variables is inferred from the model.

        At the same time, the metabolic state is also decoded from the initial state of all metabolic genes.
        If otherwise set by using the regulatory flag, the metabolic state is decoded from the boolean regulatory model
        simulation.

        The result of the boolean regulatory model updates the state of all targets, while the result of the
        metabolic state decoding updates the state of all reactions and metabolites that are also regulatory variables.

        Then, this new resulting metabolic-regulatory state is used to solve again the boolean regulatory model and
        decode again the metabolic state towards the retrieval of a new regulatory and metabolic states.

        This cycle is iterated until a given state is repeated, the so called metabolic-regulatory steady-state.
        Hence, dynamic RFBA is based on finding the cyclic attractors of the metabolic-regulatory networks
        Alternatively, an iteration limit may be reached.

        Finally, all states between the cyclic attractor are used for decoding final metabolic states using
        the resulting metabolic-regulatory state. Thus, a solution/simulation is obtained for each mid-state in the
        cyclic attractor of the metabolic-regulatory networks

        Objective and constraint-based model analysis method (pFBA, FBA, ...) can be altered here reversibly.

        :param initial_state: a dictionary of variable ids and their values to set as initial state
        :param iterations: The maximum number of iterations. Default: 10
        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False
        :param solver_kwargs: Keyword arguments to pass to the solver.
        :return: A DynamicSolution instance or a list of solver Solutions if to_solver is True.
        """
        if not initial_state:
            initial_state = {}

        if not solver_kwargs:
            solver_kwargs = {}

        # It takes the initial state from the model and then updates with the initial state provided as input
        initial_state = self.initial_state(initial_state)

        regulatory_solutions = []
        solver_solutions = []

        # solve using the initial state
        # noinspection PyTypeChecker
        state, solver_solution = self.next_state(state=initial_state, solver_kwargs=solver_kwargs)
        regulatory_solutions.append(state)
        solver_solutions.append(solver_solution)

        i = 1
        steady_state = False
        while not steady_state:
            # Updating state upon state. See next state for further detail
            state, solver_solution = self.next_state(state=state, solver_kwargs=solver_kwargs)

            for regulatory_solution in regulatory_solutions:

                is_duplicated = _find_duplicated_state(state=state, regulatory_solution=regulatory_solution,
                                                       regulatory_reactions=self._regulatory_reactions,
                                                       regulatory_metabolites=self._regulatory_metabolites)
                if not is_duplicated:
                    continue

                steady_state = True
                break

            # add the new state to the list of regulatory solutions
            regulatory_solutions.append(state)
            solver_solutions.append(solver_solution)

            # if the maximum number of iterations is reached, the simulation is stopped
            if i < iterations:
                i += 1

            else:
                def iteration_limit(message):
                    warn(message, UserWarning, stacklevel=2)

                iteration_limit("Iteration limit reached")
                steady_state = True

        if to_solver:
            return solver_solutions

        minimize = solver_kwargs.get('minimize', self._minimize)
        solutions = [ModelSolution.from_solver(method=self.method,
                                               solution=solver_solution,
                                               model=self.model,
                                               minimize=minimize)
                     for solver_solution in solver_solutions]
        return DynamicSolution(*solutions)

    def _steady_state_optimize(self,
                               initial_state: Dict[str, float] = None,
                               to_solver: bool = False,
                               solver_kwargs: Dict = None) -> Union[ModelSolution, Solution]:
        """
        RFBA model one-step simulation (pseudo steady-state).

        First, the boolean regulatory model is simulated using the initial state. If there is no initial state ,
        the state of all regulatory variables is inferred from the model.

        At the same time, the metabolic state is also decoded from the initial state of all metabolic genes.
        If otherwise set by using the regulatory flag, the metabolic state is decoded from the boolean regulatory model
        simulation.

        The result of the boolean regulatory model updates the state of all targets, while the result of the
        metabolic state decoding updates the state of all reactions and metabolites that are also regulatory variables.

        Then, this new state is used to decode the final metabolic state using the resulting metabolic-regulatory state.

        Objective and constraint-based model analysis method (pFBA, FBA, ...) can be altered here reversibly.

        :param initial_state: a dictionary of variable ids and their values to set as initial state
        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False
        :param solver_kwargs: A dictionary with the solver arguments. Default: None
        :return: A ModelSolution instance or a SolverSolution instance if to_solver is True.
        """
        if not initial_state:
            initial_state = {}

        if not solver_kwargs:
            solver_kwargs = {}

        if 'constraints' in solver_kwargs:
            constraints = solver_kwargs['constraints'].copy()
        else:
            constraints = {}

        # It takes the initial state from the model and then updates with the initial state provided as input
        initial_state = self.initial_state(initial_state)

        # Regulatory state from a synchronous boolean simulation
        # noinspection PyTypeChecker
        regulatory_state = self.decode_regulatory_state(state=initial_state)

        # After a simulation of the regulators outputs, the state of the targets are retrieved now
        metabolic_state = self.decode_metabolic_state(state={**initial_state, **regulatory_state})

        regulatory_constraints = self.decode_constraints(metabolic_state)
        metabolic_regulatory_constraints = {**constraints, **regulatory_constraints}

        solver_kwargs['constraints'] = metabolic_regulatory_constraints
        solution = self.solver.solve(**solver_kwargs)

        if solution.values:
            solution.values = {**initial_state, **regulatory_state, **solution.values}

        if to_solver:
            return solution

        minimize = solver_kwargs.get('minimize', self._minimize)
        return ModelSolution.from_solver(method=self.method, solution=solution, model=self.model, minimize=minimize)

    def _optimize(self,
                  initial_state: Dict[str, float] = None,
                  dynamic: bool = False,
                  iterations: int = 10,
                  to_solver: bool = False,
                  solver_kwargs: Dict = None,
                  **kwargs) -> Union[DynamicSolution, ModelSolution, List[Solution]]:
        """
        RFBA simulation.

        First, the boolean regulatory model is solved using the initial state. If there is no initial state ,
        the state of all regulatory variables is inferred from the model.

        At the same time, the metabolic state is also decoded from the initial state of all metabolic genes.
        If otherwise set by using the regulatory flag, the metabolic state is decoded from the boolean regulatory model
        simulation.

        The result of the boolean regulatory model updates the state of all targets, while the result of the
        metabolic state decoding updates the state of all reactions and metabolites that are also regulatory variables.

        Then, this new resulting metabolic-regulatory state is used to solve again the boolean regulatory model and
        decode again the metabolic state towards the retrieval of a new regulatory and metabolic states.

        If the dynamic flag is set to True, the simulation will continue until the metabolic state remains the same
        between two consecutive simulations, the so called metabolic-regulatory steady-state.
        Hence, dynamic RFBA is based on finding the cyclic attractors of the metabolic-regulatory networks
        Alternatively, an iteration limit may be reached.

        Finally, all states between the cyclic attractor are used for decoding final metabolic states using
        the resulting metabolic-regulatory state. Thus, a solution/simulation is obtained for each mid-state in the
        cyclic attractor of the metabolic-regulatory networks

        Objective and constraint-based model analysis method (pFBA, FBA, ...) can be altered here reversibly.

        :param initial_state: a dictionary of variable ids and their values to set as initial state
        :param dynamic: If True, the model is simulated until the metabolic-regulatory steady-state is reached.
        :param iterations: The maximum number of iterations. Default: 10
        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False
        :param solver_kwargs: Keyword arguments to pass to the solver. See LinearProblem.optimize for details.
        :return: A DynamicSolution instance or a list of solver Solutions if to_solver is True.
        """
        if dynamic:
            return self._dynamic_optimize(initial_state=initial_state, iterations=iterations,
                                          to_solver=to_solver, solver_kwargs=solver_kwargs)

        return self._steady_state_optimize(initial_state=initial_state,
                                           to_solver=to_solver, solver_kwargs=solver_kwargs)

    def optimize(self,
                 to_solver: bool = False,
                 solver_kwargs: Dict = None,
                 initial_state: Dict[str, float] = None,
                 dynamic: bool = False,
                 iterations: int = 10,
                 **kwargs) -> Union[DynamicSolution, ModelSolution, Solution, List[Solution]]:
        """
        RFBA simulation.

        First, the boolean regulatory model is solved using the initial state. If there is no initial state ,
        the state of all regulatory variables is inferred from the model.

        At the same time, the metabolic state is also decoded from the initial state of all metabolic genes.
        If otherwise set by using the regulatory flag, the metabolic state is decoded from the boolean regulatory model
        simulation.

        The result of the boolean regulatory model updates the state of all targets, while the result of the
        metabolic state decoding updates the state of all reactions and metabolites that are also regulatory variables.

        Then, this new resulting metabolic-regulatory state is used to solve again the boolean regulatory model and
        decode again the metabolic state towards the retrieval of a new regulatory and metabolic states.

        If the dynamic flag is set to True, the simulation will continue until the metabolic state remains the same
        between two consecutive simulations, the so called metabolic-regulatory steady-state.
        Hence, dynamic RFBA is based on finding the cyclic attractors of the metabolic-regulatory networks
        Alternatively, an iteration limit may be reached.

        Finally, all states between the cyclic attractor are used for decoding final metabolic states using
        the resulting metabolic-regulatory state. Thus, a solution/simulation is obtained for each mid-state in the
        cyclic attractor of the metabolic-regulatory networks

        Objective and constraint-based model analysis method (pFBA, FBA, ...) can be altered here reversibly.

        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False
        :param solver_kwargs: Keyword arguments to pass to the solver. See LinearProblem.optimize for details.
        :param initial_state: a dictionary of variable ids and their values to set as initial state
        :param dynamic: If True, the model is simulated until the metabolic-regulatory steady-state is reached.
        :param iterations: The maximum number of iterations. Default: 10
        :return: A DynamicSolution instance or a list of solver Solutions if to_solver is True.
        """
        # build solver if out of sync
        if not self.synchronized:
            self.build()

        if not solver_kwargs:
            solver_kwargs = {}

        # concrete optimize
        return self._optimize(to_solver=to_solver, solver_kwargs=solver_kwargs,
                              initial_state=initial_state, dynamic=dynamic, iterations=iterations,
                              **kwargs)
