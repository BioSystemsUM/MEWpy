from typing import Union
from warnings import warn

from mewpy.solvers.solver import Solver
from mewpy.model import Model, MetabolicModel, RegulatoryModel
from mewpy.util.constants import ModelConstants
from mewpy.lp import MetabolicLinearizer
from mewpy.solution import ModelSolution, DynamicSolution


# TODO: type hinting and documentation
class RFBA(MetabolicLinearizer):

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = True,
                 attach: bool = False):

        """
        Regulatory Flux Balance Analysis (RFBA) of a metabolic-regulatory model.
        Implementation of a steady-state and dynamic versions of RFBA for a integrated metabolic-regulatory model.

        For more details consult Covert et al. 2004 at https://doi.org/10.1038/nature02456

        :param model: a mewpy Model, MetabolicModel, RegulatoryModel or all. The model is used to retrieve
        variables and constraints to the linear problem

        :param solver: A Solver, CplexSolver, GurobiSolver or OptLangSolver instance.
        Alternatively, the name of the solver is also accepted.
        The solver interface will be used to load and solve a linear problem in a given solver.
        If none, a new solver is instantiated. An instantiated solver may be used but it will be overwritten
        if build is true.

        :param build: Whether to build the linear problem upon instantiation. Default: False
        :param attach: Whether to attach the linear problem to the model upon instantiation. Default: False
        """

        self._regulatory_reactions = []
        self._regulatory_metabolites = []

        super().__init__(model=model,
                         solver=solver,
                         build=build,
                         attach=attach)

    def build(self):

        if self._variables or self._constraints:
            self.clean()

        if self.model.is_metabolic():
            self.metabolite_reaction_lookup(reactions=self.model.yield_reactions(),
                                            metabolites=self.model.yield_metabolites())

            self._regulatory_reactions = [rxn.id
                                          for rxn in self.model.yield_reactions()
                                          if rxn.is_regulator()]

            self._regulatory_metabolites = [met.id
                                            for met in self.model.yield_metabolites()
                                            if met.is_regulator()]

            self.update()

            self.set_objective(linear={var.id: value for var, value in self.model.objective.items()},
                               minimize=False)

    def initial_state(self):

        if self.model.is_regulatory():
            return {regulator.id: regulator.coefficient.active_coefficient
                    for regulator in self.model.yield_regulators()}

        return {}

    def sim_bool(self, state):

        """
        Solves the boolean regulatory network for a given specific state.
        It also updates all targets having a valid regulatory interaction associated with it for the resulting state

        :param state: dict of regulatory variable keys (regulators) and a value of 0, 1 or float
        (reactions and metabolites predicates)

        :return: dict of target keys and a value of the resulting state
        """

        result = {}

        # solving regulatory rule by regulatory rule (synchronously though, as asynchronously would take too much time)
        # Each target has one and only one regulatory interaction
        # The regulatory model object only supports this strategy, though this can easily be altered. If so,
        # this solve method must be changed too

        if self.model.is_regulatory():

            # interaction over model.interactions.
            for interaction in self.model.yield_interactions():

                # regulatory event over interaction.regulatory events.
                # an interaction can have multiple regulatory events.
                # Regularly, there is one event for each possible coefficient of the target variable.
                # If the regulatory event is evaluated to True or 1, this is the outcome coefficient of the target
                # By default the outcome of the target is a zero coefficient

                is_empty = True
                active_coefficient = None

                for coefficient, regulatory_event in interaction.regulatory_events.items():

                    if not regulatory_event.is_none:

                        is_empty = False

                        regulatory_event_eval = regulatory_event.evaluate(values=state)

                        if regulatory_event_eval:
                            active_coefficient = coefficient

                            break

                if is_empty:
                    target_val = interaction.target.coefficient.active_coefficient

                elif active_coefficient is not None:
                    target_val = active_coefficient

                else:
                    target_val = 0.0

                result[interaction.target.id] = target_val

        return result

    def decode_metabolic_state(self, state):

        """
        Method responsible for decoding the RFBA metabolic state, namely the state of all metabolic genes associated
        at least with one reaction in the GPRs rule.

        :param state: dict of regulatory/metabolic variable keys (metabolic target) and a value of 0 or 1

        :return: dict of constraints of the resulting metabolic state
        """

        # This method retrieves the constraints associated with a given metabolic/regulatory state

        constraints = {}

        if self.model.is_metabolic():

            # gpr over model.reactions
            for rxn in self.model.yield_reactions():

                if rxn.gpr.is_none:
                    continue

                else:

                    res = rxn.gpr.evaluate(values=state)

                    if not res:
                        constraints[rxn.id] = (0.0, 0.0)

        return constraints

    def next_state(self, state, constraints):

        """
        Retrieves the next state for a given state

        Solves the boolean regulatory model using method regulatory_simulation(state) and decodes the metabolic state
        for that state or initial state (according to the flag regulatory_state)

        :param state: dict of regulatory/metabolic variable keys (regulatory and metabolic target) and a value of 0, 1
        or float (reactions and metabolites predicates)
        :param constraints:

        :return: dict of all regulatory/metabolic variables keys and a value of the resulting state
        """

        # Regulatory rules from the boolean network
        regulatory_state = self.sim_bool(state=state)

        # Reactions and metabolites inactivated by the genes, reactions and metabolites
        # in the initial state or regulatory state simulation
        metabolic_state = {gene.id: gene.coefficient.active_coefficient
                           for gene in self.model.yield_genes()}
        metabolic_regulatory_state = {**metabolic_state, **regulatory_state}

        regulatory_constraints = self.decode_metabolic_state(metabolic_regulatory_state)
        metabolic_regulatory_constraints = {**constraints, **regulatory_constraints}

        solver_solution = self._solver.solve(linear=None,
                                             quadratic=None,
                                             minimize=None,
                                             model=None,
                                             constraints=metabolic_regulatory_constraints,
                                             get_values=True,
                                             shadow_prices=False,
                                             reduced_costs=False,
                                             pool_size=0,
                                             pool_gap=None)

        if not solver_solution.values:
            solver_solution.values = {}

        next_state = {**state, **metabolic_regulatory_state}

        # update the regulatory/metabolic state
        for reaction in self._regulatory_reactions:

            next_state[reaction] = solver_solution.values.get(reaction, 0.0)

        for metabolite in self._regulatory_metabolites:

            reaction = self.model.get(metabolite).exchange_reaction

            if reaction:

                next_state[metabolite] = solver_solution.values.get(reaction.id, 0.0)

        return next_state

    def steady_state_optimize(self,
                              initial_state=None,
                              objective=None,
                              minimize=False,
                              constraints=None,
                              to_solver=False,
                              get_values=True,
                              shadow_prices=False,
                              reduced_costs=False,
                              pool_size=0,
                              pool_gap=None):

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

        :return: ModelSolution object with metabolic/regulatory solution
        """

        self.set_objective(linear=objective, minimize=minimize)

        # It takes the initial state from the model and then updates with the initial state provide as input
        if not constraints:
            constraints = {}
        else:
            constraints = constraints.copy()

        if not initial_state:
            initial_state = self.initial_state()

        else:
            initial_state = {**self.initial_state(), **initial_state}

        # Regulatory rules from the boolean network
        regulatory_state = self.sim_bool(state=initial_state)

        # Reactions and metabolites inactivated by the genes, reactions and metabolites
        # in the initial state or regulatory state simulation
        metabolic_state = {gene.id: gene.coefficient.active_coefficient
                           for gene in self.model.yield_genes()}
        metabolic_regulatory_state = {**metabolic_state, **regulatory_state}

        regulatory_constraints = self.decode_metabolic_state(metabolic_regulatory_state)
        metabolic_regulatory_constraints = {**constraints, **regulatory_constraints}

        solver_solution = self._solver.solve(linear=None,
                                             quadratic=None,
                                             minimize=None,
                                             model=None,
                                             constraints=metabolic_regulatory_constraints,
                                             get_values=get_values,
                                             shadow_prices=shadow_prices,
                                             reduced_costs=reduced_costs,
                                             pool_size=pool_size,
                                             pool_gap=pool_gap)

        if solver_solution.values:
            solver_solution.values = {**initial_state, **metabolic_regulatory_state, **solver_solution.values}

        if to_solver:
            return solver_solution

        if minimize:
            sense = 'minimize'

        else:
            sense = 'maximize'

        return ModelSolution.from_solver(method=self.method,
                                         solution=solver_solution,
                                         model=self.model,
                                         objective_direction=sense)

    def dynamic_optimize(self,
                         initial_state=None,
                         iter_limit=10,
                         objective=None,
                         minimize=False,
                         constraints=None,
                         to_solver=False,
                         get_values=True,
                         shadow_prices=False,
                         reduced_costs=False,
                         pool_size=0,
                         pool_gap=None):

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

        :return:
        """

        self.set_objective(linear=objective, minimize=minimize)

        # It takes the initial state from the model and then updates with the initial state provide as input
        if not constraints:
            constraints = {}
        else:
            constraints = constraints.copy()

        if not initial_state:
            initial_state = self.initial_state()

        else:
            initial_state = {**self.initial_state(), **initial_state}

        regulatory_solution = []

        # solve using the initial state
        state = self.next_state(state=initial_state, constraints=constraints)
        regulatory_solution.append(state)

        i = 1
        sol = 1
        steady_state = False
        while not steady_state:

            # Updating state upon state. See next state for further detail
            state = self.next_state(state=state, constraints=constraints)
            regulatory_solution.append(state)

            n_sols = len(regulatory_solution) - 1

            for sol in range(n_sols):

                # A reaction can have different flux values between two solutions, though the state remains the same.
                # Thus, non-zero flux stands for ON state while zero flux stands for OFF state

                replica_v = []

                for key, val in regulatory_solution[sol].items():

                    sol_val = regulatory_solution[n_sols][key]

                    if key in self._regulatory_reactions or key in self._regulatory_metabolites:

                        if abs(val) > ModelConstants.TOLERANCE:
                            val = True
                        else:
                            val = False

                        if abs(sol_val) > ModelConstants.TOLERANCE:
                            sol_val = True
                        else:
                            sol_val = False

                    replica_v.append(val == sol_val)

                is_replica = all(replica_v)

                if is_replica:
                    steady_state = True

                    break

            if i < iter_limit:
                i += 1
            else:

                def iteration_limit(message):
                    warn(message, UserWarning, stacklevel=2)

                iteration_limit("Iteration limit reached")
                steady_state = True

        # The solution comprehends and FBA simulation between the first and second similar regulatory and metabolic
        # states using the corresponding regulatory states to switch on or off the associated reactions
        solution = []
        for i in range(sol, len(regulatory_solution)):

            regulatory_state = regulatory_solution[i]

            regulatory_state_constraints = self.decode_metabolic_state(regulatory_state)

            state_constraints = {**constraints, **regulatory_state_constraints}

            solver_solution = self._solver.solve(linear=None,
                                                 quadratic=None,
                                                 minimize=None,
                                                 model=None,
                                                 constraints=state_constraints,
                                                 get_values=get_values,
                                                 shadow_prices=shadow_prices,
                                                 reduced_costs=reduced_costs,
                                                 pool_size=pool_size,
                                                 pool_gap=pool_gap)

            if solver_solution.values:
                solver_solution.values = {**regulatory_state, **solver_solution.values}

            if to_solver:
                solution.append(solver_solution)

                continue

            if minimize:
                sense = 'minimize'

            else:
                sense = 'maximize'

            sol = ModelSolution.from_solver(method=self.method,
                                            solution=solver_solution,
                                            model=self.model,
                                            objective_direction=sense)

            solution.append(sol)

        if to_solver:
            return solution

        return DynamicSolution(*solution)

    def optimize(self,
                 initial_state=None,
                 dynamic=False,
                 iter_limit=10,
                 objective=None,
                 minimize=False,
                 constraints=None,
                 to_solver=False,
                 get_values=True,
                 shadow_prices=False,
                 reduced_costs=False,
                 pool_size=0,
                 pool_gap=None):

        if dynamic:
            return self.dynamic_optimize(initial_state=initial_state,
                                         iter_limit=iter_limit,
                                         objective=objective,
                                         minimize=minimize,
                                         constraints=constraints,
                                         to_solver=to_solver,
                                         get_values=get_values,
                                         shadow_prices=shadow_prices,
                                         reduced_costs=reduced_costs,
                                         pool_size=pool_size,
                                         pool_gap=pool_gap)

        return self.steady_state_optimize(initial_state=initial_state,
                                          objective=objective,
                                          minimize=minimize,
                                          constraints=constraints,
                                          to_solver=to_solver,
                                          get_values=get_values,
                                          shadow_prices=shadow_prices,
                                          reduced_costs=reduced_costs,
                                          pool_size=pool_size,
                                          pool_gap=pool_gap)
