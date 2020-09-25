from warnings import warn
from mewpy.regulation import IntegratedModel
from mewpy.simulation import SStatus


class RFBAModel(IntegratedModel):

    """

    RFBAModel class object inherits everything from the IntegratedModel class.
    This object is the standard for simulating a RFBA model.
    First, a boolean regulatory model is solved synchronously
    It provides the simulate method where one can know all about the metabolic and regulatory analysis of RFBA

    """

    def __init__(self,
                 identifier,
                 name=None,
                 cbm_model=None,
                 cbm_simulation_interface=None,
                 targets=None,
                 regulators=None,
                 regulatory_interactions=None,
                 initial_state=None):

        # Integrated model
        super().__init__(identifier,
                         name,
                         cbm_model,
                         cbm_simulation_interface,
                         targets,
                         regulators,
                         regulatory_interactions,
                         initial_state)

    def get_current_state(self):

        """
        Retrieves the current metabolic and regulatory state

        :return: dict, for convenience returns dict of all regulatory variables where key is the id of the regulatory variable and value can be 0, 1 or float (reactions and metabolites predicates)
        """

        return {var.id: var.expression_coef for var in self.regulatory_variables_gen()}

    def __update_metabolic_genes_state__(self, state):

        # Updates only the state of the metabolic genes ignoring TFs and others

        for key, val in state.items():

            if key in self.metabolic_regulatory_genes:
                self.metabolic_regulatory_genes[key].expression_coef = val

    def __update_reactions_state__(self, decode):

        # Updates the regulatory variables considered to be reactions and available in the model.
        # It does so for the cbm simulation. The absolute value is considered, as the flux predicates might not
        # consider reversibility

        if decode.status != SStatus.OPTIMAL and decode.status != SStatus.SUBOPTIMAL:

            for val in self.metabolic_regulatory_reactions.values():
                val.expression_coef = float(0.0)

        else:

            for val in self.metabolic_regulatory_reactions.values():
                val.expression_coef = abs(float(decode.fluxes[val.cbm_model_id]))

    def __update_metabolites_state__(self):

        # The next comment version seems to be more robust to actually capture the presence or absence of a
        # metabolite.
        # The current version is taken from the matlab version, which just considers extracellular metabolites and
        # the associated reaction lower bound. This can also be adjusted for intracellular following the same line of
        # thought

        for val in self.metabolic_regulatory_metabolites.values():
            reaction = self.cbm_simulation_interface.get_boundary_reaction(val.cbm_model_id)

            val.expression_coef = abs(float(self.cbm_simulation_interface.get_reaction_bounds(reaction)[0]))

    # def __update_metabolites_state__(self, decode):
    #
    #     if decode.status != SStatus.OPTIMAL and decode.status != SStatus.SUBOPTIMAL:
    #
    #         for val in self.metabolic_regulatory_metabolites.values():
    #             val.expression_coef = float(0.0)
    #
    #     else:
    #
    #         for val in self.metabolic_regulatory_metabolites.values():
    #
    #             reactions = self.cbm_simulation_interface.get_metabolite_reactions(val.cbm_model_id)
    #
    #             value = 0
    #             for reaction in reactions:
    #
    #                 if (self.cbm_simulation_interface.is_product(reaction, val.cbm_model_id) and
    #                         (float(decode.fluxes[reaction]) > float(1E-10))):
    #
    #                     if abs(float(decode.fluxes[reaction])) > value:
    #                         value = abs(float(decode.fluxes[reaction]))
    #
    #                 if (self.cbm_simulation_interface.is_reactant(reaction, val.cbm_model_id) and
    #                         (float(decode.fluxes[reaction]) < float(-1E-10))):
    #
    #                     if abs(float(decode.fluxes[reaction])) > value:
    #                         value = abs(float(decode.fluxes[reaction]))
    #
    #             val.expression_coef = value

    def solve_regulatory_model(self, state, silently=False):

        """
        Solves the boolean regulatory network for a given specific state.
        It also updates all targets having a valid regulatory interaction associated with it for the resulting state


        :param state: dict, key is the id of the regulatory variable while value can be 0, 1 or float (reactions and metabolites predicates)
        :param silently: bool, if True, the state (expression_coef) of the regulatory variables is not changed
        :return: dict, for convenience returns dict where key is the id of the target and value can be 0, 1 or float (reactions and metabolites predicates)
        """

        res = {}

        # solving regulatory rule by regulatory rule (synchronously though, as asynchronously would take too much time)
        # Each target has one and only one regulatory interaction
        # The regulatory model object only supports this strategy, though this can easily be altered. If so,
        # this solve method must be changed too
        for val in self.regulatory_interactions_gen():

            state_map = {arg: state[arg] if arg in state else self.initial_state[arg]
                         for arg in val.variables}

            sol = val.evaluate(state_map)
            if sol is None:
                continue

            if not silently:
                val.target.expression_coef = sol
            res[val.id] = sol

        return res

    def next_state(self, state):

        """
        Retrieves the next state for a given state. This is mandatory for now.
        In the future, if no state is provided the current state will be used

        Solves the boolean regulatory model using method solve_regulatory_model(state) and decodes the metabolic state
        for the given state


        :param state: dict, dict, key is the id of the regulatory variable while value can be 0, 1 or float (reactions and metabolites predicates)
        :return: dict, for convenience returns dict of all regulatory variables where key is the id of the regulatory variable and value can be 0, 1 or float (reactions and metabolites predicates)
        """

        # Regulatory Rules from the boolean network
        self.solve_regulatory_model(state)

        # Reactions and metabolites from active genes, reactions and metabolites in the initial state
        decode = self.simulate_cbm_model(constraints=self.decode_metabolic_state(state))
        self.__update_reactions_state__(decode)

        # IF THE UPDATE METABOLITES STATE IS CHANGED FOR ACCOUNTING DECODED SOLUTIONS SUCH AS THE REACTIONS,
        # UNCOMMENT THIS!!! Otherwise, just use it once in the simulate
        # self.__update_metabolites_state__()

        return self.get_current_state()

    def __simulate__(self):

        self._regulatory_solution = []

        # For safety, we want to store a new dict and not the pointers
        new_state = {key: val for key, val in self.initial_state_gen()}
        self.regulatory_solution.append(new_state)

        # Perhaps this will change. See next state
        self.__update_metabolites_state__()

        # solve using the initial state
        new_state = self.next_state(new_state)
        self.regulatory_solution.append(new_state)
        return new_state

    def simulate(self, objective=None, maximize=True, method=None):

        """
        RFBA model one-step simulation (pseudo steady-state).

        First, the boolean regulatory model is solved for the initial state. If the initial state has not yet been
        set, the state of all regulatory variables is considered zero.

        At the same time, the metabolic state is also decoded from the initial state of all metabolic genes (once
        again, the metabolic genes state is zero if the initial state is not altered)

        The result of the boolean regulatory model updates the state of all targets, while the result of the
        metabolic state decoding updates the state of all reactions and metabolites that are also regulatory variables.

        Then, this new state is used to decode the final metabolic state using the resulting metabolic-regulatory state.

        Objective and constraint-based model analysis method (pFBA, FBA, ...) can be altered here reversibly.

        :param objective: None or dict, objective for the optimization. None for the current objective
        :param maximize: bool, direction for the optimization.
        :param method: None or SimulationMethod, constraint-based model analysis method
        :return: SimulationResult object, it contains the metabolic and regulatory states too.
        
        
        Solutions are also stored under solution property
        """

        self.objective = objective
        self.maximize = maximize
        self.cbm_simulation_method = method

        self.__simulate__()

        fba_sol = self.simulate_cbm_model(constraints=self.decode_metabolic_state(self.regulatory_solution[-1]))
        fba_sol.regulatory_solution = self.regulatory_solution[-1]
        self._solution = fba_sol

        self.__populate_solution__()

        return self._solution

    def dynamic_simulate(self, objective=None, maximize=True, method=None, iter_limit=1000):

        """
        RFBA model dynamic simulation (until the metabolic-regulatory state is reached).

        First, the boolean regulatory model is solved for the initial state. If the initial state has not yet been
        set, the state of all regulatory variables is considered zero.

        At the same time, the metabolic state is also decoded from the initial state of all metabolic genes (once
        again, the metabolic genes state is zero if the initial state is not altered)

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

        

        :param objective: None or dict, objective for the optimization. None for the current objective
        :param maximize: bool, direction for the optimization.
        :param method: None or SimulationMethod, constraint-based model analysis method
        :param iter_limit: int, iteration limit
        :return: list, list of SimulationResult objects containing the metabolic and regulatory states too.
        
        Solutions are also stored under solution property.
        
        """

        self.objective = objective
        self.maximize = maximize
        self.cbm_simulation_method = method

        new_state = self.__simulate__()

        i = 1
        sol = 1
        steady_state = False

        while not steady_state:

            # Updating state upon state. See next state for further detail
            new_state = self.next_state(new_state)
            self.regulatory_solution.append(new_state)

            n_sols = len(self.regulatory_solution) - 1

            for sol in range(n_sols):

                # A reaction can have different flux values between two solutions, though the state remains the same.
                # Thus, non-zero flux stands for ON state while zero flux stands for OFF state

                replica_v = []

                for key, val in self.regulatory_solution[sol].items():

                    sol_val = self.regulatory_solution[n_sols][key]

                    if key in self.metabolic_regulatory_reactions or key in self.metabolic_regulatory_metabolites:

                        if abs(val) > 1E-10:
                            val = True
                        else:
                            val = False

                        if abs(sol_val) > 1E-10:
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
        self._solution = []
        n_sols = len(self.regulatory_solution)
        for i in range(sol, n_sols):
            fba_sol = self.simulate_cbm_model(constraints=self.decode_metabolic_state(self.regulatory_solution[i]))
            fba_sol.regulatory_solution = self.regulatory_solution[i]
            self.solution.append(fba_sol)

        self.__populate_solution__()

        return self.solution
