import numpy as np
from sympy import postorder_traversal
from sympy.logic.boolalg import And, Or, Not, BooleanFalse, BooleanTrue
from sympy.core.relational import StrictGreaterThan, StrictLessThan, GreaterThan, LessThan
from sympy.core.numbers import Zero, One
from cobamp.core.linear_systems import GenericLinearSystem, VAR_CONTINUOUS, VAR_BINARY
from cobamp.core.optimization import LinearSystemOptimizer
from mewpy.simulation.simulation import SimulationResult
from mewpy.regulation import IntegratedModel

_SRFBA_TOL = 1E-10


class SRFBAModel(IntegratedModel):
    """

    SRFBAModel class object inherits everything from the IntegratedModel class.
    This object is the standard for simulating a SRFBA model by solving a boolean regulatory model together with
    the metabolic model as Mixed Integer Linear Problem (MILP) asynchronously
    It provides the simulate method where one can know all about the metabolic and regulatory analysis of SRFBA

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

        self._A = None
        self._lbs = None
        self._ubs = None
        self._a_lb = None
        self._a_ub = None
        self._variables = None
        self._variables_types = {'boolean': [],
                                 'continuous': []}
        self._problem = None

        # integrated model
        super().__init__(identifier,
                         name,
                         cbm_model,
                         cbm_simulation_interface,
                         targets,
                         regulators,
                         regulatory_interactions,
                         initial_state)

    @property
    def A(self):
        return getattr(self, '_A', None)

    @property
    def lbs(self):
        return getattr(self, '_lbs', None)

    @property
    def ubs(self):
        return getattr(self, '_ubs', None)

    @property
    def a_lb(self):
        return getattr(self, '_a_lb', None)

    @property
    def a_ub(self):
        return getattr(self, '_a_ub', None)

    @property
    def variables(self):
        return getattr(self, '_variables', None)

    @property
    def problem(self):
        return getattr(self, '_problem', None)

    def __add_column__(self, variable, lb=0, ub=1, v_type='boolean'):

        # adding column/variable to the problem (return it if already exists)
        # It returns the column index

        if variable in self.variables:
            return self.variables[variable]

        for row in self.A:
            row.append(0)
        self.lbs.append(lb)
        self.ubs.append(ub)
        self.variables[variable] = len(self.A[0]) - 1
        self._variables_types[v_type].append(len(self.A[0]) - 1)

        return len(self.A[0]) - 1

    def __add_row__(self, elements, a_lb=0, a_ub=1):

        # adding row to the problem and return its index in the A matrix

        new_row = [0] * len(self.A[0])
        for i, val in elements:
            new_row[i] = val
        self.A.append(new_row)
        self.a_lb.append(a_lb)
        self.a_ub.append(a_ub)

        return len(self.A) - 1

    def __decode_reaction_bounds__(self, reaction, reaction_bool_variable):

        # Relation between reaction boolean value and the reaction constrains
        # For that, the following reactions must be added
        # V - Y*Vmax < 0
        # V - Y*Vmin > 0
        # where V is the reaction in S matrix
        # where Y is the reaction boolean variable

        _lb, _ub = self.lbs[reaction], self.ubs[reaction]

        self.__add_row__([(reaction, 1), (reaction_bool_variable, -_ub)], -999999, 0)

        self.__add_row__([(reaction, 1), (reaction_bool_variable, -_lb)], 0, 999999)

    def __boolean_rule_decode__(self, variable, expression):

        # Results of sympy's expression parsing:
        # gene = X and 1 => X
        # gene = X or 1 => 1 (sympy's BooleanTrue)
        # gene = X and 0 => 0 (sympy's BooleanFalse)
        # gene = X or 0 => X
        # gene = 1 => 1 (sympy's one)
        # gene = 0 => 0 (sympy's zero)
        # gene = X > 10 => (X, 1 (sympy's number)), namely two args

        # if expression is atom and defines a variable always On or Off, add a On/Off row
        if expression.is_Atom:

            if expression.is_Boolean or expression.is_Number:

                if isinstance(expression, BooleanFalse) or isinstance(expression, Zero):

                    # add row and set mip bounds to 0;0
                    self.__add_row__([(variable, 1)], a_ub=0)
                    return

                elif isinstance(expression, BooleanTrue) or isinstance(expression, One):

                    # add row and set mip bounds to 1;1
                    self.__add_row__([(variable, 1)], a_lb=1)
                    return

                else:
                    print("Either the target result is undetermined or "
                          "something is wrong with expression: {}".format(str(expression)))

                    # add row and set mip bounds to 0;1
                    self.__add_row__([(variable, 1)])
                    return

            # elif expression.is_Number:
            #
            #     if isinstance(expression, Zero):
            #
            #         # add row and set mip bounds to 0;0
            #         self.__add_row__([(variable, 1)], a_ub=0)
            #         return
            #
            #     elif isinstance(expression, One):
            #
            #         # add row and set mip bounds to 1;1
            #         self.__add_row__([(variable, 1)], a_lb=1)
            #         return
            #
            #     else:
            #         print("Either the target result is undetermined or "
            #               "something is wrong with expression: {}".format(str(expression)))
            #
            #         # add row and set mip bounds to 0;1
            #         self.__add_row__([(variable, 1)])
            #         return

            else:

                # Else it must be a symbol or nothing
                # The rest of the function can handle it
                pass

        traversal_results = {}
        last_index = 0

        # postorder_traversal recursively iterates an expression in the reverse order
        # e.g. expr = A and (B or C)
        # list(postorder_traversal(expr)) # [C, B, A, B | C, A & (B | C)]
        for arg in postorder_traversal(expression):

            if isinstance(arg, Or):

                # Following Boolean algebra, an Or (a = b or c) can be translated as: a = b + c - b*c
                # Alternatively, an Or can be written as lb < b + c - a < ub

                # So, for the mid term expression a = b or c
                # We have therefore the equation: -2 <= 2*b + 2*c – 4*a <= 1

                # add Or variable / column
                # or_col_idx = self.__add_column__('mid_term_' + str(self.A.shape[1]))
                or_col_idx = self.__add_column__('mid_term_' + str(len(self.A[0])))

                # Or left and right operators
                op_l = traversal_results[arg.args[0]]
                op_r = traversal_results[arg.args[1]]

                # add Or row and set mip bounds to -2;1
                self.__add_row__([(or_col_idx, -4), (op_l, 2), (op_r, 2)], a_lb=-2, a_ub=1)
                traversal_results[arg] = or_col_idx
                last_index = or_col_idx

                # building a nested Or subexpression
                for sub_arg in range(2, len(arg.args)):
                    # add Or variable / column
                    # or_col_idx = self.__add_column__('mid_term_' + str(self.A.shape[1]))
                    or_col_idx = self.__add_column__('mid_term_' + str(len(self.A[0])))

                    # Or left (last Or) and right operators
                    op_l = traversal_results[arg]
                    op_r = traversal_results[arg.args[sub_arg]]

                    # add Or row and set mip bounds to -2;1
                    self.__add_row__([(or_col_idx, -4), (op_l, 2), (op_r, 2)], a_lb=-2, a_ub=1)
                    traversal_results[arg] = or_col_idx
                    last_index = or_col_idx

            elif isinstance(arg, And):

                # Following Boolean algebra, an And (a = b and c) can be translated as: a = b*c
                # Alternatively, an And can be written as lb < b + c - a < ub

                # So, for the mid term expression a = b and c
                # We have therefore the equation: -1 <= 2*b + 2*c – 4*a <= 3

                # add And variable / column
                # and_col_idx = self.__add_column__('mid_term_' + str(self.A.shape[1]))
                and_col_idx = self.__add_column__('mid_term_' + str(len(self.A[0])))

                # And left and right operators
                op_l = traversal_results[arg.args[0]]
                op_r = traversal_results[arg.args[1]]

                # add And row and set mip bounds to -1;3
                self.__add_row__([(and_col_idx, -4), (op_l, 2), (op_r, 2)], a_lb=-1, a_ub=3)
                traversal_results[arg] = and_col_idx
                last_index = and_col_idx

                # building a nested And subexpression
                for sub_arg in range(2, len(arg.args)):
                    # add And variable / column
                    # and_col_idx = self.__add_column__('mid_term_' + str(self.A.shape[1]))
                    and_col_idx = self.__add_column__('mid_term_' + str(len(self.A[0])))

                    # And left (last And) and right operators
                    op_l = traversal_results[arg]
                    op_r = traversal_results[arg.args[sub_arg]]

                    # add Or row and set mip bounds to -2;1
                    self.__add_row__([(and_col_idx, -4), (op_l, 2), (op_r, 2)], a_lb=-1, a_ub=3)
                    traversal_results[arg] = and_col_idx
                    last_index = and_col_idx

            elif isinstance(arg, Not):

                # Following Boolean algebra, an Not (a = not b) can be translated as: a = 1 - b
                # Alternatively, an Not can be written as a + b = 1

                # So, for the mid term expression a = not b
                # We have therefore the equation: 1 < a + b < 1

                # add Not variable / column
                # not_col_idx = self.__add_column__('mid_term_' + str(self.A.shape[1]))
                not_col_idx = self.__add_column__('mid_term_' + str(len(self.A[0])))

                # Not right operators
                op_r = traversal_results[arg.args[0]]

                # add Not row and set mip bounds to 1;1
                self.__add_row__([(not_col_idx, 1), (op_r, 1)], a_lb=1)
                traversal_results[arg] = not_col_idx
                last_index = not_col_idx

            elif isinstance(arg, StrictGreaterThan) or isinstance(arg, GreaterThan):

                # Following Propositional logic, a predicate (a => r > value) can be translated as: r - value > 0
                # Alternatively, flux predicate a => r>value can be a(value + tolerance - r_UB) + r < value + tolerance
                # & a(r_LB - value - tolerance) + r > r_LB

                # So, for the mid term expression a => r > value
                # We have therefore the equations: a(value + tolerance - r_UB) + r < value + tolerance
                # & a(r_LB - value - tolerance) + r > r_LB

                # In matlab: c * (Vmax - const) + v <= Vmax
                # In matlab: const <= c * (const - Vmin) + v

                # In matlab: compareVal = -compareVal - epsilon
                # In matlab: matrix_values = [vmax - compare_vale, 1] (-inf, vmax)
                # In matlab: matrix_values = [compare_vale - vmin, 1] (compareVal, inf)

                # add Greater variable / column
                # greater_col_idx = self.__add_column__('mid_term_' + str(self.A.shape[1]))
                greater_col_idx = self.__add_column__('mid_term_' + str(len(self.A[0])))

                # Greater left operators
                # op_l can be growth, ph, ...
                # otherwise, it can be a reaction or a boundary reaction associated with a metabolite
                # Either way, they should be already in the matrix

                if arg.args[0].is_Symbol:
                    op_l = traversal_results[arg.args[0]]
                    c_val = arg.args[1]

                else:
                    op_l = traversal_results[arg.args[1]]
                    c_val = arg.args[0]

                _lb, _ub = self.lbs[op_l], self.ubs[op_l]

                # add Greater row (a(value + tolerance - r_UB) + r < value + tolerance) and set mip bounds to
                # -inf;comparison_val
                self.__add_row__([(greater_col_idx, c_val + _SRFBA_TOL - _ub), (op_l, 1)], a_lb=-999999999,
                                 a_ub=c_val + _SRFBA_TOL)
                # self.__add_row__([(greater_col_idx, _ub + c_val), (op_l, 1)], a_lb=-999999, a_ub=_ub)

                # add Greater row (a(r_LB - value - tolerance) + r > r_LB) and set mip bounds to lb;inf
                self.__add_row__([(greater_col_idx, _lb - c_val - _SRFBA_TOL), (op_l, 1)], a_lb=_lb, a_ub=999999999)
                # self.__add_row__([(greater_col_idx, - c_val - _lb), (op_l, 1)], a_lb=-c_val, a_ub=999999)

                traversal_results[arg] = greater_col_idx
                last_index = greater_col_idx

            elif isinstance(arg, StrictLessThan) or isinstance(arg, LessThan):

                # Following Propositional logic, a predicate (a => r < value) can be translated as: r - value < 0
                # Alternatively, flux predicate a => r<value can be a(value + tolerance - r_LB) + r > value + tolerance
                # & a(r_UB - value - tolerance) + r < r_UB

                # So, for the mid term expression a => r > value
                # We have therefore the equations: a(value + tolerance - r_LB) + r > value + tolerance
                # & a(r_UB - value - tolerance) + r < r_UB

                # In matlab: c * (const - Vmax) + v <= const
                # In matlab: Vmin <= c * (Vmin - const) + v

                # In matlab: compareVal = -compareVal + epsilon
                # In matlab: matrix_values = [compareVal-vmax, 1] (-inf, compareVal)
                # In matlab: matrix_values = [vmin-compareVal, 1] (vmin, inf)

                # add Less variable / column
                # less_col_idx = self.__add_column__('mid_term_' + str(self.A.shape[1]))
                less_col_idx = self.__add_column__('mid_term_' + str(len(self.A[0])))

                # Less left operators
                # op_l can be growth, ph, ...
                # otherwise, it can be a reaction or a boundary reaction associated with a metabolite
                # Either way, they should be already in the matrix

                if arg.args[0].is_Symbol:
                    op_l = traversal_results[arg.args[0]]
                    c_val = arg.args[1]

                else:
                    op_l = traversal_results[arg.args[1]]
                    c_val = arg.args[0]

                _lb, _ub = self.lbs[op_l], self.ubs[op_l]

                # add Less row (a(value + tolerance - r_LB) + r > value + tolerance) and set mip bounds to
                # -inf;-comparison_val
                self.__add_row__([(less_col_idx, c_val + _SRFBA_TOL - _lb), (op_l, 1)], a_lb=c_val + _SRFBA_TOL,
                                 a_ub=999999999)
                # self.__add_row__([(less_col_idx, -c_val - _ub), (op_l, 1)], a_lb=-999999, a_ub=-c_val)

                # add Less row (a(r_UB - value - tolerance) + r < r_UB) and set mip bounds to lb;inf
                self.__add_row__([(less_col_idx, _ub - c_val - _SRFBA_TOL), (op_l, 1)], a_lb=-999999999, a_ub=_ub)
                # self.__add_row__([(less_col_idx, _lb + c_val), (op_l, 1)], a_lb=_lb, a_ub=999999)

                traversal_results[arg] = less_col_idx
                last_index = less_col_idx

            else:

                # If the argument is not an And, Or, Not or Relational
                # it must be a symbol or a number from the flux predicates

                # If it is a symbol, it must be added to the matrix
                if arg.is_Symbol:

                    # However, some things should be handle before adding the variable:

                    # special case of ph
                    # ph can be treated as a reaction
                    # as other reactions, the lower and upper bounds can be set to the initial state or vary between
                    # 0 and 14
                    # the row bounds (a_lb and a_ub) must be equal to zero as it is done on the reactions
                    if arg.name.lower() == 'ph':

                        if arg.name in self.initial_state:
                            arg_idx = self.__add_column__(arg.name, self.initial_state[arg.name],
                                                          self.initial_state[arg.name], v_type='continuous')
                            traversal_results[arg] = arg_idx
                            last_index = arg_idx

                        else:
                            arg_idx = self.__add_column__(arg.name, 0, 14, v_type='continuous')
                            traversal_results[arg] = arg_idx
                            last_index = arg_idx

                    # special case of growth (always positive)
                    # growth can be treated as a reaction
                    # as other reactions, the lower and upper bounds can be set to the initial state or vary between
                    # 0.1 and infinite (always positive)
                    # the row bounds (a_lb and a_ub) must be equal to zero as it is done on the reactions
                    elif arg.name.lower() == 'growth':

                        if arg.name in self.initial_state:
                            arg_idx = self.__add_column__(arg.name, self.initial_state[arg.name],
                                                          self.initial_state[arg.name], v_type='continuous')
                            traversal_results[arg] = arg_idx
                            last_index = arg_idx

                        else:
                            arg_idx = self.__add_column__(arg.name, 0 + _SRFBA_TOL, 999999, v_type='continuous')
                            traversal_results[arg] = arg_idx
                            last_index = arg_idx

                    else:

                        var_id = arg.name

                        if var_id in self._aliases_map:
                            var_id = self._aliases_map[var_id]

                        if var_id in self.metabolic_regulatory_reactions:

                            # if the arg is a reaction, the real reaction index should be returned
                            reaction = self._metabolic_regulatory_reactions[var_id].cbm_model_id

                            rxn_idx = self.__add_column__(reaction)
                            traversal_results[arg] = rxn_idx
                            last_index = rxn_idx

                        elif var_id in self.metabolic_regulatory_metabolites:

                            # if the arg is a metabolite, the associated reaction index should be returned

                            reaction = self.cbm_simulation_interface.get_boundary_reaction(
                                self._metabolic_regulatory_metabolites[arg.name].cbm_model_id)

                            if not reaction:
                                # The metabolite does not have a reaction associated with it. This is bad sign though
                                reaction = var_id

                            rxn_idx = self.__add_column__(reaction)
                            traversal_results[arg] = rxn_idx
                            last_index = rxn_idx

                        elif var_id in self.targets:
                            # if the arg is in targets only the column is added, as in the next iterations the row will be
                            # created

                            tg_idx = self.__add_column__(var_id)
                            traversal_results[arg] = tg_idx
                            last_index = tg_idx

                        else:

                            # if the args of a sympy's expression are regulators or genes,
                            # which are not defined in self.targets or self.reactions, they will never be created elsewhere
                            # they must be added to the MILP matrix as columns
                            # The state can be inferred from the initial state if the variables are in the initial state
                            # else, the variables are considered as unknown

                            if var_id in self.initial_state:

                                arg_idx = self.__add_column__(var_id, self.initial_state[var_id],
                                                              self.initial_state[var_id])
                                traversal_results[arg] = arg_idx
                                last_index = arg_idx

                            else:
                                arg_idx = self.__add_column__(var_id)
                                traversal_results[arg] = arg_idx
                                last_index = arg_idx

        # add gene row which means that the gene variable in the mip matrix is associated with the last mid term
        # expression, namely the whole expression
        # set mip bounds to 0;0
        self.__add_row__([(variable, 1), (last_index, -1)], a_ub=0)
        return

    def build(self):

        """
        Method for building SRFBA's MILP problem. It must be ran whenever a change or set of changes (adding a
        regulator, regulatory interaction, ...) is made to the model!
        Otherwise, the latest changes won't be considered for the simulation.

        Simulate method builds the problem when the flag build is used.

        The problem includes:
            - matrix A (rows are odes and columns are variables)
            - ode's lower and upper bounds
            - continuous and boolean variables
            - variables' lower and upper bounds

        A boolean variable is created for each target, regulator, regulatory condition (e.g. pH), metabolic gene,
        reaction and mid-term operators available in the regulatory interactions and GPRs.
        For instance:

        PFL: pyruvate + coa -> formate + acetyl-coa

        R_PFL, G1

        A = [[-1, 0, 0], # (0,0)
            [ -1, 0, 0], # (0,0)
            [  1, 0, 0], # (0,0)
            [  1, 0, 0], # (0,0)
            [  0, 1, -1]] # (0,0)

        lbs, ubs = (0,0,0), (1000, 1, 1)

        The MILP problem is then passed to the cobamp engine for simulating generic linear systems.

        :return:
        """

        # cbm_model S matrix (metabolites, reactions)
        self._A = self.cbm_simulation_interface.get_S()

        self._lbs, self._ubs = self.cbm_simulation_interface.get_bounds()

        self._variables = {r: i for i, r in enumerate(self.reactions)}
        self._variables_types['continuous'] = list(self.variables.values())

        # using list plus append results into significant speed improvements
        self._A = [list(row) for row in self.A]
        self._lbs, self._ubs = list(self.lbs), list(self.ubs)

        self._a_lb = [0] * len(self.A)
        self._a_ub = [0] * len(self.A)

        # regulatory rules over genes
        for reg_inter in self.regulatory_interactions_gen():

            # add target boolean variable / add column / Add MILP bool variable
            # If already exists, only the index is returned
            target_bool_idx = self.__add_column__(reg_inter.target.id)

            # If the expression is empty, it means that the target or reaction can take any boolean value (0;1)
            # add row and set mip bounds to 0;1
            if reg_inter.sympify is None:

                # The initial state can have however the state of this variable
                if reg_inter.target.id in self.initial_state:
                    self.__add_row__([(target_bool_idx, 1)], a_lb=self.initial_state[reg_inter.target.id],
                                     a_ub=self.initial_state[reg_inter.target.id])

                else:

                    self.__add_row__([(target_bool_idx, 1)])

            else:

                self.__boolean_rule_decode__(target_bool_idx, reg_inter.sympify)

        # gprs over reactions
        for i, (rxn, (_, _, sympify, _, _, _)) in enumerate(self._gprs_evaluator.items()):

            # add reaction boolean variable / add column / Add MILP bool variable
            # If already exists, only the index is returned
            reaction_bool_idx = self.__add_column__('B_' + rxn)

            # If the expression is empty, it means that the target or reaction can take any boolean value (0;1)
            # add row and set mip bounds to 0;1
            if sympify is None:

                if rxn in self._aliases_map:
                    rxn = self._aliases_map[rxn]

                # The initial state can have however the state of this variable
                if rxn in self.initial_state:
                    self.__add_row__([(reaction_bool_idx, 1)], a_lb=self.initial_state[rxn],
                                     a_ub=self.initial_state[rxn])

                else:

                    self.__add_row__([(reaction_bool_idx, 1)])

            else:

                self.__boolean_rule_decode__(reaction_bool_idx, sympify)

            self.__decode_reaction_bounds__(i, reaction_bool_idx)

        self._A = np.array(self.A)

        self._problem = GenericLinearSystem(self.A,
                                            VAR_CONTINUOUS,
                                            self.lbs,
                                            self.ubs,
                                            self.a_lb,
                                            self.a_ub,
                                            list(self.variables.keys()))
        # build cobamp problem
        self.problem.build_problem()

        # var types can be continuous or binary, but only the variables after the reactions are binary
        self.problem.set_variable_types([self.problem.model.variables[i] for i in self._variables_types['continuous']],
                                        VAR_CONTINUOUS)
        self.problem.set_variable_types([self.problem.model.variables[i] for i in self._variables_types['boolean']],
                                        VAR_BINARY)

    def simulate(self, build=False, objective=None, maximize=True):

        """

        SRFBA model simulation. Simulation is performed using cobamp's linear system optimizer.

        See build method for further insight over SRFBA's MILP problem and how it is solved.

        :param build: bool, build the milp problem?. False is the default
        :param objective: None or dict, objective for the optimization. None for the current objective
        :param maximize: bool, direction for the optimization.
        :return: SimulationResult object, it contains the metabolic and regulatory states too.
        
        Solutions are also stored under solution property
        """

        if build:
            self.build()

        if not self.A:
            self.build()

        self.objective = objective
        self.maximize = maximize

        lso = LinearSystemOptimizer(self.problem, build=False)

        # objective function
        objective_ids = self.cbm_simulation_interface.get_objective()
        c = np.zeros(self.A.shape[1])
        for var in objective_ids:
            c[self.variables[var]] = 1
        minimize = not self.maximize

        self.problem.set_objective(c, minimize)

        solution = lso.optimize()

        #building mewpy's solution object
        flxs = {}
        bools = {}
        regulatory_conditions = {}

        for key, val in solution.var_values().items():

            if key in self.reactions:
                flxs[key] = val
            elif key in self.regulatory_conditions:
                regulatory_conditions[key] = val
            else:
                bools[key] = val

        solution = SimulationResult(self.cbm_model,
                                    solution.objective_value(),
                                    fluxes=flxs,
                                    status=solution.status(),
                                    envcond=self.environmental_conditions,
                                    model_constraints=self.cbm_simulation_interface.constraints,
                                    simul_constraints={},
                                    maximize=maximize)

        solution.regulatory_solution = bools
        solution.regulatory_conditions = regulatory_conditions

        self._regulatory_solution = bools.update(regulatory_conditions)
        self._solution = solution

        return self.solution