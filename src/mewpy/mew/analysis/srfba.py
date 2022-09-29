from functools import partial
from typing import Union, Dict, TYPE_CHECKING

from mewpy.util.constants import ModelConstants

from mewpy.mew.analysis import FBA
from mewpy.mew.lp import ConstraintContainer, VariableContainer, concat_constraints, integer_coefficients
from mewpy.mew.solution import ModelSolution
from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel
from mewpy.solvers import Solution
from mewpy.solvers.solver import Solver, VarType

if TYPE_CHECKING:
    from mewpy.mew.variables import Reaction, Interaction


class SRFBA(FBA):

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = False,
                 attach: bool = False):
        """
        Steady-state Regulatory Flux Balance Analysis (SRFBA) of a metabolic-regulatory model.
        Implementation of a steady-state version of SRFBA for an integrated metabolic-regulatory model.

        For more details consult Shlomi et al. 2007 at https://dx.doi.org/10.1038%2Fmsb4100141

        :param model: a mewpy Model, MetabolicModel, RegulatoryModel or all. The model is used to retrieve
        variables and constraints to the linear problem
        :param solver: A Solver, CplexSolver, GurobiSolver or OptLangSolver instance.
        Alternatively, the name of the solver is also accepted.
        The solver interface will be used to load and solve a linear problem in a given solver.
        If none, a new solver is instantiated. An instantiated solver may be used, but it will be overwritten
        if build is true.
        :param build: Whether to build the linear problem upon instantiation. Default: False
        :param attach: Whether to attach the linear problem to the model upon instantiation. Default: False
        """
        super().__init__(model=model, solver=solver, build=build, attach=attach)
        self._model_default_lb = ModelConstants.REACTION_LOWER_BOUND
        self._model_default_ub = ModelConstants.REACTION_UPPER_BOUND

    @property
    def model_default_lb(self) -> float:
        """
        The default lower bound for the model reactions.
        :return:
        """
        if self.synchronized:
            return self._model_default_lb

        self._model_default_lb = min(reaction.lower_bound for reaction in self.model.yield_reactions())
        return self._model_default_lb

    @property
    def model_default_ub(self) -> float:
        """
        The default upper bound for the model reactions.
        :return:
        """
        if self.synchronized:
            return self._model_default_ub

        self._model_default_ub = max(reaction.upper_bound for reaction in self.model.yield_reactions())
        return self._model_default_ub

    def gpr_constraint(self, reaction: 'Reaction'):
        """
        It creates a constraint for a given GPR where variables are genes and constraints are designed for
        each reaction.
        A linear GPR follows a boolean algebra-based model.

        The GPR is translated as a set of variables and constraints using the method `linearize_expression`.

        First, it creates the relation between reaction boolean value and the reaction constrains
        For that, the following reactions are added
        V - Y*Vmax <= 0
        V - Y*Vmin => 0
        where V is the reaction in S matrix
        where Y is the reaction boolean variable
        It adds reaction boolean variable and the constraint

        Then,
        Constraints are created using the multiple methods for each operator. Consult `and_constraint`, `or_constraint`,
        `not_constraint`, `greater_constraint`, `less_constraint`, `equal_constraint`, `true_constraint`,
        `false_constraint`, `symbol_constraint` for more information.

        Variables are created using the multiple methods for each operator. Consult `and_variable`, `or_variable`,
        `not_variable`, `greater_variable`, `less_variable`, `equal_variable`, `true_variable`,
        `false_variable`, `symbol_variable` for more information.

        :param reaction: the reaction
        :return: gpr variables and constraints
        """
        boolean_variable = f'bool_{reaction.id}'

        variables = [VariableContainer(name=boolean_variable,
                                       sub_variables=[boolean_variable],
                                       lbs=[0.0],
                                       ubs=[1.0],
                                       variables_type=[VarType.INTEGER])]

        lb, ub = reaction.bounds

        coefs = [{reaction.id: 1.0, boolean_variable: -float(ub)},
                 {reaction.id: 1.0, boolean_variable: -float(lb)}]
        lbs = [self.model_default_lb - float(ub), 0.0]
        ubs = [0.0, self.model_default_ub - float(lb)]

        constraints = [ConstraintContainer(name=None,
                                           coefs=coefs,
                                           lbs=lbs,
                                           ubs=ubs)]

        expression_variables, expression_cnt = self.linearize_expression(boolean_variable=boolean_variable,
                                                                         symbolic=reaction.gpr.symbolic)

        variables.extend(expression_variables)
        constraints.extend(expression_cnt)

        constraints = [concat_constraints(constraints=constraints, name=reaction.id)]
        return variables, constraints

    def interaction_constraint(self, interaction: 'Interaction'):
        """
        It creates a constraint for a given interaction where variables are regulators and constraints are designed for
        each target.
        A linear interaction follows a boolean algebra-based model.

        The interaction is translated as a set of variables and constraints using the method `linearize_expression`.

        Constraints are created using the multiple methods for each operator. Consult `and_constraint`, `or_constraint`,
        `not_constraint`, `greater_constraint`, `less_constraint`, `equal_constraint`, `true_constraint`,
        `false_constraint`, `symbol_constraint` for more information.

        Variables are created using the multiple methods for each operator. Consult `and_variable`, `or_variable`,
        `not_variable`, `greater_variable`, `less_variable`, `equal_variable`, `true_variable`,
        `false_variable`, `symbol_variable` for more information.

        :param interaction: the interaction
        :return: interaction variables and constraints
        """
        symbolic = None
        for coefficient, expression in interaction.regulatory_events.items():

            if coefficient > 0.0:
                symbolic = expression.symbolic

        if symbolic is None:
            return [], []

        lb = float(min(interaction.target.coefficients))
        ub = float(max(interaction.target.coefficients))

        if (lb, ub) in integer_coefficients:
            v_type = VarType.INTEGER
        else:
            v_type = VarType.CONTINUOUS

        variables = [VariableContainer(name=interaction.target.id,
                                       sub_variables=[interaction.target.id],
                                       lbs=[lb],
                                       ubs=[ub],
                                       variables_type=[v_type])]
        constraints = []

        expression_variables, expression_cnt = self.linearize_expression(boolean_variable=interaction.target.id,
                                                                         symbolic=symbolic)

        variables.extend(expression_variables)
        constraints.extend(expression_cnt)

        constraints = [concat_constraints(constraints=constraints, name=interaction.target.id)]
        return variables, constraints

    @staticmethod
    def and_constraint(symbolic):
        """
        Following Boolean algebra, an And (a = b and c) can be translated as: a = b*c
        Alternatively, an And can be written as lb < b + c - a < ub

        So, for the midterm expression a = b and c
        We have therefore the equation: -1 <= 2*b + 2*c – 4*a <= 3
        :param symbolic: the symbolic expression
        :return: a constraint
        """
        name = symbolic.key()
        names = [f'{name}_{i}' for i, _ in enumerate(symbolic.variables[:-1])]

        and_op = names[0]
        op_l = symbolic.variables[0]
        op_r = symbolic.variables[1]

        _coefs = []
        _lbs = []
        _ubs = []

        if op_l.is_one or op_l.is_true:

            _coef = {and_op: -4.0, op_r.key(): 2.0}
            _state = (-3.0, 1.0)

        elif op_r.is_one or op_r.is_true:

            _coef = {and_op: -4.0, op_l.key(): 2.0}
            _state = (-3.0, 1.0)

        elif op_l.is_zero or op_l.is_false:

            _coef = {and_op: -4.0, op_r.key(): 2.0}
            _state = (-1.0, 3.0)

        elif op_r.is_zero or op_r.is_false:

            _coef = {and_op: -4.0, op_l.key(): 2.0}
            _state = (-1.0, 3.0)

        else:

            _coef = {and_op: -4.0, op_l.key(): 2.0, op_r.key(): 2.0}
            _state = (-1.0, 3.0)

        _coefs.append(_coef)
        _lbs.append(_state[0])
        _ubs.append(_state[1])

        children = []

        if len(symbolic.variables) > 2:
            children = symbolic.variables[2:]

        # building a nested And subexpression
        for i, op_r in enumerate(children):

            op_l = names[i]
            and_op = names[i + 1]

            if op_r.is_one or op_r.is_true:

                _coef = {and_op: -4.0, op_l: 2.0}
                _state = (-3.0, 1.0)

            elif op_r.is_zero or op_r.is_false:

                _coef = {and_op: -4.0, op_l: 2.0}
                _state = (-1.0, 3.0)

            else:
                _coef = {and_op: -4.0, op_l: 2.0, op_r.key(): 2.0}
                _state = (-1.0, 3.0)

            _coefs.append(_coef)
            _lbs.append(_state[0])
            _ubs.append(_state[1])

        return ConstraintContainer(name=None, coefs=_coefs, lbs=_lbs, ubs=_ubs)

    @staticmethod
    def or_constraint(symbolic):
        """
        Following Boolean algebra, an Or (a = b or c) can be translated as: a = b + c - b*c
        Alternatively, an Or can be written as lb < b + c - a < ub

        So, for the midterm expression a = b or c
        We have therefore the equation: -2 <= 2*b + 2*c – 4*a <= 1
        :param symbolic: the symbolic expression
        :return: a constraint
        """
        name = symbolic.key()
        names = [f'{name}_{i}' for i, _ in enumerate(symbolic.variables[:-1])]

        or_op = names[0]
        op_l = symbolic.variables[0]
        op_r = symbolic.variables[1]

        _coefs = []
        _lbs = []
        _ubs = []

        if op_l.is_one or op_l.is_true:

            _coef = {or_op: -4.0, op_r.key(): 2.0}
            _state = (-4.0, -1.0)

        elif op_r.is_one or op_r.is_true:

            _coef = {or_op: -4.0, op_l.key(): 2.0}
            _state = (-4.0, -1.0)

        elif op_l.is_zero or op_l.is_false:

            _coef = {or_op: -4.0, op_r.key(): 2.0}
            _state = (-2.0, 0.0)

        elif op_r.is_zero or op_r.is_false:

            _coef = {or_op: -4.0, op_l.key(): 2.0}
            _state = (-2.0, 0.0)

        else:

            _coef = {or_op: -4.0, op_l.key(): 2.0, op_r.key(): 2.0}
            _state = (-2.0, 1.0)

        _coefs.append(_coef)
        _lbs.append(_state[0])
        _ubs.append(_state[1])

        children = []

        if len(symbolic.variables) > 2:
            children = symbolic.variables[2:]

        # building a nested Or subexpression
        for i, op_r in enumerate(children):

            op_l = names[i]
            or_op = names[i + 1]

            if op_r.is_one or op_r.is_true:

                _coef = {or_op: -4.0, op_l: 2.0}
                _state = (-4.0, -1.0)

            elif op_r.is_zero or op_r.is_false:

                _coef = {or_op: -4.0, op_l: 2.0}
                _state = (-2.0, 1.0)

            else:

                _coef = {or_op: -4.0, op_l: 2.0, op_r.key(): 2.0}
                _state = (-2.0, 1.0)

            _coefs.append(_coef)
            _lbs.append(_state[0])
            _ubs.append(_state[1])

        return ConstraintContainer(name=None, coefs=_coefs, lbs=_lbs, ubs=_ubs)

    @staticmethod
    def not_constraint(symbolic):
        """
        Following Boolean algebra, a Not (a = not b) can be translated as: a = 1 - b
        Alternatively, a Not can be written as a + b = 1

        So, for the midterm expression a = not b
        We have therefore the equation: 1 < a + b < 1
        :param symbolic: the symbolic expression
        :return: a constraint
        """
        op_l = symbolic.variables[0]

        # Not right operators
        if op_l.is_numeric:

            _coef = {symbolic.key(): 1.0}
            _state = (float(op_l.value), float(op_l.value))

        else:

            # add Not row and set mip bounds to 1;1
            _coef = {symbolic.key(): 1.0, op_l.key(): 1.0}
            _state = (1.0, 1.0)

        return ConstraintContainer(name=None, coefs=[_coef], lbs=[_state[0]], ubs=[_state[1]])

    def greater_constraint(self, symbolic):
        """
        Following Propositional logic, a predicate (a => r > value) can be translated as: r - value > 0
        Alternatively, flux predicate a => r>value can be a(value + tolerance - r_UB) + r < value + tolerance
        & a(r_LB - value - tolerance) + r > r_LB

        So, for the midterm expression a => r > value
        We have therefore the equations: a(value + tolerance - r_UB) + r < value + tolerance
        & a(r_LB - value - tolerance) + r > r_LB
        :param symbolic: the symbolic expression
        :return: a constraint
        """

        greater_op = symbolic.key()
        op_l = symbolic.variables[0]
        op_r = symbolic.variables[1]

        if op_l.is_numeric:
            operand = op_r
            c_val = float(op_l.value)

        else:
            operand = op_l
            c_val = float(op_r.value)

        _lb, _ub = operand.bounds

        _lb = float(_lb)
        _ub = float(_ub)
        # add Greater row (a(value + tolerance - r_UB) + r < value + tolerance) and set mip bounds to
        # -inf;comparison_val
        # add Greater row (a(r_LB - value - tolerance) + r > r_LB) and set mip bounds to lb;inf
        _coefs = [
            {greater_op: c_val + ModelConstants.TOLERANCE - _ub, operand.key(): 1.0},
            {greater_op: _lb - c_val - ModelConstants.TOLERANCE, operand.key(): 1.0}
        ]

        _lbs = [self.model_default_lb, _lb]
        _ubs = [c_val + ModelConstants.TOLERANCE, self.model_default_ub]

        return ConstraintContainer(name=None, coefs=_coefs, lbs=_lbs, ubs=_ubs)

    def less_constraint(self, symbolic):
        """
        Following Propositional logic, a| predicate (a => r < value) can be translated as: r - value < 0
        Alternatively, flux predicate a => r<value can be a(value + tolerance - r_LB) + r > value + tolerance
        & a(r_UB - value - tolerance) + r < r_UB

        So, for the midterm expression a => r > value
        We have therefore the equations: a(value + tolerance - r_LB) + r > value + tolerance
        & a(r_UB - value - tolerance) + r < r_UB
        :param symbolic: the symbolic expression
        :return: a constraint
        """
        less_op = symbolic.key()
        op_l = symbolic.variables[0]
        op_r = symbolic.variables[1]

        if op_l.is_numeric:
            operand = op_r
            c_val = float(op_l.value)

        else:
            operand = op_l
            c_val = float(op_r.value)

        _lb, _ub = operand.bounds
        _lb = float(_lb)
        _ub = float(_ub)
        # add Less row (a(value + tolerance - r_LB) + r > value + tolerance) and set mip bounds to
        # -inf;-comparison_val
        # add Less row (a(r_UB - value - tolerance) + r < r_UB) and set mip bounds to lb;inf
        _coefs = [
            {less_op: c_val + ModelConstants.TOLERANCE - _lb, operand.key(): 1.0},
            {less_op: _ub - c_val - ModelConstants.TOLERANCE, operand.key(): 1.0}
        ]
        _lbs = [c_val + ModelConstants.TOLERANCE, self.model_default_lb]
        _ubs = [self.model_default_ub, _ub]

        return ConstraintContainer(name=None, coefs=_coefs, lbs=_lbs, ubs=_ubs)

    @staticmethod
    def equal_constraint(symbolic):
        """
        Following Propositional logic, a predicate (a => r = value) can be translated as: r - value = 0
        :param symbolic: the symbolic expression
        :return: a constraint
        """
        less_op = symbolic.key()
        op_l = symbolic.variables[0]
        op_r = symbolic.variables[1]

        if op_l.is_numeric:
            operand = op_r
            c_val = float(op_l.value)

        else:
            operand = op_l
            c_val = float(op_r.value)

        _coefs = [{less_op: - c_val, operand.key(): 1.0}]
        _lbs = [0]
        _ubs = [0]

        return ConstraintContainer(name=None, coefs=_coefs, lbs=_lbs, ubs=_ubs)

    @staticmethod
    def none_constraint(_):
        """
        The target or reaction can take any boolean value (0;1),
        so a constraint with bounds to 0;1 is added to the problem
        :param _:
        :return:
        """
        return ConstraintContainer(name=None, coefs=[{}], lbs=[0.0], ubs=[1.0])

    @staticmethod
    def false_constraint(_):
        """
        Constraint with 0;0 bounds
        :param _:
        :return:
        """
        return ConstraintContainer(name=None, coefs=[{}], lbs=[0.0], ubs=[0.0])

    @staticmethod
    def true_constraint(_):
        """
        Constraint with 1;1 bounds
        :param _:
        :return:
        """
        return ConstraintContainer(name=None, coefs=[{}], lbs=[1.0], ubs=[1.0])

    @staticmethod
    def number_constraint(symbolic):
        """
        Constraint with number;number bounds
        :param symbolic:
        :return:
        """
        return ConstraintContainer(name=None, coefs=[{}], lbs=[float(symbolic.value)], ubs=[float(symbolic.value)])

    @staticmethod
    def symbol_constraint(symbolic):
        """
        Constraint with -1*symbol and 0;0 bounds
        :param symbolic:
        :return:
        """
        return ConstraintContainer(name=None, coefs=[{symbolic.key(): -1.0}], lbs=[0.0], ubs=[0.0])

    def get_lp_constraint(self,
                          symbolic,
                          operators=True,
                          bool_atoms=True,
                          numeric_atoms=True,
                          symbolic_atoms=True,
                          empty_symbolic=True,
                          ):
        """
        Get the constraint corresponding to the symbolic expression
        :param symbolic: the symbolic expression
        :param operators: if True, the operators constraints are returned
        :param bool_atoms: if True, the boolean atoms constraints are returned
        :param numeric_atoms: if True, the numeric atoms constraints are returned
        :param symbolic_atoms: if True, the symbolic atoms constraints are returned
        :param empty_symbolic: if True, the empty symbolic constraints are returned
        :return: a constraint or None if there is no constraint for the given symbolic expression
        """
        if operators:

            if symbolic.is_and:
                return partial(self.and_constraint, symbolic)

            elif symbolic.is_or:
                return partial(self.or_constraint, symbolic)

            elif symbolic.is_not:
                return partial(self.not_constraint, symbolic)

            elif symbolic.is_greater or symbolic.is_greater_equal:
                return partial(self.greater_constraint, symbolic)

            elif symbolic.is_less or symbolic.is_less_equal:
                return partial(self.less_constraint, symbolic)

            elif symbolic.is_equal:
                return partial(self.equal_constraint, symbolic)

        if bool_atoms:
            if symbolic.is_true:

                return partial(self.true_constraint, symbolic)

            elif symbolic.is_false:

                return partial(self.false_constraint, symbolic)

        if numeric_atoms:
            if symbolic.is_numeric:
                return partial(self.number_constraint, symbolic)

        if symbolic_atoms:
            if symbolic.is_symbol:
                return partial(self.symbol_constraint, symbolic)

        if empty_symbolic:
            if symbolic.is_none:
                return partial(self.none_constraint, symbolic)

        return

    @staticmethod
    def _variable_operator(symbolic):
        name = symbolic.key()
        names = []
        lbs = []
        ubs = []
        var_types = []

        for i, _ in enumerate(symbolic.variables[:-1]):
            names.append(f'{name}_{i}')
            lbs.append(0.0)
            ubs.append(1.0)
            var_types.append(VarType.INTEGER)

        return VariableContainer(name=name, sub_variables=names, lbs=lbs, ubs=ubs, variables_type=var_types)

    def and_variable(self, symbolic):
        """
        The AND operator is translated as a variable with 0;1 bounds
        :param symbolic: the symbolic expression
        :return: a variable
        """
        return self._variable_operator(symbolic=symbolic)

    def or_variable(self, symbolic):
        """
        The OR operator is translated as a variable with 0;1 bounds
        :param symbolic: the symbolic expression
        :return: a variable
        """
        return self._variable_operator(symbolic=symbolic)

    @staticmethod
    def not_variable(symbolic):
        """
        The NOT operator is translated as a variable with 0;1 bounds
        :param symbolic: the symbolic expression
        :return: a variable
        """
        name = symbolic.key()
        sub_variable_name = f'{name}_0'

        return VariableContainer(name=symbolic.key(),
                                 sub_variables=[sub_variable_name],
                                 lbs=[0.0],
                                 ubs=[1.0],
                                 variables_type=[VarType.INTEGER])

    def greater_variable(self, symbolic):
        """
        The greater operator is translated as a variable with 0;1 bounds
        :param symbolic: the symbolic expression
        :return: a variable
        """
        return self._variable_operator(symbolic=symbolic)

    def less_variable(self, symbolic):
        """
        The less operator is translated as a variable with 0;1 bounds
        :param symbolic: the symbolic expression
        :return: a variable
        """
        return self._variable_operator(symbolic=symbolic)

    def equal_variable(self, symbolic):
        """
        The equal operator is translated as a variable with 0;1 bounds
        :param symbolic: the symbolic expression
        :return: a variable
        """
        return self._variable_operator(symbolic=symbolic)

    @staticmethod
    def symbol_variable(symbolic):
        """
        The symbol is translated as a variable with alpha < symbol < beta bounds
        :param symbolic: the symbolic expression
        :return: a variable
        """
        name = symbolic.key()
        lb, ub = symbolic.bounds
        lb = float(lb)
        ub = float(ub)

        if (lb, ub) in integer_coefficients:
            v_type = VarType.INTEGER
        else:
            v_type = VarType.CONTINUOUS

        return VariableContainer(name=name, sub_variables=[name], lbs=[lb], ubs=[ub], variables_type=[v_type])

    def get_lp_variable(self, symbolic):
        """
        Get the variable corresponding to the symbolic expression
        :param symbolic: the symbolic expression
        :return: a variable or None if there is no variable for the given symbolic expression
        """

        if symbolic.is_and:
            return partial(self.and_variable, symbolic)

        elif symbolic.is_or:
            return partial(self.or_variable, symbolic)

        elif symbolic.is_not:
            return partial(self.not_variable, symbolic)

        elif symbolic.is_greater or symbolic.is_greater_equal:
            return partial(self.greater_variable, symbolic)

        elif symbolic.is_less or symbolic.is_less_equal:
            return partial(self.less_variable, symbolic)

        elif symbolic.is_equal:
            return partial(self.equal_variable, symbolic)

        elif symbolic.is_symbol:
            return partial(self.symbol_variable, symbolic)

        return

    def linearize_atomic_expression(self, boolean_variable, symbolic):
        """
        It builds the variables and constraints corresponding to the linearization of the atomic expression
        :param boolean_variable: the boolean variable corresponding to the atomic expression
        :param symbolic: the symbolic expression
        :return: a list of variables and a list of constraints
        """
        variables = []
        constraints = []

        if symbolic.is_symbol:
            var = self.symbol_variable(symbolic=symbolic)
            variables.append(var)

        linearizer = self.get_lp_constraint(symbolic, operators=False)

        expression_cnt = linearizer()

        expression_cnt.coefs[0][boolean_variable] = 1.0

        constraints.append(expression_cnt)

        return variables, constraints

    def linearize_complex_expression(self, boolean_variable, symbolic):
        """
        It builds the variables and constraints corresponding to the linearization of the complex expression

        It iterates an expression in the reverse order
        For example, expr = A and (B or C) yields the following elements:
            - C
            - B
            - A
            - (B or C)
            - A and (B or C)
        For each element, the linearization is performed and the resulting variables and constraints
        are added to the list. To see which variables and constraints are added for each element,
        see the constraint methods:
            - `and_constraint`
            - `or_constraint`
            - `not_constraint`
            - `greater_constraint`
            - `less_constraint`
            - `equal_constraint`
            - `symbol_constraint`
            - `none_constraint`

        and the variable methods:
            - `and_variable`
            - `or_variable`
            - `not_variable`
            - `greater_variable`
            - `less_variable`
            - `equal_variable`
            - `symbol_variable`

        :param boolean_variable: the boolean variable corresponding to the complex expression
        :param symbolic: the symbolic expression
        :return: a list of variables and a list of constraints
        """
        variables = []
        constraints = []
        last_variable = None
        for atom in symbolic:

            # An operator expression will be decomposed into multiple operator expressions, namely into a nested
            # expression. For instance, an A & B & C & D will become And(D, And(C, And(A, B))). Each nested operator
            # expression will be a column in the matrix and then linked in the columns. However, this operator
            # expression will be the same for all maters and have the same column identifier regardless of the length.
            # Thus, all operator columns (variables) will be under the columns linked list engine. When retrieving
            # the indexes of the operator expression, the last index
            # of the slice should be used to get the last real column, as the result of this column is the one that
            # really matters. The hashes of the columns for the simulation engine should be the row name plus
            # str of the operator
            last_variable = atom

            lp_variable = self.get_lp_variable(symbolic=atom)

            if lp_variable is not None:
                var = lp_variable()
                variables.append(var)

            lp_constraint = self.get_lp_constraint(atom,
                                                   bool_atoms=False,
                                                   numeric_atoms=False,
                                                   symbolic_atoms=False,
                                                   empty_symbolic=False)

            if lp_constraint is not None:
                constraint = lp_constraint()
                constraints.append(constraint)

        # identifying the last index to link the outcome of this variable to the boolean variable associated to the
        # expression
        last_variable_name = last_variable.key()
        names = [f'{last_variable_name}_{i}' for i, _ in enumerate(last_variable.variables[:-1])]

        if names:
            last_variable_name = names[-1]
        else:
            last_variable_name = last_variable_name

        # add gene row which means that the gene variable in the mip matrix is associated with the last midterm
        # expression, namely the whole expression
        # set mip bounds to 0;0
        expression_cnt = ConstraintContainer(name=None,
                                             coefs=[{boolean_variable: 1.0, last_variable_name: -1.0}],
                                             lbs=[0.0],
                                             ubs=[0.0])
        constraints.append(expression_cnt)

        return variables, constraints

    def linearize_expression(self, boolean_variable, symbolic):
        """
        It builds the variables and constraints corresponding to the linearization of the expression.

        It iterates an expression in the reverse order
        For example, expr = A and (B or C) yields the following elements:
            - C
            - B
            - A
            - (B or C)
            - A and (B or C)
        For each element, the linearization is performed and the resulting variables and constraints
        are added to the list. To see which variables and constraints are added for each element,
        see the constraint methods:
            - `and_constraint`
            - `or_constraint`
            - `not_constraint`
            - `greater_constraint`
            - `less_constraint`
            - `equal_constraint`
            - `symbol_constraint`
            - `none_constraint`

        and the variable methods:
            - `and_variable`
            - `or_variable`
            - `not_variable`
            - `greater_variable`
            - `less_variable`
            - `equal_variable`
            - `symbol_variable`
        :param boolean_variable: the boolean variable corresponding to the expression
        :param symbolic: the symbolic expression
        :return: a list of variables and a list of constraints
        """
        # if expression is atom and defines a variable always On or Off, add an On/Off row
        if symbolic.is_atom:
            return self.linearize_atomic_expression(boolean_variable=boolean_variable, symbolic=symbolic)

        return self.linearize_complex_expression(boolean_variable=boolean_variable, symbolic=symbolic)

    def _build_interactions(self):
        """
        It builds the algebraic constraints for SRFBA, namely the GPR constraints and the interaction constraints.
        It is called automatically upon instantiation if build is True.
        :return:
        """
        variables = []
        constraints = []
        for interaction in self.model.yield_interactions():
            interaction_variables, interaction_constraints = self.interaction_constraint(interaction)
            variables.extend(interaction_variables)
            constraints.extend(interaction_constraints)

        self.add_variables(*variables)
        self.add_constraints(*constraints)

    def _build_gprs(self):
        """
        It builds the algebraic constraints for SRFBA, namely the GPR constraints and the interaction constraints.
        It is called automatically upon instantiation if build is True.
        :return:
        """
        variables = []
        constraints = []

        for reaction in self.model.yield_reactions():
            gpr_variables, gpr_constraints = self.gpr_constraint(reaction)
            variables.extend(gpr_variables)
            constraints.extend(gpr_constraints)

        self.add_variables(*variables)
        self.add_constraints(*constraints)

    def _build(self):
        """
        It builds the linear problem for SRFBA. It is called automatically upon instantiation if build is True.
        The SRFBA problem is a mixed-integer linear problem (MILP) with the following structure:
            - metabolic constraints
            - GPR constraints
            - interaction constraints

        :return:
        """
        if self.model.is_metabolic() and self.model.is_regulatory():
            self._build_mass_constraints()
            self._build_gprs()
            self._build_interactions()

            self._linear_objective = {var.id: value for var, value in self.model.objective.items()}
            self._minimize = False

    def _optimize(self,
                  to_solver: bool = False,
                  solver_kwargs: Dict = None,
                  initial_state: Dict[str, float] = None,
                  **kwargs) -> Union[ModelSolution, Solution]:
        """
        It solves the linear problem. The linear problem is solved using the solver interface.

        The optimize method allows setting temporary changes to the linear problem. The changes are
        applied to the linear problem reverted to the original state afterward.
        Objective, constraints and solver parameters can be set temporarily.

        The solution is returned as a ModelSolution instance, unless to_solver is True. In this case,
        the solution is returned as a SolverSolution instance.

        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False
        :param solver_kwargs: A dictionary of solver parameters to be set temporarily. Default: None
        :param initial_state: a dictionary of variable ids and their values to set as initial state
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

        constraints = {**constraints, **initial_state}
        solver_kwargs['constraints'] = constraints

        solution = self.solver.solve(**solver_kwargs)
        return solution
