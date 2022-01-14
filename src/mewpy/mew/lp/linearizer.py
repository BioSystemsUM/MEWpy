from abc import abstractmethod
from functools import partial
from typing import List, TYPE_CHECKING

from mewpy.solvers.solver import VarType
from mewpy.util.constants import ModelConstants
from mewpy.util.utilities import iterable

from .linear_problem import LinearProblem
from .linear_containers import ConstraintContainer, VariableContainer, concat_constraints
from .linear_utils import integer_coefficients, bounds_from_symbol, bounds_from_variable
from .notification import Notification

if TYPE_CHECKING:
    from mewpy.mew.algebra import Symbolic
    from mewpy.mew.variables import Reaction, Metabolite

MEWPY_LB = ModelConstants.REACTION_LOWER_BOUND
MEWPY_UB = ModelConstants.REACTION_UPPER_BOUND
MEWPY_TOL = ModelConstants.TOLERANCE


# TODO: missing documentation and typing
class MetabolicLinearizer(LinearProblem):

    @abstractmethod
    def build(self):

        # The concrete implementation is defined by each simulation method, e.g. fba, pfba, etc

        pass

    def notification(self, notification: Notification):

        if notification.content_type == 'reactions' and notification.action == 'add':

            return self.add_reactions(notification.content)

        elif notification.content_type == 'reactions' and notification.action == 'remove':

            return self.remove_reactions(notification.content)

        elif notification.content_type == 'metabolites' and notification.action == 'add':

            return self.add_metabolites(notification.content)

        elif notification.content_type == 'metabolites' and notification.action == 'remove':

            return self.remove_metabolites(notification.content)

        else:
            return super(MetabolicLinearizer, self).notification(notification)

    def metabolite_reaction_lookup(self, reactions: List['Reaction'], metabolites: List['Metabolite']):

        reactions: List['Reaction'] = iterable(reactions)
        metabolites: List['Metabolite'] = iterable(metabolites)

        constraints = {metabolite.id: ConstraintContainer(name=metabolite.id,
                                                          coefs=[{}],
                                                          lbs=[0.0],
                                                          ubs=[0.0])
                       for metabolite in metabolites}

        variables = []
        for reaction in reactions:

            lb, ub = reaction.bounds

            var = VariableContainer(name=reaction.id,
                                    sub_variables=[reaction.id],
                                    lbs=[float(lb)],
                                    ubs=[float(ub)],
                                    variables_type=[VarType.CONTINUOUS])

            variables.append(var)

            for metabolite, coef in reaction.stoichiometry.items():
                metabolite_cnt = constraints[metabolite.id]
                metabolite_cnt.coefs[0][reaction.id] = coef

        self.add(variables)
        self.add(constraints.values())

    def add_reactions(self, reactions: List['Reaction']):

        reactions: List['Reaction'] = iterable(reactions)

        variables = []
        constraints = []
        for reaction in reactions:

            lb, ub = reaction.bounds

            var = VariableContainer(name=reaction.id,
                                    sub_variables=[reaction.id],
                                    lbs=[float(lb)],
                                    ubs=[float(ub)],
                                    variables_type=[VarType.CONTINUOUS])

            variables.append(var)

            for metabolite, coef in reaction.stoichiometry.items():
                constraint = self._constraints.get(metabolite.id, ConstraintContainer(name=metabolite.id,
                                                                                      coefs=[{}],
                                                                                      lbs=[0.0],
                                                                                      ubs=[0.0]))

                constraint.coefs[0][reaction.id] = coef
                constraints.append(constraint)

        self.add(variables)
        self.add(constraints)

    def remove_reactions(self, reactions: List['Reaction']):

        reactions: List['Reaction'] = iterable(reactions)

        variables = []
        constraints = []
        for reaction in reactions:

            if reaction.id in self._variables:
                variables.append(self._variables[reaction.id])

            for metabolite in reaction.metabolites:

                if metabolite in self._constraints:

                    constraint = self._constraints[metabolite]

                    if reaction.id in constraint.coefs[0]:
                        del constraint.coefs[0][reaction.id]

                    constraints.append(constraint)

        self.remove(variables)
        self.remove(constraints)

    def add_metabolites(self, metabolites: List['Metabolite']):

        metabolites: List['Metabolite'] = iterable(metabolites)

        constraints = [ConstraintContainer(name=metabolite.id,
                                           coefs=[{}],
                                           lbs=[0.0],
                                           ubs=[0.0])
                       for metabolite in metabolites]

        self.add(constraints)

    def remove_metabolites(self, metabolites: List['Metabolite']):

        metabolites: List['Metabolite'] = iterable(metabolites)

        self.remove([self._constraints[metabolite.id] for metabolite in metabolites
                     if metabolite.id in self._constraints])


class LogicLinearizer(LinearProblem):

    @abstractmethod
    def build(self):

        # The concrete implementation is defined by each simulation method, e.g. fba, pfba, etc

        pass

    def notification(self, notification: Notification):

        # The concrete implementation is defined by each simulation method, e.g. fba, pfba, etc

        return super(LogicLinearizer, self).notification(notification)

    @staticmethod
    def linearize_and(symbolic):

        # Following Boolean algebra, an And (a = b and c) can be translated as: a = b*c
        # Alternatively, an And can be written as lb < b + c - a < ub

        # So, for the mid term expression a = b and c
        # We have therefore the equation: -1 <= 2*b + 2*c – 4*a <= 3

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
    def linearize_or(symbolic):

        # Following Boolean algebra, an Or (a = b or c) can be translated as: a = b + c - b*c
        # Alternatively, an Or can be written as lb < b + c - a < ub

        # So, for the mid term expression a = b or c
        # We have therefore the equation: -2 <= 2*b + 2*c – 4*a <= 1

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
    def linearize_not(symbolic):

        # Following Boolean algebra, an Not (a = not b) can be translated as: a = 1 - b
        # Alternatively, an Not can be written as a + b = 1

        # So, for the mid term expression a = not b
        # We have therefore the equation: 1 < a + b < 1

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

    def linearize_greater(self, symbolic):

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

        greater_op = symbolic.key()
        op_l = symbolic.variables[0]
        op_r = symbolic.variables[1]

        if op_l.is_numeric:
            operand = op_r
            c_val = float(op_l.value)

        else:
            operand = op_l
            c_val = float(op_r.value)

        _lb, _ub = bounds_from_symbol(symbol=operand, model=self.model, default=(MEWPY_LB, MEWPY_UB))

        _lb = float(_lb)
        _ub = float(_ub)
        # add Greater row (a(value + tolerance - r_UB) + r < value + tolerance) and set mip bounds to
        # -inf;comparison_val
        # add Greater row (a(r_LB - value - tolerance) + r > r_LB) and set mip bounds to lb;inf
        _coefs = [
            {greater_op: c_val + MEWPY_TOL - _ub, operand.key(): 1.0},
            {greater_op: _lb - c_val - MEWPY_TOL, operand.key(): 1.0}
        ]

        _lbs = [MEWPY_LB, _lb]
        _ubs = [c_val + MEWPY_TOL, MEWPY_UB]

        return ConstraintContainer(name=None, coefs=_coefs, lbs=_lbs, ubs=_ubs)

    def linearize_less(self, symbolic):

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

        less_op = symbolic.key()
        op_l = symbolic.variables[0]
        op_r = symbolic.variables[1]

        if op_l.is_numeric:
            operand = op_r
            c_val = float(op_l.value)

        else:
            operand = op_l
            c_val = float(op_r.value)

        _lb, _ub = bounds_from_symbol(symbol=operand, model=self.model, default=(MEWPY_LB, MEWPY_UB))
        _lb = float(_lb)
        _ub = float(_ub)
        # add Less row (a(value + tolerance - r_LB) + r > value + tolerance) and set mip bounds to
        # -inf;-comparison_val
        # add Less row (a(r_UB - value - tolerance) + r < r_UB) and set mip bounds to lb;inf
        _coefs = [
            {less_op: c_val + MEWPY_TOL - _lb, operand.key(): 1.0},
            {less_op: _ub - c_val - MEWPY_TOL, operand.key(): 1.0}
        ]
        _lbs = [c_val + MEWPY_TOL, MEWPY_LB]
        _ubs = [MEWPY_UB, _ub]

        return ConstraintContainer(name=None, coefs=_coefs, lbs=_lbs, ubs=_ubs)

    @staticmethod
    def linearize_equal(symbolic):

        # Following Propositional logic, a predicate (a => r = value) can be translated as: r - value = 0

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
    def linearize_none(symbolic):

        # If the expression is empty, it means that the target or reaction can take any boolean value (0;1)
        # add row and set mip bounds to 0;1

        return ConstraintContainer(name=None, coefs=[{}], lbs=[0.0], ubs=[1.0])

    @staticmethod
    def linearize_false(symbolic):

        return ConstraintContainer(name=None, coefs=[{}], lbs=[0.0], ubs=[0.0])

    @staticmethod
    def linearize_true(symbolic):

        return ConstraintContainer(name=None, coefs=[{}], lbs=[1.0], ubs=[1.0])

    @staticmethod
    def linearize_number(symbolic):

        return ConstraintContainer(name=None, coefs=[{}], lbs=[float(symbolic.value)], ubs=[float(symbolic.value)])

    @staticmethod
    def linearize_symbol(symbolic):

        return ConstraintContainer(name=None, coefs=[{symbolic.key(): -1.0}], lbs=[0.0], ubs=[0.0])

    def get_linearizer(self,
                       symbolic,
                       operators=True,
                       bool_atoms=True,
                       numeric_atoms=True,
                       symbolic_atoms=True,
                       empty_symbolic=True,
                       ):

        if operators:

            if symbolic.is_and:
                return partial(self.linearize_and, symbolic)

            elif symbolic.is_or:
                return partial(self.linearize_or, symbolic)

            elif symbolic.is_not:
                return partial(self.linearize_not, symbolic)

            elif symbolic.is_greater or symbolic.is_greater_equal:
                return partial(self.linearize_greater, symbolic)

            elif symbolic.is_less or symbolic.is_less_equal:
                return partial(self.linearize_less, symbolic)

            elif symbolic.is_equal:
                return partial(self.linearize_equal, symbolic)

        if bool_atoms:
            if symbolic.is_true:

                return partial(self.linearize_true, symbolic)

            elif symbolic.is_false:

                return partial(self.linearize_false, symbolic)

        if numeric_atoms:
            if symbolic.is_numeric:
                return partial(self.linearize_number, symbolic)

        if symbolic_atoms:
            if symbolic.is_symbol:
                return partial(self.linearize_symbol, symbolic)

        if empty_symbolic:
            if symbolic.is_none:
                return partial(self.linearize_none, symbolic)

        return

    @staticmethod
    def _operator_variablize(symbolic):

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

    def variablize_and(self, symbolic):

        return self._operator_variablize(symbolic=symbolic)

    def variablize_or(self, symbolic):

        return self._operator_variablize(symbolic=symbolic)

    @staticmethod
    def variablize_not(symbolic):

        name = symbolic.key()
        sub_variable_name = f'{name}_0'

        return VariableContainer(name=symbolic.key(),
                                 sub_variables=[sub_variable_name],
                                 lbs=[0.0],
                                 ubs=[1.0],
                                 variables_type=[VarType.INTEGER])

    def variablize_greater(self, symbolic):

        return self._operator_variablize(symbolic=symbolic)

    def variablize_less(self, symbolic):

        return self._operator_variablize(symbolic=symbolic)

    def variablize_equal(self, symbolic):

        return self._operator_variablize(symbolic=symbolic)

    def variablize_symbol(self, symbolic):

        name = symbolic.key()
        lb, ub = bounds_from_symbol(symbol=symbolic, model=self.model)
        lb = float(lb)
        ub = float(ub)

        if (lb, ub) in integer_coefficients:
            v_type = VarType.INTEGER
        else:
            v_type = VarType.CONTINUOUS

        return VariableContainer(name=name, sub_variables=[name], lbs=[lb], ubs=[ub], variables_type=[v_type])

    def get_variabilizer(self, symbolic):

        if symbolic.is_and:
            return partial(self.variablize_and, symbolic)

        elif symbolic.is_or:
            return partial(self.variablize_or, symbolic)

        elif symbolic.is_not:
            return partial(self.variablize_not, symbolic)

        elif symbolic.is_greater or symbolic.is_greater_equal:
            return partial(self.variablize_greater, symbolic)

        elif symbolic.is_less or symbolic.is_less_equal:
            return partial(self.variablize_less, symbolic)

        elif symbolic.is_equal:
            return partial(self.variablize_equal, symbolic)

        elif symbolic.is_symbol:
            return partial(self.variablize_symbol, symbolic)

        return

    def linearize_expression(self, boolean_variable, symbolic):

        # Results of str expression parsing:
        # gene = X and 1 => X
        # gene = X or 1 => 1 (str BooleanTrue)
        # gene = X and 0 => 0 (str BooleanFalse)
        # gene = X or 0 => X
        # gene = 1 => 1 (str one)
        # gene = 0 => 0 (str zero)
        # gene = X > 10 => (X, 1 (str number)), namely two args

        # if expression is atom and defines a variable always On or Off, add a On/Off row
        if symbolic.is_atom:
            return self.linearize_atomic_expression(boolean_variable=boolean_variable, symbolic=symbolic)

        return self.linearize_complex_expression(boolean_variable=boolean_variable, symbolic=symbolic)

    def linearize_atomic_expression(self, boolean_variable, symbolic):

        variables = []
        constraints = []

        if symbolic.is_symbol:
            var = self.variablize_symbol(symbolic=symbolic)
            variables.append(var)

        linearizer = self.get_linearizer(symbolic, operators=False)

        expression_cnt = linearizer()

        expression_cnt.coefs[0][boolean_variable] = 1.0

        constraints.append(expression_cnt)

        return variables, constraints

    def linearize_complex_expression(self, boolean_variable, symbolic):

        # iterating an expression in the reverse order
        # e.g. expr = A and (B or C)
        # list(expression) # [C, B, A, B | C, A & (B | C)]

        variables = []
        constraints = []
        last_variable = None
        # let's add the expressions
        for atom in symbolic:

            # An operator expression will be decomposed into multiple operator expressions, namely into a nested
            # expression. For instance, an A & B & C & D will become And(D, And(C, And(A, B))). Each nested operator
            # expression will be a column in the matrix and then linked in the columns. However, this operator
            # expression will be the same for all maters and have the same column identifier regardless of the length.
            # Thus, all operator columns (variables) will be under the columns linked list engine. When retrieving
            # the operator expression indexes from the cols dictionary in the other methods, the last index
            # of the slice should be used to get the last real column, as the result of this column is the one that
            # really matters. The hashes of the columns for the simulation engine should be the the row name plus
            # str of the operator

            last_variable = atom

            variablizer = self.get_variabilizer(symbolic=atom)

            if variablizer is not None:
                var = variablizer()
                variables.append(var)

            linearizer = self.get_linearizer(atom,
                                             bool_atoms=False,
                                             numeric_atoms=False,
                                             symbolic_atoms=False,
                                             empty_symbolic=False)

            if linearizer is not None:
                constraint = linearizer()
                constraints.append(constraint)

        # identifying the last index to link the outcome of this variable to the boolean variable associated to the
        # expression
        last_variable_name = last_variable.key()
        names = [f'{last_variable_name}_{i}' for i, _ in enumerate(last_variable.variables[:-1])]

        if names:
            last_variable_name = names[-1]
        else:
            last_variable_name = last_variable_name

        # add gene row which means that the gene variable in the mip matrix is associated with the last mid term
        # expression, namely the whole expression
        # set mip bounds to 0;0
        expression_cnt = ConstraintContainer(name=None,
                                             coefs=[{boolean_variable: 1.0, last_variable_name: -1.0}],
                                             lbs=[0.0],
                                             ubs=[0.0])
        constraints.append(expression_cnt)

        return variables, constraints


class GPRLinearizer(LogicLinearizer):

    @abstractmethod
    def build(self):

        # The concrete implementation is defined by each simulation method, e.g. fba, pfba, etc

        pass

    def notification(self, notification: Notification):

        if notification.content_type == 'gprs' and notification.action == 'add':

            return self.add_gprs(notification.content)

        elif notification.content_type == 'gprs' and notification.action == 'remove':

            return self.remove_gprs(notification.content)

        else:
            return super(GPRLinearizer, self).notification(notification)

    def _gpr_constraint(self, reaction: str, symbolic: 'Symbolic'):

        # Relation between reaction boolean value and the reaction constrains
        # For that, the following reactions must be added
        # V - Y*Vmax <= 0
        # V - Y*Vmin => 0
        # where V is the reaction in S matrix
        # where Y is the reaction boolean variable

        # add reaction boolean variable / add column / Add MILP bool variable
        # If already exists, only the index is returned

        boolean_variable = f'bool_{reaction}'

        variables = []
        constraints = []

        boolean_variable_linear_variable = VariableContainer(name=boolean_variable,
                                                             sub_variables=[boolean_variable],
                                                             lbs=[0.0],
                                                             ubs=[1.0],
                                                             variables_type=[VarType.INTEGER])

        variables.append(boolean_variable_linear_variable)

        lb, ub = bounds_from_variable(reaction, model=self.model, default=(MEWPY_LB, MEWPY_UB))

        coefs = [{reaction: 1.0, boolean_variable: -float(ub)},
                 {reaction: 1.0, boolean_variable: -float(lb)}]

        lbs = [MEWPY_LB - float(ub), 0.0]
        ubs = [0.0, MEWPY_UB - float(lb)]

        cnt = ConstraintContainer(name=None,
                                  coefs=coefs,
                                  lbs=lbs,
                                  ubs=ubs)
        constraints.append(cnt)

        expression_variables, expression_cnt = self.linearize_expression(boolean_variable=boolean_variable,
                                                                         symbolic=symbolic)

        variables.extend(expression_variables)
        constraints.extend(expression_cnt)

        constraints = [concat_constraints(constraints=constraints, name=reaction)]

        return variables, constraints

    def add_gprs(self, reactions: List['Reaction']):

        reactions = iterable(reactions)

        to_add = []

        # gprs over reactions
        for rxn in reactions:
            rxn_variables, rxn_cnt = self._gpr_constraint(reaction=rxn.id,
                                                          symbolic=rxn.gpr.symbolic)

            to_add.extend(rxn_variables)
            to_add.extend(rxn_cnt)

        self.add(to_add)

    def remove_gprs(self, reactions: List['Reaction']):

        reactions = iterable(reactions)

        to_remove = []

        for rxn in reactions:
            to_remove.append(self._variables[f'bool_{rxn.id}'])
            to_remove.append(self._constraints[rxn.id])

        self.remove(to_remove)


class InteractionLinearizer(LogicLinearizer):

    @abstractmethod
    def build(self):

        # The concrete implementation is defined by each simulation method, e.g. fba, pfba, etc

        pass

    def notification(self, notification: Notification):

        if notification.content_type == 'interactions' and notification.action == 'add':

            return self.add_interactions(notification.content)

        elif notification.content_type == 'interactions' and notification.action == 'remove':

            return self.remove_interactions(notification.content)

        else:
            return super(InteractionLinearizer, self).notification(notification)

    def _interaction_constraint(self, target, symbolic):

        lb, ub = bounds_from_variable(target, model=self.model, default=(0, 1))
        lb = float(lb)
        ub = float(ub)

        if (lb, ub) in integer_coefficients:
            v_type = VarType.INTEGER
        else:
            v_type = VarType.CONTINUOUS

        variables = []
        constraints = []

        target_linear_variable = VariableContainer(name=target,
                                                   sub_variables=[target],
                                                   lbs=[lb],
                                                   ubs=[ub],
                                                   variables_type=[v_type])

        variables.append(target_linear_variable)

        expression_variables, expression_cnt = self.linearize_expression(boolean_variable=target,
                                                                         symbolic=symbolic)

        variables.extend(expression_variables)
        constraints.extend(expression_cnt)

        constraints = [concat_constraints(constraints=constraints, name=target)]

        return variables, constraints

    def add_interactions(self, interactions):

        interactions = iterable(interactions)

        to_add = []

        # expressions over interactions
        for interaction in interactions:

            symbolic = None

            for coefficient, expression in interaction.regulatory_events.items():

                if coefficient > 0.0:
                    symbolic = expression.symbolic

            if symbolic is None:
                continue

            interaction_variables, interaction_cnt = self._interaction_constraint(target=interaction.target.id,
                                                                                  symbolic=symbolic)

            to_add.extend(interaction_variables)
            to_add.extend(interaction_cnt)

        self.add(to_add)

    def remove_interactions(self, interactions):

        interactions = iterable(interactions)

        to_remove = []

        for interaction in interactions:
            to_remove.append(self._variables[interaction.target.id])
            to_remove.append(self._constraints[interaction.target.id])

        self.remove(to_remove)
