import abc
from typing import List, Any, Dict, Callable, Union, Tuple

from .algebra_utils import solution_decode, _walk


class Symbolic:

    def __init__(self, value=None, variables=None):

        """
        The Symbolic object implements the base logic for symbolic algebra manipulation and evaluation.
        Symbolic is the base class that provides the abstraction and tools for a given variable or operator be handled
        as an algebra operand or operator, respectively.

        Operators and Operands that inherit from Symbolic implement python dunder methods (e.g. __add__), so that they
        can be used in python evaluation method (eval) during binary string evaluation

        Symbolic objects can represent either atoms (e.g. Symbol) or nodes (e.g. And)
        of a binary Abstract Syntax Tree (AST). Thus, each symbolic-object can hold a given value (e.g. Symbol('x'))
        or hold variables/children (e.g. And(variables=[Symbol('x'), Symbol('y')])).

        The symbolic object at the root of the algebra expression (AST) holds the variables and subsequent descents.

        A symbolic object is usually assembled during parsing of an algebra expression encoded into a string object.

        :param value: The value that the symbolic operand is representing in the algebra expression.
        It can either be a numeric or variable value
        :param variables: The variables/children that the symbolic operator is holding in the algebra expression.
        It can either be a Symbol, Numeric, Boolean, or other atoms.
        """

        if not variables:
            variables = []

        self.value = value
        self.variables = variables

        self._name = None
        self.model_variable = None

        self.is_boolean = False
        self.is_true = False
        self.is_false = False
        self.is_and = False
        self.is_or = False
        self.is_not = False

        self.is_relational = False
        self.is_equal = False
        self.is_not_equal = False
        self.is_inequality = False
        self.is_greater = False
        self.is_greater_equal = False
        self.is_less = False
        self.is_less_equal = False

        self.is_numeric = False
        self.is_integer = False
        self.is_float = False
        self.is_one = False
        self.is_zero = False

        self.is_symbol = False
        self.is_atom = False

        self.is_none = False

    @property
    def name(self) -> str:
        """
        String representation of the symbolic object value
        """

        if self._name is None:
            return str(self.value)

        return self._name

    @name.setter
    def name(self, value):

        if value is None:
            return

        self._name = value

    @property
    def bounds(self) -> Tuple[float, float]:
        """
        Returns the bounds of the symbolic object
        """
        if self.model_variable:
            if hasattr(self.model_variable, 'bounds'):
                return self.model_variable.bounds

            if hasattr(self.model_variable, 'coefficients'):
                return min(self.model_variable.coefficients), max(self.model_variable.coefficients)

        return 0, 1

    def key(self):
        return self.to_string().replace(' ', '')

    def __repr__(self):

        if self.is_atom:

            return f"{self.__class__.__name__}({self.value})"

        else:

            return f"{self.__class__.__name__}(variables={self.variables})"

    def __str__(self):

        string_repr = ''

        if self.is_and:
            string_repr += '('
            for child in self.variables:
                string_repr += child.to_string() + ' & '
            string_repr = string_repr[:-3]
            string_repr += ')'

        elif self.is_or:
            string_repr += '('
            for child in self.variables:
                string_repr += child.to_string() + ' | '
            string_repr = string_repr[:-3]
            string_repr += ')'

        elif self.is_not:
            string_repr += '('
            for child in self.variables:
                string_repr += ' ~ ' + child.to_string()
            string_repr += ')'

        elif self.is_equal:
            string_repr += '('
            for child in self.variables:
                string_repr += child.to_string() + ' = '
            string_repr = string_repr[:-3]
            string_repr += ')'

        elif self.is_not_equal:
            string_repr += '('
            for child in self.variables:
                string_repr += child.to_string() + ' != '
            string_repr = string_repr[:-3]
            string_repr += ')'

        elif self.is_inequality:
            string_repr += '('
            for child in self.variables:
                string_repr += child.to_string() + ' > '
            string_repr = string_repr[:-3]
            string_repr += ')'

        elif self.is_greater:
            string_repr += '('
            for child in self.variables:
                string_repr += child.to_string() + ' > '
            string_repr = string_repr[:-3]
            string_repr += ')'

        elif self.is_greater_equal:
            string_repr += '('
            for child in self.variables:
                string_repr += child.to_string() + ' >= '
            string_repr = string_repr[:-3]
            string_repr += ')'

        elif self.is_less:
            string_repr += '('
            for child in self.variables:
                string_repr += child.to_string() + ' < '
            string_repr = string_repr[:-3]
            string_repr += ')'

        elif self.is_less_equal:
            string_repr += '('
            for child in self.variables:
                string_repr += child.to_string() + ' <= '
            string_repr = string_repr[:-3]
            string_repr += ')'

        else:

            string_repr += self.name

        return string_repr

    def to_string(self) -> str:
        """
        String representation of this Symbolic object
        """
        return self.__str__()

    def __iter__(self):

        self.__variables_iter = list(_walk(self, reverse=True))
        self.__i = -1
        self.__i_limit = len(self.__variables_iter)

        return self

    def __next__(self):

        self.__i += 1

        if self.__i < self.__i_limit:
            return self.__variables_iter[self.__i]

        raise StopIteration

    def atoms(self, symbols_only=False) -> List['Symbolic']:
        """
        It returns all symbolic atoms associated with this symbolic object

        :param symbols_only: Whether to return only symbols or not
        :return: It returns a list of symbolic atoms
        """

        symbols = list(_walk(self))

        if symbols_only:

            return [symbol for symbol in symbols if symbol.is_symbol]

        else:

            return [symbol for symbol in symbols if symbol.is_atom]

    def _tree_evaluation(self, values, operators, default, **kwargs):

        left = None
        right = None

        if not self.variables:
            return self._evaluate(values=values, left=left, right=right, operators=operators, default=default, **kwargs)

        # noinspection PyProtectedMember
        left = self.variables[0]._tree_evaluation(values=values, operators=operators, default=default, **kwargs)

        for child in self.variables[1:]:
            # noinspection PyProtectedMember
            right = child._tree_evaluation(values=values, operators=operators, default=default, **kwargs)
            left = self._evaluate(values=values, left=left, right=right, operators=operators, default=default, **kwargs)

        # not special case
        if right is None:
            left = self._evaluate(values=values, left=left, right=right, operators=operators, default=default, **kwargs)

        return left

    def evaluate(self,
                 values: Dict[str, Any] = None,
                 operators: Dict[str, Callable] = None,
                 default=0,
                 **kwargs) -> Union[float, int, Any]:

        """
        Evaluate a Symbolic algebra expression based on the coefficients/values of the symbols or atoms contained
        in the expression.

        The symbolic expression can be evaluated according to the custom values
        assigned to the symbols in the values' dictionary. Likewise, operators evaluation can also be altered
        by passing a callable to a given operator in the operators' dictionary.

        A symbol is also a callable.

        :param values: A dictionary of values that the symbols names must take during expression evaluation
        :param operators: A dictionary of custom operators. That is, python operators-based evaluation
        (e.g. 3 > 2 yield True) can be replaced by custom callable objects such as functions. For instance,
        3 > 2 can be evaluated with max, and thus max(3, 2) yields 3 now.
        :param default: The default value of a given symbol/variable.
        The default value is used if a given variable/symbol is missing in the values dictionary.
        :param kwargs: Additional keyword arguments for Symbolic evaluate method

        :return: The solution of the Symbolic expression evaluation as int, float or Any type.
        """

        if not values:
            values = {}

        if not operators:
            operators = {}

        return self._tree_evaluation(values=values, operators=operators, default=default, **kwargs)

    @staticmethod
    @abc.abstractmethod
    def _evaluate(**kwargs):
        return

    def __call__(self, values=None, operators=None, default=0, **kwargs):

        return self.evaluate(values=values,
                             operators=operators,
                             default=default,
                             **kwargs)


class Boolean(Symbolic):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_boolean = True

    def __and__(self, *args):
        return And(variables=[self, *args])

    __rand__ = __and__

    def __or__(self, *args):
        return Or(variables=[self, *args])

    __ror__ = __or__

    def __invert__(self):
        return Not(variables=[self])

    @staticmethod
    @abc.abstractmethod
    def _evaluate(**kwargs):
        return


class BoolFalse(Boolean):
    """
    Atomic Boolean False value
    """

    def __init__(self, *args, **kwargs):
        if args:
            args = args[1:]

        super().__init__(*args, **kwargs)

        self.value = 0

        self.is_false = True
        self.is_atom = True

    def __bool__(self):
        return False

    @staticmethod
    def _evaluate(**kwargs):
        return 0


class BoolTrue(Boolean):
    """
    Atomic Boolean True value
    """

    def __init__(self, *args, **kwargs):
        if args:
            args = args[1:]

        super().__init__(*args, **kwargs)

        self.value = 1

        self.is_true = True
        self.is_atom = True

    def __bool__(self):
        return True

    @staticmethod
    def _evaluate(**kwargs):
        return 1


class And(Boolean):
    """
    Boolean Operator AND
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_and = True

    @staticmethod
    def _evaluate(left=0, right=0, operators=None, **kwargs):

        if not operators:
            operators = {}

        operator = operators.get(And, None)

        if operator is not None:
            return operator(left, right)

        return solution_decode(int(left) & int(right))

    def __and__(self, *args):
        self.variables.extend(args)
        return self

    __rand__ = __and__


class Or(Boolean):
    """
    Boolean Operator OR
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_or = True

    @staticmethod
    def _evaluate(left=0, right=0, operators=None, **kwargs):

        if not operators:
            operators = {}

        operator = operators.get(Or, None)

        if operator is not None:
            return operator(left, right)

        return solution_decode(int(left) | int(right))

    def __or__(self, *args):
        self.variables.extend(args)
        return self

    __ror__ = __or__


class Not(Boolean):
    """
    Boolean Operator NOT
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_not = True

    @staticmethod
    def _evaluate(left=0, operators=None, **kwargs):

        if not operators:
            operators = {}

        operator = operators.get(Not, None)

        if operator is not None:
            return operator(left)

        return solution_decode(~ int(left))


class Relational(Boolean):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_relational = True

    @staticmethod
    @abc.abstractmethod
    def _evaluate(**kwargs):
        return


class Equal(Relational):
    """
    Relational Operator EQUAL
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_equal = True

    @staticmethod
    def _evaluate(left=0, right=0, operators=None, **kwargs):

        if not operators:
            operators = {}

        operator = operators.get(Equal, None)

        if operator is not None:
            return operator(left, right)

        return solution_decode(left == right)


class NotEqual(Relational):
    """
    Relational Operator NOT EQUAL
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_not_equal = True

    @staticmethod
    def _evaluate(left=0, right=0, operators=None, **kwargs):

        if not operators:
            operators = {}

        operator = operators.get(NotEqual, None)

        if operator is not None:
            return operator(left, right)

        return solution_decode(left != right)


class Inequality(Relational):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_inequality = True

    @staticmethod
    def _evaluate(left=0, right=0, operators=None, **kwargs):

        if not operators:
            operators = {}

        operator = operators.get(Inequality, None)

        if operator is not None:
            return operator(left, right)

        return solution_decode(left > right)


class Greater(Relational):
    """
    Relational Operator Greater
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_greater = True

    @staticmethod
    def _evaluate(left=0, right=0, operators=None, **kwargs):

        if not operators:
            operators = {}

        operator = operators.get(Greater, None)

        if operator is not None:
            return operator(left, right)

        return solution_decode(left > right)


class GreaterEqual(Relational):
    """
    Relational Operator Greater Equal
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_greater_equal = True

    @staticmethod
    def _evaluate(left=0, right=0, operators=None, **kwargs):

        if not operators:
            operators = {}

        operator = operators.get(GreaterEqual, None)

        if operator is not None:
            return operator(left, right)

        return solution_decode(left >= right)


class Less(Relational):
    """
    Relational Operator Less
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_less = True

    @staticmethod
    def _evaluate(left=0, right=0, operators=None, **kwargs):

        if not operators:
            operators = {}

        operator = operators.get(Less, None)

        if operator is not None:
            return operator(left, right)

        return solution_decode(left < right)


class LessEqual(Relational):
    """
    Relational Operator Less Equal
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_less_equal = True

    @staticmethod
    def _evaluate(left=0, right=0, operators=None, **kwargs):

        if not operators:
            operators = {}

        operator = operators.get(LessEqual, None)

        if operator is not None:
            return operator(left, right)

        return solution_decode(left <= right)


class AtomNumber(Symbolic):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_numeric = True
        self.is_atom = True

    def __eq__(self, other):
        return Equal(variables=[self, other])

    def __ne__(self, other):
        return NotEqual(variables=[self, other])

    def __ge__(self, other):
        return GreaterEqual(variables=[self, other])

    def __le__(self, other):
        return LessEqual(variables=[self, other])

    def __gt__(self, other):
        return Greater(variables=[self, other])

    def __lt__(self, other):
        return Less(variables=[self, other])

    @staticmethod
    @abc.abstractmethod
    def _evaluate(**kwargs):
        return


class Number(AtomNumber):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __eq__(self, other):
        return Equal(variables=[self, other])

    def __ne__(self, other):
        return NotEqual(variables=[self, other])

    def __ge__(self, other):
        return GreaterEqual(variables=[self, other])

    def __le__(self, other):
        return LessEqual(variables=[self, other])

    def __gt__(self, other):
        return Greater(variables=[self, other])

    def __lt__(self, other):
        return Less(variables=[self, other])

    @staticmethod
    @abc.abstractmethod
    def _evaluate(**kwargs):
        return


class Integer(Number):
    """
    Atomic Numeric value Integer
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.value = int(self.value)

        self.is_integer = True

    def _evaluate(self, **kwargs):
        return self.value


class Float(Number):
    """
    Atomic Numeric value Float
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.value = float(self.value)

        self.is_float = True

    def _evaluate(self, **kwargs):
        return self.value


class Zero(Boolean, AtomNumber):
    """
    Atomic Numeric value Zero
    """

    def __init__(self, *args, **kwargs):
        if args:
            args = args[1:]

        super().__init__(*args, **kwargs)

        self.value = 0

        self.is_zero = True

    def __bool__(self):
        return False

    def _evaluate(self, **kwargs):
        return 0


class One(Boolean, AtomNumber):
    """
    Atomic Numeric value One
    """

    def __init__(self, *args, **kwargs):
        if args:
            args = args[1:]

        super().__init__(*args, **kwargs)

        self.value = 1

        self.is_one = True

    def __bool__(self):
        return True

    def _evaluate(self, **kwargs):
        return 1


class Symbol(AtomNumber, Boolean):
    """
    Atomic Symbolic Variable
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_boolean = False
        self.is_numeric = False

        self.is_symbol = True

    def _evaluate(self, values=None, default=0, **kwargs):
        if values is None:
            values = {}

        if self.name in values:
            value = values[self.name]

        elif self.model_variable:
            if hasattr(self.model_variable, 'bounds'):
                value = max(self.model_variable.bounds)

            elif hasattr(self.model_variable, 'coefficients'):
                value = max(self.model_variable.coefficients)

            else:
                value = default

        else:
            value = default

        if not isinstance(value, (int, float, bool, AtomNumber)):
            raise ValueError(f'cannot evaluate {self.name} for {value}')

        return solution_decode(value)


class NoneAtom(Symbol):
    """
    Atomic Symbolic Empty value
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_none = True
        self.is_symbol = False

    @property
    def name(self):
        return ''

    def _evaluate(self, **kwargs):
        return 0
