import abc

from .algebra_utils import solution_decode


def _walk(symbolic, reverse=False):
    if reverse:

        for child in symbolic.variables:
            for subtree in _walk(child, reverse):
                yield subtree

        yield symbolic

    else:

        yield symbolic

        for child in symbolic.variables:
            yield from _walk(child, reverse)


# TODO: methods stubs and type hinting
class Symbolic:

    def __init__(self, value=None, variables=None):

        if not variables:
            variables = []

        self.value = value
        self.variables = variables

        self._name = None

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
    def name(self):

        if self._name is None:

            return str(self.value)

        return self._name

    @name.setter
    def name(self, value):

        if value is None:
            return

        self._name = value

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

    def to_string(self):
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

    def atoms(self, symbols_only=False):

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

    def evaluate(self, values=None, operators=None, default=0, **kwargs):

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

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.value = int(self.value)

        self.is_integer = True

    def _evaluate(self, **kwargs):
        return self.value


class Float(Number):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.value = float(self.value)

        self.is_float = True

    def _evaluate(self, **kwargs):
        return self.value


class Zero(Boolean, AtomNumber):

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

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_boolean = False
        self.is_numeric = False

        self.is_symbol = True

    def _evaluate(self, values=None, default=0, **kwargs):

        if not values:
            values = {}

        value = values.get(self.name, default)

        if not isinstance(value, (int, float, bool, AtomNumber)):
            raise ValueError(f'cannot evaluate {self.name} for {value}')

        return solution_decode(value)


class NoneAtom(Symbol):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.is_none = True
        self.is_symbol = False

    @property
    def name(self):
        return ''

    def _evaluate(self, **kwargs):
        return 0
