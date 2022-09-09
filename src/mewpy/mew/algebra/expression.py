from itertools import product
from typing import Dict, Union, TYPE_CHECKING, Callable, Any

import pandas as pd

from .algebra_utils import solution_decode, _walk
from .parsing import tokenize
from .symbolic import NoneAtom, Symbolic, Symbol

if TYPE_CHECKING:
    from mewpy.mew.variables import Gene, Metabolite, Reaction, Regulator, Interaction, Target, Variable


class Expression:

    def __init__(self,
                 symbolic: Symbolic = None,
                 variables: Dict[str, Union['Gene',
                                            'Interaction',
                                            'Metabolite',
                                            'Reaction',
                                            'Regulator',
                                            'Target',
                                            'Variable']] = None):

        """
        Expression allows high-level programming of symbolic algebra (logic and math)
        using mewpy variables (e.g. Gene, Interaction, Metabolite, Reaction, Regulator, Target, etc).

        An expression allows symbolic algebra evaluation with mewpy variables or Symbolic-based symbols using
        regular python binary operators or custom functions.

        An expression can be tokenized, translated or iterated.
        The truth table can be computed for the variables or symbols.

        The identifiers of all variables must match the names of all symbols, otherwise an error is raised.

        :param symbolic: The symbolic object which implements an algebra expression
        containing Symbolic operators, symbols or atoms. If None (default), an empty symbolic expression,
        aka NoneAtom, is created.
        :param variables: A dictionary of mewpy variables instances. The mewpy variables that must be associated
        with this symbolic expression. Variables identifiers, and thus dictionary keys,
        must match all symbols in the Symbolic expression
        """

        if symbolic is None:
            symbolic = NoneAtom()

        if variables is None:
            variables = {}

        self._symbolic: Symbolic = symbolic
        self._variables: Dict[str, Union['Gene',
                                         'Interaction',
                                         'Metabolite',
                                         'Reaction',
                                         'Regulator',
                                         'Target']] = variables

        self._check_association()

    def _check_association(self):

        for symbol in self.symbols.values():

            if symbol.name not in self._variables:
                raise ValueError('Expression set with incorrect variables or symbolic expression')

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------
    @property
    def symbolic(self) -> Symbolic:
        """
        Symbolic property
        :return: The symbolic logic or math algebra expression.
        """
        return self._symbolic

    @property
    def variables(self) -> Dict[str, Union['Gene', 'Interaction', 'Metabolite', 'Reaction', 'Regulator', 'Target']]:
        """
        Variables dictionary property
        :return: A copy of the mewpy variables dictionary
        """
        return self._variables.copy()

    # -----------------------------------------------------------------------------
    # Static setters
    # -----------------------------------------------------------------------------
    @variables.setter
    def variables(self, value: Dict[str, Union['Gene',
                                               'Interaction',
                                               'Metabolite',
                                               'Reaction',
                                               'Regulator',
                                               'Target']]):

        """
        Variables setter
        :param value: A dictionary of mewpy variables instances. The mewpy variables that must be associated
        with this symbolic expression. Variables identifiers, and thus dictionary keys,
        must match all symbols in the Symbolic expression
        """

        self._variables = value

        self._check_association()

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------
    @property
    def symbols(self) -> Dict[str, 'Symbol']:

        """
        Symbols property.
        :return: A dictionary containing all symbols in the symbolic algebra expression
        """
        return {symbol.name: symbol for symbol in self.symbolic.atoms(symbols_only=True)}

    @property
    def is_none(self):
        """
        Is none property.
        :return: bool True whether the current expression is empty or a NoneAtom
        """
        return self.symbolic.is_none

    # -----------------------------------------------------------------------------
    # Built-in
    # -----------------------------------------------------------------------------
    def __repr__(self):
        return repr(self.symbolic)

    def __str__(self):
        return str(self.symbolic)

    def to_string(self):
        """
        String representation of the symbolic algebra expression
        :return: string with bitwise operators and symbols names
        """

        return self.symbolic.to_string()

    def to_tokenize(self):
        """
        Tokenize symbolic algebra expression
        :return: list of all tokens (no spaces) in the string representation
        """
        return tokenize(self.to_string())

    # -----------------------------------------------------------------------------
    # Iteration
    # -----------------------------------------------------------------------------
    def __iter__(self):

        return self.symbolic.__iter__()

    def __next__(self):

        return self.symbolic.__next__()

    def walk(self, reverse=False):

        """
        Iterate/walk over Symbolic algebra expression yielding Symbolic-type symbols or operators.
        Parent operators are iterated first by default.

        :param reverse: If true, Symbolic symbols are yield first
        :return: A generator object to be iterated (yield Symbolic symbols and operators)
        """

        return _walk(self.symbolic, reverse)

    def __call__(self,
                 values: Dict[str, Union[float, int]],
                 coefficient: Union[float, int] = None,
                 operators: Dict[Symbolic, Callable] = None,
                 default: Union[float, int] = 0.0,
                 decoder: dict = None,
                 **kwargs) -> Union[float, int, Any]:

        """
        Evaluate a Symbolic algebra expression based on
        the coefficients/values of the Symbolic symbols - mewpy variables.

        An expression is a also callable. Thus, the symbolic expression is evaluated according to the Symbolic
        operators and values attributed to the symbols in the values dictionary.

        :param values: A dictionary of values that the variables identifiers, aka symbols names,
        must take during expression evaluation
        :param coefficient: If a coefficient is provided the coefficient value is returned
        in case the Symbolic expression is evaluated to True.
        :param operators: A dictionary of custom operators. That is, python operators-based evaluation
        (e.g. 3 > 2 yield True) can be replaced by custom callable objects such as functions. For instance,
        3 > 2 can be evaluated with max, and thus max(3, 2) yields 3 now.
        :param default: The default value of a given symbol/variable.
        The default value is used if a given variable/symbol is missing in the values dictionary.
        :param decoder: A custom dictionary for decoding the solution (key) into a given output (value).
        Binary output is currently set to 0 or 1 by default
        :param kwargs: Additional keyword arguments for Symbolic evaluate method

        :return: The solution of the Symbolic expression evaluation as int, float or Any type.
        """
        return self.evaluate(values=values,
                             coefficient=coefficient,
                             operators=operators,
                             default=default,
                             decoder=decoder,
                             **kwargs)

    def evaluate(self,
                 values: Dict[str, Union[float, int]],
                 coefficient: Union[float, int] = None,
                 operators: Dict[Symbolic, Callable] = None,
                 default: Union[float, int] = 0.0,
                 decoder: dict = None,
                 **kwargs):

        """
        Evaluate a Symbolic algebra expression based on
        the coefficients/values of the Symbolic symbols - mewpy variables.

        The symbolic expression is evaluated according to the Symbolic operators and values attributed to the symbols
        in the values dictionary.

        :param values: A dictionary of values that the variables identifiers, aka symbols names,
        must take during expression evaluation
        :param coefficient: If a coefficient is provided the coefficient value is returned
        in case the Symbolic expression is evaluated to True.
        :param operators: A dictionary of custom operators. That is, python operators-based evaluation
        (e.g. 3 > 2 yield True) can be replaced by custom callable objects such as functions. For instance,
        3 > 2 can be evaluated with max, and thus max(3, 2) yields 3 now.
        :param default: The default value of a given symbol/variable.
        The default value is used if a given variable/symbol is missing in the values dictionary.
        :param decoder: A custom dictionary for decoding the solution (key) into a given output (value).
        Binary output is currently set to 0 or 1 by default
        :param kwargs: Additional keyword arguments for Symbolic evaluate method

        :return: The solution of the Symbolic expression evaluation as int, float or Any type.
        """

        res = self.symbolic.evaluate(values=values, operators=operators, default=default, **kwargs)

        res = solution_decode(res, decoder)

        # if a coefficient is provided and the expression is evaluated to true, the respective coefficient will
        # be returned
        if coefficient is None:
            return res

        else:

            if res:
                return coefficient

            return res

    def truth_table(self,
                    values: Dict[str, Union[float, int]] = None,
                    default_coefficients: bool = True,
                    coefficient: Union[float, int] = None,
                    operators: Dict[Symbolic, Callable] = None,
                    default: Union[float, int] = 0.0,
                    decoder: dict = None) -> pd.DataFrame:

        """
        It calculates the truth table for this expression. The truth table is composed by the combination of values
        taken by empty, numeric and symbolic variables available in the algebra expression.

        The symbolic expression is evaluated according to the Symbolic operators and values attributed to the symbols
        in the values' dictionary.

        :param values: A dictionary of values that the variables identifiers, aka symbols names,
        must take during expression evaluation. If the values dictionary is defined,
        only one solution will be calculated for the truth table
        :param default_coefficients: A boolean value (True by default) whether only the mewpy variables default
        states/coefficients should be used for the truth table calculation. If so, only one solution will be calculated
        for the truth table
        :param coefficient: If a coefficient is provided the coefficient value is returned
        in case the Symbolic expression is evaluated to True.
        :param operators: A dictionary of custom operators. That is, python operators-based evaluation
        (e.g. 3 > 2 yield True) can be replaced by custom callable objects such as functions. For instance,
        3 > 2 can be evaluated with max, and thus max(3, 2) yields 3 now.
        :param default: The default value of a given symbol/variable.
        The default value is used if a given variable/symbol is missing in the values' dictionary.
        :param decoder: A custom dictionary for decoding the solution (key) into a given output (value).
        Binary output is currently set to 0 or 1 by default

        :return: The solution of the combination of values taken by empty, numeric and symbolic variables available
        in the algebra expression as a pandas DataFrame. DataFrame columns should be composed by all
        empty, numeric and symbolic variables in addition to the result columns
        (the result of evaluating the algebra expression). DataFrame rows should stand for the cartesian product of all
        combination of values taken by empty, numeric and symbolic variables,
        unless specific parameters are taken.

        """

        truth_table = []

        if values:

            values['result'] = self.evaluate(values=values,
                                             coefficient=coefficient,
                                             operators=operators,
                                             default=default,
                                             decoder=decoder)

            truth_table.append(values)

        else:
            if default_coefficients:
                values = {key: variable.coefficient.default_coefficient
                          for key, variable in self._variables.items()}

                values['result'] = self.evaluate(values=values,
                                                 coefficient=coefficient,
                                                 operators=operators,
                                                 default=default,
                                                 decoder=decoder)

                truth_table.append(values)

            else:

                variables = list(self.variables.keys())
                coefficients = [variable.coefficient.coefficients for variable in self._variables.values()]

                for combined_coefficients in product(*coefficients):
                    values = dict(list(zip(variables, combined_coefficients)))

                    values['result'] = self.evaluate(values=values,
                                                     coefficient=coefficient,
                                                     operators=operators,
                                                     default=default,
                                                     decoder=decoder)

                    truth_table.append(values)

        return pd.DataFrame(truth_table)
