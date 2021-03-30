from itertools import product
from typing import Dict, Union, TYPE_CHECKING, Callable

# TODO: this module depends on pandas dataframes. Should it be set as package requirement?
# noinspection PyPackageRequirements
from pandas import DataFrame

from .algebra_utils import solution_decode
from .parsing import tokenize
from .symbolic import _walk, NoneAtom, Symbolic

if TYPE_CHECKING:
    from mewpy.variables import Gene, Metabolite, Reaction, Regulator


# TODO: methods stubs and type hinting
class Expression:

    def __init__(self,
                 symbolic: Symbolic = None,
                 variables: Dict[str, Union['Gene', 'Metabolite', 'Reaction', 'Regulator']] = None):

        """
        Association between a symbolic mathematical expression and variables.

        :param symbolic:
        :param variables:

        """

        if symbolic is None:
            symbolic = NoneAtom()

        if variables is None:
            variables = {}

        self._symbolic: Symbolic = symbolic
        self._variables: Dict[str, Union['Gene', 'Metabolite', 'Reaction', 'Regulator']] = variables

        self.check_association()

    def check_association(self):

        for symbol in self.symbols.values():

            if symbol.name not in self._variables:
                raise ValueError('Expression set with incorrect variables or symbolic expression')

    @property
    def symbolic(self) -> Symbolic:
        return self._symbolic

    @property
    def variables(self) -> Dict[str, Union['Gene', 'Metabolite', 'Reaction', 'Regulator']]:
        return self._variables.copy()

    @variables.setter
    def variables(self, value):

        self._variables = value

        self.check_association()

    @property
    def symbols(self):
        return {symbol.name: symbol for symbol in self.symbolic.atoms(symbols_only=True)}

    @property
    def is_none(self):
        return self.symbolic.is_none

    def __repr__(self):
        return repr(self.symbolic)

    def __str__(self):
        return str(self.symbolic)

    def to_string(self):
        return self.symbolic.to_string()

    def to_tokenize(self):
        return tokenize(self.to_string())

    def __iter__(self):

        return self.symbolic.__iter__()

    def __next__(self):

        return self.symbolic.__next__()

    def walk(self, reverse=False):

        return _walk(self.symbolic, reverse)

    def __call__(self,
                 values: Dict[str, int],
                 coefficient: Union[float, int] = None,
                 operators: Dict[Symbolic, Callable] = None,
                 default: Union[float, int] = 0.0,
                 decoder: dict = None,
                 **kwargs):

        return self.evaluate(values=values,
                             coefficient=coefficient,
                             operators=operators,
                             default=default,
                             decoder=decoder,
                             **kwargs)

    def evaluate(self,
                 values: Dict[str, int],
                 coefficient: Union[float, int] = None,
                 operators: Dict[Symbolic, Callable] = None,
                 default: Union[float, int] = 0.0,
                 decoder: dict = None,
                 **kwargs):

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
                    active_states: bool = True,
                    coefficient: Union[float, int] = None,
                    operators: Dict[Symbolic, Callable] = None,
                    default: Union[float, int] = 0.0,
                    decoder: dict = None):

        truth_table = []

        if values:

            values['result'] = self.evaluate(values=values,
                                             coefficient=coefficient,
                                             operators=operators,
                                             default=default,
                                             decoder=decoder)

            truth_table.append(values)

        else:
            if active_states:
                values = {key: variable.coefficient.active_coefficient
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

        return DataFrame(truth_table)
