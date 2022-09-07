from typing import Union, List, Tuple, TYPE_CHECKING, Set, Iterable

from mewpy.mew.lp import Notification
from mewpy.util.constants import ModelConstants
from mewpy.util.history import recorder

if TYPE_CHECKING:
    from .variable import Variable
    from .gene import Gene
    from .interaction import Interaction
    from .metabolite import Metabolite
    from .reaction import Reaction
    from .regulator import Regulator
    from .target import Target


class Coefficient:

    def __init__(self,
                 variable: 'Variable',
                 coefficients: Iterable = None,
                 default: Union[int, float] = None):

        """
        A Coefficient object represents variable coefficients or bounds.
        Some variables can vary between two states, such as a gene that can be active or inactive,
        or a reaction that can take fluxes between an upper and lower bound.
        However, some variables can take continuous values, such as a metabolite concentration
        or a target gene expression.

        A coefficient is particularly useful to store all coefficients of a variable in a single object.

        Coefficients vary between 0 and 1 by default. The default coefficient sets the default value of the variable.
        If the default coefficient is not set, the minimum coefficient is set as default.

        :param variable: variable to which the coefficient belongs to
        :param coefficients: all coefficients that the variable can take
        :param default: the default coefficient
        """

        # variable is a must
        if not variable:
            raise ValueError('A coefficient strictly requires a variable')

        # boolean behavior if an empty coefficient is built
        if not coefficients:
            coefficients = (0.0, 1.0)

        self._variable = variable

        # all coefficients are recorded under data.
        self._data = list(coefficients)

        # the default coefficient defines the current level/bound that a variable can take in a linear problem.
        # it is set to the minimum possible by default.
        if default is None:
            default = min(coefficients)

        else:

            if default not in self._data:
                raise ValueError(f'default coefficient {default} is not available in coefficients')

        self._default_coefficient = default

    # -----------------------------------------------------------------------------
    # Built-in
    # -----------------------------------------------------------------------------
    def __str__(self):
        return str(self.coefficients)

    # -----------------------------------------------------------------------------
    # Polymorphic constructors
    # -----------------------------------------------------------------------------
    @classmethod
    def from_bounds(cls,
                    variable: 'Variable',
                    lb: float = None,
                    ub: float = None):
        """
        Creates a coefficient object from a lower and upper bound.

        :param variable: variable to which the coefficient belongs to
        :param lb: lower bound
        :param ub: upper bound
        :return: a coefficient object
        """

        # useful for reactions
        if lb is None:
            lb = ModelConstants.REACTION_LOWER_BOUND

        if ub is None:
            ub = ModelConstants.REACTION_UPPER_BOUND

        return cls(variable=variable,
                   coefficients=(lb, ub))

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------

    @property
    def coefficients(self) -> Union[List[int], List[float]]:
        """
        All coefficients that the variable can take
        :return: list of coefficients
        """
        return self._data.copy()

    @property
    def default_coefficient(self) -> float:
        """
        The default coefficient defines the current level/bound that a variable can take in a linear problem.
        :return: default coefficient as a float
        """
        return float(self._default_coefficient)

    @property
    def is_active(self) -> bool:
        """
        Checks if the default coefficient is higher than the minimum coefficient
        :return: True if the default coefficient is higher than the minimum coefficient
        """
        return self.default_coefficient > self.minimum_coefficient

    @property
    def minimum_coefficient(self) -> Union[int, float]:
        """
        The minimum coefficient is the minimum value that the variable can take.
        :return: minimum coefficient as a float
        """
        return min(self._data)

    @property
    def maximum_coefficient(self) -> Union[int, float]:
        """
        The maximum coefficient is the maximum value that the variable can take.
        :return: maximum coefficient as a float
        """
        return max(self._data)

    @property
    def lower_bound(self) -> Union[int, float]:
        """
        The lower bound is the minimum coefficient
        :return: lower bound as a float
        """
        return self.minimum_coefficient

    @property
    def upper_bound(self) -> Union[int, float]:
        """
        The upper bound is the maximum coefficient
        :return: upper bound as a float
        """
        return self.maximum_coefficient

    @property
    def bounds(self) -> Tuple[Union[int, float], Union[int, float]]:
        """
        The bounds are the minimum and maximum coefficients
        :return: bounds as a tuple of floats
        """
        return self.minimum_coefficient, self.maximum_coefficient

    # -----------------------------------------------------------------------------
    # Variable-based attributes
    # -----------------------------------------------------------------------------

    @property
    def variable(self) -> Union['Variable', 'Gene', 'Interaction', 'Metabolite', 'Reaction', 'Regulator', 'Target']:
        """
        The variable to which the coefficient belongs to
        :return: variable
        """
        return self._variable

    @property
    def id(self):
        """
        The id of the coefficient aka the id of the variable
        :return: id
        """
        return self.variable.id

    @property
    def model(self):
        """
        The model to which the coefficient belongs to
        :return: model
        """
        return self.variable.model

    @property
    def history(self):
        """
        The history of the coefficient aka the history of the variable
        :return: history
        """
        return self.variable.history

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------
    @default_coefficient.setter
    @recorder
    def default_coefficient(self, value: Union[int, float]):

        if value not in self._data:
            raise ValueError(f'{value} is not available for this coefficient object')

        self._default_coefficient = value

        self.notify()

    def _setter(self,
                coefficients: Iterable,
                active=None):

        """
        It sets new coefficients and updates the default coefficient if necessary.
        Internal use only

        :param coefficients: coefficients to be set
        :param active: the default coefficient
        :return:
        """

        self._data = list(coefficients)

        if active is None:

            if self.default_coefficient not in self._data:
                self._default_coefficient = self.minimum_coefficient

        else:

            if active in self._data:
                self._default_coefficient = active

            else:
                self._default_coefficient = self.minimum_coefficient

        self.notify()

    @coefficients.setter
    @recorder
    def coefficients(self, value: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]]):
        """
        It sets new coefficients and updates the default coefficient if necessary.
        It can handle the linear problems associated with the variable.

        :param value: coefficients to be set
        :return:
        """
        self._setter(value)

    def __getitem__(self, item: int) -> Union[int, float]:
        return self._data[item]

    def __setitem__(self, key: int, value: Union[int, float], history: bool = True):
        old_coef = self.__getitem__(key)

        self._data[value] = value

        self.notify()

        if history:
            self.history.queue_command(undo_func=self.__setitem__,
                                       undo_kwargs={'key': key,
                                                    'value': old_coef,
                                                    'history': False},
                                       func=self.__setitem__,
                                       kwargs={'key': key,
                                               'value': value,
                                               'history': False})

    def get(self, index: int, default=None) -> Union[int, float]:

        try:
            return self.__getitem__(index)

        except IndexError:

            return default

    def index(self, value: Union[int, float], default=None) -> int:

        try:
            return self._data.index(value)

        except ValueError:

            return default

    def _add(self, coefficient: Union[int, float], default: bool = False):
        """
        It adds a new coefficient to the list of coefficients.
        Internal use only

        :param coefficient: coefficient to be added
        :param default: if True, the coefficient will be set as the default coefficient
        :return:
        """

        self._data.append(coefficient)

        if default:
            self._default_coefficient = coefficient

    def add(self, coefficient: Union[int, float], default: bool = False, history: bool = True):
        """
        It adds a new coefficient to the list of coefficients.
        It can handle the linear problems associated with the variable.

        :param coefficient: coefficient to be added
        :param default: if True, the coefficient will be set as the default coefficient
        :param history: if True, the command will be queued to the history
        :return:
        """

        if history:
            self.history.queue_command(undo_func=self.remove,
                                       undo_kwargs={'coefficient': coefficient,
                                                    'history': False},
                                       func=self.add,
                                       kwargs={'coefficient': coefficient,
                                               'default': default,
                                               'history': False})

        self._add(coefficient=coefficient, default=default)

        self.notify()

    def remove(self, coefficient: Union[int, float], history: bool = True) -> Union[int, float]:
        """
        It removes a coefficient from the list of coefficients.
        It can handle the linear problems associated with the variable.

        :param coefficient: coefficient to be removed
        :param history: if True, the command will be queued to the history
        :return: the removed coefficient
        """
        removed = self._data.remove(coefficient)

        is_active = False

        if coefficient == self._default_coefficient:
            is_active = True

        if is_active:
            self._default_coefficient = self.minimum_coefficient

        if history:
            self.history.queue_command(undo_func=self.add,
                                       undo_kwargs={'coefficient': removed,
                                                    'default': is_active,
                                                    'history': False},
                                       func=self.remove,
                                       kwargs={'coefficient': coefficient,
                                               'history': False})

        self.notify()

        return removed

    def extend(self,
               coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]],
               default: float = None,
               history: bool = True):
        """
        It extends the list of coefficients with a new list of coefficients.
        It can handle the linear problems associated with the variable.

        :param coefficients: coefficients to be added
        :param default: if not None, the coefficient will be set as the default coefficient
        :param history: if True, the command will be queued to the history
        :return:
        """
        old_coefficients = self.coefficients
        old_default = self.default_coefficient

        for coefficient in coefficients:

            is_default = False

            if coefficient == default:
                is_default = True

            self._add(coefficient, default=is_default)

        if history:
            self.history.queue_command(undo_func=self._setter,
                                       undo_kwargs={'coefficients': old_coefficients,
                                                    'default': old_default},
                                       func=self.extend,
                                       kwargs={'coefficients': coefficients,
                                               'history': False})

        self.notify()

    def ko(self, history: bool = True):
        """
        It sets the minimum coefficient as the only possible coefficient that the variable can take.
        It can handle the linear problems associated with the variable.

        :param history: if True, the command will be queued to the history
        :return:
        """
        coefficient = self.minimum_coefficient

        old_coefficients = self.coefficients
        old_default = self.default_coefficient

        self._data = [coefficient]
        self._default_coefficient = coefficient

        if history:
            self.history.queue_command(undo_func=self._setter,
                                       undo_kwargs={'coefficients': old_coefficients,
                                                    'default': old_default},
                                       func=self.ko,
                                       kwargs={'coefficient': coefficient,
                                               'history': False})

        self.notify()

    def activate(self, history: bool = True):
        """
        It sets the maximum coefficient as the only possible coefficient that the variable can take.
        It can handle the linear problems associated with the variable.

        :param history: if True, the command will be queued to the history
        :return:
        """
        maximum_coefficient = self.maximum_coefficient

        old_coefficients = self.coefficients
        old_default = self.default_coefficient

        self._data = [maximum_coefficient]
        self._default_coefficient = maximum_coefficient

        if history:
            self.history.queue_command(undo_func=self._setter,
                                       undo_kwargs={'coefficients': old_coefficients,
                                                    'default': old_default,
                                                    'history': False},
                                       func=self.activate,
                                       kwargs={'maximum_coefficient': maximum_coefficient,
                                               'history': False})

        self.notify()

    # -----------------------------------------------------------------------------
    # Simulators interface
    # -----------------------------------------------------------------------------
    def notify(self):
        """
        It notifies the variable coefficient change to the simulators.
        :return:
        """

        if self.model:

            notification = Notification(content=self,
                                        content_type='coefficients',
                                        action='set')

            self.model.notify(notification)
