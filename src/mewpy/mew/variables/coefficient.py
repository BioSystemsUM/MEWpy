from typing import Union, List, Tuple, TYPE_CHECKING, Set

from mewpy.util.history import recorder
from mewpy.util.constants import ModelConstants
from mewpy.mew.algebra import Symbolic
from mewpy.mew.lp import Notification

if TYPE_CHECKING:
    from .variable import Variable
    from .gene import Gene
    from .interaction import Interaction
    from .metabolite import Metabolite
    from .reaction import Reaction
    from .regulator import Regulator
    from .target import Target


# TODO: methods stubs and type hinting
# TODO: Coefficient can be further integrated with the variables that use this
#  object. For that, coefficient composition pattern may be applied
class Coefficient:

    def __init__(self,
                 variable: 'Variable',
                 coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]] = None,
                 active: Union[int, float] = None):

        """

        A coefficient handles multiple coefficients and bounds for a given variable of any type.
        Coefficient object holds many attributes regularly associated with variables,
        such as bounds, coefficients, active coefficient, among others. These attributes are dynamically generated based
        on the coefficients provided as input.

        A coefficient is particularly useful for variables that can take multiple states

        :param variable: variable to which the coefficient belongs to
        :param coefficients: all coefficients that the variable can take
        :param active: the active coefficient
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

        # the active coefficient defines the current level/bound that a variable can take in a linear problem.
        # by default it is set to the minimum possible
        if active is None:
            active = min(coefficients)

        else:

            if active not in self._data:
                raise ValueError(f'active coefficient {active} is not available in coefficients')

        self._active_coefficient = active

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
                    lb: Union[int, float] = None,
                    ub: Union[int, float] = None):

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
    def coefficients(self) -> List[Union[int, float]]:

        return self._data.copy()

    @property
    def active_coefficient(self) -> float:
        return float(self._active_coefficient)

    @property
    def is_active(self) -> bool:
        return self.active_coefficient > ModelConstants.TOLERANCE

    @property
    def minimum_coefficient(self) -> Union[int, float]:

        return min(self._data)

    @property
    def maximum_coefficient(self) -> Union[int, float]:

        return max(self._data)

    @property
    def lower_bound(self) -> Union[int, float]:

        return self.minimum_coefficient

    @property
    def upper_bound(self) -> Union[int, float]:

        return self.maximum_coefficient

    @property
    def bounds(self) -> Tuple[Union[int, float], Union[int, float]]:

        return self.minimum_coefficient, self.maximum_coefficient

    # -----------------------------------------------------------------------------
    # Variable-based attributes
    # -----------------------------------------------------------------------------

    @property
    def variable(self) -> Union['Variable', 'Gene', 'Interaction', 'Metabolite', 'Reaction', 'Regulator', 'Target']:
        return self._variable

    @property
    def id(self):
        return self.variable.id

    @property
    def model(self):
        return self.variable.model

    @property
    def history(self):
        return self.variable.history

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------

    @active_coefficient.setter
    @recorder
    def active_coefficient(self, value: Union[int, float]):

        if value not in self._data:
            raise ValueError(f'{value} is not available for this coefficient object')

        self._active_coefficient = value

        self.notify()

    def _setter(self,
                coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]],
                active=None):

        """
        Internal use only

        :param coefficients: coefficients to be set
        :param active: the active coefficient
        :return:
        """

        self._data = list(coefficients)

        if active is None:

            if self.active_coefficient not in self._data:
                self._active_coefficient = self.minimum_coefficient

        else:

            if active in self._data:
                self._active_coefficient = active

            else:
                self._active_coefficient = self.minimum_coefficient

        self.notify()

    @coefficients.setter
    @recorder
    def coefficients(self, value: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]]):

        self._setter(value)

    def __getitem__(self, item: int) -> Union[int, float]:

        return self._data[item]

    def __setitem__(self, key: int, value: Union[int, float], history=True):

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

    def _add(self, coefficient: Union[int, float], active: bool = False):

        """
        Internal use only

        :param coefficient:
        :param active:
        :return:
        """

        self._data.append(coefficient)

        if active:
            self._active_coefficient = coefficient

    def add(self, coefficient: Union[int, float], active: bool = False, history=True):

        if history:
            self.history.queue_command(undo_func=self.remove,
                                       undo_kwargs={'coefficient': coefficient,
                                                    'history': False},
                                       func=self.add,
                                       kwargs={'coefficient': coefficient,
                                               'active': active,
                                               'history': False})

        self._add(coefficient=coefficient, active=active)

        self.notify()

    def remove(self, coefficient: Union[int, float], history=True) -> Union[int, float]:

        removed = self._data.remove(coefficient)

        is_active = False

        if coefficient == self._active_coefficient:
            is_active = True

        if is_active:
            self._active_coefficient = self.minimum_coefficient

        if history:
            self.history.queue_command(undo_func=self.add,
                                       undo_kwargs={'coefficient': removed,
                                                    'active': is_active,
                                                    'history': False},
                                       func=self.remove,
                                       kwargs={'coefficient': coefficient,
                                               'history': False})

        self.notify()

        return removed

    def extend(self,
               coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]],
               active_coefficient: float = None,
               history=True):

        old_coefficients = self.coefficients
        old_active = self.active_coefficient

        for coefficient in coefficients:

            active = False

            if coefficient == active_coefficient:
                active = True

            self._add(coefficient, active=active)

        if history:
            self.history.queue_command(undo_func=self._setter,
                                       undo_kwargs={'coefficients': old_coefficients,
                                                    'active': old_active},
                                       func=self.extend,
                                       kwargs={'coefficients': coefficients,
                                               'history': False})

        self.notify()

    def ko(self, minimum_coefficient: Union[int, float] = None, history=True):

        if not minimum_coefficient:
            minimum_coefficient = self.minimum_coefficient

        old_coefficients = self.coefficients
        old_active = self.active_coefficient

        self._data = [minimum_coefficient]
        self._active_coefficient = minimum_coefficient

        if history:
            self.history.queue_command(undo_func=self._setter,
                                       undo_kwargs={'coefficients': old_coefficients,
                                                    'active': old_active},
                                       func=self.ko,
                                       kwargs={'minimum_coefficient': minimum_coefficient,
                                               'history': False})

        self.notify()

    def activate(self, maximum_coefficient: Union[int, float] = None, history=True):

        if not maximum_coefficient:
            maximum_coefficient = self.maximum_coefficient

        old_coefficients = self.coefficients
        old_active = self.active_coefficient

        self._data = [maximum_coefficient]
        self._active_coefficient = maximum_coefficient

        if history:
            self.history.queue_command(undo_func=self._setter,
                                       undo_kwargs={'coefficients': old_coefficients,
                                                    'active': old_active,
                                                    'history': False},
                                       func=self.activate,
                                       kwargs={'maximum_coefficient': maximum_coefficient,
                                               'history': False})

        self.notify()

    # -----------------------------------------------------------------------------
    # Simulators interface
    # -----------------------------------------------------------------------------

    def notify(self):

        if self.model:

            notification = Notification(content=self,
                                        content_type='coefficients',
                                        action='set')

            self.model.notify(notification)
