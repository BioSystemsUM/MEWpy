from typing import Any, TYPE_CHECKING, Set, Union, Dict, Generator, Tuple, List

from mewpy.util.utilities import generator
from mewpy.util.serialization import serialize
from mewpy.util.history import recorder
from .coefficient import Coefficient
from .variable import Variable

if TYPE_CHECKING:
    from .interaction import Interaction
    from .regulator import Regulator


# TODO: methods stubs
class Target(Variable, variable_type='target', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier: Any,
                 coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]] = None,
                 active_coefficient: Union[int, float] = None,
                 interaction: 'Interaction' = None,
                 **kwargs):

        """
        A target is commonly associated with one and only one interaction but
        can usually be available as gene and regulator too.

        It holds information regarding the coefficients that can take and the interaction that is modelling these
        coefficients

        :param identifier: identifier, e.g. b0001
        :param coefficients: the set of coefficients that this target can take.
        These coefficients can be expanded later. 0 and 1 are added by default
        :param active_coefficient: the active coefficient
        :param interactions: the  interaction to which the target is associated with
        """

        # the coefficient initializer sets minimum and maximum coefficients of 0.0 and 1.0
        if not coefficients:
            coefficients = (0.0, 1.0)

        if not active_coefficient:
            active_coefficient = max(coefficients)

        if not interaction:
            interaction = None

        self._coefficient = Coefficient(variable=self, coefficients=coefficients, active=active_coefficient)
        self._interaction = interaction

        super().__init__(identifier,
                         **kwargs)

    # -----------------------------------------------------------------------------
    # Variable type manager
    # -----------------------------------------------------------------------------

    @property
    def types(self):

        # noinspection PyUnresolvedReferences
        _types = {Target.variable_type}

        _types.update(super(Target, self).types)

        return _types

    # -----------------------------------------------------------------------------
    # Built-in
    # -----------------------------------------------------------------------------

    def __str__(self):

        if self.interaction:

            return f'{self.interaction}'

        return super(Target, self).__str__()

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------

    @serialize('coefficient', 'coefficients', '_coefficient')
    @property
    def coefficient(self) -> Coefficient:

        if hasattr(self, '_bounds'):

            # if it is a reaction, the bounds coefficient must be returned
            return self._bounds

        elif hasattr(self, 'exchange_reaction'):

            # if it is a metabolite, the bounds coefficient of the exchange reaction must be returned
            if hasattr(self.exchange_reaction, '_bounds'):
                # noinspection PyProtectedMember
                return self.exchange_reaction._bounds

        return self._coefficient

    @serialize('interaction', 'interaction', '_interaction')
    @property
    def interaction(self) -> 'Interaction':
        return self._interaction

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------

    @interaction.setter
    @recorder
    def interaction(self, value: 'Interaction'):

        if not value:
            value = None

        value.add_target(self, history=False)

        self._interaction = value

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------

    @property
    def regulators(self) -> Dict[str, 'Regulator']:

        if self.interaction:
            return self.interaction.regulators.copy()

        return {}

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------

    def yield_regulators(self) -> Generator['Regulator', None, None]:
        return generator(self.regulators)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------

    def ko(self, minimum_coefficient: Union[int, float] = 0.0, history=True):

        return self.coefficient.ko(minimum_coefficient=minimum_coefficient, history=history)

    def update(self,
               coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]] = None,
               active_coefficient: Union[int, float] = None,
               interaction: 'Interaction' = None,
               **kwargs):

        super(Target, self).update(**kwargs)

        if coefficients is not None:
            self.coefficient.coefficients = coefficients

        if active_coefficient is not None:
            self.coefficient.active_coefficient = active_coefficient

        if interaction is not None:
            self.interaction = interaction
