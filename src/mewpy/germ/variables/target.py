from typing import Any, TYPE_CHECKING, Dict, Generator, Tuple, Sequence, Union

from mewpy.util.constants import ModelConstants
from mewpy.util.history import recorder
from mewpy.germ.models.serialization import serialize
from mewpy.util.utilities import generator
from .variable import Variable
from .variables_utils import coefficients_setter

if TYPE_CHECKING:
    from .interaction import Interaction
    from .regulator import Regulator


class Target(Variable, variable_type='target', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier: Any,
                 coefficients: Sequence[float] = None,
                 interaction: 'Interaction' = None,
                 **kwargs):

        """
        A target gene is associated with a single interaction but can be regulated by multiple regulators.
        The regulatory interaction establishes a relationship between a target gene and regulators.
        A target gene holds the coefficients that it can take. While 0 and 1 are added by default,
        these coefficients can be changed later.

        A target gene can also be represented as a metabolic gene in Metabolic-Regulatory models, as it
        can be inferred from the metabolic rules defining metabolic interactions.
        In regulatory models, a target gene can also be represented as a regulator gene, as it can be inferred
        from the regulatory rules defining regulatory interactions.

        :param identifier: identifier, e.g. b0001
        :param coefficients: the set of coefficients that this target can take.
        These coefficients can be expanded later. 0 and 1 are added by default
        :param interactions: the  interaction to which the target is associated with
        """

        # the coefficient initializer sets minimum and maximum coefficients of 0.0 and 1.0
        if not coefficients:
            coefficients = (0.0, 1.0)
        else:
            coefficients = tuple(coefficients)

        if not interaction:
            interaction = None

        self._coefficients = coefficients
        self._interaction = interaction

        super().__init__(identifier, **kwargs)

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

            return str(self.interaction)

        return f'{self.id} || {self.coefficients}'

    def _target_to_html(self):
        """
        It returns a html representation.
        """
        html_dict = {'Coefficients': self.coefficients,
                     'Active': self.is_active,
                     'Interaction': self.interaction,
                     'Regulators': ', '.join(self.regulators)}
        return html_dict

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------

    @serialize('coefficients', 'coefficients', '_coefficients')
    @property
    def coefficients(self) -> Tuple[float, ...]:
        """
        The coefficients of the target gene.
        :return: the coefficients
        """
        if hasattr(self, '_bounds'):

            # if it is a reaction, bounds must be returned
            return self._bounds

        # if it is a metabolite, the bounds coefficient of the exchange reaction must be returned
        elif hasattr(self, 'exchange_reaction'):

            if hasattr(self.exchange_reaction, '_bounds'):
                # noinspection PyProtectedMember
                return self.exchange_reaction._bounds

        return self._coefficients

    @serialize('interaction', 'interaction', '_interaction')
    @property
    def interaction(self) -> 'Interaction':
        """
        The interaction to which the target is associated with.
        :return: the interaction
        """
        return self._interaction

    @property
    def is_active(self):
        """
        It checks whether the gene is active or not
        :return: True if the gene is active, False otherwise
        """
        return max(self.coefficients) > ModelConstants.TOLERANCE

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------
    @coefficients.setter
    @recorder
    def coefficients(self, value: Union[float, Sequence[float]]):
        """
        The target coefficients setter
        :param value: The target coefficients
        :return:
        """
        coefficients_setter(self, value)

    @interaction.setter
    @recorder
    def interaction(self, value: 'Interaction'):
        """
        The interaction setter.
        It sets the interaction and adds the target to the interaction.

        It also handles the linear problems associated with the interaction.
        :param value: the interaction to which the target is associated with
        :return:
        """
        if not value:
            value = None

        value.add_target(self, history=False)

        self._interaction = value

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------
    @property
    def regulators(self) -> Dict[str, 'Regulator']:
        """
        The regulators that regulate the target gene.
        :return: the regulators as a dictionary
        """
        if self.interaction:
            return self.interaction.regulators.copy()

        return {}

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------
    def yield_regulators(self) -> Generator['Regulator', None, None]:
        """
        A generator that yields the regulators that regulate the target gene.
        :return:
        """
        return generator(self.regulators)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------
    def ko(self, minimum_coefficient: float = 0.0, history=True):
        """
        Knock-out the target gene.
        :param minimum_coefficient: the minimum coefficient
        :param history: whether to record the operation
        :return:
        """
        old_coef = tuple(self.coefficients)

        coefficients_setter(self, minimum_coefficient)

        if history:
            self.history.queue_command(undo_func=coefficients_setter,
                                       undo_kwargs={'instance': self,
                                                    'value': old_coef},
                                       func=self.ko,
                                       kwargs={'minimum_coefficient': minimum_coefficient,
                                               'history': False})

    def update(self,
               coefficients: Sequence[float] = None,
               interaction: 'Interaction' = None,
               **kwargs):
        """
        It updates the target gene with relevant information.

        It also handles the linear problems associated to the interaction.

        Note that, some update operations are not registered in history.
        It is strongly advisable to use update outside history context manager

        :param coefficients: the set of coefficients that this target can take.
        :param interaction: the  interaction to which the target is associated with
        :param kwargs: additional arguments
        :return:
        """
        super(Target, self).update(**kwargs)

        if coefficients is not None:
            self.coefficients = coefficients

        if interaction is not None:
            self.interaction = interaction
