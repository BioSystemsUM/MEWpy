from typing import Any, Dict, TYPE_CHECKING, Generator, Tuple, Sequence, Union

from mewpy.util.constants import ModelConstants
from mewpy.util.history import recorder
from mewpy.mew.models.serialization import serialize
from mewpy.util.utilities import generator
from .variable import Variable
from .variables_utils import coefficients_setter

if TYPE_CHECKING:
    from .interaction import Interaction
    from .target import Target


class Regulator(Variable, variable_type='regulator', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier: Any,
                 coefficients: Sequence[float] = None,
                 interactions: Dict[str, 'Interaction'] = None,
                 **kwargs):
        """
        A regulator is commonly associated with interactions and
        can usually be available as metabolite or reaction or target too.

        Regulators can control the gene expression of several target genes.

        It holds information regarding the coefficients that can take and the interactions to which is associated

        :param identifier: identifier, e.g. b0001
        :param coefficients: the set of coefficients that this regulator can take.
        These coefficients can be expanded later. 0 and 1 are added by default
        :param interactions: the dictionary of interactions to which the regulator is associated with
        """

        # the coefficient initializer sets minimum and maximum coefficients of 0.0 and 1.0
        if not coefficients:
            coefficients = (0.0, 1.0)
        else:
            coefficients = tuple(coefficients)

        if not interactions:
            interactions = {}

        self._coefficients = coefficients
        self._interactions = interactions

        super().__init__(identifier, **kwargs)

    # -----------------------------------------------------------------------------
    # Variable type manager
    # -----------------------------------------------------------------------------

    @property
    def types(self):

        # noinspection PyUnresolvedReferences
        _types = {Regulator.variable_type}

        _types.update(super(Regulator, self).types)

        return _types

    def __str__(self):

        if self.is_reaction() or self.is_metabolite():
            return super(Regulator, self).__str__()

        return f'{self.id} || {self.coefficients}'

    def _regulator_to_html(self):
        """
        It returns a html representation.
        """
        html_dict = {'Coefficients': self.coefficients,
                     'Active': self.is_active,
                     'Interactions': ', '.join(self.interactions),
                     'Targets': ', '.join(self.targets),
                     'Environmental stimulus': self.environmental_stimulus}
        return html_dict

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------

    @serialize('coefficients', 'coefficients', '_coefficients')
    @property
    def coefficients(self) -> Tuple[float, ...]:
        """
        The coefficient of the regulator
        :return: the coefficient
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

    @serialize('interactions', 'interactions', '_interactions')
    @property
    def interactions(self) -> Dict[str, 'Interaction']:
        """
        The dictionary of interactions to which the regulator is associated with
        :return: the dictionary of interactions
        """
        return self._interactions.copy()

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

    @interactions.setter
    def interactions(self, value: Dict[str, 'Interaction']):
        """
        The dictionary of interactions to set to the regulator

        It does not perform any check or action on the interactions and targets.
        :param value: the dictionary of interactions
        :return:
        """
        if not value:
            value = {}

        self._interactions = value

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------

    @property
    def targets(self) -> Dict[str, 'Target']:
        """
        The dictionary of targets to which the regulator is associated with
        :return: the dictionary of targets
        """
        return {interaction.target.id: interaction.target for interaction in self.yield_interactions()
                if interaction.target is not None}

    @property
    def environmental_stimulus(self) -> bool:
        """
        True if the regulator is an environmental stimulus
        :return: True if the regulator is an environmental stimulus
        """
        if self.types == {'regulator'}:
            return True

        return False

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------

    def yield_interactions(self) -> Generator['Interaction', None, None]:
        """
        Yields the interactions to which the regulator is associated with
        :return: the generator of interactions
        """
        return generator(self._interactions)

    def yield_targets(self) -> Generator['Target', None, None]:
        """
        Yields the targets to which the regulator is associated with
        :return: the generator of targets
        """
        return generator(self.targets)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------

    def ko(self, minimum_coefficient: float = 0.0, history=True):
        """
        Knock out the regulator by setting the minimum coefficient
        :param minimum_coefficient: the minimum coefficient
        :param history: if True, the previous coefficient is stored in the history
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
               interactions: Dict[str, 'Interaction'] = None,
               **kwargs):
        """
        It updates the regulator
        Note that, some update operations are not registered in history.
        It is strongly advisable to use update outside history context manager
        :param coefficients: the set of coefficients that this regulator can take.
        These coefficients can be expanded later. 0 and 1 are added by default
        :param interactions: the dictionary of interactions to which the regulator is associated with
        :param kwargs: the keyword arguments to pass to the super class
        """
        super().update(**kwargs)

        if coefficients is not None:
            self.coefficients = coefficients

        if interactions is not None:
            self._interactions.update(interactions)
