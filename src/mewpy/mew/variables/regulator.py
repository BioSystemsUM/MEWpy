from typing import Any, Set, Union, Dict, TYPE_CHECKING, Generator, Tuple, List

from mewpy.util.utilities import generator
from mewpy.util.serialization import serialize
from .coefficient import Coefficient
from .variable import Variable

if TYPE_CHECKING:
    from .interaction import Interaction
    from .target import Target


# TODO: methods stubs
class Regulator(Variable, variable_type='regulator', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier: Any,
                 coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]] = None,
                 active_coefficient: Union[int, float] = None,
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
        :param active_coefficient: the default coefficient
        :param interactions: the dictionary of interactions to which the regulator is associated with
        """

        # the coefficient initializer sets minimum and maximum coefficients of 0.0 and 1.0
        if not coefficients:
            coefficients = (0.0, 1.0)

        if not active_coefficient:
            active_coefficient = min(coefficients)

        if not interactions:
            interactions = {}

        self._coefficient = Coefficient(variable=self, coefficients=coefficients, default=active_coefficient)
        self._interactions = interactions

        super().__init__(identifier,
                         **kwargs)

    # -----------------------------------------------------------------------------
    # Variable type manager
    # -----------------------------------------------------------------------------

    @property
    def types(self):

        # noinspection PyUnresolvedReferences
        _types = {Regulator.variable_type}

        _types.update(super(Regulator, self).types)

        return _types

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------

    @serialize('coefficient', 'coefficients', '_coefficient')
    @property
    def coefficient(self) -> Coefficient:
        """
        The coefficient of the regulator
        :return: the coefficient
        """
        if hasattr(self, '_bounds'):

            # if it is a reaction, the bounds coefficient must be returned
            return self._bounds

        elif hasattr(self, 'exchange_reaction'):

            # if it is a metabolite, the bounds coefficient of the exchange reaction must be returned
            if hasattr(self.exchange_reaction, '_bounds'):
                # noinspection PyProtectedMember
                return self.exchange_reaction._bounds

        return self._coefficient

    @serialize('interactions', 'interactions', '_interactions')
    @property
    def interactions(self) -> Dict[str, 'Interaction']:
        """
        The dictionary of interactions to which the regulator is associated with
        :return: the dictionary of interactions
        """
        return self._interactions.copy()

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------
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

        elif self.is_regulator() and self.is_metabolite():
            return True

        elif self.is_regulator() and self.is_reaction():
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

    def ko(self, minimum_coefficient: Union[int, float] = 0.0, history=True):
        """
        Knock out the regulator by setting the minimum coefficient
        :param minimum_coefficient: the minimum coefficient
        :param history: if True, the previous coefficient is stored in the history
        """
        return self.coefficient.ko(coefficient=minimum_coefficient, history=history)

    def update(self,
               coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]] = None,
               active_coefficient: Union[int, float] = None,
               interactions: Dict[str, 'Interaction'] = None,
               **kwargs):
        """
        It updates the regulator
        Note that, some update operations are not registered in history.
        It is strongly advisable to use update outside history context manager
        :param coefficients: the set of coefficients that this regulator can take.
        These coefficients can be expanded later. 0 and 1 are added by default
        :param active_coefficient: the default coefficient
        :param interactions: the dictionary of interactions to which the regulator is associated with
        :param kwargs: the keyword arguments to pass to the super class
        """
        super().update(**kwargs)

        if coefficients is not None:
            self.coefficient.coefficients = coefficients

        if active_coefficient is not None:
            self.coefficient.default_coefficient = active_coefficient

        if interactions is not None:
            self._interactions.update(interactions)
