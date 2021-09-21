from typing import Any, Dict, TYPE_CHECKING, Set, Union, List, Tuple, Generator

from mewpy.util.utilities import generator
from mewpy.util.serialization import serialize
from .coefficient import Coefficient
from .variable import Variable

if TYPE_CHECKING:
    from .reaction import Reaction


# TODO: methods stubs
class Gene(Variable, variable_type='gene', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier: Any,
                 coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]] = None,
                 active_coefficient: Union[int, float] = None,
                 reactions: Dict[str, 'Reaction'] = None,
                 **kwargs):

        """
        A metabolic gene is regularly associated with reactions and can usually be available as
        target too.
        It holds information regarding the coefficients that can take and the reactions to which is associated


        :param identifier: identifier, e.g. b0001
        :param coefficients: the set of coefficients that this gene can take.
        These coefficients can be expanded later. 0 and 1 are added by default
        :param active_coefficient: the active coefficient
        :param reactions: the dictionary of reactions to which the gene is associated with

        """

        # the coefficient initializer sets minimum and maximum coefficients of 0.0 and 1.0
        if not coefficients:
            coefficients = (0.0, 1.0)

        if not active_coefficient:
            active_coefficient = max(coefficients)

        if not reactions:
            reactions = {}

        self._coefficient = Coefficient(variable=self, coefficients=coefficients, active=active_coefficient)
        self._reactions = reactions

        super().__init__(identifier,
                         **kwargs)

    # -----------------------------------------------------------------------------
    # Built-in
    # -----------------------------------------------------------------------------

    @property
    def types(self):

        # noinspection PyUnresolvedReferences
        _types = {Gene.variable_type}

        _types.update(super(Gene, self).types)

        return _types

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

    @serialize('reactions', 'reactions', '_reactions')
    @property
    def reactions(self) -> Dict[str, 'Reaction']:
        return self._reactions.copy()

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------

    @reactions.setter
    def reactions(self, value: Dict[str, 'Reaction']):

        if not value:
            value = {}

        self._reactions = value

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------

    def yield_reactions(self) -> Generator['Reaction', None, None]:

        return generator(self._reactions)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------

    def ko(self, minimum_coefficient: Union[int, float] = 0.0, history=True):

        return self.coefficient.ko(minimum_coefficient=minimum_coefficient, history=history)

    def update(self,
               coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]] = None,
               active_coefficient: Union[int, float] = None,
               reactions: Dict[str, 'Reaction'] = None,
               **kwargs):

        super(Gene, self).update(**kwargs)

        if coefficients is not None:
            self.coefficient.coefficients = coefficients

        if active_coefficient is not None:
            self.coefficient.active_coefficient = active_coefficient

        if reactions is not None:
            self._reactions.update(reactions)
