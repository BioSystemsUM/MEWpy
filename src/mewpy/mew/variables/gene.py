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
        A metabolic gene is regularly associated with metabolic reactions. A metabolic gene is inferred from a metabolic
        model by parsing the Gene-Protein-Reactions associations usually encoded as boolean rules.

        A gene usually stores the gene coefficients that it gene can take during GPR rule evaluation. It usually takes
        either 0 (inactive) or 1 (active)
        A gene object also holds information regarding the reactions that it is associated with.

        In addition, a metabolic gene can also be represented as a target gene in Metabolic-Regulatory models, as it
        can be inferred from the regulatory rules defining regulatory interactions.

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
        """
        The gene coefficients
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

    @serialize('reactions', 'reactions', '_reactions')
    @property
    def reactions(self) -> Dict[str, 'Reaction']:
        """
        The reactions associated with this gene
        """
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
        """
        It yields all reactions
        """
        return generator(self._reactions)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------
    def ko(self, minimum_coefficient: Union[int, float] = 0.0, history=True):
        """
        It performs a gene knock-out. This is accomplished by setting the coefficient object to zero.

        This operation can be reverted using the model history

        :param minimum_coefficient: The minimum coefficient that represents a gene KO. 0.0 is used as default
        :param history: Whether to register this operation in the model history
        :return: None
        """
        return self.coefficient.ko(minimum_coefficient=minimum_coefficient, history=history)

    def update(self,
               coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]] = None,
               active_coefficient: Union[int, float] = None,
               reactions: Dict[str, 'Reaction'] = None,
               **kwargs):
        """
        It performs an update operation to this gene.
        The update operation is similar to a dictionary update.

        Note that, some update operations are not registered in history.
        It is strongly advisable to use update outside history context manager

        :param coefficients: The gene coefficients
        :param active_coefficient: The active gene coefficient
        :param reactions: The reactions associated with this gene
        :param kwargs: Other arguments for the base variable, such as identifier, name, etc
        :return:
        """

        super(Gene, self).update(**kwargs)

        if coefficients is not None:
            self.coefficient.coefficients = coefficients

        if active_coefficient is not None:
            self.coefficient.active_coefficient = active_coefficient

        if reactions is not None:
            self._reactions.update(reactions)
