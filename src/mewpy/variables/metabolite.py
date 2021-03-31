from re import findall
from typing import Any, Dict, Generator, Union, TYPE_CHECKING

from mewpy.util.utilities import generator, chemical_formula_re
from mewpy.util.serialization import serialize
from mewpy.util.history import recorder
from mewpy.util.constants import atomic_weights
from .variable import Variable

if TYPE_CHECKING:
    from .reaction import Reaction


# TODO: methods stubs
class Metabolite(Variable, variable_type='metabolite', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier: Any,
                 charge: int = None,
                 compartment: str = None,
                 formula: str = None,
                 reactions: Dict[str, 'Reaction'] = None,
                 **kwargs):

        """
        A metabolite is regularly associated with reactions and
        can usually be available as regulator too.
        It holds information regarding the charge, compartment, formula and reactions to which is associated
        Some dynamic information is inferred from the formula such as molecular weight and atoms.
        Other information is inferred from the reactions associated, such as the exchange reaction

        :param identifier: identifier, e.g. h2o_e
        :param charge: the charge of the metabolite
        :param compartment: the compartment of this respective metabolite
        :param formula: a string-like representation of the chemical formula
        :param reactions: the dictionary of reactions to which the metabolite is associated with
        """

        if not charge and charge is not 0:
            charge = None

        if not compartment:
            compartment = None

        if not formula:
            formula = ''

        if not reactions:
            reactions = {}

        self._charge = charge
        self._compartment = compartment
        self._formula = formula
        self._reactions = reactions

        super().__init__(identifier,
                         **kwargs)

    # -----------------------------------------------------------------------------
    # Variable type manager
    # -----------------------------------------------------------------------------

    @property
    def types(self):

        # noinspection PyUnresolvedReferences
        _types = {Metabolite.variable_type}

        _types.update(super(Metabolite, self).types)

        return _types

    # -----------------------------------------------------------------------------
    # Built-in
    # -----------------------------------------------------------------------------

    def __str__(self):

        return f'{self.id}: {self.name} {self.formula}'

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------

    @serialize('charge', 'charge', '_charge')
    @property
    def charge(self) -> int:

        if self._charge is None:
            return 0

        return self._charge

    @serialize('compartment', 'compartment', '_compartment')
    @property
    def compartment(self) -> str:
        return self._compartment

    @serialize('formula', 'formula', '_formula')
    @property
    def formula(self) -> str:
        return self._formula

    @serialize('reactions', 'reactions', '_reactions')
    @property
    def reactions(self) -> Dict[str, 'Reaction']:
        return self._reactions.copy()

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------

    @charge.setter
    @recorder
    def charge(self, value):

        if not value and value is not 0:
            value = None

        self._charge = value

    @formula.setter
    @recorder
    def formula(self, value):

        if not value:
            value = ''

        self._formula = value

    @compartment.setter
    @recorder
    def compartment(self, value):

        if not value:
            value = None

        self._compartment = value

    @reactions.setter
    def reactions(self, value):

        if not value:
            value = {}

        self._reactions = value

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------

    @property
    def atoms(self) -> Dict[str, int]:

        all_elements = findall(chemical_formula_re, self.formula)

        atoms = {}
        for atom, count in all_elements:

            if not count:
                count = '1'

            atoms[atom] = atoms.get(atom, 0) + int(count)

        return atoms

    @property
    def molecular_weight(self) -> Union[float, int]:

        return sum([atomic_weights[atom] * count for atom, count in self.atoms.items()])

    @property
    def exchange_reaction(self) -> 'Reaction':

        for reaction in self.yield_reactions():
            if reaction.boundary:
                return reaction

    @property
    def exchange_reactions(self) -> Dict[str, 'Reaction']:

        exchanges = {}

        for reaction in self.yield_reactions():
            if reaction.boundary:
                exchanges[reaction.id] = reaction

        return exchanges

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------

    def yield_reactions(self) -> Generator['Reaction', None, None]:

        return generator(self._reactions)

    def yield_exchange_reactions(self) -> Generator['Reaction', None, None]:

        return generator(self.exchange_reactions)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------

    def update(self,
               charge: int = None,
               compartment: str = None,
               formula: str = None,
               reactions: Dict[str, 'Reaction'] = None,
               **kwargs):

        super(Metabolite, self).update(**kwargs)

        if charge is not None:
            self.charge = charge

        if compartment is not None:
            self.compartment = compartment

        if formula is not None:
            self.formula = formula

        if reactions is not None:
            self._reactions.update(reactions)
