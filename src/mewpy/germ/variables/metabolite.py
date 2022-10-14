from re import findall
from typing import Any, Dict, Generator, Union, TYPE_CHECKING

from mewpy.util.utilities import generator, chemical_formula_re
from mewpy.germ.models.serialization import serialize
from mewpy.util.history import recorder
from mewpy.util.constants import atomic_weights
from .variable import Variable

if TYPE_CHECKING:
    from .reaction import Reaction


class Metabolite(Variable, variable_type='metabolite', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier: Any,
                 charge: int = None,
                 compartment: str = None,
                 formula: str = None,
                 reactions: Dict[str, 'Reaction'] = None,
                 **kwargs):

        """
        A metabolite is regularly associated with reactions.
        In metabolic-regulatory models, metabolites can be associated with regulators too.

        It holds information regarding the charge, compartment, formula and reactions to which is associated
        Some dynamic information is inferred from the formula such as molecular weight and atoms.
        Other information is inferred from the reactions associated, such as the exchange reaction

        :param identifier: identifier, e.g. h2o_e
        :param charge: the charge of the metabolite
        :param compartment: the compartment of this respective metabolite
        :param formula: a string-like representation of the chemical formula
        :param reactions: the dictionary of reactions to which the metabolite is associated with
        """

        if not charge and charge != 0:
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

        return f'{self.id} || {self.name} || {self.formula}'

    def _metabolite_to_html(self):
        """
        It returns a html dict representation.
        """
        html_dict = {'Compartment': self.compartment,
                     'Formula': self.formula,
                     'Molecular weight': self.molecular_weight,
                     'Charge': self.charge,
                     'Reactions': ', '.join(self.reactions)}
        return html_dict

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------
    @serialize('charge', 'charge', '_charge')
    @property
    def charge(self) -> int:
        """
        The charge of the metabolite
        :return: charge as int
        """

        if self._charge is None:
            return 0

        return self._charge

    @serialize('compartment', 'compartment', '_compartment')
    @property
    def compartment(self) -> str:
        """
        The compartment of the metabolite
        :return: compartment as str
        """
        return self._compartment

    @serialize('formula', 'formula', '_formula')
    @property
    def formula(self) -> str:
        """
        The chemical formula of the metabolite
        :return: formula as str
        """
        return self._formula

    @serialize('reactions', 'reactions', '_reactions')
    @property
    def reactions(self) -> Dict[str, 'Reaction']:
        """
        The reactions to which the metabolite is associated with
        :return: reactions as dict
        """
        return self._reactions.copy()

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------
    @charge.setter
    @recorder
    def charge(self, value):
        """
        Sets the charge of the metabolite
        :param value: charge as int
        :return:
        """

        if not value and value != 0:
            value = None

        self._charge = value

    @formula.setter
    @recorder
    def formula(self, value):
        """
        Sets the chemical formula of the metabolite
        :param value: formula as str
        :return:
        """

        if not value:
            value = ''

        self._formula = value

    @compartment.setter
    @recorder
    def compartment(self, value):
        """
        Sets the compartment of the metabolite
        :param value: compartment as str
        :return:
        """

        if not value:
            value = None

        self._compartment = value

    @reactions.setter
    def reactions(self, value):
        """
        Sets the reactions to which the metabolite is associated with
        It does not perform additional operations in the associated reactions and models
        :param value: reactions as dict
        :return:
        """

        if not value:
            value = {}

        self._reactions = value

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------
    @property
    def atoms(self) -> Dict[str, int]:
        """
        Atoms and frequency of the metabolite
        :return: atoms as dict
        """

        all_elements = findall(chemical_formula_re, self.formula)

        atoms = {}
        for atom, count in all_elements:

            if not count:
                count = '1'

            atoms[atom] = atoms.get(atom, 0) + int(count)

        return atoms

    @property
    def molecular_weight(self) -> Union[float, int]:
        """
        The molecular weight of the metabolite
        :return: molecular weight as float
        """

        return sum([atomic_weights[atom] * count for atom, count in self.atoms.items()])

    @property
    def exchange_reaction(self) -> 'Reaction':
        """
        The exchange reaction of the metabolite.
        It finds the first boundary reaction in which the metabolite is involved
        :return: exchange reaction as Reaction
        """

        for reaction in self.yield_reactions():
            if reaction.boundary:
                return reaction

    @property
    def exchange_reactions(self) -> Dict[str, 'Reaction']:
        """
        The exchange reactions of the metabolite.
        It finds all boundary reactions in which the metabolite is involved
        :return: exchange reactions as dict
        """

        exchanges = {}

        for reaction in self.yield_reactions():
            if reaction.boundary:
                exchanges[reaction.id] = reaction

        return exchanges

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------
    def yield_reactions(self) -> Generator['Reaction', None, None]:
        """
        Yields the reactions to which the metabolite is associated with
        :return: reactions as generator
        """

        return generator(self._reactions)

    def yield_exchange_reactions(self) -> Generator['Reaction', None, None]:
        """
        Yields the exchange reactions to which the metabolite is associated with
        :return: exchange reactions as generator
        """
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

        """
        It updates the metabolite with the provided information

        Note that, some update operations are not registered in history.
        It is strongly advisable to use update outside history context manager

        :param charge: the charge of the metabolite
        :param compartment: the compartment of this respective metabolite
        :param formula: a string-like representation of the chemical formula
        :param reactions: the dictionary of reactions to which the metabolite is associated with
        :param kwargs: additional attributes
        :return:
        """

        super(Metabolite, self).update(**kwargs)

        if charge is not None:
            self.charge = charge

        if compartment is not None:
            self.compartment = compartment

        if formula is not None:
            self.formula = formula

        if reactions is not None:
            self._reactions.update(reactions)
