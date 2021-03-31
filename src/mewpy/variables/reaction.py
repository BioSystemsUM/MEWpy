from typing import Any, Dict, Union, TYPE_CHECKING, Tuple, Generator, List

from mewpy.algebra import Expression, parse_expression
from mewpy.lp import Notification
from mewpy.util.utilities import generator
from mewpy.util.serialization import serialize
from mewpy.util.history import recorder
from mewpy.util.constants import ModelConstants
from mewpy.io.engines.engines_utils import expression_warning
from .coefficient import Coefficient
from .variable import Variable, variables_from_symbolic

if TYPE_CHECKING:
    from .metabolite import Metabolite
    from .gene import Gene
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


# TODO: methods stubs
class Reaction(Variable, variable_type='reaction', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier: Any,
                 bounds: Tuple[Union[float, int], Union[float, int]] = None,
                 stoichiometry: Dict['Metabolite', Union[float, int]] = None,
                 gpr: Expression = None,
                 **kwargs):

        """
        A set of reactions is the basis for all metabolic models.
        This reaction object holds information for bounds, metabolites, stoichiometry and
        gene-protein-reaction rules (the associated between reactions, enzymes/proteins and genes)

        A reaction usually contains the stoichiometry of the metabolites that take part in the reaction.
        A reaction usually contains the GPR rule/logic of the associated genes
        A reaction usually contains the bounds for the fluxes values that can take during simulation

        Dynamic information can be obtained from the previous attributes, such as metabolites, reactants, products,
        genes, boundary, reversibility, charge and mass balance, equation, compartments, etc.

        A reaction is the main object for metabolic models and thus the object that offers more manipulation solutions.

        :param identifier: identifier, e.g. PYK
        :param bounds: tuple having both lower and upper bounds, respectively the minimum and maximum coefficient
        that a reaction can have
        :param stoichiometry: dictionary of metabolites objects and their corresponding stoichiometric value
        :param gpr: an expression object with the boolean logic of the associated genes.
        Genes objects must be associated with the expression object.
        """

        # the coefficient bounds initializer sets minimum and maximum bounds of MEWPY_LB and MEWPY_UB
        if not bounds:
            lb, ub = ModelConstants.REACTION_LOWER_BOUND, ModelConstants.REACTION_UPPER_BOUND

        else:
            lb, ub = bounds

        if not gpr:
            gpr = Expression()

        if not stoichiometry:
            stoichiometry = {}

        self._bounds = Coefficient.from_bounds(variable=self, lb=lb, ub=ub)
        self._gpr = Expression()
        self._stoichiometry = {}

        super().__init__(identifier,
                         **kwargs)

        self.replace_stoichiometry(stoichiometry, history=False)

        # The gpr must be an expression object. The expression object properties are used to check whether the genes
        # container is correct.
        self.add_gpr(gpr, history=False)

        # There is an alternative constructor for building a reaction with a stringify-like gpr rule.
        # Parsing takes the hard job to infer the correct relation between genes

    # -----------------------------------------------------------------------------
    # Variable type manager
    # -----------------------------------------------------------------------------

    @property
    def types(self):

        # noinspection PyUnresolvedReferences
        _types = {Reaction.variable_type}

        _types.update(super(Reaction, self).types)

        return _types

    # -----------------------------------------------------------------------------
    # Built-in
    # -----------------------------------------------------------------------------

    def __str__(self):
        return f'{self.id}: {self.equation}'

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------

    @serialize('stoichiometry', 'stoichiometry', '_stoichiometry')
    @property
    def stoichiometry(self) -> Dict['Metabolite', Union[int, float]]:
        return self._stoichiometry.copy()

    @serialize('coefficient', 'bounds', '_bounds')
    @property
    def coefficient(self) -> Coefficient:
        return self._bounds

    @serialize('gpr', 'gpr', '_gpr')
    @property
    def gpr(self) -> Expression:
        return self._gpr

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------

    @stoichiometry.setter
    @recorder
    def stoichiometry(self, value: Dict['Metabolite', Union[int, float]]):

        self.replace_stoichiometry(value, history=False)

    @gpr.setter
    @recorder
    def gpr(self, value: Expression):

        self.remove_gpr(history=False)
        self.add_gpr(value, history=False)

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------

    @property
    def metabolites(self) -> Dict[str, 'Metabolite']:
        return {met.id: met for met in self._stoichiometry}

    @property
    def products(self) -> Dict[str, 'Metabolite']:
        return {met.id: met for met, st in self._stoichiometry.items()
                if st > 0.0}

    @property
    def reactants(self) -> Dict[str, 'Metabolite']:
        return {met.id: met for met, st in self._stoichiometry.items()
                if st < 0.0}

    @property
    def compartments(self) -> set:

        return {met.compartment
                for met in self.yield_metabolites()
                if met.compartment is not None}

    @property
    def bounds(self) -> Tuple[Union[int, float], Union[int, float]]:
        return self.coefficient.bounds

    @property
    def lower_bound(self) -> Union[int, float]:
        return self.coefficient.lower_bound

    @property
    def upper_bound(self) -> Union[int, float]:
        return self.coefficient.upper_bound

    @property
    def reversibility(self) -> bool:
        return self.lower_bound < - ModelConstants.TOLERANCE and self.upper_bound > ModelConstants.TOLERANCE

    @property
    def equation(self) -> str:

        equation = ''

        for _id, met in self.reactants.items():
            equation += f'{-self._stoichiometry[met]} {_id} + '

        if self.reversibility:
            equation = equation[:-2] + '<-> '
        else:
            equation = equation[:-2] + '-> '

        for _id, met in self.products.items():
            equation += f'{self._stoichiometry[met]} {_id} + '

        if not equation:
            return ''

        elif self.products:
            return equation[:-2]

        else:
            return equation[:-1]

    @property
    def genes(self) -> Dict[str, 'Gene']:

        return self.gpr.variables.copy()

    @property
    def gene_protein_reaction_rule(self) -> str:
        return self.gpr.to_string()

    @property
    def boundary(self) -> bool:
        return not (self.reactants and self.products)

    @property
    def charge_balance(self) -> Dict[str, Union[int, float]]:
        return {'reactants': sum([self._stoichiometry[met] * met.charge for met in self.yield_reactants()]),
                'products': sum([self._stoichiometry[met] * met.charge for met in self.yield_products()])}

    @property
    def mass_balance(self) -> Dict[str, Union[int, float]]:

        balance = {}

        for metabolite, st in self._stoichiometry.items():

            for atom, count in metabolite.atoms.items():
                amount = count * st

                balance[atom] = balance.get(atom, 0) + amount

        return balance

    # -----------------------------------------------------------------------------
    # Dynamic attributes setters
    # -----------------------------------------------------------------------------

    @bounds.setter
    @recorder
    def bounds(self, value: Tuple[Union[float, int], Union[float, int]]):

        if not value:
            value = (ModelConstants.REACTION_LOWER_BOUND, ModelConstants.REACTION_UPPER_BOUND)

        self.coefficient.coefficients = value

    @lower_bound.setter
    @recorder
    def lower_bound(self, value: Union[float, int]):

        if not value:
            value = ModelConstants.REACTION_LOWER_BOUND

        value = (value,  self.upper_bound)

        self.coefficient.coefficients = value

    @upper_bound.setter
    @recorder
    def upper_bound(self, value: Union[float, int]):

        if not value:
            value = ModelConstants.REACTION_UPPER_BOUND

        value = (self.lower_bound, value)

        self.coefficient.coefficients = value

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------

    def yield_genes(self) -> Generator['Gene', None, None]:

        return generator(self.genes)

    def yield_metabolites(self) -> Generator['Metabolite', None, None]:

        return generator(self.metabolites)

    def yield_reactants(self) -> Generator['Metabolite', None, None]:

        return generator(self.reactants)

    def yield_products(self) -> Generator['Metabolite', None, None]:

        return generator(self.products)

    # -----------------------------------------------------------------------------
    # Polymorphic constructors
    # -----------------------------------------------------------------------------

    @classmethod
    def from_stoichiometry(cls,
                           stoichiometry: Dict['Metabolite', Union[float, int]],
                           identifier: Any,
                           bounds: Tuple[Union[float, int], Union[float, int]] = None,
                           gpr: Expression = None,
                           model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                           **kwargs) -> 'Reaction':

        instance = cls(identifier=identifier,
                       bounds=bounds,
                       stoichiometry=stoichiometry,
                       gpr=gpr,
                       model=model,
                       **kwargs)

        return instance

    @classmethod
    def from_gpr_string(cls,
                        gpr_string: str,
                        identifier: Any,
                        bounds: Tuple[Union[float, int], Union[float, int]] = None,
                        stoichiometry: Dict['Metabolite', Union[float, int]] = None,
                        model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                        **kwargs) -> 'Reaction':

        try:

            symbolic = parse_expression(gpr_string)

            genes = variables_from_symbolic(symbolic=symbolic, types=('gene', ), model=model)

            expression = Expression(symbolic=symbolic, variables=genes)

        except SyntaxError as exc:

            expression_warning(f'{gpr_string} cannot be parsed')

            raise exc

        instance = cls(identifier=identifier,
                       bounds=bounds,
                       stoichiometry=stoichiometry,
                       gpr=expression,
                       model=model,
                       **kwargs)

        return instance

    @classmethod
    def from_gpr_expression(cls,
                            gpr: Expression,
                            identifier: Any,
                            bounds: Tuple[Union[float, int], Union[float, int]] = None,
                            stoichiometry: Dict['Metabolite', Union[float, int]] = None,
                            model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                            **kwargs) -> 'Reaction':

        if not isinstance(gpr, Expression):
            raise TypeError(f'expression must be an {Expression} object')

        instance = cls(identifier=identifier,
                       bounds=bounds,
                       stoichiometry=stoichiometry,
                       gpr=gpr,
                       model=model,
                       **kwargs)

        return instance

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------

    def ko(self, minimum_coefficient: Union[int, float] = 0.0, history=True):

        return self.coefficient.ko(minimum_coefficient=minimum_coefficient, history=history)

    def update(self,
               bounds: Tuple[Union[float, int], Union[float, int]] = None,
               stoichiometry: Dict['Metabolite', Union[float, int]] = None,
               gpr: Expression = None,
               **kwargs):

        super(Reaction, self).update(**kwargs)

        if bounds is not None:
            self.bounds = bounds

        if gpr is not None:
            self.gpr = gpr

        if stoichiometry is not None:
            self.stoichiometry = stoichiometry

    def replace_stoichiometry(self,
                              stoichiometry: Dict['Metabolite', Union[float, int]],
                              remove_orphans_from_model: bool = True,
                              history=True):

        if history:
            self.history.queue_command(undo_func=self.replace_stoichiometry,
                                       undo_kwargs={'stoichiometry': self.stoichiometry,
                                                    'remove_orphans_from_model': remove_orphans_from_model,
                                                    'history': False},
                                       func=self.replace_stoichiometry,
                                       kwargs={'stoichiometry': stoichiometry,
                                               'remove_orphans_from_model': remove_orphans_from_model,
                                               'history': history})

        if not stoichiometry and not self._stoichiometry:
            self._stoichiometry = {}

            return

        self.remove_metabolites(self.yield_metabolites(),
                                remove_orphans_from_model=remove_orphans_from_model,
                                history=False)

        self.add_metabolites(stoichiometry,
                             history=False)

    def add_metabolites(self,
                        stoichiometry: Dict['Metabolite', Union[float, int]],
                        history=True):

        to_add = []

        for metabolite, coefficient in stoichiometry.items():

            # noinspection PyProtectedMember
            metabolite._reactions[self.id] = self

            self._stoichiometry[metabolite] = coefficient

            to_add.append(metabolite)

        if self.model:
            # the add interface will add all metabolites to the simulators. The simulators will retrieve the new
            # constraints for each metabolite and replace it in the LP
            self.model.add(to_add, 'metabolite', comprehensive=False, history=False)

            notification = Notification(content=(self, ),
                                        content_type='reactions',
                                        action='add')

            self.model.notify(notification)

        if history:
            self.history.queue_command(undo_func=self.remove_metabolites,
                                       undo_kwargs={'metabolites': list(stoichiometry.keys()),
                                                    'remove_orphans_from_model': True,
                                                    'history': False},
                                       func=self.add_metabolites,
                                       kwargs={'stoichiometry': stoichiometry,
                                               'history': history})

    def remove_metabolites(self,
                           metabolites: List['Metabolite'],
                           remove_orphans_from_model: bool = True,
                           history=True):

        if isinstance(metabolites, dict):
            metabolites = list(metabolites.values())

        orphan_mets = []
        old_stoichiometry = {}

        for metabolite in metabolites:

            # noinspection PyProtectedMember
            del metabolite._reactions[self.id]

            # noinspection PyProtectedMember
            if not metabolite._reactions:
                orphan_mets.append(metabolite)

            old_stoichiometry[metabolite] = self._stoichiometry[metabolite]

            del self._stoichiometry[metabolite]

        if self.model:

            # the add interface will add all metabolites to the simulators. The simulators will retrieve the new
            # constraints for each metabolite and replace it in the LP

            if orphan_mets and remove_orphans_from_model:
                self.model.remove(orphan_mets, 'metabolite', remove_orphans=False, history=False)

            notification = Notification(content=(self, ),
                                        content_type='reactions',
                                        action='add')

            self.model.notify(notification)

        if history:
            self.history.queue_command(undo_func=self.add_metabolites,
                                       undo_kwargs={'stoichiometry': old_stoichiometry,
                                                    'history': False},
                                       func=self.remove_metabolites,
                                       kwargs={'metabolites': metabolites,
                                               'remove_orphans_from_model': remove_orphans_from_model,
                                               'history': history})

    def add_gpr(self, gpr: Expression, history=True):

        if not isinstance(gpr, Expression):
            raise TypeError(f'expression must be an {Expression} object. '
                            f'To set None, provide an empty Expression()')

        if history:
            self.history.queue_command(undo_func=self.remove_gpr,
                                       undo_kwargs={'remove_orphans': True,
                                                    'history': False},
                                       func=self.add_gpr,
                                       kwargs={'gpr': gpr,
                                               'history': history})

        if self.gpr:
            self.remove_gpr(history=False)

        to_add = []

        for gene in gpr.variables.values():
            gene.update(reactions={self.id: self},
                        model=self.model)

            to_add.append(gene)

        self._gpr = gpr

        if self.model:
            self.model.add(to_add, 'gene', comprehensive=False, history=False)

            notification = Notification(content=(self,),
                                        content_type='gprs',
                                        action='add')

            self.model.notify(notification)

    def remove_gpr(self, remove_orphans: bool = True, history=True):

        if history:
            self.history.queue_command(undo_func=self.add_gpr,
                                       undo_kwargs={'gpr': self.gpr,
                                                    'history': False},
                                       func=self.remove_gpr,
                                       kwargs={'remove_orphans': remove_orphans,
                                               'history': history})

        to_remove = []

        for gene in self.yield_genes():

            # noinspection PyProtectedMember
            del gene._reactions[self.id]

            # noinspection PyProtectedMember
            if not gene._reactions:
                to_remove.append(gene)

        self._gpr = Expression()

        if self.model:

            if remove_orphans:
                self.model.remove(to_remove, 'gene', remove_orphans=False, history=False)

            notification = Notification(content=(self,),
                                        content_type='gprs',
                                        action='add')

            self.model.notify(notification)
