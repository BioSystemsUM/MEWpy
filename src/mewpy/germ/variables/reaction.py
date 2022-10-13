from typing import Any, Dict, Union, TYPE_CHECKING, Tuple, Generator, List

from mewpy.germ.algebra import Expression, parse_expression
from mewpy.util.utilities import generator
from mewpy.germ.models.serialization import serialize
from mewpy.util.history import recorder
from mewpy.util.constants import ModelConstants
from mewpy.io.engines.engines_utils import expression_warning
from .variable import Variable, variables_from_symbolic

if TYPE_CHECKING:
    from .metabolite import Metabolite
    from .gene import Gene
    from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel


def _bounds_setter(instance, old_bounds):
    instance._bounds = tuple(old_bounds)

    if instance.model:
        instance.model.notify()


class Reaction(Variable, variable_type='reaction', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier: Any,
                 bounds: Tuple[float, float] = None,
                 stoichiometry: Dict['Metabolite', float] = None,
                 gpr: Expression = None,
                 **kwargs):

        """
        Reactions are the primary variables of metabolic models.
        A reaction holds information for bounds, metabolites, stoichiometry and
        gene-protein-reaction rules (the associated between reactions, enzymes/proteins and genes)

        A reaction usually contains the stoichiometry of the metabolites that take part in the reaction.
        A reaction usually contains the GPR rule/logic of the associated genes
        A reaction usually contains the bounds for the fluxes values that can take during simulation

        Dynamic information can be obtained from the previous attributes, such as metabolites, reactants, products,
        genes, boundary, reversibility, charge and mass balance, equation, compartments, etc.

        :param identifier: identifier, e.g. PYK
        :param bounds: tuple having both lower and upper bounds, respectively the minimum and maximum coefficient
        that a reaction can have
        :param stoichiometry: dictionary of metabolites objects and their corresponding stoichiometric value
        :param gpr: an expression object with the boolean logic of the associated genes.
        Genes objects must be associated with the expression object.
        """
        # the coefficient bounds initializer sets minimum and maximum bounds of MEWPY_LB and MEWPY_UB
        if not bounds:
            bounds = ModelConstants.REACTION_LOWER_BOUND, ModelConstants.REACTION_UPPER_BOUND

        elif len(bounds) == 1:
            bounds = (bounds[0], bounds[0])

        elif len(bounds) == 2:
            bounds = tuple(bounds)

        else:
            raise ValueError('Bounds must be a tuple of length 1 or 2')

        if not gpr:
            gpr = Expression()

        if not stoichiometry:
            stoichiometry = {}

        self._bounds = bounds
        self._gpr = Expression()
        self._stoichiometry = {}

        super().__init__(identifier,
                         **kwargs)

        self.replace_stoichiometry(stoichiometry, history=False)

        # The gpr must be an expression object. The expression object properties are used to check whether the genes
        # container is correct.
        self.add_gpr(gpr, history=False)

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
        return f'{self.id} || {self.equation}'

    def _reaction_to_html(self):
        """
        It returns a html representation.
        """
        html_dict = {'Equation': self.equation,
                     'Bounds': self.bounds,
                     'Reversibility': self.reversibility,
                     'Metabolites': ', '.join(self.metabolites),
                     'Boundary': self.boundary,
                     'GPR': self.gene_protein_reaction_rule,
                     'Genes': ', '.join(self.genes),
                     'Compartments': ', '.join(self.compartments),
                     'Charge balance': self.charge_balance,
                     'Mass balance': self.mass_balance}

        return html_dict

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------
    @serialize('bounds', 'bounds', '_bounds')
    @property
    def bounds(self) -> Tuple[float, float]:
        """
        Bounds of the reaction. The lower and upper bounds of the reaction
        :return: bounds tuple
        """
        return self._bounds

    @serialize('stoichiometry', 'stoichiometry', '_stoichiometry')
    @property
    def stoichiometry(self) -> Dict['Metabolite', Union[int, float]]:
        """
        Stoichiometry of the reaction. A dictionary of metabolites and their corresponding stoichiometric value
        :return: stoichiometry as a dictionary
        """
        return self._stoichiometry.copy()

    @serialize('gpr', 'gpr', '_gpr')
    @property
    def gpr(self) -> Expression:
        """
        Gene-protein-reaction rule of the reaction. An expression object that holds the boolean logic of the associated
        genes. Genes objects must be associated with the expression object.
        :return: expression object
        """
        return self._gpr

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------
    @bounds.setter
    @recorder
    def bounds(self, value: Union[Tuple[float, float], Tuple[float], float]):
        """
        The setter for the bounds of the reaction.
        It accepts a tuple of lower and upper bounds.
        This setter replaces the current bounds with a new one.
        It also handles the linear problems associated to this reaction.
        :param value: tuple of lower and upper bounds
        :return:
        """
        if value is None:
            value = ModelConstants.REACTION_LOWER_BOUND, ModelConstants.REACTION_UPPER_BOUND

        elif isinstance(value, (int, float)):
            value = (value, value)

        elif len(value) == 1:
            value = (value[0], value[0])

        elif len(value) == 2:
            value = tuple(value)

        else:
            raise ValueError('Invalid value for bounds')

        self._bounds = value

        if self.model:
            self.model.notify()

    @stoichiometry.setter
    @recorder
    def stoichiometry(self, value: Dict['Metabolite', Union[int, float]]):
        """
        The setter for the stoichiometry of the reaction.
        It accepts a dictionary of metabolites and their corresponding stoichiometric value.
        This setter replaces the current stoichiometry with a new one.
        It also handles the linear problems associated to this reaction.
        :param value: stoichiometry as a dictionary
        :return: stoichiometry as a dictionary
        """

        self.replace_stoichiometry(value, history=False)

    @gpr.setter
    @recorder
    def gpr(self, value: Expression):
        """
        The setter for the gpr of the reaction.
        It accepts an expression object that holds the boolean logic of the associated genes.
        Genes objects must be associated with the expression object.
        This setter replaces the current gpr with a new one.
        It also handles the linear problems associated to this reaction.
        :param value: expression object
        :return:
        """

        self.remove_gpr(history=False)
        self.add_gpr(value, history=False)

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------
    @property
    def metabolites(self) -> Dict[str, 'Metabolite']:
        """
        Metabolites of the reaction
        :return: dictionary of metabolites
        """
        return {met.id: met for met in self._stoichiometry}

    @property
    def products(self) -> Dict[str, 'Metabolite']:
        """
        Products of the reaction
        :return: dictionary of products
        """
        return {met.id: met for met, st in self._stoichiometry.items()
                if st > 0.0}

    @property
    def reactants(self) -> Dict[str, 'Metabolite']:
        """
        Reactants of the reaction
        :return: dictionary of reactants
        """
        return {met.id: met for met, st in self._stoichiometry.items()
                if st < 0.0}

    @property
    def compartments(self) -> set:
        """
        Compartments of the reaction
        :return: set of compartments
        """
        return {met.compartment
                for met in self.yield_metabolites()
                if met.compartment is not None}

    @property
    def lower_bound(self) -> Union[int, float]:
        """
        Lower bound of the reaction
        :return: lower bound as a float
        """
        return self.bounds[0]

    @property
    def upper_bound(self) -> Union[int, float]:
        """
        Upper bound of the reaction
        :return: upper bound as a float
        """
        return self.bounds[1]

    @property
    def reversibility(self) -> bool:
        """
        Reversibility of the reaction
        :return: reversibility as a boolean
        """
        return self.lower_bound < - ModelConstants.TOLERANCE and self.upper_bound > ModelConstants.TOLERANCE

    @property
    def equation(self) -> str:
        """
        Equation of the reaction
        :return: equation as a string
        """

        reactants = ' + '.join([f'{abs(st)} {met.id}' for met, st in self._stoichiometry.items() if st < 0.0])
        products = ' + '.join([f'{abs(st)} {met.id}' for met, st in self._stoichiometry.items() if st > 0.0])

        if self.reversibility:
            return f'{reactants} <-> {products}'
        else:
            return f'{reactants} -> {products}'

    @property
    def genes(self) -> Dict[str, 'Gene']:
        """
        Genes of the reaction
        :return: dictionary of genes
        """
        return self.gpr.variables.copy()

    @property
    def gene_protein_reaction_rule(self) -> str:
        """
        Gene-protein-reaction rule of the reaction
        :return: gene-protein-reaction rule as a string
        """
        return self.gpr.to_string()

    @property
    def boundary(self) -> bool:
        """
        Boundary of the reaction
        :return: boundary as a boolean
        """
        return not (self.reactants and self.products)

    @property
    def charge_balance(self) -> Dict[str, Union[int, float]]:
        """
        Charge balance of the reaction
        :return: charge balance as a dictionary
        """
        return {'reactants': sum([self._stoichiometry[met] * met.charge for met in self.yield_reactants()]),
                'products': sum([self._stoichiometry[met] * met.charge for met in self.yield_products()])}

    @property
    def mass_balance(self) -> Dict[str, Union[int, float]]:
        """
        Mass balance of the reaction
        :return: mass balance as a dictionary
        """

        balance = {}

        for metabolite, st in self._stoichiometry.items():

            for atom, count in metabolite.atoms.items():
                amount = count * st

                balance[atom] = balance.get(atom, 0) + amount

        return balance

    # -----------------------------------------------------------------------------
    # Dynamic attributes setters
    # -----------------------------------------------------------------------------
    @lower_bound.setter
    @recorder
    def lower_bound(self, value: Union[float, int]):
        """
        The setter for the lower bound of the reaction.
        It accepts a float or an integer.
        This setter replaces the current lower bound with a new one.
        It also handles the linear problems associated to this reaction.
        :param value: lower bound
        :return:
        """

        if not value:
            value = ModelConstants.REACTION_LOWER_BOUND

        value = (value, self.upper_bound)

        self._bounds = value

        if self.model:
            self.model.notify()

    @upper_bound.setter
    @recorder
    def upper_bound(self, value: Union[float, int]):
        """
        The setter for the upper bound of the reaction.
        It accepts a float or an integer.
        This setter replaces the current upper bound with a new one.
        It also handles the linear problems associated to this reaction.
        :param value: upper bound
        :return:
        """

        if not value:
            value = ModelConstants.REACTION_UPPER_BOUND

        self._bounds = (self.lower_bound, value)

        if self.model:
            self.model.notify()

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------
    def yield_genes(self) -> Generator['Gene', None, None]:
        """
        Generator of genes of the reaction
        :return: generator of genes
        """
        return generator(self.genes)

    def yield_metabolites(self) -> Generator['Metabolite', None, None]:
        """
        Generator of metabolites of the reaction
        :return: generator of metabolites
        """
        return generator(self.metabolites)

    def yield_reactants(self) -> Generator['Metabolite', None, None]:
        """
        Generator of reactants of the reaction
        :return: generator of reactants
        """
        return generator(self.reactants)

    def yield_products(self) -> Generator['Metabolite', None, None]:
        """
        Generator of products of the reaction
        :return: generator of products
        """
        return generator(self.products)

    # -----------------------------------------------------------------------------
    # Polymorphic constructors
    # -----------------------------------------------------------------------------
    @classmethod
    def from_stoichiometry(cls,
                           identifier: Any,
                           stoichiometry: Dict['Metabolite', Union[float, int]],
                           bounds: Tuple[Union[float, int], Union[float, int]] = None,
                           gpr: Expression = None,
                           model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                           **kwargs) -> 'Reaction':
        """
        A reaction is defined by its stoichiometry, establishing the relationship between the metabolites and the
        reaction. The stoichiometry is a dictionary of metabolites and their stoichiometric coefficients.

        A reaction can be created from a stoichiometry dictionary.

        :param identifier: identifier of the reaction
        :param stoichiometry: stoichiometry dictionary
        :param bounds: tuple of lower and upper bounds
        :param gpr: gene-protein-reaction rule
        :param model: model of the reaction
        :param kwargs: additional arguments
        :return: a new reaction object
        """
        instance = cls(identifier=identifier,
                       bounds=bounds,
                       stoichiometry=stoichiometry,
                       gpr=gpr,
                       model=model,
                       **kwargs)

        return instance

    @classmethod
    def from_gpr_string(cls,
                        identifier: Any,
                        rule: str,
                        bounds: Tuple[Union[float, int], Union[float, int]] = None,
                        stoichiometry: Dict['Metabolite', Union[float, int]] = None,
                        model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                        **kwargs) -> 'Reaction':
        """
        A reaction is regularly associated with genes. The GPR rule is a boolean expression that relates the genes
        to the reaction. The GPR rule is a string that can be parsed into a boolean expression.

        A reaction can be created from a gene-protein-reaction rule string.

        :param identifier: identifier of the reaction
        :param rule: gene-protein-reaction rule string
        :param bounds: tuple of lower and upper bounds
        :param stoichiometry: stoichiometry dictionary
        :param model: model of the reaction
        :param kwargs: additional arguments
        :return: a new reaction object
        """
        try:

            symbolic = parse_expression(rule)

            genes = variables_from_symbolic(symbolic=symbolic, types=('gene',), model=model)

            expression = Expression(symbolic=symbolic, variables=genes)

        except SyntaxError as exc:

            expression_warning(f'{rule} cannot be parsed')

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
                            identifier: Any,
                            gpr: Expression,
                            bounds: Tuple[Union[float, int], Union[float, int]] = None,
                            stoichiometry: Dict['Metabolite', Union[float, int]] = None,
                            model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                            **kwargs) -> 'Reaction':
        """
        A reaction is regularly associated with genes. The GPR rule is a boolean expression that relates the genes
        to the reaction. The GPR rule is a string that can be parsed into a boolean expression.

        A reaction can be created from a gene-protein-reaction expression.

        :param identifier: identifier of the reaction
        :param gpr: gene-protein-reaction expression
        :param bounds: tuple of lower and upper bounds
        :param stoichiometry: stoichiometry dictionary
        :param model: model of the reaction
        :param kwargs: additional arguments
        :return: a new reaction object
        """
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
        """
        Knockout the reaction by setting the lower and upper bounds to zero.
        It also handles the linear problems associated to this reaction.

        :param minimum_coefficient: minimum coefficient of the reaction
        :param history: record the operation
        :return:
        """
        old_bounds = tuple(self.bounds)

        self._bounds = (minimum_coefficient, minimum_coefficient)

        if self.model:
            self.model.notify()

        if history:
            self.history.queue_command(undo_func=_bounds_setter,
                                       undo_kwargs={'instance': self,
                                                    'old_bounds': old_bounds},
                                       func=self.ko,
                                       kwargs={'minimum_coefficient': minimum_coefficient,
                                               'history': False})
        return

    def update(self,
               bounds: Tuple[Union[float, int], Union[float, int]] = None,
               stoichiometry: Dict['Metabolite', Union[float, int]] = None,
               gpr: Expression = None,
               **kwargs):
        """
        Update the reaction with new bounds, stoichiometry and gene-protein-reaction rule.

        It also handles the linear problems associated to this reaction.

        Note that, some update operations are not registered in history.
        It is strongly advisable to use update outside history context manager

        :param bounds: tuple of lower and upper bounds
        :param stoichiometry: stoichiometry dictionary
        :param gpr: gene-protein-reaction rule
        :param kwargs: additional arguments
        :return:
        """

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
        """
        Replace the stoichiometry of the reaction with a new one.

        It also handles the linear problems associated to this reaction.

        :param stoichiometry: stoichiometry dictionary
        :param remove_orphans_from_model: whether to remove orphan metabolites from the model
        :param history: whether to record the operation
        :return:
        """
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

        self.remove_metabolites(list(self.yield_metabolites()),
                                remove_orphans_from_model=remove_orphans_from_model,
                                history=False)

        self.add_metabolites(stoichiometry,
                             history=False)

    def add_metabolites(self,
                        stoichiometry: Dict['Metabolite', Union[float, int]],
                        history=True):
        """
        Add metabolites to the reaction.
        The stoichiometry dictionary is updated with the new metabolites and coefficients.

        It also handles the linear problems associated to this reaction.
        :param stoichiometry: stoichiometry dictionary
        :param history: whether to record the operation
        :return:
        """
        to_add = []

        for metabolite, coefficient in stoichiometry.items():
            # noinspection PyProtectedMember
            metabolite._reactions[self.id] = self

            self._stoichiometry[metabolite] = coefficient

            to_add.append(metabolite)

        if self.model:
            # the add interface will add all metabolites to the simulators. The simulators will retrieve the new
            # constraints for each metabolite and replace it in the LP
            self.model.add(*to_add, comprehensive=False, history=False)

            self.model.notify()

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
        """
        Remove metabolites from the reaction.
        Metabolites and their coefficients are removed from the stoichiometry dictionary.

        It also handles the linear problems associated to this reaction.
        :param metabolites: list of metabolites
        :param remove_orphans_from_model: whether to remove orphan metabolites from the model
        :param history: whether to record the operation
        :return:
        """
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
                self.model.remove(*orphan_mets, remove_orphans=False, history=False)

            self.model.notify()

        if history:
            self.history.queue_command(undo_func=self.add_metabolites,
                                       undo_kwargs={'stoichiometry': old_stoichiometry,
                                                    'history': False},
                                       func=self.remove_metabolites,
                                       kwargs={'metabolites': metabolites,
                                               'remove_orphans_from_model': remove_orphans_from_model,
                                               'history': history})

    def add_gpr(self, gpr: Expression, history=True):
        """
        Add a gene-protein-reaction rule to the reaction.

        It also handles the linear problems associated to this reaction.
        :param gpr: gene-protein-reaction expression object
        :param history: whether to record the operation
        :return:
        """
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
            self.model.add(*to_add, comprehensive=False, history=False)

            self.model.notify()

    def remove_gpr(self, remove_orphans: bool = True, history=True):
        """
        Remove the gene-protein-reaction rule from the reaction.

        It also handles the linear problems associated to this reaction.
        :param remove_orphans: whether to remove orphan genes from the model
        :param history: whether to record the operation
        :return:
        """
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
                self.model.remove(*to_remove, remove_orphans=False, history=False)

            self.model.notify()
