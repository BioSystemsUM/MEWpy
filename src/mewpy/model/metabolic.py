from typing import TYPE_CHECKING, Any, Union, Generator, Dict, List, Tuple, Set

from mewpy.model.model import Model
from mewpy.mew.lp import Notification
from mewpy.util.history import recorder
from mewpy.util.serialization import serialize
from mewpy.util.utilities import iterable, generator

if TYPE_CHECKING:
    from mewpy.mew.algebra import Expression
    from mewpy.mew.variables import Gene, Metabolite, Reaction


# TODO: methods stubs
class MetabolicModel(Model, model_type='metabolic', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier: Any,
                 compartments: Dict[str, str] = None,
                 genes: Dict[str, 'Gene'] = None,
                 metabolites: Dict[str, 'Metabolite'] = None,
                 objective: Dict['Reaction', Union[float, int]] = None,
                 reactions: Dict[str, 'Reaction'] = None,
                 **kwargs):

        """
        A metabolic model contains essentially reactions, metabolites and genes.
        The metabolic model can be also also associated with a given objective for the simulation/analysis.
        The metabolic model can be loaded with compartments, although these can be inferred from the available
        metabolites.

        Other information is retrieved from these attributes, namely demand, exchange and sink reactions as well as gprs

        The metabolic model, as with other models, provides a clean interface for manipulation with the add, remove and
        update methods.

        :param identifier: identifier, e.g. iMC1010
        :param compartments: a dictionary with additional compartments not encoded in the metabolites
        :param genes: a dictionary with Gene objects. See variables.Gene for more info
        :param metabolites: a dictionary with Metabolite objects. See variables.Metabolite for more info
        :param objective: a dictionary with the Reaction objects that must be considered objective functions of
        the simulations together with the respective coefficients
        :param reactions: a dictionary with Reaction objects. See variables.Reaction for more info
        """

        # compartments attribute can be shared across the children, thus name mangling
        self.__compartments = {}
        self._genes = {}
        self._metabolites = {}
        self._objective = {}
        self._reactions = {}

        super().__init__(identifier,
                         **kwargs)

        # the setters will handle adding and removing variables to the correct containers
        self.compartments = compartments
        self.genes = genes
        self.metabolites = metabolites
        self.objective = objective
        self.reactions = reactions

    # -----------------------------------------------------------------------------
    # Model type manager
    # -----------------------------------------------------------------------------

    @serialize('types', None)
    @property
    def types(self):

        # noinspection PyUnresolvedReferences
        _types = {MetabolicModel.model_type}

        _types.update(super(MetabolicModel, self).types)

        return _types

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------

    @serialize('genes', 'genes', '_genes')
    @property
    def genes(self) -> Dict[str, 'Gene']:
        return self._genes.copy()

    @serialize('metabolites', 'metabolites', '_metabolites')
    @property
    def metabolites(self) -> Dict[str, 'Metabolite']:
        return self._metabolites.copy()

    @serialize('objective', 'objective', '_objective')
    @property
    def objective(self) -> Dict['Reaction', Union[float, int]]:
        return self._objective.copy()

    @serialize('reactions', 'reactions', '_reactions')
    @property
    def reactions(self) -> Dict[str, 'Reaction']:
        return self._reactions.copy()

    @property
    def compartments(self) -> Dict[str, str]:

        compartments = {met.compartment: self.__compartments.get(met.compartment, '')
                        for met in self.yield_metabolites()
                        if met.compartment is not None}

        compartments.update(self.__compartments)

        compartments.update(super(MetabolicModel, self).compartments)

        return compartments

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------

    @compartments.setter
    @recorder
    def compartments(self, value: Dict[str, str]):

        if not value:
            value = {}

        self.__compartments.update(value)

    @genes.setter
    @recorder
    def genes(self, value: Dict[str, 'Gene']):

        if not value:
            value = {}

        self.remove(list(self.yield_genes()), 'gene', history=False)
        self.add(list(value.values()), 'gene', history=False)

    @metabolites.setter
    @recorder
    def metabolites(self, value: Dict[str, 'Metabolite']):

        if not value:
            value = {}

        self.remove(list(self.yield_metabolites()), 'metabolite', history=False)
        self.add(list(value.values()), 'metabolite', history=False)

    @objective.setter
    @recorder
    def objective(self, value: Dict['Reaction', Union[float, int]]):

        if not value:
            value = {}

        if isinstance(value, str):

            value = {self.get(value): 1}

        elif hasattr(value, 'types'):

            value = {value: 1}

        elif isinstance(value, dict):

            value = {self.get(var, var): val for var, val in value.items()}

        else:
            raise ValueError(f'{value} is not a valid objective')

        self._objective = value

        linear_obj = {var.id: coef for var, coef in self._objective.items()}

        content = {'linear': linear_obj,
                   'minimize': False}

        notification = Notification(content=content,
                                    content_type='objectives',
                                    action='set')

        self.notify(notification)

    @reactions.setter
    @recorder
    def reactions(self, value: Dict[str, 'Reaction']):

        if not value:
            value = {}

        self.remove(list(self.yield_reactions()), 'reaction', history=False)
        self.add(list(value.values()), 'reaction', history=False)

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------

    @property
    def external_compartment(self) -> Union[str, None]:

        if not self.compartments:
            return

        if not self._reactions:
            return

        boundary_compartments = {}

        for rxn in self.yield_reactions():

            if rxn.boundary:

                for compartment in rxn.compartments:
                    boundary_compartments[compartment] = boundary_compartments.get(compartment, 0) + 1

        external_compartment = None
        highest_count = 0

        for compartment, count in boundary_compartments.items():

            if count > highest_count:
                external_compartment = compartment
                highest_count = count

        return external_compartment

    def _get_boundaries(self):

        external_compartment = self.external_compartment

        if external_compartment is None:
            return {}, {}, {}

        all_boundaries = [rxn for rxn_id, rxn in self._reactions.items()
                          if rxn.boundary]

        exchanges = {}
        sinks = {}
        demands = {}

        for variable in all_boundaries:

            if variable.types == {'reaction'}:

                if external_compartment in variable.compartments:
                    exchanges[variable.id] = variable

                else:
                    if variable.reversibility:
                        sinks[variable.id] = variable

                    else:
                        demands[variable.id] = variable

        return exchanges, sinks, demands

    @property
    def demands(self) -> Dict[str, 'Reaction']:

        _, _, demands = self._get_boundaries()

        return demands

    @property
    def exchanges(self) -> Dict[str, 'Reaction']:

        exchanges, _, _ = self._get_boundaries()

        return exchanges

    @property
    def sinks(self) -> Dict[str, 'Reaction']:

        _, sinks, _ = self._get_boundaries()

        return sinks

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------

    def yield_compartments(self) -> Generator[str, None, None]:

        return generator(self.compartments)

    def yield_demands(self) -> Generator['Reaction', None, None]:

        return generator(self.demands)

    def yield_exchanges(self) -> Generator['Reaction', None, None]:

        return generator(self.exchanges)

    def yield_genes(self) -> Generator['Gene', None, None]:

        return generator(self._genes)

    def yield_gprs(self) -> Generator['Expression', None, None]:

        return (value.gpr for value in self._reactions.values())

    def yield_metabolites(self) -> Generator['Metabolite', None, None]:

        return generator(self._metabolites)

    def yield_reactions(self) -> Generator['Reaction', None, None]:

        return generator(self._reactions)

    def yield_sinks(self) -> Generator['Reaction', None, None]:

        return generator(self.sinks)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------

    def get(self, identifier: Any, default=None) -> Union['Gene', 'Metabolite', 'Reaction']:

        if identifier in self._metabolites:
            return self._metabolites[identifier]

        elif identifier in self._reactions:
            return self._reactions[identifier]

        elif identifier in self._genes:
            return self._genes[identifier]

        else:
            return super(MetabolicModel, self).get(identifier=identifier, default=default)

    def add(self,
            variables: Union[List[Union['Gene', 'Metabolite', 'Reaction']],
                             Tuple[Union['Gene', 'Metabolite', 'Reaction']],
                             Set[Union['Gene', 'Metabolite', 'Reaction']]],
            *types: str,
            comprehensive: bool = True,
            history=True):

        variables = iterable(variables)

        if not types:
            types = [var.types for var in variables]

        elif len(types) == len(variables):
            types = [{_type} if isinstance(_type, str) else set(_type) for _type in types]

        elif len(types) != len(variables):
            types = [set(types) for var in variables]

        reactions = []
        new_variables = []
        new_types = []

        for var_types, var in zip(types, variables):

            if 'gene' in var_types:
                self._add_gene(var)
                var_types.remove('gene')

            elif 'metabolite' in var_types:
                self._add_metabolite(var)
                var_types.remove('metabolite')

            elif 'reaction' in var_types:
                self._add_reaction(var, comprehensive=comprehensive)
                reactions.append(var)
                var_types.remove('reaction')

            if var_types:
                new_types.append(var_types)
                new_variables.append(var)

        if reactions:

            notification = Notification(content=reactions,
                                        content_type='reactions',
                                        action='add')

            self.notify(notification)

            notification = Notification(content=reactions,
                                        content_type='gprs',
                                        action='add')

            self.notify(notification)

        if history:
            self.history.queue_command(undo_func=self.remove,
                                       undo_kwargs={'variables': variables,
                                                    'remove_orphans': True,
                                                    'history': False},
                                       func=self.add,
                                       kwargs={'variables': variables,
                                               'comprehensive': comprehensive,
                                               'history': history})

        super(MetabolicModel, self).add(new_variables,
                                        *new_types,
                                        comprehensive=comprehensive,
                                        history=False)

    def remove(self,
               variables: Union[List[Union['Gene', 'Metabolite', 'Reaction']],
                                Tuple[Union['Gene', 'Metabolite', 'Reaction']],
                                Set[Union['Gene', 'Metabolite', 'Reaction']]],
               *types: str,
               remove_orphans: bool = False, history=True):

        variables = iterable(variables)

        if not types:
            types = [var.types for var in variables]

        elif len(types) == len(variables):
            types = [{_type} if isinstance(_type, str) else set(_type) for _type in types]

        elif len(types) != len(variables):
            types = [set(types) for var in variables]

        reactions = []
        metabolites = []
        new_variables = []
        new_types = []

        for var_types, var in zip(types, variables):

            if 'gene' in var_types:
                self._remove_gene(var)
                var_types.remove('gene')

            elif 'metabolite' in var_types:
                metabolites.append(var)
                self._remove_metabolite(var)
                var_types.remove('metabolite')

            elif 'reaction' in var_types:
                reactions.append(var)
                self._remove_reaction(var)
                var_types.remove('reaction')

            if var_types:
                new_types.append(var_types)
                new_variables.append(var)

        if reactions:

            notification = Notification(content=reactions,
                                        content_type='reactions',
                                        action='remove')

            self.notify(notification)

            notification = Notification(content=reactions,
                                        content_type='gprs',
                                        action='remove')

            self.notify(notification)

            if remove_orphans:
                orphan_mets, _ = self._remove_metabolic_orphans(reactions)

                notification = Notification(content=orphan_mets,
                                            content_type='metabolites',
                                            action='remove')

                self.notify(notification)

        # metabolites take precedence since linear coefficients are added to linear problems
        # by the reactions stoichiometry dictionary.
        # So, if metabolites are being removed as a result of removing reactions and its metabolites
        # all together from the model, the linear coefficients obtained from the reactions should be the last thing to
        # persist in the lp.
        if metabolites:
            notification = Notification(content=metabolites,
                                        content_type='metabolites',
                                        action='remove')

            self.notify(notification)

        if history:
            self.history.queue_command(undo_func=self.add,
                                       undo_kwargs={'variables': variables,
                                                    'comprehensive': True,
                                                    'history': False},
                                       func=self.remove,
                                       kwargs={'variables': variables,
                                               'remove_orphans': remove_orphans,
                                               'history': history})

        super(MetabolicModel, self).remove(new_variables,
                                           *new_types,
                                           remove_orphans=remove_orphans,
                                           history=False)

    def update(self,
               compartments: Dict[str, str] = None,
               objective: Dict['Reaction', Union[float, int]] = None,
               variables: Union[List[Union['Gene', 'Metabolite', 'Reaction']],
                                Tuple[Union['Gene', 'Metabolite', 'Reaction']],
                                Set[Union['Gene', 'Metabolite', 'Reaction']]] = None,
               **kwargs):

        if compartments is not None:
            self.compartments = compartments

        if variables is not None:
            self.add(variables=variables)

        if objective is not None:
            self.objective = objective

        super(MetabolicModel, self).update(**kwargs)

    # -----------------------------------------------------------------------------
    # Helper functions for the operations/Manipulations
    # -----------------------------------------------------------------------------

    def _add_reaction(self, reaction, comprehensive=True):

        # adding a reaction is regularly a in-depth append method, as both metabolites and genes associated with the
        # reaction are also regularly added to the model. Although, this behaviour can be avoided passing
        # comprehensive=False

        if reaction.id not in self._reactions:

            if comprehensive:

                for metabolite in reaction.yield_metabolites():
                    self._add_metabolite(metabolite)

                for gene in reaction.yield_genes():
                    self._add_gene(gene)

            self._add_variable_to_container(reaction, self._reactions)

    def _add_metabolite(self, metabolite):

        if metabolite.id not in self._metabolites:
            self._add_variable_to_container(metabolite, self._metabolites)

    def _add_gene(self, gene):

        if gene.id not in self._genes:
            self._add_variable_to_container(gene, self._genes)

    def _remove_reaction(self, reaction, remove_orphans=False):

        if reaction.id in self._reactions:
            self._remove_variable_from_container(reaction, self._reactions)

    def _remove_metabolite(self, metabolite):

        if metabolite.id in self._metabolites:
            self._remove_variable_from_container(metabolite, self._metabolites)

    def _remove_gene(self, gene):

        if gene.id in self._genes:
            self._remove_variable_from_container(gene, self._genes)

    def _remove_metabolic_orphans(self, reactions):

        orphan_mets = self._get_orphans(to_remove=reactions,
                                        first_container='metabolites',
                                        second_container='reactions')

        orphan_genes = self._get_orphans(to_remove=reactions,
                                         first_container='genes',
                                         second_container='reactions')

        if orphan_mets:

            for met in orphan_mets:
                self._remove_metabolite(met)

        if orphan_genes:
            for gene in orphan_genes:
                self._remove_gene(gene)

        return orphan_mets, orphan_genes
