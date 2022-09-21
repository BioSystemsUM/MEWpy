from typing import TYPE_CHECKING, Any, Union, Generator, Dict, List, Tuple, Set

from .model import Model
from mewpy.util.history import recorder
from mewpy.mew.models.serialization import serialize
from mewpy.util.utilities import generator

if TYPE_CHECKING:
    from mewpy.mew.algebra import Expression
    from mewpy.mew.variables import Gene, Metabolite, Reaction


class MetabolicModel(Model, model_type='metabolic', register=True, constructor=True, checker=True):
    """
    A mew metabolic model consists of a classic Genome-Scale Metabolic (GEM) model,
    containing reactions, metabolites and genes.

    GEM models are systems biology tools used to predict the phenotype of an organism or cellular community
    in range of environmental and genetic conditions.
    To perform phenotype prediction, a Mew metabolic model can be attached to several simulation methods:
        - FBA
        - pFBA
        - FVA
        - ...
    Thus, a mew metabolic model can be associated with a given objective function for the analysis of the model.

    The metabolic model can be loaded with compartments, although these can be inferred from the available
    metabolites.

    A mew metabolic model can hold additional information as follows:
        - demand reactions
        - exchange reactions
        - sink reactions
        - GPRs
        - External compartment

    The metabolic model, as with other models, provides a clean interface for manipulation with the add, remove and
    update methods. One can perform the following operations:
        - Add reactions, metabolites and genes
        - Remove reactions, metabolites and genes
        - Update the objective function
    """
    def __init__(self,
                 identifier: Any,
                 compartments: Dict[str, str] = None,
                 genes: Dict[str, 'Gene'] = None,
                 metabolites: Dict[str, 'Metabolite'] = None,
                 objective: Dict['Reaction', Union[float, int]] = None,
                 reactions: Dict[str, 'Reaction'] = None,
                 **kwargs):

        """
        A mew metabolic model consists of a classic Genome-Scale Metabolic (GEM) model,
        containing reactions, metabolites and genes.

        GEM models are systems biology tools used to predict the phenotype of an organism or cellular community
        in range of environmental and genetic conditions.
        To perform phenotype prediction, a Mew metabolic model can be attached to several simulation methods:
            - FBA
            - pFBA
            - FVA
            - ...
        Thus, a mew metabolic model can be associated with a given objective function for the analysis of the model.

        The metabolic model can be loaded with compartments, although these can be inferred from the available
        metabolites.

        A mew metabolic model can hold additional information as follows:
            - demand reactions
            - exchange reactions
            - sink reactions
            - GPRs
            - External compartment

        The metabolic model, as with other models, provides a clean interface for manipulation with the add, remove and
        update methods. One can perform the following operations:
            - Add reactions, metabolites and genes
            - Remove reactions, metabolites and genes
            - Update the objective function

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
        """
        Returns the types of the model
        :return: a set with the types of the model
        """
        _types = {MetabolicModel.model_type}

        _types.update(super(MetabolicModel, self).types)

        return _types

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------
    @serialize('genes', 'genes', '_genes')
    @property
    def genes(self) -> Dict[str, 'Gene']:
        """
        It returns a dictionary with the genes of the model. The key is the gene identifier and the value is the
        `Gene` object. To retrieve an iterator with the genes use `yield_genes` method.
        Note that the genes attribute retrieves a copy of the genes' container.
        To update the genes container set new `genes` or use `add` and `remove` methods.
        :return: a dictionary with the genes of the model
        """
        return self._genes.copy()

    @serialize('metabolites', 'metabolites', '_metabolites')
    @property
    def metabolites(self) -> Dict[str, 'Metabolite']:
        """
        It returns a dictionary with the metabolites of the model. The key is the metabolite identifier and the value is
        the `Metabolite` object. To retrieve an iterator with the metabolites use `yield_metabolites` method.
        Note that the metabolites attribute retrieves a copy of the metabolites' container.
        To update the metabolites container set new `metabolites` or use `add` and `remove` methods.
        :return: a dictionary with the metabolites of the model
        """
        return self._metabolites.copy()

    @serialize('objective', 'objective', '_objective')
    @property
    def objective(self) -> Dict['Reaction', Union[float, int]]:
        """
        It returns a dictionary with the objective functions of the model.
        The key is the `Reaction` object and the value is the respective coefficient.
        Note that the objective attribute retrieves a copy of the objective's container.
        To update the objective container set a new `objective` or use `update` method.
        :return: a dictionary with the objective functions of the model
        """
        return self._objective.copy()

    @serialize('reactions', 'reactions', '_reactions')
    @property
    def reactions(self) -> Dict[str, 'Reaction']:
        """
        It returns a dictionary with the reactions of the model. The key is the reaction identifier and the value is
        the `Reaction` object. To retrieve an iterator with the reactions use `yield_reactions` method.
        Note that the reactions attribute retrieves a copy of the reactions' container.
        To update the reactions container set new `reactions` or use `add` and `remove` methods.
        :return: a dictionary with the reactions of the model
        """
        return self._reactions.copy()

    @property
    def compartments(self) -> Dict[str, str]:
        """
        It returns a dictionary with the compartments of the model. The key is the compartment identifier a
        nd the value is the compartment name.
        To retrieve an iterator with the compartments use `yield_compartments` method.
        Note that the compartments attribute retrieves a copy of the compartments' container.
        To update the compartments container set new `compartments`.
        :return:
        """

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
        """
        It sets the compartments of the model. The key is the compartment identifier a
        nd the value is the compartment name.
        :param value: a dictionary with the compartments of the model
        :return:
        """
        if not value:
            value = {}

        self.__compartments.update(value)

    @genes.setter
    @recorder
    def genes(self, value: Dict[str, 'Gene']):
        """
        It sets the genes of the model. The key is the gene identifier and the value is the `Gene` object.
        :param value: a dictionary with the genes of the model
        :return:
        """
        if not value:
            value = {}

        self.remove(*self.yield_genes(), history=False)
        self.add(*value.values(), history=False)

    @metabolites.setter
    @recorder
    def metabolites(self, value: Dict[str, 'Metabolite']):
        """
        It sets the metabolites of the model. The key is the metabolite identifier and the value is the `Metabolite`
        object.
        :param value: a dictionary with the metabolites of the model
        :return:
        """
        if not value:
            value = {}

        self.remove(*self.yield_metabolites(), history=False)
        self.add(*value.values(), history=False)

    @objective.setter
    @recorder
    def objective(self, value: Dict['Reaction', Union[float, int]]):
        """
        It sets the objective functions of the model. The key is the `Reaction` object and the value is the respective
        coefficient.
        :param value: a dictionary with the objective functions of the model
        :return:
        """
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

        for simulator in self.simulators:
            simulator.set_objective(linear=linear_obj, minimize=False)

    @reactions.setter
    @recorder
    def reactions(self, value: Dict[str, 'Reaction']):
        """
        It sets the reactions of the model. The key is the reaction identifier and the value is the `Reaction` object.
        :param value: a dictionary with the reactions of the model
        :return:
        """
        if not value:
            value = {}

        self.remove(*self.yield_reactions(), history=False)
        self.add(*value.values(), history=False)

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------

    @property
    def external_compartment(self) -> Union[str, None]:
        """
        It returns the external compartment of the model. This compartment usually corresponds
        to the extracellular space.
        The external compartment is identified as the compartment having more external metabolites.
        External metabolites are often associated with boundary reactions.
        :return: the external compartment of the model
        """
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

    def _get_boundaries(self) -> Tuple[Dict[str, 'Reaction'], Dict[str, 'Reaction'], Dict[str, 'Reaction']]:
        """
        It returns the boundary reactions of the model.
        :return: a tuple with exchanges, sinks, demands reactions of the model
        """
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
        """
        It returns the demand reactions of the model.
        Demand reactions are reactions that consume a metabolite from its compartment.
        :return: a dictionary with the demand reactions of the model
        """
        _, _, demands = self._get_boundaries()
        return demands

    @property
    def exchanges(self) -> Dict[str, 'Reaction']:
        """
        It returns the exchange reactions of the model.
        Exchange reactions are reactions define the environmental conditions of a metabolic model.
        These reactions can provide or consume a metabolite from the extracellular space.
        :return: a dictionary with the exchange reactions of the model
        """
        exchanges, _, _ = self._get_boundaries()
        return exchanges

    @property
    def sinks(self) -> Dict[str, 'Reaction']:
        """
        It returns the sink reactions of the model.
        Sink reactions are reactions that either consume or produce a metabolite in its compartment.
        :return: a dictionary with the sink reactions of the model
        """
        _, sinks, _ = self._get_boundaries()
        return sinks

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------
    def yield_compartments(self) -> Generator[str, None, None]:
        """
        It yields the compartments of the model.
        :return: a generator with the compartments of the model
        """
        return generator(self.compartments)

    def yield_demands(self) -> Generator['Reaction', None, None]:
        """
        It yields the demand reactions of the model.
        :return: a generator with the demand reactions of the model
        """
        return generator(self.demands)

    def yield_exchanges(self) -> Generator['Reaction', None, None]:
        """
        It yields the exchange reactions of the model.
        :return: a generator with the exchange reactions of the model
        """
        return generator(self.exchanges)

    def yield_genes(self) -> Generator['Gene', None, None]:
        """
        It yields the genes of the model.
        :return: a generator with the genes of the model
        """
        return generator(self._genes)

    def yield_gprs(self) -> Generator['Expression', None, None]:
        """
        It yields the GPRs of the model.
        :return: a generator with the GPRs of the model
        """
        return (value.gpr for value in self._reactions.values())

    def yield_metabolites(self) -> Generator['Metabolite', None, None]:
        """
        It yields the metabolites of the model.
        :return: a generator with the metabolites of the model
        """
        return generator(self._metabolites)

    def yield_reactions(self) -> Generator['Reaction', None, None]:
        """
        It yields the reactions of the model.
        :return: a generator with the reactions of the model
        """
        return generator(self._reactions)

    def yield_sinks(self) -> Generator['Reaction', None, None]:
        """
        It yields the sink reactions of the model.
        :return: a generator with the sink reactions of the model
        """
        return generator(self.sinks)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------
    def get(self, identifier: Any, default=None) -> Union['Gene', 'Metabolite', 'Reaction']:
        """
        It returns the object associated with the identifier.
        In case the identifier is not found, it returns the default value.
        For metabolic models, the identifier can be a gene, a metabolite or a reaction.
        :param identifier: the identifier of the object
        :param default: the default value to return in case the identifier is not found
        :return: the object associated with the identifier
        """
        if identifier in self._metabolites:
            return self._metabolites[identifier]

        elif identifier in self._reactions:
            return self._reactions[identifier]

        elif identifier in self._genes:
            return self._genes[identifier]

        else:
            return super(MetabolicModel, self).get(identifier=identifier, default=default)

    def update(self,
               compartments: Dict[str, str] = None,
               objective: Dict['Reaction', Union[float, int]] = None,
               variables: Union[List[Union['Gene', 'Metabolite', 'Reaction']],
                                Tuple[Union['Gene', 'Metabolite', 'Reaction']],
                                Set[Union['Gene', 'Metabolite', 'Reaction']]] = None,
               **kwargs):
        """
        It updates the model with relevant information, namely the compartments, objective and variables.

        :param compartments: the compartments to be updated
        :param objective: the objective to be updated
        :param variables: the variables to be updated
        :param kwargs: additional arguments
        :return:
        """
        if compartments is not None:
            self.compartments = compartments

        if variables is not None:
            self.add(*variables)

        if objective is not None:
            self.objective = objective

        super(MetabolicModel, self).update(**kwargs)
