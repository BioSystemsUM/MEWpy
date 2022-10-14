from typing import Any, TYPE_CHECKING, Union, Generator, Dict, List, Tuple, Set

from .model import Model
from mewpy.util.history import recorder
from mewpy.germ.models.serialization import serialize
from mewpy.util.utilities import generator

if TYPE_CHECKING:
    from mewpy.germ.variables import Interaction, Regulator, Target, Metabolite, Reaction


class RegulatoryModel(Model, model_type='regulatory', register=True, constructor=True, checker=True):
    """
    A germ regulatory model can represent a Transcriptional Regulatory Network (TRN), containing
    interactions between regulators and targets.
    However, this regulatory model is not a standard directed TRN graph, but a set of
    interactions, each one containing a single target and a set of regulators.
    Each interaction contains a set of regulatory events that determine the state or coefficient of the target gene.
    A regulatory event consists of a boolean algebra expression and the corresponding coefficient.
    If the boolean expression is evaluated to True, the resulting coefficient is applied to the target gene.
    Otherwise, the coefficient is ignored.

    The model also contains a set of environmental stimuli, which are associated with
    interactions. The environmental stimuli can be used to define the initial state of the model.

    The regulatory model can be loaded with compartments, although these can be inferred from the available
    regulators.

    The regulatory model, as with other models, provides a clean interface for manipulation with the add, remove and
    update methods. One can perform the following operations:
        - Add a new interaction, regulator or target
        - Remove an interaction, regulator or target
        - Update the compartments of the model
    """
    def __init__(self,
                 identifier: Any,
                 compartments: Dict[str, str] = None,
                 interactions: Dict[str, 'Interaction'] = None,
                 regulators: Dict[str, 'Regulator'] = None,
                 targets: Dict[str, 'Target'] = None,
                 **kwargs):

        """
        A germ regulatory model can represent a Transcriptional Regulatory Network (TRN), containing
        interactions between regulators and targets.
        However, this regulatory model is not a standard directed TRN graph, but a set of
        interactions, each one containing a single target and a set of regulators.
        Each interaction contains a set of regulatory events that determine the state or coefficient of the target gene.
        A regulatory event consists of a boolean algebra expression and the corresponding coefficient.
        If the boolean expression is evaluated to True, the resulting coefficient is applied to the target gene.
        Otherwise, the coefficient is ignored.

        The model also contains a set of environmental stimuli, which are associated with
        interactions. The environmental stimuli can be used to define the initial state of the model.

        The regulatory model can be loaded with compartments, although these can be inferred from the available
        regulators.

        The regulatory model, as with other models, provides a clean interface for manipulation with the add, remove and
        update methods. One can perform the following operations:
            - Add a new interaction, regulator or target
            - Remove an interaction, regulator or target
            - Update the compartments of the model

        :param identifier: identifier, e.g. lac_operon_regulation
        :param compartments: a dictionary with additional compartments not encoded in the regulators
        :param interactions: a dictionary with Interaction objects. See variables.Interaction for more info
        :param regulators: a dictionary with Regulator objects. See variables.Regulator for more info
        :param targets: a dictionary with Target objects. See variables.Target for more info
        """
        # compartments attribute can be shared across the children, thus name mangling
        self.__compartments = {}
        self._interactions = {}
        self._regulators = {}
        self._targets = {}

        super().__init__(identifier,
                         **kwargs)

        # the setters will handle adding and removing variables to the correct containers
        self.compartments = compartments
        self.interactions = interactions
        self.regulators = regulators
        self.targets = targets

    # -----------------------------------------------------------------------------
    # Model type manager
    # -----------------------------------------------------------------------------
    @serialize('types', None)
    @property
    def types(self):
        """
        Returns the model types
        :return: a list with the model types
        """
        _types = {RegulatoryModel.model_type}

        _types.update(super(RegulatoryModel, self).types)

        return _types

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------

    @serialize('interactions', 'interactions', '_interactions')
    @property
    def interactions(self) -> Dict[str, 'Interaction']:
        """
        It returns a dictionary with the interactions of the model. The keys are the identifiers of the interactions
        and the values are the `Interaction` objects. To retrieve an iterator with the interactions, use the
        `yield_interactions` method.
        Note that the interactions attribute retrieves a copy of the interactions' container.
        To update the interactions container set new `interactions` or use `add` and `remove` methods.
        :return: a dictionary with the interactions of the model
        """
        return self._interactions.copy()

    @serialize('regulators', 'regulators', '_regulators')
    @property
    def regulators(self) -> Dict[str, 'Regulator']:
        """
        It returns a dictionary with the regulators of the model. The keys are the identifiers of the regulators
        and the values are the `Regulator` objects. To retrieve an iterator with the regulators, use the
        `yield_regulators` method.
        Note that the regulators attribute retrieves a copy of the regulators' container.
        To update the regulators container set new `regulators` or use `add` and `remove` methods.
        :return: a dictionary with the regulators of the model
        """
        return self._regulators.copy()

    @serialize('targets', 'targets', '_targets')
    @property
    def targets(self) -> Dict[str, 'Target']:
        """
        It returns a dictionary with the targets of the model. The keys are the identifiers of the targets
        and the values are the `Target` objects. To retrieve an iterator with the targets, use the
        `yield_targets` method.
        Note that the targets attribute retrieves a copy of the targets' container.
        To update the targets container set new `targets` or use `add` and `remove` methods.
        :return: a dictionary with the targets of the model
        """
        return self._targets.copy()

    @property
    def compartments(self) -> Dict[str, str]:
        """
        It returns a dictionary with the compartments of the model. The keys are the identifiers of the compartments
        and the values are the compartment names. To retrieve an iterator with the compartments, use the
        `yield_compartments` method.
        Note that the compartments attribute retrieves a copy of the compartments' container.
        To update the compartments container set new `compartments`.
        :return: a dictionary with the compartments of the model
        """
        compartments = {regulator.compartment: self.__compartments.get(regulator.compartment, '')
                        for regulator in self.yield_regulators()
                        if regulator.is_metabolite() and regulator.compartment is not None}

        compartments.update(self.__compartments)

        compartments.update(super(RegulatoryModel, self).compartments)

        return compartments

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------
    @compartments.setter
    @recorder
    def compartments(self, value: Dict[str, str]):
        """
        It sets the compartments of the model. The keys are the identifiers of the compartments
        and the values are the compartment names.
        :param value: a dictionary with the compartments of the model
        :return:
        """
        if not value:
            value = {}

        self.__compartments.update(value)

    @interactions.setter
    @recorder
    def interactions(self, value: Dict[str, 'Interaction']):
        """
        It sets the interactions of the model. The keys are the identifiers of the interactions
        and the values are the `Interaction` objects.
        :param value: a dictionary with the interactions of the model
        :return:
        """

        if not value:
            value = {}

        self.remove(*self.yield_interactions(), history=False)
        self.add(*value.values(), history=False)

    @regulators.setter
    @recorder
    def regulators(self, value: Dict[str, 'Regulator']):
        """
        It sets the regulators of the model. The keys are the identifiers of the regulators
        and the values are the `Regulator` objects.
        :param value: a dictionary with the regulators of the model
        :return:
        """
        if not value:
            value = {}

        self.remove(*self.yield_regulators(), history=False)
        self.add(*value.values(), history=False)

    @targets.setter
    @recorder
    def targets(self, value: Dict[str, 'Target']):
        """
        It sets the targets of the model. The keys are the identifiers of the targets
        and the values are the `Target` objects.
        :param value: a dictionary with the targets of the model
        :return:
        """
        if not value:
            value = {}

        self.remove(*self.yield_targets(), history=False)
        self.add(*value.values(), history=False)

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------
    @property
    def environmental_stimuli(self) -> Dict[str, 'Regulator']:
        """
        It returns a dictionary with the environmental stimuli of the model. The keys are the identifiers of the
        environmental stimuli and the values are the `Regulator` objects. To retrieve an iterator with the environmental
        stimuli, use the `yield_environmental_stimuli` method.
        :return: a dictionary with the environmental stimuli of the model
        """
        return {reg_id: regulator for reg_id, regulator in self.regulators.items()
                if regulator.environmental_stimulus}

    @property
    def regulatory_reactions(self) -> Dict[str, 'Regulator']:
        """
        It returns a dictionary with the regulatory reactions of the model. The keys are the identifiers of the
        regulatory reactions and the values are the `Regulator` objects. To retrieve an iterator with the regulatory
        reactions, use the `yield_regulatory_reactions` method.
        :return: a dictionary with the regulatory reactions of the model
        """
        return {reg_id: regulator for reg_id, regulator in self.regulators.items()
                if regulator.is_reaction()}

    @property
    def regulatory_metabolites(self) -> Dict[str, 'Regulator']:
        """
        It returns a dictionary with the regulatory metabolites of the model. The keys are the identifiers of the
        regulatory metabolites and the values are the `Regulator` objects. To retrieve an iterator with the regulatory
        metabolites, use the `yield_regulatory_metabolites` method.
        :return: a dictionary with the regulatory metabolites of the model
        """
        return {reg_id: regulator for reg_id, regulator in self.regulators.items()
                if regulator.is_metabolite()}

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------
    def yield_environmental_stimuli(self) -> Generator['Regulator', None, None]:
        """
        It returns an iterator with the environmental stimuli of the model.
        :return: a generator with the environmental stimuli of the model
        """
        return generator(self.environmental_stimuli)

    def yield_regulatory_reactions(self) -> Generator['Regulator', None, None]:
        """
        It returns an iterator with the regulatory reactions of the model.
        :return: a generator with the regulatory reactions of the model
        """
        return generator(self.regulatory_reactions)

    def yield_regulatory_metabolites(self) -> Generator['Regulator', None, None]:
        """
        It returns an iterator with the regulatory metabolites of the model.
        :return: a generator with the regulatory metabolites of the model
        """
        return generator(self.regulatory_metabolites)

    def yield_interactions(self) -> Generator['Interaction', None, None]:
        """
        It returns an iterator with the interactions of the model.
        :return: a generator with the interactions of the model
        """
        return generator(self._interactions)

    def yield_regulators(self) -> Generator[Union['Regulator', 'Metabolite', 'Reaction'], None, None]:
        """
        It returns an iterator with the regulators of the model.
        :return: a generator with the regulators of the model
        """
        return generator(self._regulators)

    def yield_targets(self) -> Generator['Target', None, None]:
        """
        It returns an iterator with the targets of the model.
        :return: a generator with the targets of the model
        """
        return generator(self._targets)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------
    def get(self, identifier: Any, default=None) -> Union['Interaction', 'Regulator', 'Target']:
        """
        It returns the object with the given identifier. If the object is not found, it returns the default value.
        For regulatory models, the identifier can be a gene, an interaction or a regulator.
        :param identifier: the identifier of the object
        :param default: the default value to return if the object is not found
        :return: the object associated the given identifier
        """
        if identifier in self._targets:
            return self._targets[identifier]

        elif identifier in self._regulators:
            return self._regulators[identifier]

        elif identifier in self._interactions:
            return self._interactions[identifier]

        else:
            return super(RegulatoryModel, self).get(identifier=identifier, default=default)

    def add(self,
            *variables: Union['Interaction', 'Regulator', 'Target'],
            comprehensive: bool = True,
            history: bool = True):
        """
        It adds the given variables to the model.
        This method accepts a single variable or a list of variables to be added to specific containers in the model.
        The containers to which the variables will be added are specified by the types.

        For instance, if a variable is simultaneously a metabolite and regulator,
        it will be added to the metabolites and regulators containers.

        If comprehensive is True, the variables and their related variables will be added to the model too.
        If history is True, the changes will be recorded in the history.

        This method notifies all simulators with the recent changes.

        :param variables: the variables to be added to the model
        :param comprehensive: if True, the variables and their related variables will be added to the model too
        :param history: if True, the changes will be recorded in the history
        :return:
        """
        if self.is_a('regulatory'):

            for variable in variables:

                if 'target' in variable.types:
                    self._add_variable_to_container(variable, '_targets')

                if 'regulator' in variable.types:
                    self._add_variable_to_container(variable, '_regulators')

                if 'interaction' in variable.types:
                    if comprehensive:

                        if variable.target is not None:
                            self._add_variable_to_container(variable.target, '_targets')

                        for regulator in variable.yield_regulators():
                            self._add_variable_to_container(regulator, '_regulators')

                    self._add_variable_to_container(variable, '_interactions')

        return super(RegulatoryModel, self).add(*variables, comprehensive=comprehensive, history=history)

    def remove(self,
               *variables: Union['Interaction', 'Regulator', 'Target'],
               remove_orphans: bool = False,
               history: bool = True):
        """
        It removes the given variables from the model.
        This method accepts a single variable or a list of variables to be removed from specific containers
        in the model.
        The containers from which the variables will be removed are specified by the types.

        For instance, if a variable is simultaneously a metabolite and regulator,
        it will be removed from the metabolites and regulators containers.

        If remove_orphans is True, the variables and their related variables will be removed from the model too.
        If history is True, the changes will be recorded in the history.

        This method notifies all simulators with the recent changes.

        :param variables: the variables to be removed from the model
        :param remove_orphans: if True, the variables and their related variables will be removed from the model too
        :param history: if True, the changes will be recorded in the history
        :return:
        """
        if self.is_a('regulatory'):

            interactions = set()

            for variable in variables:

                if 'target' in variable.types:
                    self._remove_variable_from_container(variable, '_targets')

                if 'regulator' in variable.types:
                    self._remove_variable_from_container(variable, '_regulators')

                if 'interaction' in variable.types:
                    self._remove_variable_from_container(variable, '_interactions')
                    interactions.add(variable)

            if remove_orphans:
                for interaction in interactions:
                    if interaction.target:
                        self._remove_variable_from_container(interaction.target, '_targets')

                orphan_regulators = self._get_orphans(to_remove=interactions,
                                                      first_container='regulators',
                                                      second_container='interactions')

                for regulator in orphan_regulators:
                    self._remove_variable_from_container(regulator, '_regulators')

        return super(RegulatoryModel, self).remove(*variables, remove_orphans=remove_orphans, history=history)

    def update(self,
               compartments: Dict[str, str] = None,
               variables: Union[List[Union['Interaction', 'Regulator', 'Target']],
                                Tuple[Union['Interaction', 'Regulator', 'Target']],
                                Set[Union['Interaction', 'Regulator', 'Target']]] = None,
               **kwargs):
        """
        It updates the model with relevant information, namely the compartments and the variables.

        :param compartments: the compartments to be added to the model
        :param variables: the variables to be added to the model
        :param kwargs: additional arguments
        :return:
        """
        if compartments is not None:
            self.compartments = compartments

        if variables is not None:
            self.add(*variables)

        super(RegulatoryModel, self).update(**kwargs)
