from typing import Any, TYPE_CHECKING, Union, Generator, Dict, List, Tuple, Set

from mewpy.model.model import Model
from mewpy.mew.lp import Notification
from mewpy.util.history import recorder
from mewpy.util.serialization import serialize
from mewpy.util.utilities import iterable, generator

if TYPE_CHECKING:
    from mewpy.mew.variables import Interaction, Regulator, Target, Metabolite, Reaction


# TODO: methods stubs
class RegulatoryModel(Model, model_type='regulatory', register=True, constructor=True, checker=True):
    """
    A mew regulatory model can represent a Transcriptional Regulatory Network (TRN), containing
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
        A mew regulatory model can represent a Transcriptional Regulatory Network (TRN), containing
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

        self.remove(list(self.yield_interactions()), 'interaction', history=False)
        self.add(list(value.values()), 'interaction', history=False)

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

        self.remove(list(self.yield_regulators()), 'regulator', history=False)
        self.add(list(value.values()), 'regulator', history=False)

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

        self.remove(list(self.yield_targets()), 'target', history=False)
        self.add(list(value.values()), 'target', history=False)

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

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------

    def yield_environmental_stimuli(self) -> Generator['Regulator', None, None]:
        """
        It returns an iterator with the environmental stimuli of the model.
        :return: a generator with the environmental stimuli of the model
        """
        return generator(self.environmental_stimuli)

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
            variables: Union[List[Union['Interaction', 'Regulator', 'Target']],
                             Tuple[Union['Interaction', 'Regulator', 'Target']],
                             Set[Union['Interaction', 'Regulator', 'Target']]],
            *types: str,
            comprehensive: bool = True,
            history=True):
        """
        It adds the variables to the model.

        This method accepts a single variable or a list of variables to be added to specific containers in the model.
        The containers to which the variables will be added are specified by the types.

        For instance, if a variable is simultaneously a targets and regulators,
        it will be added to the targets and regulators containers.

        If comprehensive is True, the variables and their related variables will be added to the model too.
        If history is True, the changes will be recorded in the history.

        This method notifies all simulators with the recent changes.

        :param variables: the variables to be added to the model
        :param types: the types of the variables setting which containers will store the variables.
        If none, the variables will be added to the containers according to their type.
        :param comprehensive: if True, the variables and their related variables will be added to the model too
        :param history: if True, the changes will be recorded in the history
        :return:
        """
        variables = iterable(variables)

        if not types:
            types = [var.types for var in variables]

        elif len(types) == len(variables):
            types = [{_type} if isinstance(_type, str) else set(_type) for _type in types]

        elif len(types) != len(variables):
            types = [set(types) for var in variables]

        interactions = []
        new_variables = []
        new_types = []

        for var_types, var in zip(types, variables):

            if 'target' in var_types:
                self._add_target(var)
                var_types.remove('target')

            elif 'regulator' in var_types:
                self._add_regulator(var)
                var_types.remove('regulator')

            elif 'interaction' in var_types:
                self._add_interaction(var, comprehensive=comprehensive)
                interactions.append(var)
                var_types.remove('interaction')

            if var_types:
                new_types.append(var_types)
                new_variables.append(var)

        if interactions:
            notification = Notification(content=interactions,
                                        content_type='interactions',
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

        super(RegulatoryModel, self).add(new_variables,
                                         *new_types,
                                         comprehensive=comprehensive,
                                         history=False)

    def remove(self,
               variables: Union[List[Union['Interaction', 'Regulator', 'Target']],
                                Tuple[Union['Interaction', 'Regulator', 'Target']],
                                Set[Union['Interaction', 'Regulator', 'Target']]],
               *types: str,
               remove_orphans: bool = False,
               history=True):
        """
        It removes the variables from the model.

        This method accepts a single variable or a list of variables to be removed from specific containers in the model.
        The containers from which the variables will be removed are specified by the types.

        For instance, if a variable is simultaneously a targets and regulators,
        it will be removed from the targets and regulators containers.

        If remove_orphans is True, the variables and their related variables will be removed from the model too.
        If history is True, the changes will be recorded in the history.

        This method notifies all simulators with the recent changes.

        :param variables: the variables to be removed from the model
        :param types: the types of the variables setting which containers will remove the variables.
        If none, the variables will be removed from the containers according to their type.
        :param remove_orphans: if True, the variables and their related variables will be removed from the model too
        :param history: if True, the changes will be recorded in the history
        :return:
        """
        variables = iterable(variables)

        if not types:
            types = [var.types for var in variables]

        elif len(types) == len(variables):
            types = [{_type} if isinstance(_type, str) else set(_type) for _type in types]

        elif len(types) != len(variables):
            types = [set(types) for var in variables]

        interactions = []
        new_variables = []
        new_types = []

        for var_types, var in zip(types, variables):

            if 'target' in var_types:
                self._remove_target(var)
                var_types.remove('target')

            elif 'regulator' in var_types:
                self._remove_regulator(var)
                var_types.remove('regulator')

            elif 'interaction' in var_types:
                self._remove_interaction(var, remove_orphans=remove_orphans)
                interactions.append(var)
                var_types.remove('interaction')

            if var_types:
                new_types.append(var_types)
                new_variables.append(var)

        if interactions:

            notification = Notification(content=interactions,
                                        content_type='interactions',
                                        action='remove')

            self.notify(notification)

            if remove_orphans:
                self._remove_regulatory_orphans(interactions)

        if history:
            self.history.queue_command(undo_func=self.add,
                                       undo_kwargs={'variables': variables,
                                                    'comprehensive': True,
                                                    'history': False},
                                       func=self.remove,
                                       kwargs={'variables': variables,
                                               'remove_orphans': remove_orphans,
                                               'history': history})

        super(RegulatoryModel, self).remove(new_variables,
                                            *new_types,
                                            remove_orphans=remove_orphans,
                                            history=False)

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
            self.add(variables=variables)

        super(RegulatoryModel, self).update(**kwargs)

    def _add_interaction(self, interaction, comprehensive=True):

        # adding a interaction is regularly a in-depth append method, as both regulators and target associated with
        # the interaction are also regularly added to the model. Although, this behaviour can be avoided passing
        # comprehensive=True

        if interaction.id not in self._interactions:

            if comprehensive:

                for regulator in interaction.yield_regulators():
                    self._add_regulator(regulator)

                if interaction.target is not None:
                    self._add_target(interaction.target)

            self._add_variable_to_container(interaction, self._interactions)

    def _add_regulator(self, regulator):

        if regulator.id not in self._regulators:
            self._add_variable_to_container(regulator, self._regulators)

    def _add_target(self, target):

        if target.id not in self._targets:
            self._add_variable_to_container(target, self._targets)

    def _remove_interaction(self, interaction, remove_orphans=False):

        if interaction.id in self._interactions:
            self._remove_variable_from_container(interaction, self._interactions)

    def _remove_regulator(self, regulator):

        if regulator.id in self._regulators:
            self._remove_variable_from_container(regulator, self._regulators)

    def _remove_target(self, target):

        if target.id in self._targets:
            self._remove_variable_from_container(target, self._targets)

    def _remove_regulatory_orphans(self, interactions):

        orphan_regs = self._get_orphans(to_remove=interactions,
                                        first_container='regulators',
                                        second_container='interactions')

        orphan_targets = self._get_orphans(to_remove=interactions,
                                           first_container='target',
                                           second_container='interaction')

        if orphan_regs:

            for regulator in orphan_regs:
                self._remove_regulator(regulator)

        if orphan_targets:

            for target in orphan_targets:
                self._remove_target(target)
