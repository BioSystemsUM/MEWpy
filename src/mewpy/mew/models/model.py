from typing import Any, Union, Type, TYPE_CHECKING, List, Set, Dict, Iterable, Tuple

from mewpy.util.history import HistoryManager, recorder
from mewpy.mew.models.serialization import serialize, Serializer

# Preventing circular dependencies that only happen due to type checking
if TYPE_CHECKING:
    from mewpy.mew.models import MetabolicModel, RegulatoryModel
    from mewpy.mew.variables import Gene, Interaction, Metabolite, Reaction, Regulator, Target, Variable
    from mewpy.mew.lp import LinearProblem


class MetaModel(type):
    """
    The MetaModel class is used to dynamically create Model classes.
    Models instances are created from static classes representing the different information for a given model.

    However, integrated Metabolic-Regulatory models must combine different types of variables, such as genes, reactions,
    metabolites, but also regulators, targets, interactions. This is not possible with static classes, as they miss
    the containers to store these variables.
    Thus, a Model factory is used to dynamically create Model classes
    based on the multiple types of a single integrated model.

    The Model factory is the interface to create these multi-type models.
    It also manages all models by implementing base attributes and methods.

    The MetaModel class creates a dynamic Model class as follows:
        1. It collects all static base classes to be included into the dynamic Model class
        2. It collects the model type from the base classes or from the model_type attribute
        3. It sets the name of the dynamic Model class using the name of the base classes
        4. It updates the model_type attribute with all types of the base classes
        5. It creates the dynamic Model class having all attributes and methods of the base classes
        6. It adds polymorphic constructors to the dynamic Model class based on the types of the base classes
        7. It adds type checkers to the dynamic Model class based on the types of the base classes
    """
    factories = {}

    def __new__(mcs, name, bases, attrs, **kwargs):
        # if it is the model factory, only registration is done
        factory = kwargs.get('factory', False)

        if factory:
            cls = super(MetaModel, mcs).__new__(mcs, name, bases, attrs)

            MetaModel.factories[name] = cls

            return cls

        # Dynamic typing being used. In this case, a proper name and model type must be provided
        dynamic = kwargs.get('dynamic', False)
        if dynamic:
            names = [base.model_type for base in bases]

            name = ''.join([name.title() for name in names])
            name += 'Model'

            kwargs['model_type'] = '-'.join(names)

        # The model type is always added to the subclasses. If it is not given upon subclass creation,
        # the subclass name is to be used
        model_type = kwargs.get('model_type', name.lower())
        attrs['model_type'] = model_type

        return super(MetaModel, mcs).__new__(mcs, name, bases, attrs)

    def __init__(cls, name, bases, attrs, **kwargs):
        super().__init__(name, bases, attrs)

        factory = kwargs.get('factory', False)
        if factory:
            # collection of all containers that must be serialized for a particular class or subclass
            containers = cls.get_serializable_containers(attrs)

            attrs['_containers_registry']['model'] = containers

            # Skip further building of the Model factory
            return

        # Dynamic typing being used. In this case, all children have already been constructed, so everything can be
        # skipped
        dynamic = kwargs.get('dynamic', False)
        if dynamic:
            return

        # Several attributes and methods must be added automatically to the Model factory children based on the
        # model type.
        model_type = attrs['model_type']

        for base in bases:

            factory = MetaModel.factories.get(base.__name__)

            if factory and hasattr(cls, 'model_type'):

                # if some class inherits from the Model factory, it must be registered in the factory for type
                # checking
                if cls.model_type not in factory.get_registry():
                    raise TypeError(f'{cls.model_type} does not inherit from Model')

                # If set otherwise upon class creation, all subclasses of the Model factory will contribute to
                # their parent factory with an alternative initializer. The polymorphic constructor, e.g. from_{
                # model_type}, is added to the Model factory
                constructor = kwargs.get('constructor', True)

                if constructor:
                    # noinspection PyProtectedMember
                    model_type_initializer = Model._from(cls)
                    setattr(factory, f'from_{model_type}', model_type_initializer)

                # If set otherwise upon class creation, all subclasses of the Model factory will contribute to
                # their parent factory with a declared type checker. The type checker, e.g. is_{model_type},
                # is added to the Model factory
                checker = kwargs.get('checker', True)

                if checker:
                    # noinspection PyProtectedMember
                    model_type_checker = Model._is(model_type)
                    setattr(factory, f'is_{model_type}', model_type_checker)

                # collection of all containers that must be serialized for a particular class or subclass
                containers = cls.get_serializable_containers(attrs)

                Model.register_containers(containers, cls)

    @staticmethod
    def get_serializable_containers(attrs: Dict[str, Any]) -> Dict[str, Tuple[str, str, str]]:
        """
        Collects all containers that must be serialized for a particular class or subclass

         All containers have property descriptors marked with the @serialize decorator. This decorator marks attributes
        for serialization, deserialization and pickle serialization/deserialization. The decorator also provides
        information about the attribute name.

        If serialization, deserialization or pickle serialization/deserialization is not desired, the attribute must
        be marked with the @serialize decorator with the corresponding argument set to None.
        :param attrs: the attributes of the class or subclass
        :return: a dictionary with the container name as key and a tuple with the attribute name for serialization,
        deserialization and pickle serialization/deserialization
        """
        containers = {}

        for name, method in attrs.items():

            if hasattr(method, 'fget'):

                if hasattr(method.fget, 'serialize') and hasattr(method.fget, 'deserialize') and hasattr(method.fget,
                                                                                                         'pickle'):
                    containers[name] = (method.fget.serialize, method.fget.deserialize, method.fget.pickle)

        return containers


class Model(Serializer, metaclass=MetaModel, factory=True):
    """
    The Model class is the base class and factory for all models.
    It is the interface to create multi-type models providing several base attributes and methods.

    It implements the basic attributes:
        - types: the types of the model
        - id: the id of the model
        - name: the name of the model
        - simulators: the simulation methods attached to the model. Linear problems such as FBA, RFBA, etc.
        - compartments: the compartments of the model
        - contexts: the contexts associated with the model.
        - history: the history of the model.

    It implements the basic methods:
        - factory: the factory for dynamic model types
        - from_types: creates a model from a list of types
        - copy: creates a copy of the model
        - deepcopy: creates a deep copy of the model
        - get: returns a variable in the model
        - add: adds a variable to the model
        - remove: removes a variable from the model
        - update: updates the model with a new information
        - attach: attaches a simulation method to the model
        - detach: detaches a simulation method from the model
        - notify: notifies all simulation methods attached to the model regarding a change in the model
        - undo: undoes the last change of the model
        - redo: redoes the last change of the model
        - reset: resets the model to its initial state
        - restore: restores the model to a previous state
        - from polymorphic constructors: creates a model from a specific type
        - is polymorphic type checkers: checks if a model is of a specific type

    It implements the basic serialization methods:
        - to_dict: serializes the model to a dictionary
        - from_dict: deserializes the model from a dictionary

    """
    # -----------------------------------------------------------------------------
    # Factory management
    # -----------------------------------------------------------------------------
    _registry = {}
    _containers_registry = {'model': {}}

    def __init_subclass__(cls, **kwargs):
        """
        This method is called when a subclass of Model is created. It is used to register the subclass in the
        Model factory.

        Internal use only.
        :param kwargs: the keyword arguments
        """
        super(Model, cls).__init_subclass__(**kwargs)

        # the child type
        model_type = getattr(cls, 'model_type', cls.__name__.lower())

        cls.register_type(model_type, cls)

    @staticmethod
    def get_registry() -> Dict[str, Type['Model']]:
        """
        Returns the registry of the Model factory.

        Internal use only.
        :return: the registry of the Model factory
        """
        return Model._registry.copy()

    @staticmethod
    def get_containers_registry() -> Dict[str, Dict[str, Tuple[str, str, str]]]:
        """
        Returns the containers registry of the Model factory.

        Internal use only.
        :return: the containers registry of the Model factory
        """
        return Model._containers_registry.copy()

    @staticmethod
    def register_type(model_type: str, child: Type['Model']):
        """
        Registers a type in the Model factory.

        Internal use only.
        :param model_type: the type to register
        :param child: the child class
        :return:
        """
        Model._registry[model_type] = child

    # -----------------------------------------------------------------------------
    # Serialization
    # -----------------------------------------------------------------------------
    @staticmethod
    def register_containers(containers, child):
        """
        Registers the containers of a child class in the containers registry of the Model factory. This is useful for
        serialization and deserialization.

        Internal use only.
        :param containers: the containers to register
        :param child: the child class
        :return:
        """
        if hasattr(child, 'model_type'):

            Model._containers_registry[child.model_type] = containers

        elif child is Model:
            Model._containers_registry['model'] = containers

    @property
    def containers(self) -> Dict[str, Tuple[str, str, str]]:
        """
        Returns the containers of the model. This is useful for serialization and deserialization.

        Internal use only.
        :return: the containers of the model as a dictionary with the container name as key and a tuple with the
        attribute name for serialization, deserialization and pickle serialization/deserialization
        """
        class_containers = self.get_containers_registry()

        containers = class_containers['model'].copy()

        for model_type in self.types:

            for name, (serialize_name, deserialize_name, pickle_name) in class_containers[model_type].items():
                containers[name] = (serialize_name, deserialize_name, pickle_name)

        return containers

    # -----------------------------------------------------------------------------
    # Factory polymorphic constructor
    # -----------------------------------------------------------------------------
    @classmethod
    def factory(cls, *args: str) -> Union[Type['Model'],
                                          Type['MetabolicModel'],
                                          Type['RegulatoryModel']]:
        """
        It creates a dynamic Model class from a list of types. The types must be registered in the Model factory.

        Example:
            >>> Model.factory('metabolic', 'regulatory')
            <class 'mewpy.model.MetabolicRegulatoryModel'>

        :param args: the types of the model
        :return: the model type
        """
        if not args:
            args = ()

        registry = cls.get_registry()

        types = tuple([registry[name] for name in args])

        if len(types) == 1:
            return types[0]

        _Model = MetaModel('Model', types, {}, dynamic=True)

        return _Model

    # -----------------------------------------------------------------------------
    # Factory polymorphic initializer
    # -----------------------------------------------------------------------------
    @classmethod
    def from_types(cls, types: Iterable[str], **kwargs) -> Union['Model',
                                                                 'MetabolicModel',
                                                                 'RegulatoryModel']:
        """
        It creates a model instance from a list of types and a dictionary of containers and attributes.
        The types must be registered in the Model factory.

        :param types: the types of the model
        :param kwargs: the containers and attributes to initilize the model
        :return: the model instance
        """
        ModelType = cls.factory(*types)
        return ModelType(**kwargs)

    # -----------------------------------------------------------------------------
    # Factory helper methods
    # -----------------------------------------------------------------------------
    @staticmethod
    def _from(cls):
        """
        Method to be added to the subclasses of Model to create a model from a specific type.

        Internal use only.
        :param cls: the class of the model
        :return: the model to be added to the subclass of Model
        """

        def from_(identifier, **kwargs):
            return cls(identifier, **kwargs)

        return from_

    @staticmethod
    def _is(model_type):
        """
        Method to be added to the subclasses of Model to check if a model is of a specific type.

        Internal use only.
        :param model_type: the type to check
        :return: the method to be added to the subclass of Model
        """

        def is_(self):
            return self.is_a(model_type)

        return is_

    # -----------------------------------------------------------------------------
    # Base initializer
    # -----------------------------------------------------------------------------
    def __init__(self,
                 identifier: Any,
                 name: str = None):

        """
        The model is the base class for all models, such as metabolic model or regulatory model.
        See also model.MetabolicModel and model.RegulatoryModel for concrete implementations of a model type.

        The model object is the factory for several model types.
        This model object provides all batteries to create new model types and to manage the model types' inheritance.

        The factory type assists with systematic tasks such as attributes and containers serialization,
        history management, simulator observation, among others.

        For that, attributes, containers and inheritance are registered in the factory. The MetaModel manages
        both the factory and the model types derived from the factory.

        :param identifier: identifier,  e.g. iMC1010
        :param name: the name of the model
        """
        self._check_inheritance()

        if not identifier:
            identifier = ''

        if not name:
            name = identifier

        self._id = ''
        self._name = ''
        self._simulators = []
        self._types = set()

        # History: reversible changes to the model
        self._contexts = []
        self._history = HistoryManager()

        self._id = identifier
        self.name = name

        # models share different compartments. Although it is not the most subtle solution, name mangling is used to
        # store the specific compartments of each child. An alternative would be to check for shared attributes
        # between all children and perform this solution for them
        self.__compartments = {}

    def _check_inheritance(self):
        """
        It checks if the class is a subclass of Model and if it is not, it raises an error.
        :return:
        """
        registry = self.get_registry()

        for model_type in self.types:
            if model_type not in registry:
                raise ValueError(f'{model_type} is not registered as subclass of {self.__class__.__name__}')

    # -----------------------------------------------------------------------------
    # Built-in
    # -----------------------------------------------------------------------------

    def __str__(self):
        return f'Model {self.id} - {self.name}'

    def __repr__(self):
        return self.__str__()

    # noinspection PyUnresolvedReferences
    def _repr_html_(self):
        """
        It returns a html representation of the gene.
        """

        objective = getattr(self, 'objective', None)
        if objective:
            objective = next(iter(objective)).id
        else:
            objective = None

        if self.is_metabolic() and self.is_regulatory():
            return f"""
            <table>
                <tr>
                    <th><b>Model</b></th>
                    <td>{self.id}</td>
                </tr>
                <tr>
                    <th>Name</th>
                    <td>{self.name}</td>
                </tr>
                <tr>
                    <th>Types</th>
                    <td>{', '.join(self.types)}</td>
                </tr>
                <tr>
                    <th>Compartments</th>
                    <td>{', '.join(self.compartments)}</td>
                </tr>
                <tr>
                    <th>Reactions</th>
                    <td>{len(self.reactions)}</td>
                </tr>
                <tr>
                    <th>Metabolites</th>
                    <td>{len(self.metabolites)}</td>
                </tr>
                <tr>
                    <th>Genes</th>
                    <td>{len(self.genes)}</td>
                </tr>
                <tr>
                    <th>Exchanges</th>
                    <td>{len(self.exchanges)}</td>
                </tr>
                <tr>
                    <th>Demands</th>
                    <td>{len(self.demands)}</td>
                </tr>
                <tr>
                    <th>Sinks</th>
                    <td>{len(self.sinks)}</td>
                </tr>
                <tr>
                    <th>Objective</th>
                    <td>{objective}</td>
                </tr>
                <tr>
                    <th>Regulatory interactions</th>
                    <td>{len(self.interactions)}</td>
                </tr>
                <tr>
                    <th>Targets</th>
                    <td>{len(self.targets)}</td>
                </tr>
                <tr>
                    <th>Regulators</th>
                    <td>{len(self.regulators)}</td>
                </tr>
                <tr>
                    <th>Regulatory reactions</th>
                    <td>{len(self.regulatory_reactions)}</td>
                </tr>
                <tr>
                    <th>Regulatory metabolites</th>
                    <td>{len(self.regulatory_metabolites)}</td>
                </tr>
                <tr>
                    <th>Environmental stimuli</th>
                    <td>{len(self.environmental_stimuli)}</td>
                </tr>
            </table>
            """
        elif self.is_metabolic():
            return f"""
            <table>
                <tr>
                    <th><b>Model</b></th>
                    <td>{self.id}</td>
                </tr>
                <tr>
                    <th>Name</th>
                    <td>{self.name}</td>
                </tr>
                <tr>
                    <th>Types</th>
                    <td>{', '.join(self.types)}</td>
                </tr>
                <tr>
                    <th>Compartments</th>
                    <td>{', '.join(self.compartments)}</td>
                </tr>
                <tr>
                    <th>Reactions</th>
                    <td>{len(self.reactions)}</td>
                </tr>
                <tr>
                    <th>Metabolites</th>
                    <td>{len(self.metabolites)}</td>
                </tr>
                <tr>
                    <th>Genes</th>
                    <td>{len(self.genes)}</td>
                </tr>
                <tr>
                    <th>Exchanges</th>
                    <td>{len(self.exchanges)}</td>
                </tr>
                <tr>
                    <th>Demands</th>
                    <td>{len(self.demands)}</td>
                </tr>
                <tr>
                    <th>Sinks</th>
                    <td>{len(self.sinks)}</td>
                </tr>
                <tr>
                    <th>Objective</th>
                    <td>{objective}</td>
                </tr>
            </table>
            """
        elif self.is_regulatory():
            return f"""
            <table>
                <tr>
                    <th><b>Model</b></th>
                    <td>{self.id}</td>
                </tr>
                <tr>
                    <th>Name</th>
                    <td>{self.name}</td>
                </tr>
                <tr>
                    <th>Types</th>
                    <td>{', '.join(self.types)}</td>
                </tr>
                <tr>
                    <th>Compartments</th>
                    <td>{', '.join(self.compartments)}</td>
                </tr>
                <tr>
                    <th>Regulatory interactions</th>
                    <td>{len(self.interactions)}</td>
                </tr>
                <tr>
                    <th>Targets</th>
                    <td>{len(self.targets)}</td>
                </tr>
                <tr>
                    <th>Regulators</th>
                    <td>{len(self.regulators)}</td>
                </tr>
                <tr>
                    <th>Regulatory reactions</th>
                    <td>{len(self.regulatory_reactions)}</td>
                </tr>
                <tr>
                    <th>Regulatory metabolites</th>
                    <td>{len(self.regulatory_metabolites)}</td>
                </tr>
                <tr>
                    <th>Environmental stimuli</th>
                    <td>{len(self.environmental_stimuli)}</td>
                </tr>
            </table>
            """
        return f"""
                <table>
                    <tr>
                        <th><b>Model</b></th>
                        <td>{self.id}</td>
                    </tr>
                    <tr>
                        <th>Name</th>
                        <td>{self.name}</td>
                    </tr>
                    <tr>
                        <th>Types</th>
                        <td>{', '.join(self.types)}</td>
                    </tr>
                </table>
                """

    # -----------------------------------------------------------------------------
    # Model type manager
    # -----------------------------------------------------------------------------
    @serialize('types', None, None)
    @property
    def types(self) -> Set[str]:
        """
        It returns the types of the model.
        :return: the types of the model as a set of strings
        """
        return set()

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------
    @serialize('id', None, '_id')
    @property
    def id(self) -> Any:
        """
        It returns the identifier of the model.
        :return: the identifier of the model
        """
        return self._id

    @serialize('name', 'name', '_name')
    @property
    def name(self) -> str:
        """
        It returns the name of the model.
        :return: the name of the model
        """
        return self._name

    @property
    def simulators(self) -> List['LinearProblem']:
        """
        It returns the list of simulation methods associated with the model.
        :return: the list of simulation methods
        """
        return self._simulators

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------
    @name.setter
    @recorder
    def name(self, value: str):
        """
        It sets the name of the model.
        :param value: the name of the model
        :return:
        """
        if not value:
            value = ''

        self._name = value

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------

    def get(self, identifier: Any, default=None) -> Union['Gene',
                                                          'Interaction',
                                                          'Metabolite',
                                                          'Reaction',
                                                          'Regulator',
                                                          'Target']:
        """
        It returns the variable with the given identifier.
        :param identifier: the identifier of the variable
        :param default: the default value to return if the variable is not found
        :return: the variable with the given identifier
        """
        return default

    def add(self,
            *variables: Union['Gene',
                              'Interaction',
                              'Metabolite',
                              'Reaction',
                              'Regulator',
                              'Target'],
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
        if self.is_a('metabolic'):

            for variable in variables:

                if 'gene' in variable.types:
                    self._add_variable_to_container(variable, '_genes')

                if 'metabolite' in variable.types:
                    self._add_variable_to_container(variable, '_metabolites')

                if 'reaction' in variable.types:
                    if comprehensive:

                        for metabolite in variable.yield_metabolites():
                            self._add_variable_to_container(metabolite, '_metabolites')

                        for gene in variable.yield_genes():
                            self._add_variable_to_container(gene, '_genes')

                    self._add_variable_to_container(variable, '_reactions')

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

        if history:
            self.history.queue_command(undo_func=self.remove,
                                       undo_args=variables,
                                       undo_kwargs={'remove_orphans': True,
                                                    'history': False},
                                       func=self.add,
                                       args=variables,
                                       kwargs={'comprehensive': comprehensive,
                                               'history': history})

        self.notify()

    def remove(self,
               *variables: Union['Gene',
                                 'Interaction',
                                 'Metabolite',
                                 'Reaction',
                                 'Regulator',
                                 'Target',
                                 'Variable'],
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
        if self.is_a('metabolic'):

            reactions = set()

            for variable in variables:

                if 'gene' in variable.types:
                    self._remove_variable_from_container(variable, '_genes')

                if 'metabolite' in variable.types:
                    self._remove_variable_from_container(variable, '_metabolites')

                if 'reaction' in variable.types:
                    self._remove_variable_from_container(variable, '_reactions')
                    reactions.add(variable)

            if remove_orphans:
                orphan_metabolites = self._get_orphans(to_remove=reactions,
                                                       first_container='metabolites',
                                                       second_container='reactions')

                for metabolite in orphan_metabolites:
                    self._remove_variable_from_container(metabolite, '_metabolites')

                orphan_genes = self._get_orphans(to_remove=reactions,
                                                 first_container='genes',
                                                 second_container='reactions')

                for gene in orphan_genes:
                    self._remove_variable_from_container(gene, '_genes')

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

        if history:
            self.history.queue_command(undo_func=self.add,
                                       undo_args=variables,
                                       undo_kwargs={'comprehensive': True,
                                                    'history': False},
                                       func=self.remove,
                                       args=variables,
                                       kwargs={'remove_orphans': remove_orphans,
                                               'history': history})

        self.notify()

    def update(self, name: str = None):
        """
        It updates the model with relevant information.
        :param name: the name of the model
        :return:
        """
        if name:
            self.name = name

    def _add_variable_to_container(self, variable, container):
        """
        It adds the given variable to the given container.

        Helper method to be used internally.
        :param variable: the variable to be added to the container
        :param container: the container to which the variable will be added
        :return:
        """

        if isinstance(container, str):
            container = getattr(self, container, container)

        if variable.id not in container:
            container[variable.id] = variable

            if variable.model is not self:
                variable.model = self
                # noinspection PyProtectedMember
                variable._model_ref = 1

            elif variable.model is self:
                # noinspection PyProtectedMember
                variable._model_ref += 1

    def _remove_variable_from_container(self, variable, container):
        """
        It removes the given variable from the given container.

        Helper method to be used internally.
        :param variable: the variable to be removed from the container
        :param container: the container from which the variable will be removed
        :return:
        """
        if isinstance(container, str):
            container = getattr(self, container, container)

        if variable.id in container:

            del container[variable.id]

            # noinspection PyProtectedMember
            if variable._model_ref < 2:
                variable.model = None

            elif variable._model_ref > 1:
                # noinspection PyProtectedMember
                variable._model_ref -= 1

    # -----------------------------------------------------------------------------
    # Simulators observer pattern
    # -----------------------------------------------------------------------------
    def attach(self, simulator: 'LinearProblem'):
        """
        It attaches the given simulation method (simulator) to the model.
        Once a simulator is attached to the model, it will be notified with model changes.
        Hence, the simulator is synchronized with the model and following optimizations reflect
        the changes in the model.

        :param simulator: the simulator to be attached to the model
        :return:
        """
        self.simulators.append(simulator)

    def detach(self, simulator):
        """
        It detaches the given simulation method (simulator) from the model.
        Once a simulator is detached from the model, it will not be notified with model changes.
        Hence, the simulator is not synchronized with the model and following optimizations do not reflect
        the changes in the model.

        :param simulator: the simulator to be detached from the model
        :return:
        """
        self.simulators.remove(simulator)

    def notify(self):
        """
        It notifies all simulators with the recent changes in the model.

        :return:
        """
        for simulator in self.simulators:
            simulator.update()

    # -----------------------------------------------------------------------------
    # History manager command pattern
    # -----------------------------------------------------------------------------
    @property
    def contexts(self) -> List[HistoryManager]:
        """
        It returns the list of contexts the model is currently in.
        :return: the list of contexts the model is currently in
        """
        return self._contexts

    @property
    def history(self) -> HistoryManager:
        """
        It returns the history manager of the model.
        If the model is in a context, the history manager of the last context is returned.
        :return: the history manager of the model
        """
        if self.contexts:
            return self.contexts[-1]

        return self._history

    def undo(self):
        """
        It undoes the last action performed on the model.
        :return:
        """
        self.history.undo()

    def redo(self):
        """
        It redoes the last action performed on the model.
        :return:
        """
        self.history.redo()

    def clean_history(self):
        """
        It cleans the history of the model.
        :return:
        """
        self._history = HistoryManager()

    def clean_context(self, index=None):
        """
        It cleans the context of the model.
        :param index: the index of the context to be cleaned
        :return:
        """
        if not index:
            index = -1

        self.contexts[index] = HistoryManager()

    def reset(self):
        """
        It resets the model to its initial state.
        :return:
        """
        self.history.reset()

    def restore(self):
        """
        It restores the model to its last state.
        :return:
        """
        self.history.restore()

    def __enter__(self):
        """
        It enters the context of the model.
        It creates a new context and adds it to the list of contexts.
        All changes performed within this context will be reset when the context is exited.

        :return: the model itself
        """
        self.contexts.append(HistoryManager())

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        It exits the context of the model performing a history reset.

        :param exc_type:
        :param exc_val:
        :param exc_tb:
        :return:
        """
        context = self.contexts.pop()
        context.reset()

    # -----------------------------------------------------------------------------
    # Abstract implementations
    # -----------------------------------------------------------------------------
    # The following polymorphic initializers are just registered here to avoid type checking errors
    @classmethod
    def from_metabolic(cls,
                       identifier: Any,
                       name: str = None,
                       compartments: Dict[str, str] = None,
                       genes: Dict[str, 'Gene'] = None,
                       metabolites: Dict[str, 'Metabolite'] = None,
                       objective: Dict['Reaction', Union[float, int]] = None,
                       reactions: Dict[str, 'Reaction'] = None) -> 'MetabolicModel':
        """
        It creates a metabolic model from the given information.
        :param identifier: the identifier of the model
        :param name: the name of the model
        :param compartments: the compartments of the model
        :param genes: the genes of the model
        :param metabolites: the metabolites of the model
        :param objective: the objective of the model
        :param reactions: the reactions of the model
        :return: a metabolic model
        """
        ...

    @classmethod
    def from_regulatory(cls,
                        identifier: Any,
                        name: str = None,
                        compartments: Dict[str, str] = None,
                        interactions: Dict[str, 'Interaction'] = None,
                        regulators: Dict[str, 'Regulator'] = None,
                        targets: Dict[str, 'Target'] = None) -> 'RegulatoryModel':
        """
        It creates a regulatory model from the given information.
        :param identifier: the identifier of the model
        :param name: the name of the model
        :param compartments: the compartments of the model
        :param interactions: the interactions of the model
        :param regulators: the regulators of the model
        :param targets: the targets of the model
        :return: a regulatory model
        """
        ...

    # -----------------------------------------------------------------------------
    # Type checker
    # -----------------------------------------------------------------------------
    # The following type checkers are just registered here to avoid type checking errors
    def is_metabolic(self) -> bool:
        """
        It checks whether the model is a metabolic model.
        :return: True if the model is a metabolic model, False otherwise
        """
        ...

    def is_regulatory(self) -> bool:
        """
        It checks whether the model is a regulatory model.
        :return: True if the model is a regulatory model, False otherwise
        """
        ...

    def is_a(self, model_type) -> bool:
        """
        It checks whether the model is of the given type.
        :param model_type: the type of the model
        :return: True if the model is of the given type, False otherwise
        """
        if model_type in self.types:
            return True

        return False

    # -----------------------------------------------------------------------------
    # Common attributes - Name mangling
    # -----------------------------------------------------------------------------
    @property
    def compartments(self) -> Dict[str, str]:
        """
        It returns the compartments of the model.
        :return: the compartments of the model as a dictionary
        """
        return self.__compartments

    # -----------------------------------------------------------------------------
    # Helper methods
    # -----------------------------------------------------------------------------

    @staticmethod
    def _get_orphans(to_remove, first_container, second_container):
        """
        It returns the orphans of the given containers.

        Internal use only.
        :param to_remove: the variables to be removed
        :param first_container: the first container
        :param second_container: the second container
        :return:
        """
        orphans = set()

        for variable in to_remove:

            container_iter = getattr(variable,
                                     f'yield_{first_container}',
                                     lambda: [getattr(variable, f'{first_container}')])

            for variable_2 in container_iter():

                remove_variable = True

                container_iter2 = getattr(variable_2,
                                          f'yield_{second_container}',
                                          lambda: [getattr(variable_2, f'{second_container}')])

                for variable_3 in container_iter2():

                    if variable_3 not in to_remove:
                        remove_variable = False

                if remove_variable:
                    orphans.add(variable_2)

        return orphans


def build_model(types: Iterable[str], kwargs: Dict[str, Any]) -> Union['Model', 'MetabolicModel', 'RegulatoryModel']:
    """
    It builds a model from the given types and arguments. Check the `Model.from_types()` method for details.
    :param types: the types of the model
    :param kwargs: the arguments of the model
    :return: a new model instance for the given types and arguments
    """
    return Model.from_types(types, **kwargs)
