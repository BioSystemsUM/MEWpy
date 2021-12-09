from typing import Any, Union, Type, TYPE_CHECKING, List, Set, Dict

from mewpy.util.history import HistoryManager, recorder
from mewpy.util.serialization import serialize, Serializer
from mewpy.mew.lp import Notification

# Preventing circular dependencies that only happen due to type checking
if TYPE_CHECKING:
    from mewpy.model import MetabolicModel, RegulatoryModel
    from mewpy.mew.variables import Gene, Interaction, Metabolite, Reaction, Regulator, Target
    from mewpy.mew.lp import LinearProblem


class MetaModel(type):
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
    def get_serializable_containers(attrs):

        # property getter is marked with serialize and deserialize attribute by the serialize decorator. If
        # deserialization is not to be used, deserialize should be marked with None

        containers = {}

        for name, method in attrs.items():

            if hasattr(method, 'fget'):

                if hasattr(method.fget, 'serialize') and hasattr(method.fget, 'deserialize') and hasattr(method.fget,
                                                                                                         'pickle'):
                    containers[name] = (method.fget.serialize, method.fget.deserialize, method.fget.pickle)

        return containers


# TODO: methods stubs and type hinting
class Model(Serializer, metaclass=MetaModel, factory=True):
    # -----------------------------------------------------------------------------
    # Factory management
    # -----------------------------------------------------------------------------

    _registry = {}
    _containers_registry = {'model': {}}

    def __init_subclass__(cls, **kwargs):

        super(Model, cls).__init_subclass__(**kwargs)

        # the child type
        model_type = getattr(cls, 'model_type', cls.__name__.lower())

        cls.register_type(model_type, cls)

    @staticmethod
    def get_registry():
        return Model._registry.copy()

    @staticmethod
    def get_containers_registry():
        return Model._containers_registry.copy()

    @staticmethod
    def register_type(model_type, child):

        Model._registry[model_type] = child

    # -----------------------------------------------------------------------------
    # Serialization
    # -----------------------------------------------------------------------------

    @staticmethod
    def register_containers(containers, child):

        if hasattr(child, 'model_type'):

            Model._containers_registry[child.model_type] = containers

        elif child is Model:
            Model._containers_registry['model'] = containers

    @property
    def containers(self):

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
    def from_types(cls, types: List[str], **kwargs) -> Union['Model',
                                                             'MetabolicModel',
                                                             'RegulatoryModel']:
        ModelType = cls.factory(*types)
        return ModelType(**kwargs)

    # -----------------------------------------------------------------------------
    # Factory helper methods
    # -----------------------------------------------------------------------------

    @staticmethod
    def _from(cls):

        def from_(identifier, **kwargs):
            return cls(identifier, **kwargs)

        return from_

    @staticmethod
    def _is(model_type):

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

        The model is the most base type of a given model type, such as metabolic model or regulatory model.
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

        registry = self.get_registry()

        for model_type in self.types:
            if model_type not in registry:
                raise ValueError(f'{model_type} is not registered as subclass of {self.__class__.__name__}')

    # -----------------------------------------------------------------------------
    # Built-in
    # -----------------------------------------------------------------------------

    def __str__(self):
        return self.name

    def __repr__(self):
        if self.types:
            types = [v_type.title() for v_type in self.types]

            return ', '.join(types) + f': {self.id}'

        return f'Model: {self.id}'

    # -----------------------------------------------------------------------------
    # Variable type manager
    # -----------------------------------------------------------------------------

    @serialize('types', None, None)
    @property
    def types(self) -> Set[str]:
        return set()

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------

    @serialize('id', None, '_id')
    @property
    def id(self) -> Any:
        return self._id

    @serialize('name', 'name', '_name')
    @property
    def name(self) -> str:
        return self._name

    @property
    def simulators(self) -> List['LinearProblem']:
        return self._simulators

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------

    @name.setter
    @recorder
    def name(self, value: str):

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
        return default

    # add and remove are just registered here to avoid type checking errors

    def add(self, variables, *types, comprehensive=True, history=True):

        pass

    def remove(self, variables, *types, remove_orphans=False, history=True):

        pass

    def update(self, name: str = None):

        if name:
            self.name = name

    def _add_variable_to_container(self, variable, container):

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

        self.simulators.append(simulator)

    def detach(self, simulator):

        self.simulators.remove(simulator)

    def notify(self, notification: Notification):

        for simulator in self.simulators:
            simulator.notification(notification)

    # -----------------------------------------------------------------------------
    # History manager
    # -----------------------------------------------------------------------------

    # The history manager follows the command pattern
    @property
    def contexts(self) -> List[HistoryManager]:
        return self._contexts

    @property
    def history(self) -> HistoryManager:

        if self.contexts:
            return self.contexts[-1]

        return self._history

    def undo(self):

        self.history.undo()

    def redo(self):

        self.history.redo()

    def clean_history(self):

        self._history = HistoryManager()

    def clean_context(self, index=None):

        if not index:
            index = -1

        self.contexts[index] = HistoryManager()

    def reset(self):

        self.history.reset()

    def restore(self):

        self.history.restore()

    def __enter__(self):

        self.contexts.append(HistoryManager())

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):

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
        ...

    @classmethod
    def from_regulatory(cls,
                        identifier: Any,
                        name: str = None,
                        compartments: Dict[str, str] = None,
                        interactions: Dict[str, 'Interaction'] = None,
                        regulators: Dict[str, 'Regulator'] = None,
                        targets: Dict[str, 'Target'] = None) -> 'RegulatoryModel':
        ...

    # The following type checkers are just registered here to avoid type checking errors
    def is_metabolic(self) -> bool:
        ...

    def is_regulatory(self) -> bool:
        ...

    # -----------------------------------------------------------------------------
    # Type checker
    # -----------------------------------------------------------------------------

    def is_a(self, model_type) -> bool:
        if model_type in self.types:
            return True

        return False

    # -----------------------------------------------------------------------------
    # Common attributes - Name mangling
    # -----------------------------------------------------------------------------

    # TODO: perhaps name mangling should be managed by the metaclass, as more common attributes can occur
    @property
    def compartments(self):

        return self.__compartments

    # -----------------------------------------------------------------------------
    # Helper methods
    # -----------------------------------------------------------------------------

    @staticmethod
    def _get_orphans(to_remove, first_container, second_container):

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


def build_model(types, kwargs):
    return Model.from_types(types, **kwargs)
