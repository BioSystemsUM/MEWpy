from typing import Any, Union, Type, TYPE_CHECKING, List, Set, Tuple, Dict

from mewpy.util.history import HistoryManager, recorder
from mewpy.util.serialization import serialize, Serializer

# Preventing circular dependencies that only happen due to type checking
if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel
    from mewpy.mew.algebra import Expression, Symbolic

    from .gene import Gene
    from .interaction import Interaction
    from .metabolite import Metabolite
    from .reaction import Reaction
    from .regulator import Regulator
    from .target import Target


class MetaVariable(type):
    factories = {}

    def __new__(mcs, name, bases, attrs, **kwargs):

        # if it is the variable factory, only registration is done
        factory = kwargs.get('factory', False)

        if factory:
            cls = super(MetaVariable, mcs).__new__(mcs, name, bases, attrs)

            MetaVariable.factories[name] = cls

            return cls

        # Dynamic typing being used. In this case, a proper name and variable type must be provided
        dynamic = kwargs.get('dynamic', False)
        if dynamic:
            names = [base.variable_type for base in bases]

            name = ''.join([name.title() for name in names])
            name += 'Variable'

            kwargs['variable_type'] = '-'.join(names)

        # The variable type is always added to the subclasses. If it is not given upon subclass creation,
        # the subclass name is to be used
        variable_type = kwargs.get('variable_type', name.lower())
        attrs['variable_type'] = variable_type

        return super(MetaVariable, mcs).__new__(mcs, name, bases, attrs)

    def __init__(cls, name, bases, attrs, **kwargs):

        super().__init__(name, bases, attrs)

        factory = kwargs.get('factory', False)
        if factory:
            # collection of all attributes that must be serialized for a particular class or subclass
            attributes = cls.get_serializable_attributes(attrs)

            attrs['_attributes_registry']['variable'] = attributes

            # Skip further building of the Variable factory
            return

        # Dynamic typing being used. In this case, all children have already been constructed, so everything can be
        # skipped
        dynamic = kwargs.get('dynamic', False)
        if dynamic:
            return

        # Several attributes and methods must be added automatically to the Variable factory children based on the
        # variable type.
        variable_type = attrs['variable_type']

        for base in bases:

            factory = MetaVariable.factories.get(base.__name__)

            if factory and hasattr(cls, 'variable_type'):

                # if some class inherits from the Variable factory, it must be registered in the factory for type
                # checking
                if cls.variable_type not in factory.get_registry():
                    raise TypeError(f'{cls.variable_type} does not inherit from Variable')

                # If set otherwise upon class creation, all subclasses of the Variable factory will contribute to
                # their parent factory with an alternative initializer. The polymorphic constructor, e.g. from_{
                # variable_type}, is added to the Variable factory
                constructor = kwargs.get('constructor', True)

                if constructor:
                    # noinspection PyProtectedMember
                    var_type_initializer = Variable._from(cls)
                    setattr(factory, f'from_{variable_type}', var_type_initializer)

                # If set otherwise upon class creation, all subclasses of the Variable factory will contribute to
                # their parent factory with a declared type checker. The type checker, e.g. is_{variable_type},
                # is added to the Variable factory
                checker = kwargs.get('checker', True)

                if checker:
                    # noinspection PyProtectedMember
                    var_type_checker = Variable._is(variable_type)
                    setattr(factory, f'is_{variable_type}', var_type_checker)

                # collection of all attributes that must be serialized for a particular class or sub class
                attributes = cls.get_serializable_attributes(attrs)

                Variable.register_attributes(attributes, cls)

    @staticmethod
    def get_serializable_attributes(attrs):

        # property getter is marked with serialize and deserialize attribute by the serialize decorator. If
        # deserialization is not to be used, deserialize should be marked with None

        attributes = {}

        for name, method in attrs.items():

            if hasattr(method, 'fget'):

                if hasattr(method.fget, 'serialize') and hasattr(method.fget, 'deserialize') and hasattr(method.fget,
                                                                                                         'pickle'):
                    attributes[name] = (method.fget.serialize, method.fget.deserialize, method.fget.pickle)

        return attributes


# TODO: methods stubs
class Variable(Serializer, metaclass=MetaVariable, factory=True):
    # -----------------------------------------------------------------------------
    # Factory management
    # -----------------------------------------------------------------------------

    _registry = {}
    _attributes_registry = {'variable': {}}

    def __init_subclass__(cls, **kwargs):

        super(Variable, cls).__init_subclass__(**kwargs)

        # the child type
        variable_type = getattr(cls, 'variable_type', cls.__name__.lower())

        cls.register_type(variable_type, cls)

    @staticmethod
    def get_registry():
        return Variable._registry.copy()

    @staticmethod
    def get_attributes_registry():
        return Variable._attributes_registry.copy()

    @staticmethod
    def register_type(variable_type, child):

        Variable._registry[variable_type] = child

    # -----------------------------------------------------------------------------
    # Serialization
    # -----------------------------------------------------------------------------

    @staticmethod
    def register_attributes(attributes, child):

        if hasattr(child, 'variable_type'):
            Variable._attributes_registry[child.variable_type] = attributes

        elif child is Variable:
            Variable._attributes_registry['variable'] = attributes

    @property
    def attributes(self):

        class_attributes = self.get_attributes_registry()

        attributes = class_attributes['variable'].copy()

        for variable_type in self.types:

            for name, (serialize_name, deserialize_name, pickle_name) in class_attributes[variable_type].items():
                attributes[name] = (serialize_name, deserialize_name, pickle_name)

        return attributes

    # -----------------------------------------------------------------------------
    # Factory polymorphic constructor
    # -----------------------------------------------------------------------------

    @classmethod
    def factory(cls, *args: str) -> Union[Type['Variable'],
                                          Type['Gene'],
                                          Type['Interaction'],
                                          Type['Metabolite'],
                                          Type['Reaction'],
                                          Type['Regulator'],
                                          Type['Target']]:

        if not args:
            args = ()

        registry = cls.get_registry()

        types = tuple([registry[name] for name in args])

        if len(types) == 1:
            return types[0]

        _Variable = MetaVariable('Variable', types, {}, dynamic=True)

        return _Variable

    # -----------------------------------------------------------------------------
    # Factory polymorphic initializer
    # -----------------------------------------------------------------------------

    @classmethod
    def from_types(cls, types: List[str], **kwargs) -> Union['Variable',
                                                             'Gene',
                                                             'Interaction',
                                                             'Metabolite',
                                                             'Reaction',
                                                             'Regulator',
                                                             'Target']:
        VariableType = cls.factory(*types)
        return VariableType(**kwargs)

    # -----------------------------------------------------------------------------
    # Factory helper methods
    # -----------------------------------------------------------------------------

    @staticmethod
    def _from(cls):

        def from_(identifier, **kwargs):
            return cls(identifier, **kwargs)

        return from_

    @staticmethod
    def _is(variable_type):

        def is_(self):
            return self.is_a(variable_type)

        return is_

    # -----------------------------------------------------------------------------
    # Base initializer
    # -----------------------------------------------------------------------------

    def __init__(self,
                 identifier: Any,
                 name: str = None,
                 aliases: Set[str] = None,
                 model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None):

        """

        Variable: The variable is the most base type of a given variable type,
        such as gene, interaction, metabolite, reaction, regulator and target, among others.
        See also variables.Gene, variables.Interaction, variables.Metabolite, variables.Reaction,
        variables.Regulator and variables.Target for concrete implementations of a variable type.

        The variable object is actually the factory for several variable types.
        This variable object provides all batteries to create new variable types
        and to manage the variables types' inheritance.

        The factory type assists with systematic tasks such as attributes and containers serialization,
        history management, among others.

        For that, attributes, containers and inheritance are registered in the factory. The MetaVariable manages
        both the factory and the variable types derived from the factory.

        :param identifier: identifier, e.g. b0001
        :param name: the name of the model
        :param aliases: aliases for this variable
        :param model: the model to which this variable belongs. Since the default is None, variables do not need to
        be directly associated with a model
        """

        self._check_inheritance()

        if not aliases:
            aliases = set()

        if not model:
            model = None

        if not name:
            name = ''

        self._id = identifier
        self._aliases = aliases
        self._model = model
        self._name = name

        self._contexts = []
        self._history = HistoryManager()

        self._model_ref = 0

    def _check_inheritance(self):

        registry = self.get_registry()

        for variable_type in self.types:
            if variable_type not in registry:
                raise ValueError(f'{variable_type} is not registered as subclass of {self.__class__.__name__}')

    # -----------------------------------------------------------------------------
    # Built-in
    # -----------------------------------------------------------------------------

    def __str__(self):
        return self.name

    def __repr__(self):
        if self.types:
            types = [v_type.title() for v_type in self.types]

            return ', '.join(types) + f': {self.id}'

        return f'Variable: {self.id}'

    # -----------------------------------------------------------------------------
    # Variable type manager
    # -----------------------------------------------------------------------------

    @serialize('types', None)
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

    @serialize('aliases', 'aliases', '_aliases')
    @property
    def aliases(self) -> Set[str]:
        return self._aliases

    @property
    def model(self) -> Union['Model', 'MetabolicModel', 'RegulatoryModel']:
        return self._model

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------

    @name.setter
    @recorder
    def name(self, value: str):

        if not value:
            value = ''

        self._name = value

    @aliases.setter
    @recorder
    def aliases(self, value: Set[str]):

        if not value:
            value = set()

        self._aliases = value

    @model.setter
    def model(self, value: Union['Model', 'MetabolicModel', 'RegulatoryModel']):

        if not value:
            value = None

        self._model = value

    # -----------------------------------------------------------------------------
    # History manager
    # -----------------------------------------------------------------------------

    # The history manager follows the command pattern
    @property
    def contexts(self) -> List[HistoryManager]:

        if self.model:
            return self.model.contexts

        return self._contexts

    @property
    def history(self) -> HistoryManager:

        if self.contexts:
            return self.contexts[-1]

        if self.model:
            return self.model.history

        return self._history

    def undo(self):

        self.history.undo()

    def redo(self):

        self.history.redo()

    def clean_history(self):

        self._history = HistoryManager()

    def clean_context(self):

        self.contexts[:-1] = HistoryManager()

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
    # Operations/Manipulations
    # -----------------------------------------------------------------------------

    def update(self,
               name: str = None,
               aliases: Set[str] = None,
               model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None):

        if name is not None:
            self.name = name

        if aliases is not None:
            self.aliases.update(set(aliases))

        if model is not None:
            self.model = model

    # -----------------------------------------------------------------------------
    # Abstract implementations
    # -----------------------------------------------------------------------------

    # The following polymorphic initializers are just registered here to avoid type checking errors
    @classmethod
    def from_gene(cls,
                  identifier: Any,
                  name: str = None,
                  aliases: set = None,
                  model: 'Model' = None,
                  coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]] = None,
                  active_coefficient: Union[int, float] = None,
                  reactions: Dict[str, 'Reaction'] = None) -> 'Gene':
        ...

    @classmethod
    def from_interaction(cls,
                         identifier: Any,
                         name: str = None,
                         aliases: set = None,
                         model: 'Model' = None,
                         target: 'Target' = None,
                         regulatory_events: Dict[Union[float, int], 'Expression'] = None) -> 'Interaction':
        ...

    @classmethod
    def from_metabolite(cls,
                        identifier: Any,
                        name: str = None,
                        aliases: set = None,
                        model: 'Model' = None,
                        charge: int = None,
                        compartment: str = None,
                        formula: str = None,
                        reactions: Dict[str, 'Reaction'] = None) -> 'Metabolite':
        ...

    @classmethod
    def from_reaction(cls,
                      identifier: Any,
                      name: str = None,
                      aliases: set = None,
                      model: 'Model' = None,
                      bounds: Tuple[Union[float, int], Union[float, int]] = None,
                      stoichiometry: Dict['Metabolite', Union[float, int]] = None,
                      gpr: 'Expression' = None) -> 'Reaction':
        ...

    @classmethod
    def from_regulator(cls,
                       identifier: Any,
                       name: str = None,
                       aliases: set = None,
                       model: 'Model' = None,
                       coefficients: Set[Union[float, int]] = None,
                       active_coefficient: Union[float, int] = None,
                       interactions: Dict[str, 'Interaction'] = None) -> 'Regulator':
        ...

    @classmethod
    def from_target(cls,
                    identifier: Any,
                    name: str = None,
                    aliases: set = None,
                    model: 'Model' = None,
                    coefficients: Set[Union[float, int]] = None,
                    active_coefficient: Union[float, int] = None,
                    interaction: 'Interaction' = None) -> 'Target':
        ...

    # The following type checkers are just registered here to avoid type checking errors
    def is_gene(self) -> bool:

        ...

    def is_interaction(self) -> bool:

        ...

    def is_metabolite(self) -> bool:

        ...

    def is_reaction(self) -> bool:

        ...

    def is_regulator(self) -> bool:

        ...

    def is_target(self) -> bool:

        ...

    # -----------------------------------------------------------------------------
    # Type checker
    # -----------------------------------------------------------------------------

    def is_a(self, variable_type) -> bool:

        if variable_type in self.types:
            return True

        return False


# -----------------------------------------------------------------------------
# Helper methods
# -----------------------------------------------------------------------------
def variables_from_symbolic(symbolic: 'Symbolic',
                            types: Union[Set[str], List[str], Tuple[str]],
                            model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None):
    variables = {}

    for symbol in symbolic.atoms(symbols_only=True):

        if symbol.name in variables:
            continue

        if model:

            variable = model.get(symbol.name)

            if variable:

                for variable_type in types:

                    if not variable.is_a(variable_type):
                        raise TypeError(f'{symbol.name} is not a {variable_type} in model {model.id}')

        else:
            variable = Variable.from_types(identifier=symbol.name,
                                           types=types,
                                           model=model)

        variables[variable.id] = variable

    return variables


def build_variable(types, kwargs):
    return Variable.from_types(types, **kwargs)
