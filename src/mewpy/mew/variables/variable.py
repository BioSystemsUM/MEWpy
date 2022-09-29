from typing import Any, Union, Type, TYPE_CHECKING, List, Set, Tuple, Dict, Iterable

from mewpy.util.history import HistoryManager, recorder
from mewpy.mew.models.serialization import serialize, Serializer

# Preventing circular dependencies that only happen due to type checking
if TYPE_CHECKING:
    from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel
    from mewpy.mew.algebra import Expression, Symbolic

    from .gene import Gene
    from .interaction import Interaction
    from .metabolite import Metabolite
    from .reaction import Reaction
    from .regulator import Regulator
    from .target import Target


class MetaVariable(type):
    """
    The MetaVariable class is used to dynamically create Variables classes.
    Variables instances are created from static classes representing the different information for a given variable.

    However, integrated Metabolic-Regulatory models may have different types of variables, such as genes, reactions,
    metabolites, regulators, targets, interactions, but also metabolite-regulator, gene-target, etc.
    Therefore, to avoid the duplication of certain variables, a Variable factory is used to dynamically create
    Variable classes based on the multiple types of a single variable.

    The Variable factory is a class that is the interface to create these multi-type variables.
    It also manages all variables by implementing base attributes and methods.

    The MetaVariable class creates a dynamic Variable class as follows:
        1. It collects all static base classes to be included into the dynamic Variable class
        2. It collects the variable type from the base classes or from the variable_type attribute
        3. It sets the name of the dynamic Variable class using the name of the base classes
        4. It updates the variable_type attribute of the dynamic Variable class with all types of the base classes
        5. It creates the dynamic Variable class having all attributes and methods of the base classes
        6. It adds polymorphic constructors to the dynamic Variable class based on the types of the base classes
        7. It adds type checkers to the dynamic Variable class based on the types of the base classes
    """
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
    def get_serializable_attributes(attrs: Dict[str, Any]) -> Dict[str, Tuple[str, str, str]]:
        """
        Collects all attributes that must be serialized for a particular class or subclass.

        All variables have property descriptors marked with the @serialize decorator. This decorator marks attributes
        for serialization, deserialization and pickle serialization/deserialization. The decorator also provides
        information about the attribute name.

        If serialization, deserialization or pickle serialization/deserialization is not desired, the attribute must
        be marked with the @serialize decorator with the corresponding argument set to None.

        :param attrs: the attributes of the class or subclass
        :return: a dictionary with the attribute name as key and a tuple with the attribute name for serialization,
        deserialization and pickle
        """

        attributes = {}

        for name, method in attrs.items():

            if hasattr(method, 'fget'):

                if hasattr(method.fget, 'serialize') and hasattr(method.fget, 'deserialize') and hasattr(method.fget,
                                                                                                         'pickle'):
                    attributes[name] = (method.fget.serialize, method.fget.deserialize, method.fget.pickle)

        return attributes


class Variable(Serializer, metaclass=MetaVariable, factory=True):
    """
    The Variable class is the base class for all variables and the factory for all variables.
    It is the interface to create variables of different types providing several base attributes and methods.

    Base attributes:
        - types: the types of the variable
        - id: the identifier of the variable
        - name: the name of the variable
        - aliases: the aliases of the variable
        - model: the model to which the variable belongs
        - contexts: the contexts in which the variable is used
        - history: the history of the variable

    Base methods:
        - factory: the factory for dynamic variable types
        - from_types: creates a variable from a list of types
        - copy: creates a copy of the variable
        - deepcopy: creates a deep copy of the variable
        - undo: undoes the last change of the variable
        - redo: redoes the last change of the variable
        - reset: resets the variable to its initial state
        - restore: restores the variable to a previous state
        - update: updates the variable with a new information
        - from polymorphic constructors: creates a variable from a specific type
        - is polymorphic type checkers: checks if a variable is of a specific type

    Serialization:
        - to_dict: serializes the variable to a dictionary
        - from_dict: deserializes the variable from a dictionary

    The Variable class holds the registry of all potential variables.
    """

    # -----------------------------------------------------------------------------
    # Factory management
    # -----------------------------------------------------------------------------
    _registry = {}
    _attributes_registry = {'variable': {}}

    def __init_subclass__(cls, **kwargs):
        """
        This method is called when a subclass is created. It is used to register the subclass in the registry of the
        Variable factory.

        Internal use only.
        :param kwargs:
        :return:
        """

        super(Variable, cls).__init_subclass__(**kwargs)

        # the child type
        variable_type = getattr(cls, 'variable_type', cls.__name__.lower())

        cls.register_type(variable_type, cls)

    @staticmethod
    def get_registry() -> Dict[str, Type['Variable']]:
        """
        Returns the registry of the Variable factory.

        Internal use only.
        :return: the registry of the Variable factory as a dictionary
        """
        return Variable._registry.copy()

    @staticmethod
    def get_attributes_registry() -> Dict[str, Dict[str, Tuple[str, str, str]]]:
        """
        Returns the registry of the attributes of the Variable factory.

        Internal use only.
        :return: the registry of the attributes of the Variable factory as a dictionary
        """
        return Variable._attributes_registry.copy()

    @staticmethod
    def register_type(variable_type: str, child: Type['Variable']):
        """
        Registers a child in the registry of the Variable factory.

        Internal use only.
        :param variable_type: the type of the child
        :param child: the child to be registered
        :return:
        """
        Variable._registry[variable_type] = child

    # -----------------------------------------------------------------------------
    # Serialization
    # -----------------------------------------------------------------------------
    @staticmethod
    def register_attributes(attributes, child):
        """
        Registers the attributes of a child class in the attributes registry of the Variable factory. This is useful for
        serialization and deserialization.

        Internal use only.
        :param attributes: the attributes of the child class
        :param child: the child class
        :return:
        """

        if hasattr(child, 'variable_type'):
            Variable._attributes_registry[child.variable_type] = attributes

        elif child is Variable:
            Variable._attributes_registry['variable'] = attributes

    @property
    def attributes(self) -> Dict[str, Tuple[str, str, str]]:
        """
        Returns the attributes of the variable. This is useful for serialization and deserialization.

        Internal use only.
        :return: the attributes of the variable as a dictionary
        """
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
        """
        It creates a dynamic Variable class from a list of types.
        All types must be registered in the Variable factory.

        Example:
        >>> Variable.factory('gene', 'target')
        <class 'mewpy.mew.variables.GeneTargetVariable'>

        :param args: the types of the variable
        :return: the variable type
        """

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
    def from_types(cls, types: Iterable[str], **kwargs) -> Union['Variable',
                                                                 'Gene',
                                                                 'Interaction',
                                                                 'Metabolite',
                                                                 'Reaction',
                                                                 'Regulator',
                                                                 'Target']:
        """
        It creates a variable instance from a list of types and a dictionary of attributes.
        All types must be registered in the Variable factory so the factory can create a dynamic Variable class.
        Then, the variable is initialized with the provided attributes.

        Example:
        >>> Variable.from_types(['gene', 'target'], id='GENE_1', name='Gene 1')
        <mewpy.mew.variables.GeneTargetVariable object at 0x7f8b8c0b7a90>

        :param types: the types of the variable
        :param kwargs: the attributes to initialize the variable
        :return: the variable instance
        """
        VariableType = cls.factory(*types)
        return VariableType(**kwargs)

    # -----------------------------------------------------------------------------
    # Factory helper methods
    # -----------------------------------------------------------------------------
    @staticmethod
    def _from(cls):
        """
        Method to be added to the subclasses of Variable to create a variable from a specific type.

        Internal use only.
        :param cls: the class of the variable
        :return: the method to be added to the subclasses of Variable
        """

        def from_(identifier, **kwargs):
            return cls(identifier, **kwargs)

        return from_

    @staticmethod
    def _is(variable_type):
        """
        Method to be added to the subclasses of Variable to check if a variable is of a specific type.

        Internal use only.
        :param variable_type: the type of the variable
        :return: the method to be added to the subclasses of Variable
        """

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
        Variable is the most base type of all variables in MEWpy,
        such as gene, interaction, metabolite, reaction, regulator and target, among others.
        See also variables.Gene, variables.Interaction, variables.Metabolite, variables.Reaction,
        variables.Regulator and variables.Target for concrete implementations of a variable type.

        The variable is also the factory for several variable types.
        The Variable factory provides all batteries to create new variable types
        and to manage the variables types' inheritance.

        The factory type assists with systematic tasks such as attributes and containers serialization,
        history management, among others.

        For that, attributes, containers and inheritance are registered in the factory. The MetaVariable manages
        both the factory and the variable types derived from the factory.

        :param identifier: identifier, e.g. b0001
        :param name: the name of the variable
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

        # history
        self._contexts = []
        self._history = HistoryManager()

        self._model_ref = 0

    def _check_inheritance(self):
        """
        It checks if the class is a subclass of Variable and if it is not, it raises an error.
        :return:
        """

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
        return self.__str__()

    def _variable_to_html(self):
        """
        It returns an HTML dict representation of the variable.
        """
        html_dict = {'Identifier': self.id,
                     'Name': self.name,
                     'Aliases': ', '.join(self.aliases),
                     'Model': self.model.id if self.model else None,
                     'Types': ', '.join(self.types)}
        return html_dict

    def _repr_html_(self):
        """
        It returns a html representation.
        """
        html_dict = self._variable_to_html()

        for type_ in self.types:
            method = getattr(self, f'_{type_}_to_html', lambda: {})
            html_dict.update(method())

        html_representation = ''
        for key, value in html_dict.items():
            html_representation += f'<tr><th>{key}</th><td>{value}</td></tr>'

        return f"""
            <table>
                {html_representation}
            </table>
        """

    # -----------------------------------------------------------------------------
    # Variable type manager
    # -----------------------------------------------------------------------------
    @serialize('types', None)
    @property
    def types(self) -> Set[str]:
        """
        It returns the types of the variable.
        :return: the types of the variable as a set of strings
        """
        return set()

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------

    @serialize('id', None, '_id')
    @property
    def id(self) -> Any:
        """
        It returns the identifier of the variable.
        :return: the identifier of the variable
        """
        return self._id

    @serialize('name', 'name', '_name')
    @property
    def name(self) -> str:
        """
        It returns the name of the variable.
        :return: the name of the variable as a string
        """
        return self._name

    @serialize('aliases', 'aliases', '_aliases')
    @property
    def aliases(self) -> Set[str]:
        """
        It returns the aliases of the variable.
        :return: the aliases of the variable as a set of strings
        """
        return self._aliases

    @property
    def model(self) -> Union['Model', 'MetabolicModel', 'RegulatoryModel']:
        """
        It returns the model to which the variable belongs.
        Note that variables do not need to be directly associated with a model.
        :return: the model to which the variable belongs
        """
        return self._model

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------
    @name.setter
    @recorder
    def name(self, value: str):
        """
        It sets the name of the variable.
        :param value: the name of the variable
        :return:
        """
        if not value:
            value = ''

        self._name = value

    @aliases.setter
    @recorder
    def aliases(self, value: Set[str]):
        """
        It sets the aliases of the variable.
        :param value: the aliases of the variable
        :return:
        """
        if not value:
            value = set()

        self._aliases = value

    @model.setter
    def model(self, value: Union['Model', 'MetabolicModel', 'RegulatoryModel']):
        """
        It sets the model to which the variable belongs.
        This setter does not perform any check or action in the model.

        :param value: the model to which the variable belongs
        :return:
        """
        if not value:
            value = None

        self._model = value

    # -----------------------------------------------------------------------------
    # History manager
    # -----------------------------------------------------------------------------

    # The history manager follows the command pattern
    @property
    def contexts(self) -> List[HistoryManager]:
        """
        It returns the contexts of the variable.
        :return: the contexts of the variable as a list of HistoryManager
        """

        if self.model:
            return self.model.contexts

        return self._contexts

    @property
    def history(self) -> HistoryManager:
        """
        It returns the history manager of the variable.
        If the variable is within a context, the history manager of the last context is returned.
        :return: the history manager of the variable
        """
        if self.contexts:
            return self.contexts[-1]

        if self.model:
            return self.model.history

        return self._history

    def undo(self):
        """
        It undoes the last action performed on the variable.
        :return:
        """
        self.history.undo()

    def redo(self):
        """
        It redoes the last action performed on the variable.
        :return:
        """
        self.history.redo()

    def clean_history(self):
        """
        It cleans the history of the variable.
        :return:
        """
        self._history = HistoryManager()

    def clean_context(self):
        """
        It cleans the context of the variable.
        :return:
        """
        self.contexts[:-1] = HistoryManager()

    def reset(self):
        """
        It resets the variable to its initial state.
        :return:
        """
        self.history.reset()

    def restore(self):
        """
        It restores the variable to its last state.
        :return:
        """
        self.history.restore()

    def __enter__(self):
        """
        It enters the context of the variable.
        It creates a new context and adds it to the list of contexts.
        All changes performed within this context will be reset when the context is exited.

        :return: the variable itself
        """
        self.contexts.append(HistoryManager())

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        It exits the context of the variable performing a history reset.

        :param exc_type:
        :param exc_val:
        :param exc_tb:
        :return:
        """
        context = self.contexts.pop()
        context.reset()

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------
    def update(self,
               name: str = None,
               aliases: Set[str] = None,
               model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None):
        """
        It updates the attributes of the variable.
        :param name: the name of the variable
        :param aliases: the aliases of the variable
        :param model: the model to which the variable belongs
        :return:
        """
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
        """
        It creates a gene from the given parameters.
        :param identifier: the identifier of the gene
        :param name: the name of the gene
        :param aliases: the aliases of the gene
        :param model: the model to which the gene belongs
        :param coefficients: the coefficients of the gene
        :param active_coefficient: the default coefficient of the gene
        :param reactions: the reactions associated with the gene
        :return: a gene object
        """
        ...

    @classmethod
    def from_interaction(cls,
                         identifier: Any,
                         name: str = None,
                         aliases: set = None,
                         model: 'Model' = None,
                         target: 'Target' = None,
                         regulatory_events: Dict[Union[float, int], 'Expression'] = None) -> 'Interaction':
        """
        It creates an interaction from the given parameters.
        :param identifier: the identifier of the interaction
        :param name: the name of the interaction
        :param aliases: the aliases of the interaction
        :param model: the model to which the interaction belongs
        :param target: the target of the interaction
        :param regulatory_events: the regulatory events associated with the interaction
        :return: an interaction object
        """
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
        """
        It creates a metabolite from the given parameters.
        :param identifier: the identifier of the metabolite
        :param name: the name of the metabolite
        :param aliases: the aliases of the metabolite
        :param model: the model to which the metabolite belongs
        :param charge: the charge of the metabolite
        :param compartment: the compartment of the metabolite
        :param formula: the formula of the metabolite
        :param reactions: the reactions associated with the metabolite
        :return: a metabolite object
        """
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
        """
        It creates a reaction from the given parameters.
        :param identifier: the identifier of the reaction
        :param name: the name of the reaction
        :param aliases: the aliases of the reaction
        :param model: the model to which the reaction belongs
        :param bounds: the bounds of the reaction
        :param stoichiometry: the stoichiometry of the reaction
        :param gpr: the gene-protein-reaction association of the reaction
        :return: a reaction object
        """
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
        """
        It creates a regulator from the given parameters.
        :param identifier: the identifier of the regulator
        :param name: the name of the regulator
        :param aliases: the aliases of the regulator
        :param model: the model to which the regulator belongs
        :param coefficients: the coefficients of the regulator
        :param active_coefficient: the default coefficient of the regulator
        :param interactions: the interactions associated with the regulator
        :return: a regulator object
        """
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
        """
        It creates a target from the given parameters.
        :param identifier: the identifier of the target
        :param name: the name of the target
        :param aliases: the aliases of the target
        :param model: the model to which the target belongs
        :param coefficients: the coefficients of the target
        :param active_coefficient: the default coefficient of the target
        :param interaction: the interaction associated with the target
        :return: a target object
        """
        ...

    # The following type checkers are just registered here to avoid type checking errors
    def is_gene(self) -> bool:
        """
        It checks whether the variable is a gene.
        :return: True if the variable is a gene, False otherwise
        """
        ...

    def is_interaction(self) -> bool:
        """
        It checks whether the variable is an interaction.
        :return: True if the variable is an interaction, False otherwise
        """
        ...

    def is_metabolite(self) -> bool:
        """
        It checks whether the variable is a metabolite.
        :return: True if the variable is a metabolite, False otherwise
        """
        ...

    def is_reaction(self) -> bool:
        """
        It checks whether the variable is a reaction.
        :return: True if the variable is a reaction, False otherwise
        """
        ...

    def is_regulator(self) -> bool:
        """
        It checks whether the variable is a regulator.
        :return: True if the variable is a regulator, False otherwise
        """
        ...

    def is_target(self) -> bool:
        """
        It checks whether the variable is a target.
        :return: True if the variable is a target, False otherwise
        """
        ...

    # -----------------------------------------------------------------------------
    # Type checker
    # -----------------------------------------------------------------------------

    def is_a(self, variable_type) -> bool:
        """
        It checks whether the variable is of the given type.
        :param variable_type: the type to check
        :return: True if the variable is of the given type, False otherwise
        """
        if variable_type in self.types:
            return True

        return False


# -----------------------------------------------------------------------------
# Helper methods
# -----------------------------------------------------------------------------
def variables_from_symbolic(symbolic: 'Symbolic',
                            types: Union[Set[str], List[str], Tuple[str]],
                            model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None) -> Dict[str, 'Variable']:
    """
    It returns the variables of the given type from the given symbolic expression.
    :param symbolic: the symbolic expression to parse and extract the variables from
    :param types: the types of the variables to return
    :param model: the model to which the variables belong
    :return: the variables of the given type from the given symbolic expression
    """
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


def build_variable(types: Iterable[str], kwargs: Dict[str, Any]) -> Union['Variable',
                                                                          'Gene',
                                                                          'Interaction',
                                                                          'Metabolite',
                                                                          'Reaction',
                                                                          'Regulator',
                                                                          'Target']:
    """
    It builds a variable from the given types and keyword arguments. Check the `Variable.from_types()` method for more
    details.
    :param types: the types of the variable
    :param kwargs: the keyword arguments to build the variable
    :return: a new variable of the given types
    """
    return Variable.from_types(types, **kwargs)
