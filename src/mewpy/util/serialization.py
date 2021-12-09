from typing import Union, TYPE_CHECKING, Type, Dict
import sys

from mewpy.mew.algebra import Expression
from mewpy.mew.algebra import parse_expression

if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel
    from mewpy.mew.variables import Variable, Gene, Interaction, Metabolite, Reaction, Regulator, Target

sys.setrecursionlimit(50000)


def serialize(serialization_name, deserialization_name=None, pickle_name=None):
    def wrapper(attr):
        attr.fget.serialize = serialization_name
        attr.fget.deserialize = deserialization_name
        attr.fget.pickle = pickle_name

        return attr

    return wrapper


class Serializer:

    # -----------------------------------------------------------------------------
    # Serializing to dictionary, tuple, string and int or float
    # -----------------------------------------------------------------------------

    @staticmethod
    def _regulatory_events_serializer(regulatory_events):

        return {key: expression.to_string() for key, expression in regulatory_events.items()}

    @staticmethod
    def _expression_serializer(expression):

        return expression.to_string()

    @staticmethod
    def _obj_serializer(obj):

        if isinstance(obj, set):

            return tuple(obj)

        elif hasattr(obj, 'id'):

            return obj.id

        return obj

    @staticmethod
    def _coefficient_serializer(coefficient):

        return tuple(coefficient.coefficients)

    @staticmethod
    def _variable_container_serializer(container):

        return tuple(container.keys())

    @staticmethod
    def _key_container_serializer(container):

        return {key.id: value for key, value in container.items()}

    @staticmethod
    def _model_container_serializer(container):

        return {variable.id: variable.to_dict(serialization_format='json') for variable in container.values()}

    def _get_attribute_serializer(self, attr):

        if attr in ('id', 'name', 'types', 'aliases', 'target', 'charge', 'compartment', 'formula', 'compartments',
                    'interaction'):

            return self._obj_serializer

        elif attr == 'coefficient':

            return self._coefficient_serializer

        elif attr in ('reactions', 'regulators', 'genes', 'interactions', 'targets'):

            return self._variable_container_serializer

        elif attr == 'regulatory_events':

            return self._regulatory_events_serializer

        elif attr == 'gpr':

            return self._expression_serializer

        elif attr == 'stoichiometry':

            return self._key_container_serializer

        else:
            return lambda *args, **kwargs: {}

    def _variable_serializer(self: Union['Serializer', 'Variable', 'Model']):

        variable = {}

        for _, (attr, _, _) in self.attributes.items():
            serializer = self._get_attribute_serializer(attr=attr)

            attribute = getattr(self, attr)

            variable[attr] = serializer(attribute)

        return variable

    def _get_container_serializer(self, attr, variables=True):

        if attr in ('id', 'name', 'types'):

            return self._obj_serializer

        elif attr == 'objective':

            return self._key_container_serializer

        elif attr in ('genes', 'metabolites', 'reactions', 'interactions', 'regulators', 'targets'):

            if variables:
                return self._model_container_serializer

            return dict

        else:
            return lambda *args, **kwargs: {}

    def _model_serializer(self: Union['Serializer', 'Variable', 'Model'],
                          variables=True):

        model = {}

        for _, (attr, _, _) in self.containers.items():
            model[attr] = {}

            serializer = self._get_container_serializer(attr=attr, variables=variables)

            container = getattr(self, attr, {})

            model[attr] = serializer(container)

        return model

    # -----------------------------------------------------------------------------
    # Deserializing from dictionary, tuple, string and int or float
    # -----------------------------------------------------------------------------

    @classmethod
    def _model_deserializer(cls: Union[Type['Serializer'], Type['Variable'], Type['Model']],
                            obj,
                            variables=False):

        identifier = obj.get('id')
        types = obj.get('types')

        model = cls.from_types(types=types, identifier=identifier)

        if variables:
            children = cls._build_children(obj=obj, model=model)

        else:
            children = cls._get_children(obj=obj, model=model)

        update_attributes = {}
        for attr_name, (serialize_name, deserialize_name, _) in model.containers.items():

            if deserialize_name is not None:
                deserializer = cls._get_container_deserializer(attr=deserialize_name)

                container = obj[serialize_name]
                container = deserializer(container, children=children)

                if serialize_name in ('genes', 'metabolites', 'reactions', 'interactions', 'regulators', 'targets'):
                    setattr(model, f'_{serialize_name}', container)

                else:
                    update_attributes[deserialize_name] = container

        model.update(**update_attributes)

        return model

    @staticmethod
    def _build_children(obj, model):

        from mewpy.mew.variables import Variable

        children = {}

        for attr_name, (serialize_name, deserialize_name, _) in model.containers.items():

            if deserialize_name is not None and serialize_name in ('genes', 'metabolites', 'reactions',
                                                                   'interactions', 'regulators', 'targets'):

                container = obj[serialize_name]

                for var_id, variable in container.items():

                    if var_id in children:

                        continue

                    else:

                        children[var_id] = Variable.from_types(variable['types'],
                                                               identifier=variable['id'],
                                                               model=model)

        return children

    @staticmethod
    def _get_children(obj, model):

        children = {}

        for attr_name, (serialize_name, deserialize_name, _) in model.containers.items():

            if deserialize_name is not None and serialize_name in ('genes', 'metabolites', 'reactions',
                                                                   'interactions', 'regulators', 'targets'):

                container = obj[serialize_name]

                for var_id, variable in container.items():

                    if var_id in children:

                        continue

                    else:

                        children[var_id] = variable

        return children

    @classmethod
    def _get_container_deserializer(cls, attr):

        if attr == 'name':

            return cls._obj_deserializer

        elif attr in ('genes', 'metabolites', 'reactions',
                      'interactions', 'regulators', 'targets'):

            return cls._model_container_deserializer

        elif attr == 'objective':

            return cls._key_container_deserializer

        else:
            return lambda *args, **kwargs: {}

    @staticmethod
    def _obj_deserializer(obj, *args, **kwargs):

        return obj

    @classmethod
    def _model_container_deserializer(cls, obj, children):

        from mewpy.mew.variables import Variable

        container = {}

        for var_id, variable in obj.items():

            if isinstance(variable, Variable):

                container[var_id] = variable

            else:

                variable_dict = variable

                variable_obj = children[var_id]

                variable_attributes = cls._variable_attributes(variable_dict, variable=variable_obj, children=children)

                variable_obj.update(**variable_attributes)

                container[var_id] = variable_obj

        return container

    @staticmethod
    def _key_container_deserializer(obj, children=None, children_types=None):

        if children is False:
            return {}

        elif children is None:

            if children_types:
                from mewpy.mew.variables import Variable
                return {Variable.from_types(children_types, identifier=variable): val
                        for variable, val in obj.items()}

            else:
                from mewpy.mew.variables import Metabolite
                return {Metabolite(key): val for key, val in obj.items()}

        else:
            return {children[variable]: val for variable, val in obj.items()}

    @classmethod
    def _variable_deserializer(cls: Union[Type['Serializer'], Type['Variable'], Type['Model']],
                               obj):

        identifier = obj.get('id')
        types = obj.get('types')

        variable = cls.from_types(types=types, identifier=identifier)

        variable_attributes = cls._variable_attributes(obj, variable=variable)

        variable.update(**variable_attributes)

        return variable

    @classmethod
    def _variable_attributes(cls, obj, variable, children=None):

        variable_attributes = {}

        # Filtering attributes that must be used by the update method. Child variables associated with this variable
        # are also collected for further building
        for attr_name, (serialize_name, deserialize_name, _) in variable.attributes.items():

            if deserialize_name is not None:
                deserializer, children_types = cls._get_attribute_deserializer(attr=deserialize_name)

                attribute = obj[serialize_name]
                variable_attributes[deserialize_name] = deserializer(attribute,
                                                                     children=children,
                                                                     children_types=children_types)

        return variable_attributes

    @classmethod
    def _get_attribute_deserializer(cls, attr):

        if attr in ('name', 'aliases', 'coefficients', 'bounds', 'charge', 'compartment', 'formula'):

            return cls._obj_deserializer, None

        elif attr == 'regulatory_events':

            return cls._regulatory_events_deserializer, ['regulator']

        elif attr == 'gpr':

            return cls._expression_deserializer, ['gene']

        elif attr == 'reactions':

            return cls._variable_container_deserializer, ['reaction']

        elif attr == 'interactions':

            return cls._variable_container_deserializer, ['interaction']

        elif attr == 'stoichiometry':

            return cls._key_container_deserializer, ['metabolite']

        elif attr == 'target':

            return cls._variable_attribute_deserializer, ['target']

        elif attr == 'interaction':

            return cls._variable_attribute_deserializer, ['interaction']

        else:
            return lambda *args, **kwargs: {}, None

    @staticmethod
    def _expression_deserializer(obj, children=None, children_types=None):

        if children is None:

            symbolic = parse_expression(obj)

            from mewpy.mew.variables import Variable

            variables = {symbol.name: Variable.from_types(children_types, identifier=symbol.name)
                         for symbol in symbolic.atoms(symbols_only=True)}

            return Expression(symbolic=symbolic, variables=variables)

        else:

            symbolic = parse_expression(obj)

            variables = {symbol.name: children[symbol.name]
                         for symbol in symbolic.atoms(symbols_only=True)}

            return Expression(symbolic=symbolic, variables=variables)

    @staticmethod
    def _regulatory_events_deserializer(obj, children=None, children_types=None):

        return {state: Serializer._expression_deserializer(expression, children, children_types)
                for state, expression in obj.items()}

    @staticmethod
    def _variable_container_deserializer(obj, children=None, children_types=None):

        if children is None:

            from mewpy.mew.variables import Variable
            return {key: Variable.from_types(children_types, identifier=key) for key in obj}

        else:
            return {key: children[key] for key in obj}

    @staticmethod
    def _variable_attribute_deserializer(obj, children=None, children_types=None):

        if children is None:

            from mewpy.mew.variables import Variable
            return Variable.from_types(children_types, identifier=obj)

        else:
            return children[obj]

    # -----------------------------------------------------------------------------
    # reduce for pickle serialization
    # -----------------------------------------------------------------------------
    def __reduce__(self: Union['Serializer', 'Model', 'Variable']):

        # for further detail: https://docs.python.org/3/library/pickle.html#object.__reduce__

        from mewpy.model import Model, build_model
        from mewpy.mew.variables import Variable, build_variable

        if isinstance(self, Model):
            return build_model, (tuple(self.types), {'identifier': self.id}), self._dict_to_pickle()

        if isinstance(self, Variable):
            return build_variable, (tuple(self.types), {'identifier': self.id}), self._dict_to_pickle()

        return super(Serializer, self).__reduce__()
    # -----------------------------------------------------------------------------
    # State for pickle serialization
    # -----------------------------------------------------------------------------

    def __getstate__(self):
        return self._dict_to_pickle()

    def __setstate__(self, state):
        self.__dict__.update(state)

    # -----------------------------------------------------------------------------
    # Pickle-like serialization
    # -----------------------------------------------------------------------------

    def _pickle_variable_serializer(self: Union['Serializer', 'Variable'],
                                    to_state=True):

        attributes = {}

        for attribute, (_, _, pickle_name) in self.attributes.items():

            if pickle_name is not None:
                attribute_obj = getattr(self, pickle_name)

                if to_state:
                    attributes[pickle_name] = attribute_obj
                else:
                    attributes[attribute] = attribute_obj

        return attributes

    def _pickle_model_serializer(self: Union['Serializer', 'Model'],
                                 to_state=True):

        containers = {}

        for container, (_, _, pickle_name) in self.containers.items():

            if pickle_name is not None:
                container_obj = getattr(self, pickle_name)

                if to_state:
                    containers[pickle_name] = container_obj
                else:
                    containers[container] = container_obj

        return containers

    def _dict_to_pickle(self: Union['Serializer', 'Model', 'Variable']):

        if hasattr(self, 'containers'):
            return self._pickle_model_serializer(to_state=True)

        if hasattr(self, 'attributes'):
            return self._pickle_variable_serializer(to_state=True)

        return {}

    @classmethod
    def _pickle_variable_deserializer(cls: Union[Type['Variable']],
                                      obj):

        identifier = obj.get('id')
        types = obj.get('types')

        variable = cls.from_types(types=types, identifier=identifier)

        new_dict = {}

        for attribute, (_, _, pickle_name) in variable.attributes.items():

            if pickle_name is not None:
                new_dict[pickle_name] = obj.get(attribute)

        variable.__dict__.update(new_dict)

        return variable

    @classmethod
    def _pickle_model_deserializer(cls: Union[Type['Model']],
                                   obj):

        identifier = obj.get('id')
        types = obj.get('types')

        model = cls.from_types(types=types, identifier=identifier)

        new_dict = {}

        for container, (_, _, pickle_name) in model.containers.items():

            if pickle_name is not None:
                new_dict[pickle_name] = obj.get(container)

        model.__dict__.update(new_dict)

        return model

    # FIXME: make sure variables point to the correct model
    def to_dict(self: Union['Serializer', 'Variable', 'Model'],
                serialization_format: str = 'json',
                variables: bool = False) -> Dict[str, Union[dict,
                                                            'Gene',
                                                            'Interaction',
                                                            'Metabolite',
                                                            'Reaction',
                                                            'Regulator',
                                                            'Target']]:

        is_model = False
        if hasattr(self, 'containers'):
            is_model = True

        if serialization_format == 'pickle':

            if is_model:
                dict_obj = self._pickle_model_serializer(to_state=False)
            else:
                dict_obj = self._pickle_variable_serializer(to_state=False)

            dict_obj['types'] = self.types
            return dict_obj

        elif serialization_format == 'json':
            if is_model:
                return self._model_serializer(variables)

            return self._variable_serializer()

        return {}

    @classmethod
    def from_dict(cls: Union[Type['Serializer'], Type['Variable'], Type['Model']],
                  obj: Dict[str, Union[dict,
                                       'Gene',
                                       'Interaction',
                                       'Metabolite',
                                       'Reaction',
                                       'Regulator',
                                       'Target']],
                  serialization_format: str = 'json',
                  variables: bool = False) -> Union['Gene',
                                                    'Interaction',
                                                    'Metabolite',
                                                    'Reaction',
                                                    'Regulator',
                                                    'Target',
                                                    'MetabolicModel',
                                                    'RegulatoryModel']:

        is_model = False
        if hasattr(cls, 'containers'):
            is_model = True

        if serialization_format == 'pickle':

            if is_model:
                return cls._pickle_model_deserializer(obj)

            return cls._pickle_variable_deserializer(obj)

        elif serialization_format == 'json':

            if is_model:
                return cls._model_deserializer(obj, variables)

            return cls._variable_deserializer(obj)

        return

    # -----------------------------------------------------------------------------
    # Copy leveraging serialization
    # -----------------------------------------------------------------------------

    def __copy__(self) -> Union['Gene',
                                'Interaction',
                                'Metabolite',
                                'Reaction',
                                'Regulator',
                                'Target',
                                'MetabolicModel',
                                'RegulatoryModel']:

        obj_dict = self.to_dict(serialization_format='pickle')

        return self.from_dict(obj_dict, serialization_format='pickle')

    def copy(self) -> Union['Gene',
                            'Interaction',
                            'Metabolite',
                            'Reaction',
                            'Regulator',
                            'Target',
                            'MetabolicModel',
                            'RegulatoryModel']:

        return self.__copy__()

    def __deepcopy__(self) -> Union['Gene',
                                    'Interaction',
                                    'Metabolite',
                                    'Reaction',
                                    'Regulator',
                                    'Target',
                                    'MetabolicModel',
                                    'RegulatoryModel']:

        obj_dict = self.to_dict(variables=True)

        return self.from_dict(obj_dict, variables=True)

    def deepcopy(self) -> Union['Gene',
                                'Interaction',
                                'Metabolite',
                                'Reaction',
                                'Regulator',
                                'Target',
                                'MetabolicModel',
                                'RegulatoryModel']:

        return self.__deepcopy__()
