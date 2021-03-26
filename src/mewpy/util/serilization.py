from typing import Union, TYPE_CHECKING, Type

from mewpy.algebra import Expression
from mewpy.algebra import parse_expression

if TYPE_CHECKING:
    from mewpy.model import Model
    from mewpy.variables import Variable


def serialize(serialization_name, deserialization_name):
    def wrapper(attr):
        attr.fget.serialize = serialization_name
        attr.fget.deserialize = deserialization_name

        return attr

    return wrapper


class Serializer:

    @staticmethod
    def _regulatory_events_serializer(expressions):

        return {key: expression.to_string() for key, expression in expressions.items()}

    @staticmethod
    def _expression_serializer(expression):

        return str(expression)

    @staticmethod
    def _python_native_obj_serializer(obj):

        if isinstance(obj, set):
            return tuple(obj)
        return obj

    @staticmethod
    def _variable_obj_serializer(variable):

        if hasattr(variable, 'id'):
            return variable.id

        return variable

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

        # noinspection PyProtectedMember
        return {variable.id: variable._variable_serializer() for variable in container.values()}

    def _get_model_serializer(self, attr, serialize_variables=True):

        if attr == 'id':

            return self._python_native_obj_serializer

        elif attr == 'name':

            return self._python_native_obj_serializer

        if attr == 'types':

            return self._python_native_obj_serializer

        elif attr == 'objective':

            return self._key_container_serializer

        elif attr in ('genes', 'metabolites', 'reactions', 'interactions', 'regulators', 'targets'):

            if serialize_variables:
                return self._model_container_serializer

            return dict

        else:
            return lambda *args, **kwargs: {}

    def _get_variable_serializer(self, attr):

        if attr == 'id':

            return self._python_native_obj_serializer

        elif attr == 'name':

            return self._python_native_obj_serializer

        elif attr == 'types':

            return self._python_native_obj_serializer

        elif attr == 'aliases':

            return self._python_native_obj_serializer

        elif attr == 'coefficient':

            return self._coefficient_serializer

        elif attr == 'reactions':

            return self._variable_container_serializer

        elif attr == 'regulatory_events':

            return self._regulatory_events_serializer

        elif attr == 'regulators':

            return self._variable_container_serializer

        elif attr == 'target':

            return self._variable_obj_serializer

        elif attr == 'charge':

            return self._python_native_obj_serializer

        elif attr == 'compartment':

            return self._python_native_obj_serializer

        elif attr == 'formula':

            return self._python_native_obj_serializer

        elif attr == 'bounds':

            return self._python_native_obj_serializer

        elif attr == 'compartments':

            return self._python_native_obj_serializer

        elif attr == 'equation':

            return self._python_native_obj_serializer

        elif attr == 'genes':

            return self._variable_container_serializer

        elif attr == 'gpr':

            return self._expression_serializer

        elif attr == 'stoichiometry':

            return self._key_container_serializer

        elif attr == 'metabolites':

            return self._variable_container_serializer

        elif attr == 'boundary':

            return self._python_native_obj_serializer

        elif attr == 'reversibility':

            return self._python_native_obj_serializer

        elif attr == 'interactions':

            return self._variable_container_serializer

        elif attr == 'targets':

            return self._variable_container_serializer

        elif attr == 'interaction':

            return self._variable_obj_serializer

        else:
            return lambda *args, **kwargs: {}

    def _variable_serializer(self: Union['Serializer', 'Variable', 'Model']):

        variable = {}

        for _, (attr, _) in self.attributes.items():
            serializer = self._get_variable_serializer(attr)

            attribute = getattr(self, attr)

            variable[attr] = serializer(attribute)

        return variable

    def _model_serializer(self: Union['Serializer', 'Variable', 'Model'], variables=True):

        model = {}

        for _, (attr, _) in self.containers.items():
            model[attr] = {}

            serializer = self._get_model_serializer(attr, variables)

            container = getattr(self, attr, {})

            model[attr] = serializer(container)

        return model

    def to_dict(self: Union['Serializer', 'Variable', 'Model'], variables=True):

        if hasattr(self, 'containers'):
            return self._model_serializer(variables)

        if hasattr(self, 'attributes'):
            return self._variable_serializer()

        return {}

    @staticmethod
    def _get_attribute_variables(name, attribute):

        if attribute is None:
            return [], []

        elif name == 'reactions':

            return attribute, ['reaction'] * len(attribute)

        elif name == 'regulators':

            return attribute, ['regulator'] * len(attribute)

        elif name == 'target':

            return (attribute,), ['target'] * 1

        elif name == 'genes':

            return attribute, ['gene'] * len(attribute)

        elif name == 'metabolites':

            return attribute, ['metabolite'] * len(attribute)

        elif name == 'interactions':

            return attribute, ['interaction'] * len(attribute)

        elif name == 'targets':

            return attribute, ['target'] * len(attribute)

        elif name == 'interaction':

            return (attribute,), ['interaction'] * 1

        else:
            return [], []

    @staticmethod
    def _deserialize_stoichiometry(attribute, variables):

        return {variables.get(variable): val for variable, val in attribute.items()}

    @staticmethod
    def _deserialize_variable_container(attribute, variables):

        return {variables.get(variable).id: variables.get(variable) for variable in attribute}

    @staticmethod
    def _deserialize_variable(attribute, variables):

        return variables.get(attribute)

    @staticmethod
    def _deserialize_expression(attribute, variables):

        symbolic = parse_expression(attribute)

        variables = {symbol.name: variables[symbol.name]
                     for symbol in symbolic.atoms(symbols_only=True)}

        return Expression(symbolic=symbolic, variables=variables)

    @staticmethod
    def _deserialize_regulatory_events(attribute, variables):

        return {state: Serializer._deserialize_expression(expression, variables)
                for state, expression in attribute.items()}

    @classmethod
    def _deserialize_attribute(cls, name, attribute, variables):

        if attribute is None:
            return None

        elif name == 'gpr':

            return cls._deserialize_expression(attribute, variables)

        elif name == 'regulatory_events':

            return cls._deserialize_regulatory_events(attribute, variables)

        elif name == 'stoichiometry':

            return cls._deserialize_stoichiometry(attribute, variables)

        elif name == 'reactions':

            return cls._deserialize_variable_container(attribute, variables)

        elif name == 'regulators':

            return cls._deserialize_variable_container(attribute, variables)

        elif name == 'target':

            return cls._deserialize_variable(attribute, variables)

        elif name == 'genes':

            return cls._deserialize_variable_container(attribute, variables)

        elif name == 'metabolites':

            return cls._deserialize_variable_container(attribute, variables)

        elif name == 'interactions':

            return cls._deserialize_variable_container(attribute, variables)

        elif name == 'targets':

            return cls._deserialize_variable_container(attribute, variables)

        elif name == 'interaction':

            return cls._deserialize_variable(attribute, variables)

        else:
            return attribute

    @classmethod
    def _variable_deserializer(cls: Union[Type['Serializer'], Type['Variable'], Type['Model']],
                               obj):

        identifier = obj.pop('id')
        types = obj.pop('types')

        variable = cls.from_types(types=types, identifier=identifier)
        variable_attributes = {}
        children = {variable.id: set(types)}

        # Filtering attributes that must be used by the update method. Child variables associated with this variable
        # are also collected for further building
        for attr_name, (serialize_name, deserialize_name) in variable.attributes.items():

            attribute = obj.get(serialize_name)

            if deserialize_name is not None:

                variable_attributes[deserialize_name] = attribute

            # Identifying all types for the child variables. If it is not a variable container/attribute,
            # empty lists are returned
            child_variables, child_variables_types = cls._get_attribute_variables(attr_name, attribute)

            for child_variable, child_variable_type in zip(child_variables, child_variables_types):
                previous_types = children.get(child_variable, set())

                previous_types.add(child_variable_type)

                children[child_variable] = previous_types

        children = {_id: cls.from_types(types=_types, identifier=_id)
                    for _id, _types in children.items()}

        children[variable.id] = variable

        transformed_attributes = {}

        for attr, value in variable_attributes.items():
            attribute = cls._deserialize_attribute(name=attr, attribute=value, variables=children)
            transformed_attributes[attr] = attribute

        variable.update(**transformed_attributes)

        return variable

    @staticmethod
    def _get_variables_container(name, attribute):

        if attribute is None:

            return {}

        elif name in ('genes', 'metabolites', 'reactions', 'interactions', 'regulators', 'targets'):

            return attribute

        else:
            return {}

    @classmethod
    def _model_deserializer(cls: Union[Type['Serializer'], Type['Variable'], Type['Model']],
                            obj):

        from mewpy.variables import Variable

        # id and types keys from serialization
        identifier = obj.pop('id')
        types = obj.pop('types')

        model = cls.from_types(types=types, identifier=identifier)

        model_attributes = {}
        variables_containers = {}
        variables = {}

        # Filtering containers that must be used by the update method. Child variables associated with this variable
        # are also collected for further building
        # noinspection PyProtectedMember
        for attr_name, (serialize_name, deserialize_name) in model.containers.items():

            if deserialize_name is None:
                continue

            # some attributes can be set to None, so using _skip_ flag
            attribute = obj.get(serialize_name)

            # Identifying if it is a container. If so, all variables will be build. Otherwise, the attribute is set
            # to the model attributes
            variables_container = cls._get_variables_container(attr_name, attribute)

            if variables_container:

                for variable_id, variable in variables_container.items():

                    if variable_id not in variables:

                        # sometimes full serialization is not performed. And the model dict object only has
                        # containers with variables already built
                        if not isinstance(variable, Variable):

                            # otherwise, build the variables
                            variable = Variable.from_types(variable.get('types'), identifier=variable_id)

                        variables[variable_id] = variable
                        variables_containers[variable_id] = attr_name

            else:
                model_attributes[deserialize_name] = attribute

        for variable_id, variable in variables.items():

            variable_container = variables_containers[variable_id]

            variable_dict_attributes = obj[variable_container][variable_id]

            if not isinstance(variable_dict_attributes, dict):
                variable_dict_attributes = variable_dict_attributes.to_dict()

            variable_attributes = {}

            # noinspection PyProtectedMember
            for attr_name, (serialize_name, deserialize_name) in variable.attributes.items():

                if deserialize_name is None:
                    continue

                attribute = variable_dict_attributes.get(serialize_name)

                attribute = cls._deserialize_attribute(name=attr_name, attribute=attribute, variables=variables)

                variable_attributes[deserialize_name] = attribute

            variable.update(**variable_attributes)

        model_attributes['variables'] = list(variables.values())

        model.update(**model_attributes)

        return model

    @classmethod
    def from_dict(cls: Union[Type['Serializer'], Type['Variable'], Type['Model']], obj):

        if hasattr(cls, 'containers'):
            return cls._model_deserializer(obj)

        elif hasattr(cls, 'attributes'):
            return cls._variable_deserializer(obj)

        else:
            return
