from warnings import warn
from mewpy.io import read_tabular_aliases
from mewpy.regulation.regulatory_interaction import RegulatoryInteraction
from mewpy.regulation.regulatory_variable import RegulatoryVariable
from sympy.logic.boolalg import BooleanTrue, BooleanFalse
from sympy.core.numbers import Zero, One


class RegulatoryModel(object):

    def __init__(self, identifier,
                 name=None,
                 targets=None,
                 regulators=None,
                 regulatory_interactions=None):

        if not identifier:
            raise BaseException("An ID must be provided")

        if not name:
            name = identifier

        if not targets:
            targets = {}

        if not regulators:
            regulators = {}

        if not regulatory_interactions:
            regulatory_interactions = {}

        self._id = None

        self._name = None

        self._targets = None

        self._regulators = None

        # TODO: add regulator and target should change the property state
        self._regulatory_variables = {}

        self._regulatory_interactions = None

        self._regulatory_solution = []

        self.__build_regulatory_model__(identifier=identifier,
                                        name=name,
                                        targets=targets,
                                        regulators=regulators,
                                        regulatory_interactions=regulatory_interactions)

    def __str__(self):
        return self._name

    def __repr__(self):

        return '{}({}, {}) at {}'.format(self.__class__.__name__, self._id, self._name, id(self))

    def __build_regulatory_model__(self, identifier,
                                   name,
                                   targets,
                                   regulators,
                                   regulatory_interactions):

        if targets:

            if not isinstance(targets, dict):
                raise TypeError(
                    "Targets must be a dict of RegulatoryVariables objects")

            if not isinstance(targets[list(targets.keys())[0]], RegulatoryVariable):
                raise TypeError(
                    "Targets must be a dict of RegulatoryVariables objects")

            self._targets = targets

        else:
            self._targets = {}

        if regulators:

            if not isinstance(regulators, dict):
                raise TypeError("Regulators must be a dict of RegulatoryVariables objects")

            if not isinstance(regulators[list(regulators.keys())[0]], RegulatoryVariable):
                raise TypeError("Regulators must be a dict of RegulatoryVariables objects")

            self._regulators = regulators

        else:
            self._regulators = {}

        if regulatory_interactions:

            if not isinstance(regulatory_interactions, dict):
                raise TypeError(
                    "Regulatory interactions must be a dict of RegulatoryInteraction objects")

            if not isinstance(regulatory_interactions[list(regulatory_interactions.keys())[0]], RegulatoryInteraction):
                raise TypeError("Regulatory interactions must be a dict of RegulatoryInteraction objects. "
                                "Use RegulatoryInteraction() to convert")

            self._regulatory_interactions = regulatory_interactions

        else:
            self._regulatory_interactions = {}

        self.id = identifier
        self.name = name

    @classmethod
    def from_dict(cls, identifier,
                  name=None,
                  targets=None,
                  regulators=None,
                  regulatory_interactions=None):
        # simulation_interface=None):

        return cls(identifier, name, targets, regulators, regulatory_interactions)  # , simulation_interface)

    @property
    def id(self):
        return getattr(self, '_id', None)

    @property
    def name(self):
        return getattr(self, '_name', None)

    @property
    def targets(self):
        return getattr(self, '_targets', None)

    @property
    def regulators(self):
        return getattr(self, '_regulators', None)

    @property
    def regulatory_variables(self):
        return getattr(self, '_regulatory_variables', None)

    @property
    def regulatory_interactions(self):
        return getattr(self, '_regulatory_interactions', None)

    @property
    def regulatory_solution(self):
        return getattr(self, '_regulatory_solution', None)

    @id.setter
    def id(self, value):

        if value == self._id:
            pass

        elif not isinstance(value, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        else:
            self._id = value

    @name.setter
    def name(self, value):

        if value == self._name:
            pass

        elif not isinstance(value, str):
            raise TypeError("The name must be a string. Use str() to convert")

        else:
            self._name = value

    @targets.setter
    def targets(self, value):

        if value == self._targets:
            pass

        else:
            raise BaseException("Targets property cannot be changed. "
                                "To add new targets to the model, please use the add_new_target method.")

    @regulators.setter
    def regulators(self, value):

        if value == self._regulators:
            pass

        else:
            raise BaseException("Regulators property cannot be changed. "
                                "To add new regulators to the model, please use the add_new_regulator method.")

    @regulatory_variables.setter
    def regulatory_variables(self, value):

        if value == self._regulatory_variables:
            pass

        else:
            raise BaseException("regulatory variables property cannot be changed. "
                                "To add new targets or regulators to the model, please use the add_new_target method.")

    @regulatory_interactions.setter
    def regulatory_interactions(self, value):

        if value == self._regulatory_interactions:
            pass

        elif not value:
            self._regulatory_interactions = {}

        else:

            if not isinstance(value, dict):
                raise TypeError(
                    "Regulatory interactions must be a dict of RegulatoryInteraction objects")

            if not isinstance(value[list(value.keys())[0]], RegulatoryInteraction):
                raise TypeError("Regulatory interactions must be a dict of RegulatoryInteraction objects. "
                                "Use RegulatoryInteraction() to convert")

            def regulatory_interactions(message):
                warn(message, UserWarning, stacklevel=2)

            regulatory_interactions("Setting/altering the regulatory interactions associated with a model replaces the "
                                    "whole model. Regulators and targets will be ignored and rebuild. "
                                    "To add new regulatory interactions to the stack use the "
                                    "add_new_regulatory_interaction method")

            self._regulatory_interactions = value

            for reg_interaction in value.values():

                for reg_val in reg_interaction.regulators_gen():

                    if reg_val.id not in self._regulators:
                        self._regulators[reg_val.id] = reg_val

                if reg_interaction.target.id not in self._targets:
                    self._targets[reg_interaction.target.id] = reg_interaction.target

    def get_regulator(self, identifier):

        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        return self._regulators[identifier]

    def get_target(self, identifier):

        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        return self._targets[identifier]

    def get_regulatory_variable(self, identifier):

        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        if identifier in self.regulators:
            return self._regulators[identifier]
        elif identifier in self.targets:
            return self._targets[identifier]
        else:
            return None

    def get_regulatory_interaction(self, identifier):

        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        return self._regulatory_interactions[identifier]

    def regulators_gen(self):

        return (value for value in self._regulators.values())

    def targets_gen(self):

        return (value for value in self._targets.values())

    def regulatory_variables_gen(self):

        return (value for value in self._regulatory_variables.values())

    def regulatory_interactions_gen(self):

        return (value for value in self._regulatory_interactions.values())

    def simulate(self):

        raise NotImplementedError("This is the base classe for type-specific regulatory models. Instantiate these "
                                  "instead to perform simulations")

    def add_regulator(self, identifier,
                      name=None,
                      expression_coef=None,
                      aliases=None,
                      regulatory_interactions=None,
                      targets=None,
                      regulators=None):

        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        if identifier in self.regulators:
            def add_regulator_warn(message):
                warn(message, UserWarning, stacklevel=2)

            add_regulator_warn("The ID already exists in the dictionary of regulators of the regulatory "
                               "model. Ignoring the alteration")

            return self.regulators[identifier]

        if identifier in self.targets:
            self.targets[identifier].is_regulator = True

            return self.targets[identifier]

        regulator = RegulatoryVariable(
            identifier=identifier,
            name=name,
            expression_coef=expression_coef,
            aliases=aliases,
            regulatory_interactions=regulatory_interactions,
            targets=targets,
            regulators=regulators,
            model=self,
            is_regulator=True,
            is_target=False)

        self.regulatory_variables[identifier] = regulator

        return regulator

    def add_regulator_object(self, regulator):

        if not isinstance(regulator, RegulatoryVariable):
            raise TypeError("Regulator must be RegulatoryVariable object")

        if regulator.id in self.regulators:
            def add_regulator_warn(message):
                warn(message, UserWarning, stacklevel=2)

            add_regulator_warn("The ID already exists in the dictionary of regulators of the regulatory "
                               "model. Ignoring the alteration")
            return

        if regulator.id in self.targets:
            regulator.is_target = True

        for reg_inter in regulator.regulatory_interactions_gen():
            reg_inter.model = self

        for reg in regulator.regulators_gen():
            reg.model = self

        for tar in regulator.targets_gen():
            tar.model = self

        regulator.model = self

    def __add_regulator__(self, identifier,
                          name=None,
                          expression_coef=None,
                          aliases=None,
                          regulatory_interactions=None,
                          targets=None,
                          regulators=None):

        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        if identifier in self.regulators:

            regulator = self.regulators[identifier]

            if name is not None:
                regulator.name = name
            if expression_coef is not None:
                regulator.expression_coef = expression_coef
            if aliases is not None:
                for alias in aliases:
                    if alias not in regulator.aliases:
                        regulator.aliases.append(alias)
            if regulatory_interactions is not None:
                regulator.regulatory_interactions.update(regulatory_interactions)
            if targets is not None:
                regulator.outputs.update(targets)
            if regulators is not None:
                regulator.inputs.update(regulators)

            return regulator

        if identifier in self.targets:

            regulator = self.targets[identifier]

            if name is not None:
                regulator.name = name
            if expression_coef is not None:
                regulator.expression_coef = expression_coef
            if aliases is not None:
                for alias in aliases:
                    if alias not in regulator.aliases:
                        regulator.aliases.append(alias)
            if regulatory_interactions is not None:
                regulator.regulatory_interactions.update(regulatory_interactions)
            if targets is not None:
                regulator.outputs.update(targets)
            if regulators is not None:
                regulator.inputs.update(regulators)
            regulator.is_regulator = True

            return regulator

        regulator = RegulatoryVariable(
            identifier=identifier,
            name=name,
            expression_coef=expression_coef,
            aliases=aliases,
            regulatory_interactions=regulatory_interactions,
            targets=targets,
            regulators=regulators,
            model=self,
            is_regulator=True,
            is_target=False)

        self.regulatory_variables[identifier] = regulator

        return regulator

    def remove_regulator(self, identifier):

        if identifier not in self.regulators:
            def remove_regulator_id(message):
                warn(message, UserWarning, stacklevel=2)

            remove_regulator_id("The ID does not exist in the dictionary of regulators of the regulatory "
                                "model. Ignoring the alteration")

            return

        for key in self.regulators[identifier].regulatory_interactions.keys():
            self.remove_regulatory_interaction(key)

        if identifier in self.regulators:
            del self.regulators[identifier]

        if identifier in self.targets:
            del self.targets[identifier]

    def add_target(self, identifier,
                   name=None,
                   expression_coef=None,
                   aliases=None,
                   regulatory_interactions=None,
                   targets=None,
                   regulators=None):

        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        if identifier in self.targets:
            def add_target_warn(message):
                warn(message, UserWarning, stacklevel=2)

            add_target_warn("The ID already exists in the dictionary of targets of the regulatory "
                            "model. Ignoring the alteration")

            return self.targets[identifier]

        if identifier in self.regulators:
            self.regulators[identifier].is_target = True

            return self.regulators[identifier]

        target = RegulatoryVariable(
            identifier=identifier,
            name=name,
            expression_coef=expression_coef,
            aliases=aliases,
            regulatory_interactions=regulatory_interactions,
            targets=targets,
            regulators=regulators,
            model=self,
            is_regulator=False,
            is_target=True)

        self.regulatory_variables[identifier] = target

        return target

    def add_target_object(self, target):

        if not isinstance(target, RegulatoryVariable):
            raise TypeError("Target must be RegulatoryVariable object")

        if target.id in self.targets:
            def add_regulator_warn(message):
                warn(message, UserWarning, stacklevel=2)

            add_regulator_warn("The ID already exists in the dictionary of targets of the regulatory "
                               "model. Ignoring the alteration")
            return

        if target.id in self.regulators:
            target.is_regulator = True

        for reg_inter in target.regulatory_interactions_gen():
            reg_inter.model = self

        for reg in target.regulators_gen():
            reg.model = self

        for tar in target.targets_gen():
            tar.model = self

        target.model = self

    def __add_target__(self, identifier,
                       name=None,
                       expression_coef=None,
                       aliases=None,
                       regulatory_interactions=None,
                       targets=None,
                       regulators=None):

        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        if identifier in self.targets:

            target = self.targets[identifier]

            if name is not None:
                target.name = name
            if expression_coef is not None:
                target.expression_coef = expression_coef
            if aliases is not None:
                for alias in aliases:
                    if alias not in target.aliases:
                        target.aliases.append(alias)
            if regulatory_interactions is not None:
                target.regulatory_interactions.update(regulatory_interactions)
            if targets is not None:
                target.outputs.update(targets)
            if regulators is not None:
                target.inputs.update(regulators)

            return target

        if identifier in self.regulators:

            target = self.regulators[identifier]

            if name is not None:
                target.name = name
            if expression_coef is not None:
                target.expression_coef = expression_coef
            if aliases is not None:
                for alias in aliases:
                    if alias not in target.aliases:
                        target.aliases.append(alias)
            if regulatory_interactions is not None:
                target.regulatory_interactions.update(regulatory_interactions)
            if targets is not None:
                target.outputs.update(targets)
            if regulators is not None:
                target.inputs.update(regulators)

            target.is_target = True

            return target

        target = RegulatoryVariable(
            identifier=identifier,
            name=name,
            expression_coef=expression_coef,
            aliases=aliases,
            regulatory_interactions=regulatory_interactions,
            targets=targets,
            regulators=regulators,
            model=self,
            is_regulator=False,
            is_target=True)

        self.regulatory_variables[identifier] = target

        return target

    def remove_target(self, identifier):

        if identifier not in self.targets:
            def remove_target_id(message):
                warn(message, UserWarning, stacklevel=2)

            remove_target_id("The ID does not exist in the dictionary of targets of the regulatory "
                             "model. Ignoring the alteration")

            return

        for key in self.targets[identifier].regulatory_interactions.keys():
            self.remove_regulatory_interaction(key)

        if identifier in self.targets:
            del self.targets[identifier]

        if identifier in self.regulators:
            del self.regulators[identifier]

    def add_regulatory_interaction(self, identifier,
                                   rule,
                                   target,
                                   target_aliases=None,
                                   regulators=None):

        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        if not isinstance(rule, str):
            raise TypeError(
                "The regulatory rule must be a string. Use str() to convert")

        if not target_aliases:
            target_aliases = []

        if not isinstance(target, str) and not isinstance(target, RegulatoryVariable):
            raise TypeError(
                "Target must be a string or a RegulatoryVariable object")

        if identifier in self.regulatory_interactions:
            def add_reg_inter_warn(message):
                warn(message, UserWarning, stacklevel=2)

            add_reg_inter_warn("The ID already exists in the dictionary of regulatory interactions of the regulatory "
                               "model. Ignoring the alteration")
            return

        return RegulatoryInteraction.from_str(identifier=identifier,
                                              rule=rule,
                                              target=target,
                                              target_aliases=target_aliases,
                                              regulators=regulators,
                                              model=self)

    def add_regulatory_interaction_object(self, regulatory_interaction):

        if not isinstance(regulatory_interaction, RegulatoryInteraction):
            raise TypeError(
                "The regulatory interaction must be RegulatoryInteraction object")

        if regulatory_interaction.id in self.regulatory_interactions:
            def add_reg_inter_warn(message):
                warn(message, UserWarning, stacklevel=2)

            add_reg_inter_warn("The ID already exists in the dictionary of regulatory interactions of the regulatory "
                               "model. Ignoring the alteration")
            return

        for reg in regulatory_interaction.regulators_gen():
            reg.model = self

        regulatory_interaction.target.model = self

        regulatory_interaction.model = self

    def remove_regulatory_interaction(self, identifier):

        # TODO: Remove regulator, target and regulatoy_interaction don't follow a destructive approach.
        #  Even though their direct elements (regulatory interactions, targets and regulators) are removed,
        #  the regulatory interactions of these direct elements are nonetheless kept in the regulatory model.
        #  Implement a destructive flag, that destroys any relation with a regulatoy interaction, regulator and target.

        if identifier not in self.regulatory_interactions:
            def remove_reg_inter_id(message):
                warn(message, UserWarning, stacklevel=2)

            remove_reg_inter_id("The ID does not exist in the dictionary of regulatory interactions of the regulatory "
                                "model. Ignoring the alteration")

            return

        target = self.regulatory_interactions[identifier].target

        for tar in target.targets_gen():
            if tar.id in self.targets:
                del self.targets[tar.id]
            if tar.id in self.regulators:
                del self.regulators[tar.id]

        for reg in target.regulators_gen():
            if reg.id in self.regulators:
                del self.regulators[reg.id]
            if reg.id in self.targets:
                del self.targets[reg.id]

        for regulator in self.regulatory_interactions[identifier].regulators_gen():

            for tar in regulator.targets_gen():
                if tar.id in self.targets:
                    del self.targets[tar.id]
                if tar.id in self.regulators:
                    del self.regulators[tar.id]

            for reg in regulator.regulators_gen():
                if reg.id in self.regulators:
                    del self.regulators[reg.id]
                if reg.id in self.targets:
                    del self.targets[reg.id]

        if identifier in self.regulatory_interactions:
            del self.regulatory_interactions[identifier]

    def remove_regulatory_variable(self, identifier):

        if identifier not in self.targets and identifier not in self.regulators:
            def remove_reg_var_id(message):
                warn(message, UserWarning, stacklevel=2)

            remove_reg_var_id("The ID does not exist in the dictionary of targets and regulators of the regulatory "
                              "model. Ignoring the alteration")

        else:

            if identifier in self.targets:

                for key in self.targets[identifier].regulatory_interactions.keys():
                    self.remove_regulatory_interaction(key)

                if identifier in self.targets:
                    del self.targets[identifier]

            if identifier in self.regulators:

                for key in self.regulators[identifier].regulatory_interactions.keys():
                    self.remove_regulatory_interaction(key)

                if identifier in self.regulators:
                    del self.regulators[identifier]

    def update_aliases_from_tabular_format_file(self, aliases_file,
                                                sep=None,
                                                id_col=None,
                                                aliases_cols=None,
                                                header=None):

        if not isinstance(aliases_file, str):
            raise TypeError("File or file path must be a string")

        df = read_tabular_aliases(aliases_file, sep=sep, id_col=id_col, aliases_cols=aliases_cols, header=header)

        aliases_cols = [col for col in df.columns if col != 'ids']

        for row in df.index:

            new_aliases = list(df.loc[row, aliases_cols])

            if row in self.regulatory_variables:
                self.regulatory_variables[row].aliases = new_aliases

    # TODO: Context manager for storing history. Useful for temporary changes such as KO

    # TODO: the dictionaries used in the regulatory model and variables must have some of the methods re-written,
    #  namely the set attribute


def solution_decode(solution):

    decode = {
        True: 1,
        False: 0,
        1: 1,
        0: 0,
        Zero: 0,
        One: 1,
        BooleanFalse: 0,
        BooleanTrue: 1,
        None: None
    }

    if solution in decode:
        return decode[solution]

    else:
        return solution
