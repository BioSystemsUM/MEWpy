from warnings import warn
from mewpy.utils.parsing import generic_regulatory_rule_parser
from mewpy.regulation.regulatory_variable import RegulatoryVariable


class RegulatoryInteraction(object):

    def __init__(self, identifier,
                 rule,
                 target,
                 regulators=None,
                 model=None,
                 aliases=None):

        if not aliases:
            aliases = []

        self._id = None
        self._aliases = None
        self.old_rule = None
        self._rule = None
        self._elements = None
        self._target = None
        self._regulators = None
        self._model = None

        self.__build_regulatory_interaction__(identifier,
                                              aliases,
                                              rule,
                                              target,
                                              regulators,
                                              model)

    def __str__(self):
        return self.target.id + ' ' + self.rule

    def __repr__(self):

        return self.id

    # def __repr__(self):
    #
    #     return '{}({}, {}, {}) at {}'.format(self.__class__.__name__, self._id, self._regulatory_interaction_rule,
    #     self._target, id(self))

    def __delete__(self, instance):
        del self.id
        del self.rule
        del self.elements
        del self.old_rule
        del self.target
        del self.regulators
        del self.model

    def __build_regulatory_interaction__(self,
                                         identifier,
                                         aliases,
                                         rule,
                                         target,
                                         regulators,
                                         model):

        if not identifier:
            raise BaseException("An ID must be provided")

        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        if not rule:
            def reg_rule_msg(message):
                warn(message, UserWarning)

            reg_rule_msg(
                "The regulatory interaction rule is empty for the target {}".format(target.id))

            rule = ''

        if regulators:

            if not isinstance(regulators, dict) or \
                    not isinstance(regulators[list(regulators.keys())[0]], RegulatoryVariable):
                raise TypeError(
                    "Regulators must be a dict of RegulatoryVariables objects")

            self._regulators = regulators

        else:
            self._regulators = {}

        self.id = identifier
        self.aliases = aliases

        self.target = target

        self.model = model

        self.rule = rule

    @classmethod
    def from_str(cls, identifier,
                 rule,
                 target,
                 target_aliases=None,
                 regulators=None,
                 model=None):

        if not isinstance(identifier, str):
            raise TypeError("IDs must be a list of strings")

        if not isinstance(rule, str):
            raise TypeError(
                "Rules must be a dict object of list of strings. Use list() to convert")

        if isinstance(target, str):

            if model is not None:

                tar_var = model.__add_target__(identifier=target,
                                               name=target,
                                               aliases=target_aliases,
                                               regulators=regulators)

            else:

                tar_var = RegulatoryVariable(identifier=target,
                                             name=target,
                                             aliases=target_aliases,
                                             is_target=True)

            return cls(identifier=identifier,
                       rule=rule,
                       target=tar_var,
                       model=model)

        elif isinstance(target, RegulatoryVariable):

            return cls(identifier=identifier,
                       rule=rule,
                       target=target,
                       model=model)

        else:
            raise TypeError(
                "Target must be a string or a RegulatoryVariable object")

    @classmethod
    def from_dict(cls, identifiers,
                  rules,
                  targets_aliases=None,
                  model=None):

        if not isinstance(identifiers, list):
            raise TypeError("IDs must be a list of strings")

        if not isinstance(identifiers[0], str):
            raise TypeError("IDs must be a list of strings")

        if not isinstance(rules, dict):
            raise TypeError(
                "Rules must be a dict object of list of strings. Use dict() to convert")

        if not isinstance(rules[list(rules.keys())[0]], str):
            raise TypeError(
                "Rules must be a dict object of list of strings. Use list() to convert")

        if not targets_aliases:

            if not isinstance(targets_aliases, dict):
                raise TypeError(
                    "The target_aliases must be a dict object. Use dict() to convert")

            if not isinstance(targets_aliases[list(targets_aliases.keys())[0]], list):
                raise TypeError(
                    "The target_aliases must be a dict object of lists of strings. Use list() to convert")

            if not isinstance(targets_aliases[list(targets_aliases.keys())[0]][0], str):
                raise TypeError(
                    "The target_aliases must be a dict object of lists of strings. Use str() to convert")

        else:
            targets_aliases = {}

        reg_inters = []

        for i, (key, val) in enumerate(rules.items()):

            if key in targets_aliases:
                aliases = targets_aliases[key]
            else:
                aliases = []

            reg_inter = cls.from_str(identifiers[i],
                                     rule=val,
                                     target=key,
                                     target_aliases=aliases,
                                     regulators=None,
                                     model=model)

            reg_inters.append(reg_inter)

        return reg_inters

    @property
    def id(self):
        return getattr(self, '_id', None)

    @property
    def aliases(self):
        return getattr(self, '_aliases', None)

    @property
    def rule(self):
        return getattr(self, '_rule', None)

    @property
    def elements(self):
        return getattr(self, '_elements', None)

    @property
    def target(self):
        return getattr(self, '_target', None)

    @property
    def regulators(self):
        return getattr(self, '_regulators', None)

    @property
    def model(self):
        return getattr(self, '_model', None)

    @id.setter
    def id(self, value):

        if value == self._id:
            pass

        elif not isinstance(value, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        else:

            if self.model:

                if value in self.model.regulatory_interactions:
                    def reg_id(message):
                        warn(message, UserWarning, stacklevel=2)

                    reg_id("The ID already exists in the dictionary of regulatory interactions of the regulatory "
                           "model. Ignoring the alteration")

                else:

                    if self._id in self.model.regulatory_interactions:
                        del self.model.regulatory_interactions[self._id]
                        self.model.regulatory_interactions[value] = self

                    if self.regulators:
                        for reg in self.regulators_gen():
                            if self._id in reg.regulatory_interactions:
                                del reg.regulatory_interactions[self._id]
                                reg.regulatory_interactions[value] = self

                    if self.target:
                        if self._id in self.target.regulatory_interactions:
                            del self.target.regulatory_interactions[self._id]
                            self.target.regulatory_interactions.regulatory_interactions[value] = self

                    self._id = value

            else:

                if self.regulators:
                    for reg in self.regulators_gen():
                        if self._id in reg.regulatory_interactions:
                            del reg.regulatory_interactions[self._id]
                            reg.regulatory_interactions[value] = self

                if self.target:
                    if self._id in self.target.regulatory_interactions:
                        del self.target.regulatory_interactions[self._id]
                        self.target.regulatory_interactions.regulatory_interactions[value] = self

                self._id = value

    @aliases.setter
    def aliases(self, value):

        if value == self._aliases:
            pass

        elif not value:
            pass

        else:

            if not isinstance(value, list):
                raise TypeError("The aliases must be a list. Use list() to convert")

            if value and not isinstance(value[0], str):
                raise TypeError("The aliases must be a list. Use list() to convert")

            if not self._aliases:
                self._aliases = []

            for val in value:
                if val not in self._aliases:
                    self._aliases = self._aliases + [val]

    @rule.setter
    def rule(self, value):

        if value == self._rule:
            pass

        elif not isinstance(value, str):
            raise TypeError(
                "Regulatory interaction rule must be a string. Use str() to convert")

        else:

            parser = generic_regulatory_rule_parser(value)

            # return parsed_str (str), 1
            # exp_elements (list), 2

            self.old_rule = value

            self._rule = parser[0]
            self._elements = parser[1]

            if self.regulators:

                for key in self.regulators:
                    self.regulators[key].regulatory_interactions[self.id] = self
                    self.regulators[key].is_regulator = True
                    self.regulators[key].model = self.model
                    self.regulators[key].targets[self.target.id] = self.target
                    self.target.regulators[key] = self.regulators[key]

    @elements.setter
    def elements(self, value):

        raise BaseException("regulatory_interaction_elements property cannot be changed. "
                            "please create/set a new regulatory interaction rule.")

    @target.setter
    def target(self, value):

        if value == self._target:
            pass

        elif not value:
            raise BaseException("A target must provided as a RegulatoryVariable object")

        elif not isinstance(value, RegulatoryVariable):
            raise TypeError("Target must be a RegulatoryVariable object")

        else:

            self._target = value

            self._target.is_target = True
            self._target.regulatory_interactions[self.id] = self
            self._target.model = self.model

    @regulators.setter
    def regulators(self, value):

        raise BaseException("Regulators property cannot be changed. "
                            "To change the regulators of this regulatory interacion, "
                            "please change the regulatory interaction rule.")

    @model.setter
    def model(self, value):

        # TODO: The setter should also remove the regulatory variable from the old model if exists

        from mewpy.regulation.regulatory_model import RegulatoryModel

        if value == self._model:
            pass

        elif not value:
            self._model = None

        elif not isinstance(value, RegulatoryModel):
            raise TypeError("Regulatory model must be a Regmodel object or None")

        else:

            # TODO: This should point all variables and interactions related with this one to the same model, but
            #  this might be a recursion problem

            if self.id not in value.regulatory_interactions:
                value.regulatory_interactions[self.id] = self
                self._model = value

                if self.target:
                    if self.target.id not in value.targets:
                        value.targets[self.target.id] = self.target
                    if self.target.id in value.regulators:
                        self.target.is_regulator = True

                if self.regulators:
                    for reg in self.regulators_gen():
                        if reg.id not in value.regulators:
                            value.regulators[reg.id] = reg
                        if reg.id in value.targets:
                            reg.is_target = True

            else:
                def model_id(message):
                    warn(message, UserWarning, stacklevel=2)

                model_id("The ID already exists in the dictionary of regulatory interactions of the regulatory "
                         "model. Ignoring the alteration")

    @id.deleter
    def id(self):
        del self._id

    @rule.deleter
    def rule(self):
        del self._rule

    @elements.deleter
    def elements(self):
        del self._elements

    @target.deleter
    def target(self):
        del self._target

    @regulators.deleter
    def regulators(self):
        del self._regulators

    @model.deleter
    def model(self):
        del self._model

    def get_target(self):

        return self.target

    def get_regulator(self, identifier):

        return self.regulators[identifier]

    def regulators_gen(self):

        return (value for value in self._regulators.values())

    def remove_from_model(self):

        if self.model is not None:

            self.model.remove_regulatory_interaction(self.id)

        else:
            def remove_from_model_reg(message):
                warn(message, UserWarning, stacklevel=2)

            remove_from_model_reg("The regulatory interaction is not bounded to any regulatory model")
