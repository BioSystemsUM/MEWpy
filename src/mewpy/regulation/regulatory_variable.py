from warnings import warn


class RegulatoryVariable(object):

    def __init__(self, identifier,
                 name=None,
                 expression_coef=None,
                 aliases=None,
                 regulatory_interactions=None,
                 targets=None,
                 regulators=None,
                 model=None,
                 is_regulator=False,
                 is_target=False):

        if not name:
            name = ''

        if not aliases:
            aliases = []

        if not regulatory_interactions:
            regulatory_interactions = {}

        if not targets:
            targets = {}

        if not regulators:
            regulators = {}

        self._id = None
        self._name = None
        self._expression_coef = None
        self._aliases = None
        self._regulatory_interactions = None
        self._targets = None
        self._regulators = None
        self._model = None
        self._is_regulator = None
        self._is_target = None
        self._solution = None

        self.__build_regulatory_variable__(identifier,
                                           name,
                                           expression_coef,
                                           aliases,
                                           regulatory_interactions,
                                           targets,
                                           regulators,
                                           model,
                                           is_regulator,
                                           is_target)

    def __str__(self):
        return self.name

    def __repr__(self):

        return self.id

    # def __repr__(self):
    #
    #     return '{}({}, {}, {}) at {}'.format(self.__class__.__name__, self._id, self._name, self._aliases, id(self))

    def __delete__(self, instance):
        del self.id
        del self.name
        del self.expression_coef
        del self.aliases
        del self.regulatory_interactions
        del self.targets
        del self.regulators
        del self.model
        del self.is_regulator
        del self.is_target
        del self.solution

    def __build_regulatory_variable__(self,
                                      identifier,
                                      name,
                                      expression_coef,
                                      aliases,
                                      regulatory_interactions,
                                      targets,
                                      regulators,
                                      model,
                                      is_regulator,
                                      is_target):

        if not identifier:
            raise TypeError("An ID must be provided")

        if not isinstance(regulatory_interactions, dict):
            raise TypeError(
                "regulatory interactions must be a dict of RegulatoryInteraction objects")

        from mewpy.regulation.regulatory_interaction import RegulatoryInteraction

        if regulatory_interactions and not isinstance(regulatory_interactions[list(regulatory_interactions.keys())[0]],
                                                      RegulatoryInteraction):
            raise TypeError(
                "regulatory interactions must be a dict of RegulatoryInteraction objects")

        if not isinstance(regulators, dict):
            raise TypeError(
                "Regulators must be a dict of RegulatoryVariables objects")

        if regulators and not isinstance(regulators[list(regulators.keys())[0]], RegulatoryVariable):
            raise TypeError(
                "Regulators must be a dict of RegulatoryVariables objects")

        if not isinstance(targets, dict):
            raise TypeError(
                "Targets must be a dict of RegulatoryVariables objects")

        if targets and not isinstance(targets[list(targets.keys())[0]], RegulatoryVariable):
            raise TypeError(
                "Targets must be a dict of RegulatoryVariables objects")

        self.id = identifier
        self.aliases = aliases
        self.name = name
        self.expression_coef = expression_coef

        self._regulatory_interactions = regulatory_interactions
        self._targets = targets
        self._regulators = regulators

        self.is_regulator = is_regulator
        self.is_target = is_target

        self.model = model

    @property
    def id(self):
        return getattr(self, '_id', None)

    @property
    def name(self):
        return getattr(self, '_name', None)

    @property
    def expression_coef(self):
        return getattr(self, '_expression_coef', None)

    @property
    def aliases(self):
        return getattr(self, '_aliases', None)

    @property
    def regulatory_interactions(self):
        return getattr(self, '_regulatory_interactions', None)

    @property
    def targets(self):
        return getattr(self, '_targets', None)

    @property
    def regulators(self):
        return getattr(self, '_regulators', None)

    @property
    def conditions(self):
        return getattr(self, '_conditions', None)

    @property
    def model(self):
        return getattr(self, '_model', None)

    @property
    def is_regulator(self):
        return getattr(self, '_is_regulator', None)

    @property
    def is_target(self):
        return getattr(self, '_is_target', None)

    @property
    def solution(self):
        return getattr(self, '_solution', None)

    @id.setter
    def id(self, value):

        if value == self._id:
            pass

        elif not isinstance(value, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        else:

            if self.model:

                if value in self.model.regulators and value in self.model.targets:
                    def reg_id(message):
                        warn(message, UserWarning, stacklevel=2)

                    reg_id("The ID already exists in the dictionary of regulators and targets of the regulatory "
                           "model. Ignoring the alteration")

                else:

                    if self._id in self.model.regulators:
                        del self.model.regulators[self._id]
                        self.model.regulators[value] = self

                    if self._id in self.model.targets:
                        del self.model.targets[self._id]
                        self.model.targets[value] = self

                    if self.regulatory_interactions:
                        for reg_inter in self.regulatory_interactions_gen():

                            if self._id in reg_inter.regulators:
                                del reg_inter.regulators[self._id]
                                reg_inter.regulators[value] = self

                    if self.regulators:
                        for reg in self.regulators.values():
                            del reg.targets[self._id]
                            reg.targets[value] = self

                    if self.targets:
                        for tar in self.targets.values():
                            del tar.regulators[self._id]
                            tar.regulators[value] = self

                    self._id = value

            else:

                if self.regulatory_interactions:

                    for reg_inter in self.regulatory_interactions_gen():

                        if self._id in reg_inter.regulators:
                            del reg_inter.regulators[self._id]
                            reg_inter.regulators[value] = self

                self._id = value

    @name.setter
    def name(self, value):

        if value == self._name:
            pass

        elif not isinstance(value, str):
            raise TypeError("The name must be a string. Use str() to convert")

        else:
            self._name = value

    @expression_coef.setter
    def expression_coef(self, value):

        if value == self._expression_coef:
            pass

        elif not isinstance(value, int) and not isinstance(value, float) and value is not None:
            raise TypeError("The expression coefficient must be a number (int or float) or None. "
                            "Use float() to convert")

        else:
            self._expression_coef = value

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

    @regulatory_interactions.setter
    def regulatory_interactions(self, value):

        if value == self._regulatory_interactions:
            pass

        else:
            raise BaseException("Regulatory interactions property cannot be changed. "
                                "To add new regulatory interactions to the model, "
                                "please use the add_new_regulatory_interaction method.")

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

    @model.setter
    def model(self, value):

        # TODO: The setter should also remove the regulatory variable from the old model if exists

        from mewpy.regulation.regulatory_model import RegulatoryModel

        if value == self._model:
            pass

        elif not value:
            self._model = None

        elif not isinstance(value, RegulatoryModel):
            raise TypeError(
                "Regulatory model must be a Regmodel object or None")

        else:

            # TODO: This should point all variables and interactions related with this one to the same model, but
            #             #  this might be a recursion problem

            if not self.is_target and not self.is_regulator:
                def model_id(message):
                    warn(message, UserWarning, stacklevel=2)

                model_id("The regulatory variable must be regulator or target to be added to the regulatory model. "
                         "Ignoring the alteration")

            elif self.id not in value.regulators and self.is_regulator:
                value.regulators[self.id] = self
                self._model = value

            elif self.id not in value.targets and self.is_target:
                value.targets[self.id] = self
                self._model = value

            else:
                def model_id(message):
                    warn(message, UserWarning, stacklevel=2)

                model_id("The ID already exists in the dictionary of regulators and targets of the regulatory model. "
                         "Ignoring the alteration")

    @is_regulator.setter
    def is_regulator(self, value):

        if value == self._is_regulator:
            pass

        elif not isinstance(value, bool):
            raise TypeError(
                "is_regulator must be a boolean. Use bool() to convert")

        else:
            self._is_regulator = value

            if self.model is not None:

                if self.id not in self.model.regulators:
                    self.model.regulators[self.id] = self

    @is_target.setter
    def is_target(self, value):

        if value == self._is_target:
            pass

        elif not isinstance(value, bool):
            raise TypeError(
                "is_target must be a boolean. Use bool() to convert")

        else:
            self._is_target = value

            if self.model is not None:

                if self.id not in self.model.targets:
                    self.model.targets[self.id] = self

    @solution.setter
    def solution(self, value):

        if value == self._solution:
            pass

        else:
            raise BaseException("Solution property cannot be changed. "
                                "This property is fulfilled when the regulatory model is simulated.")

    @id.deleter
    def id(self):
        del self._id

    @name.deleter
    def name(self):
        del self._name

    @aliases.deleter
    def aliases(self):
        del self._aliases

    @regulatory_interactions.deleter
    def regulatory_interactions(self):
        del self._regulatory_interactions

    @expression_coef.deleter
    def expression_coef(self):
        del self._expression_coef

    @is_regulator.deleter
    def is_regulator(self):
        del self._is_regulator

    @is_target.deleter
    def is_target(self):
        del self._is_target

    @model.deleter
    def model(self):
        del self._model

    @regulators.deleter
    def regulators(self):
        del self._regulators

    @targets.deleter
    def targets(self):
        del self._targets

    @solution.deleter
    def solution(self):
        del self._solution

    def regulators_gen(self):

        return (value for value in self._regulators.values())

    def targets_gen(self):

        return (value for value in self._targets.values())

    def regulatory_interactions_gen(self):

        return (value for value in self._regulatory_interactions.values())

    def knock_out(self):

        self.expression_coef = 0

    def remove_from_model(self):

        if self.model is not None:

            self.model.remove_regulatory_variable(self.id)

        else:
            def remove_from_model_reg(message):
                warn(message, UserWarning, stacklevel=2)

            remove_from_model_reg("The regulatory variable is not bounded to any regulatory model")
