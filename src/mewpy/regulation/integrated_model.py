import os
from warnings import warn
from mewpy.regulation import RegulatoryModel, RegulatoryVariable, RegulatoryInteraction
from mewpy.io import read_tabular_aliases, read_tabular_regulatory_model
from mewpy.regulation.regulatory_model import solution_decode
from mewpy.simulation import get_simulator, SStatus
from mewpy.simulation import SimulationMethod
from mewpy.simulation.cobra import Simulation as Cobra_model
from mewpy.simulation.reframed import Simulation as Reframed_model
from mewpy.simulation.simulation import SimulationResult
from mewpy.utils.parsing import boolean_rule_from_str, BooleanEvaluator


class IntegratedRegulatoryVariable(RegulatoryVariable):
    """
    Integrated Regulatory Variable

    """

    def __init__(self, identifier,
                 name=None,
                 expression_coef=None,
                 aliases=None,
                 regulatory_interactions=None,
                 targets=None,
                 regulators=None,
                 model=None,
                 cbm_model=None,
                 symbol=None,
                 is_regulator=False,
                 is_target=False):

        self._symbol = None
        self._cbm_model = None
        self._cbm_model_id = None

        super().__init__(identifier,
                         name=name,
                         expression_coef=expression_coef,
                         aliases=aliases,
                         regulatory_interactions=regulatory_interactions,
                         targets=targets,
                         regulators=regulators,
                         model=model,
                         is_regulator=is_regulator,
                         is_target=is_target)

        self._build_integrated_regulatory_variable(cbm_model, symbol)

    def _build_integrated_regulatory_variable(self, cbm_model, symbol):

        if symbol is not None:
            self.symbol = symbol
        self.cbm_model = cbm_model

    @property
    def id(self):
        return getattr(self, '_id', None)

    @property
    def aliases(self):
        return getattr(self, '_aliases', None)

    @property
    def symbol(self):
        return getattr(self, '_symbol', None)

    @property
    def cbm_model(self):
        return getattr(self, '_cbm_model', None)

    @property
    def cbm_model_id(self):
        return getattr(self, '_cbm_model_id', None)

    @id.setter
    def id(self, value):

        if value == self._id:
            pass

        elif not isinstance(value, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        else:

            # symbol, alias = variable_from_str(value)
            #
            # self.symbol = symbol
            # self.aliases = list(alias.values())

            # Uncomment this if we want to keep available the option for altering an id of a given variable

            # if self.model:
            #
            #     if value in self.model.regulators and value in self.model.targets:
            #         def reg_id(message):
            #             warn(message, UserWarning, stacklevel=2)
            #
            #         reg_id("The ID already exists in the dictionary of regulators and targets of the regulatory "
            #                "model. Ignoring the alteration")
            #
            #     else:
            #
            #         if self._id in self.model.regulators:
            #             del self.model.regulators[self._id]
            #             self.model.regulators[value] = self
            #
            #         if self._id in self.model.targets:
            #             del self.model.targets[self._id]
            #             self.model.targets[value] = self
            #
            #         if self.regulatory_interactions:
            #             for reg_inter in self.regulatory_interactions_gen():
            #
            #                 if self._id in reg_inter.regulators:
            #                     del reg_inter.regulators[self._id]
            #                     reg_inter.regulators[value] = self
            #
            #         if self.regulators:
            #             for reg in self.regulators.values():
            #                 del reg.targets[self._id]
            #                 reg.targets[value] = self
            #
            #         if self.targets:
            #             for tar in self.targets.values():
            #                 del tar.regulators[self._id]
            #                 tar.regulators[value] = self
            #
            #         self._id = value
            #
            # else:
            #
            #     if self.regulatory_interactions:
            #
            #         for reg_inter in self.regulatory_interactions_gen():
            #
            #             if self._id in reg_inter.regulators:
            #                 del reg_inter.regulators[self._id]
            #                 reg_inter.regulators[value] = self

            self._id = value

    @aliases.setter
    def aliases(self, value):

        if value == self._aliases:
            pass

        elif not value:
            self._aliases = []

        else:

            if not isinstance(value, list):
                raise TypeError("The aliases must be a list. Use list() to convert")

            if value and not isinstance(value[0], str):
                raise TypeError("The aliases must be a list. Use list() to convert")

            if not self._aliases:
                self._aliases = []

            for val in value:
                if val not in self._aliases:
                    # Reframed keeps the G_, M_ and R_ in sbml files
                    self._aliases = self._aliases + [val, 'G_' + val, 'M_' + val, 'R_' + val]

    @symbol.setter
    def symbol(self, value):

        from sympy.core.symbol import Symbol

        if value == self._symbol:

            pass

        elif not value:

            self._symbol = None

        elif not isinstance(value, Symbol):

            raise TypeError("The symbol must be a sympy.core.symbol.Symbol")

        else:

            self._symbol = value

    @cbm_model.setter
    def cbm_model(self, value):

        if value == self._cbm_model:
            pass

        elif not value:
            self._cbm_model = None

        else:

            if not isinstance(value, Reframed_model) and not isinstance(value, Cobra_model):

                self._cbm_model = get_simulator(value)

            else:

                self._cbm_model = value

            if self.id in self._cbm_model.genes:
                self._cbm_model_id = self.id

            elif self.id in self._cbm_model.reactions:
                self._cbm_model_id = self.id

            elif self.id in self._cbm_model.metabolites:
                self._cbm_model_id = self.id

            else:

                for alias in self.aliases:
                    if alias in self._cbm_model.genes:
                        self._cbm_model_id = alias
                        break
                    if alias in self._cbm_model.reactions:
                        self._cbm_model_id = alias
                        break
                    if alias in self._cbm_model.metabolites:
                        self._cbm_model_id = alias
                        break

    def ko(self):

        """
        Knock-out a given regulatory variable means setting the expression coefficient to zero
        :return:
        """

        self.expression_coef = 0


class IntegratedRegulatoryInteraction(RegulatoryInteraction):
    """
    Integrated Regulatory Interaction

    """

    def __init__(self, identifier,
                 rule,
                 target,
                 regulators=None,
                 model=None,
                 aliases=None):

        self._conditionals = []
        self._variables = None
        self._tree = None

        super().__init__(identifier,
                         rule,
                         target,
                         regulators=regulators,
                         model=model,
                         aliases=aliases)

    @classmethod
    def from_str(cls, identifier,
                 rule,
                 target,
                 target_aliases=None,
                 regulators=None,
                 model=None,
                 cbm_model=None):

        if not isinstance(identifier, str):
            raise TypeError("IDs must be a str object")

        if not isinstance(rule, str):
            raise TypeError("rule must be a string object")

        if isinstance(target, str):

            if model is not None:

                tar_var = model.add_target(identifier=target,
                                           name=target,
                                           aliases=target_aliases,
                                           regulators=regulators)

            else:

                tar_var = IntegratedRegulatoryVariable(identifier=target,
                                                       name=target,
                                                       aliases=target_aliases,
                                                       cbm_model=cbm_model)

            return cls(identifier=identifier,
                       rule=rule,
                       target=tar_var,
                       model=model)

        elif isinstance(target, IntegratedRegulatoryVariable):

            return cls(identifier=identifier,
                       rule=rule,
                       target=target,
                       model=model)

        else:
            raise TypeError("Target must be a string or a RegulatoryVariable object")

    @property
    def id(self):
        return getattr(self, '_id', None)

    @property
    def rule(self):
        return getattr(self, '_rule', None)

    # @property
    # def sympify(self):
    #     return getattr(self, '_sympify', None)
    #
    # @property
    # def lambdify(self):
    #     return getattr(self, '_lambdify', None)

    @property
    def variables(self):
        return getattr(self, '_variables', None)

    @property
    def tree(self):
        return getattr(self, '_tree', None)

    @id.setter
    def id(self, value):

        if value == self._id:
            pass

        elif not isinstance(value, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        else:

            # Uncomment this if we want to keep available the option for altering an id of a given variable
            # change for the regulatory variables property

            # if self.model:
            #
            #     if value in self.model.regulators and value in self.model.targets:
            #         def reg_id(message):
            #             warn(message, UserWarning, stacklevel=2)
            #
            #         reg_id("The ID already exists in the dictionary of regulators and targets of the regulatory "
            #                "model. Ignoring the alteration")
            #
            #     else:
            #
            #         if self._id in self.model.regulatory_interactions:
            #             del self.model.regulatory_interactions[self._id]
            #             self.model.regulatory_interactions[value] = self
            #
            #         if self.regulators:
            #             for reg in self.regulators_gen():
            #                 if self._id in reg.regulatory_interactions:
            #                     del reg.regulatory_interactions[self._id]
            #                     reg.regulatory_interactions[value] = self
            #
            #         if self.target:
            #             if self._id in self.target.regulatory_interactions:
            #                 del self.target.regulatory_interactions[self._id]
            #                 self.target.regulatory_interactions.regulatory_interactions[value] = self
            #
            #         self._id = value
            #
            # else:
            #
            #     if self.regulators:
            #         for reg in self.regulators_gen():
            #             if self._id in reg.regulatory_interactions:
            #                 del reg.regulatory_interactions[self._id]
            #                 reg.regulatory_interactions[value] = self
            #
            #     if self.target:
            #         if self._id in self.target.regulatory_interactions:
            #             del self.target.regulatory_interactions[self._id]
            #             self.target.regulatory_interactions.regulatory_interactions[value] = self

            self._id = value

    @rule.setter
    def rule(self, value):

        if value == self._rule:
            pass

        elif not isinstance(value, str):
            raise TypeError("Regulatory interaction rule must be a string. Use str() to convert")

        else:

            self.old_rule = value

            self._rule, self._elements, self._tree, self._variables, aliases, self._conditionals = \
                boolean_rule_from_str(value)

            if self.regulators:

                for key in self.regulators:
                    if key not in self.variables:
                        raise KeyError("The regulators property doesn't match all regulators in the regulatory "
                                       "interaction rule")

                    self.regulators[key].regulatory_interactions[self.id] = self
                    self.regulators[key].is_regulator = True
                    self.regulators[key].model = self.model
                    self.regulators[key].targets[self.target.id] = self.target
                    self.regulators[key].aliases = aliases[key]
                    # self.regulators[key].symbol = symbols[key]

                    self.target.regulators[key] = self.regulators[key]

            else:

                for reg in self.variables:

                    if self.model is not None:

                        reg_var = self.model.add_regulator(identifier=reg,
                                                           name=reg,
                                                           aliases=aliases[reg],
                                                           regulatory_interactions={self.id: self},
                                                           targets={self.target.id: self.target})

                    else:

                        reg_var = IntegratedRegulatoryVariable(identifier=reg,
                                                               name=reg,
                                                               aliases=aliases[reg],
                                                               regulatory_interactions={self.id: self},
                                                               targets={self.target.id: self.target})

                    self.regulators[reg] = reg_var
                    self.target.regulators[reg] = reg_var

    # @sympify.setter
    # def sympify(self, value):
    #
    #     raise BaseException("lambdify property cannot be changed. "
    #                         "please create/set a new regulatory interaction rule.")
    #
    # @lambdify.setter
    # def lambdify(self, value):
    #
    #     raise BaseException("lambdify property cannot be changed. "
    #                         "please create/set a new regulatory interaction rule.")

    @variables.setter
    def variables(self, value):

        raise BaseException("variables property cannot be changed. "
                            "please create/set a new regulatory interaction rule.")

    # @tree.setter
    # def tree(self, value):
    #
    #     from mewpy.utils.parsing import Node
    #
    #     if value == self._tree:
    #         pass
    #
    #     if not isinstance(value, Node):
    #         raise TypeError("The tree must be a Node object. Use Node() to convert")
    #
    #     else:
    #         self._tree = value

    @rule.deleter
    def rule(self):
        del self._rule
        del self._variables
        del self._tree

    def evaluate(self, state_map=None):

        if self.tree.value is None:
            return

        if not state_map:
            state_map = {}

        true_list = []
        _state_map = {reg.id: reg.expression_coef for reg in self.regulators.values()}

        new_state_map = {}
        for key, value in _state_map.items():

            if key in state_map:
                value = state_map[key]

            if key not in self._conditionals:

                if value:
                    true_list.append(key)

            else:
                new_state_map[key] = value

        evaluator = BooleanEvaluator(true_list, new_state_map)
        res = self.tree.evaluate(evaluator.f_operand, evaluator.f_operator)
        res = solution_decode(res)

        # print("{} evaluation result: {}".format(self.id, str(res)))

        return res

    def print_tree(self):

        self.tree.print_node()


class IntegratedModel(RegulatoryModel):
    """
    Integrated Model base class.
    Do not instantiate this class. This works as bridge/interface for linking metabolic and regulatory entities found
    in regulatory and cbm models. This is crucial for the integrated methods analysis available at mewpy (e.g. RFBA and
    SRFBA) to work properly.

    Notes:
        -Only supports boolean regulatory models (by now)

    """

    def __init__(self,
                 identifier,
                 name=None,
                 cbm_model=None,
                 cbm_simulation_interface=None,
                 targets=None,
                 regulators=None,
                 regulatory_interactions=None,
                 initial_state=None):

        # CBM model
        self._cbm_model = None
        self._cbm_simulation_interface = None

        # New from CBM model
        self._cbm_simulation_method = SimulationMethod.pFBA
        self._maximize = True
        self._gprs = {}
        self._gprs_evaluator = {}
        self._metabolic_regulatory_genes = {}
        self._metabolic_regulatory_metabolites = {}
        self._metabolic_regulatory_reactions = {}
        self._regulatory_conditions = {}
        self._aliases_map = {}
        # A self._metabolic_regulatory_gprs can be implemented by replacing a given gene into its gpr rules by the
        # regulatory rule associated

        # Regulatory model
        super().__init__(identifier,
                         name,
                         targets,
                         regulators,
                         regulatory_interactions)

        self._initial_state = {}
        self._solution = []

        # Build
        self.__build_integrated_model__(cbm_model, cbm_simulation_interface, initial_state)

    def __build_integrated_model__(self, cbm_model, cbm_simulation_interface, initial_state):

        if cbm_simulation_interface is not None:

            if cbm_model is None:
                raise BaseException("A REFRAMED or COBRApy Constraint-based model must be provided")

            else:
                self._cbm_model = cbm_model

            self.cbm_simulation_interface = cbm_simulation_interface

        else:

            self.cbm_model = cbm_model

        self.initial_state = initial_state

        self.__infer_metabolic_regulatory_variables__()

    @classmethod
    def from_tabular_format(cls, regulatory_file,
                            cbm_model,
                            cbm_simulation_interface=None,
                            identifier=None,
                            initial_state=None,
                            filter_nan=False,
                            sep=None,
                            id_col=None,
                            rule_col=None,
                            aliases_cols=None,
                            header=None,
                            special_chars=None,
                            compartments_ids=None):

        """

        The standard class method for loading an Integrated Model of any type from a tabular format file (txt, csv, etc)
        This method will build an integrated model from scratch using a cobra CBM model, reframed CBM model,
        or mewpy simulation object, and the regulatory rules written in tabular-like format file.

        NOTE THAT:
            -Only one rule is allowed per target.
            -So for only boolean rules are supported
            -Special chars present in the regulatory variables and regulatory rules will be filter and replaced by char
            description plus double underscore (e.g. / -> _slash_).
            -Regulatory variables that start with digits will also be filter and replaced

        Rules can be set as such:

            G1, FruE, G2 and (G3 or M_acetate>0)
            G2, FruR, ON
            G3, AcR, M_02>5

        :param regulatory_file: str, path to the regulatory model
        :param cbm_model: cobra or reframed model objects
        :param cbm_simulation_interface: mewpy simulation object
        :param identifier: str, integrated model identifier
        :param initial_state: dict, initial state of the regulatory variables. If none, the state of all regulatory variables is set to off (zero)
        :param filter_nan: bool, filter targets with empty rules
        :param sep: str, separator
        :param id_col: int, index of the column of the regulatory variables identifiers (those that appear in the regulatory rules)
        :param rule_col: int, index of the column of the regulatory rules
        :param aliases_cols: int, index of the column of the regulatory variables aliases
        :param header: int or None, index of the header row or None if there is none
        :param special_chars: dict, key is the special char, value is its replacement
        :param compartments_ids: dict, key is the compartment id, value is its replacement
        :return:
        """

        if not isinstance(regulatory_file, str):
            raise TypeError("File or file path must be a string")

        if not identifier:
            _, identifier = os.path.split(regulatory_file)
            identifier = os.path.splitext(identifier)[0]

        model = cls(identifier=identifier,
                    name=identifier,
                    cbm_model=cbm_model,
                    cbm_simulation_interface=cbm_simulation_interface,
                    initial_state=None)

        df = read_tabular_regulatory_model(regulatory_file, sep=sep, id_col=id_col, rule_col=rule_col,
                                           aliases_cols=aliases_cols, header=header, special_chars=special_chars,
                                           compartments_ids=compartments_ids)

        aliases_cols = [col for col in df.columns if col != 'ids' and col != 'rules']

        for gene in df.index:

            reg_inter_id = gene + '_regulation'
            reg_inter_rule = df.loc[gene, 'rules']
            target_aliases = list(df.loc[gene, aliases_cols])

            if filter_nan and str(reg_inter_rule) == '':
                continue

            model.add_regulatory_interaction(reg_inter_id, reg_inter_rule, gene, target_aliases=target_aliases)

        # last part of build integrated model, as when the initializer is ran, the regulatory variables are not yet
        # there
        model.initial_state = initial_state

        model.__infer_metabolic_regulatory_variables__()

        return model

    @property
    def cbm_model(self):
        return getattr(self, '_cbm_model', None)

    @property
    def cbm_simulation_interface(self):
        return getattr(self, '_cbm_simulation_interface', None)

    @property
    def reactions(self):
        return self.cbm_simulation_interface.reactions

    @property
    def genes(self):
        return self.cbm_simulation_interface.genes

    @property
    def metabolites(self):
        return self.cbm_simulation_interface.metabolites

    @property
    def medium(self):
        return self.cbm_simulation_interface.medium

    @property
    def objective(self):
        return self.cbm_simulation_interface.objective

    @property
    def environmental_conditions(self):
        return self.cbm_simulation_interface.environmental_conditions

    @property
    def constraints(self):
        return self.cbm_simulation_interface.constraints

    @property
    def solver(self):
        return self.cbm_simulation_interface.solver

    @property
    def reference(self):
        return self.cbm_simulation_interface.reference

    @property
    def cbm_simulation_method(self):
        return getattr(self, '_cbm_simulation_method', None)

    @property
    def maximize(self):
        return getattr(self, '_maximize', None)

    @property
    def gprs(self):
        return getattr(self, '_gprs', None)

    @property
    def metabolic_regulatory_genes(self):
        return getattr(self, '_metabolic_regulatory_genes', None)

    @property
    def metabolic_regulatory_metabolites(self):
        return getattr(self, '_metabolic_regulatory_metabolites', None)

    @property
    def metabolic_regulatory_reactions(self):
        return getattr(self, '_metabolic_regulatory_reactions', None)

    @property
    def regulatory_conditions(self):
        return getattr(self, '_regulatory_conditions', None)

    @property
    def initial_state(self):
        return getattr(self, '_initial_state', None)

    @property
    def solution(self):
        return getattr(self, '_solution', None)

    # @property
    # def metabolic_regulatory_gprs(self):
    #     return getattr(self, '_metabolic_regulatory_gprs', None)

    @cbm_model.setter
    def cbm_model(self, value):

        if value == self._cbm_model:
            pass

        elif not value:
            raise BaseException("A REFRAMED or COBRApy Constraint-based model must be provided")

        else:

            self.cbm_simulation_interface = get_simulator(value)

            self._cbm_model = value

    @cbm_simulation_interface.setter
    def cbm_simulation_interface(self, value):

        if value == self._cbm_simulation_interface:
            pass

        elif not value:
            raise BaseException("A simulation object must be provided")

        elif not isinstance(value, Reframed_model) and not isinstance(value, Cobra_model):
            raise TypeError("A simulation object must be provided")

        else:

            self._cbm_simulation_interface = value

            # Building utils
            self.__build_gprs__()
            self.__infer_metabolic_regulatory_variables__()

    @reactions.setter
    def reactions(self, value):

        raise BaseException("reactions property cannot be changed")

    @genes.setter
    def genes(self, value):

        raise BaseException("genes property cannot be changed")

    @metabolites.setter
    def metabolites(self, value):

        raise BaseException("metabolites property cannot be changed")

    @medium.setter
    def medium(self, value):

        raise BaseException("medium property cannot be changed")

    @objective.setter
    def objective(self, value):

        if value == self._cbm_simulation_interface.objective:
            pass

        elif not value:
            pass

        elif not isinstance(value, dict):
            raise TypeError("The objective must be a dict object")

        else:

            self._cbm_simulation_interface.objective = value

    @environmental_conditions.setter
    def environmental_conditions(self, value):

        if value == self._cbm_simulation_interface.environmental_conditions:
            pass

        elif not value:
            self._cbm_simulation_interface.environmental_conditions = {}

        elif not isinstance(value, dict):
            raise TypeError("The environmental conditions must be a dict object")

        else:

            self._cbm_simulation_interface.environmental_conditions = value

    @constraints.setter
    def constraints(self, value):

        if value == self._cbm_simulation_interface.constraints:
            pass

        elif not value:
            self._cbm_simulation_interface.constraints = {}

        elif not isinstance(value, dict):
            raise TypeError("The constraints must be a dict object")

        else:

            self._cbm_simulation_interface.constraints = value

    @solver.setter
    def solver(self, value):

        if value == self._cbm_simulation_interface.solver:
            pass

        elif not value:
            pass

        elif not isinstance(value, str):
            raise TypeError("The solver must be a str object")

        else:

            self._cbm_simulation_interface.solver = value

    @reference.setter
    def reference(self, value):

        raise BaseException("Reference property cannot be changed")

    @cbm_simulation_method.setter
    def cbm_simulation_method(self, value):

        if value == self._cbm_simulation_method:
            pass

        elif not value:
            pass

        else:
            self._cbm_simulation_method = value

    @maximize.setter
    def maximize(self, value):

        if value == self._maximize:
            pass

        elif value is None:
            pass

        elif isinstance(value, bool):
            raise TypeError("Maximize property must be a boolean object")

        else:
            self._maximize = value

    @gprs.setter
    def gprs(self, value):

        if value == self._gprs:
            pass

        else:
            raise BaseException("gprs property cannot be changed")

    @metabolic_regulatory_genes.setter
    def metabolic_regulatory_genes(self, value):

        if value == self._metabolic_regulatory_genes:
            pass

        else:
            raise BaseException("metabolic_regulatory_genes property cannot be changed")

    @metabolic_regulatory_metabolites.setter
    def metabolic_regulatory_metabolites(self, value):

        if value == self._metabolic_regulatory_metabolites:
            pass

        else:
            raise BaseException("metabolic_regulatory_metabolites property cannot be changed")

    @metabolic_regulatory_reactions.setter
    def metabolic_regulatory_reactions(self, value):

        if value == self._metabolic_regulatory_reactions:
            pass

        else:
            raise BaseException("metabolic_regulatory_reactions property cannot be changed")

    @regulatory_conditions.setter
    def regulatory_conditions(self, value):

        if value == self._regulatory_conditions:
            pass

        else:
            raise BaseException("regulators_environmental_conditions property cannot be changed")

    @initial_state.setter
    def initial_state(self, value):

        if not value:

            # If not defined,
            # the initial state of all variables is always zero

            self._initial_state = {}

            for variable in self.regulatory_variables_gen():
                variable.expression_coef = 0
                self._initial_state[variable.id] = 0

        elif not isinstance(value, dict):
            raise TypeError("The initial state must be a dict object. Use dict() to convert")

        else:

            # initial state comprehends the state of all regulatory variables in the regulatory model
            # partial setting is accepted and only the respective variables are altered

            # State or current state is an additional property to be added. Alternatively, this property keeps track
            # of the current state.

            # Either way, they both state and initial_state would change the state of all regulatory variables when
            # using the setter

            if not self._initial_state:

                self._initial_state = {}

                for variable in self.regulatory_variables_gen():

                    if variable.id in value:
                        variable.expression_coef = value[variable.id]
                        self._initial_state[variable.id] = value[variable.id]

                    else:
                        variable.expression_coef = 0
                        self._initial_state[variable.id] = 0

            else:

                for key, val in value.items():

                    if key in self.regulatory_variables:
                        self.regulatory_variables[key].expression_coef = val
                        self._initial_state[key] = val

    def initial_state_gen(self):

        """
        Get initial state
        :return: generator object of initial state items
        """

        return ((key, val) for key, val in self.initial_state.items())

    def metabolic_regulatory_genes_gen(self):

        """
        Get metabolic regulatory genes
        :return: generator object of metabolic regulatory genes values
        """

        return (val for val in self.metabolic_regulatory_genes.values())

    def metabolic_regulatory_reactions_gen(self):

        """
        Get metabolic regulatory genes
        :return: generator object of metabolic regulatory reactions values
        """

        return (val for val in self.metabolic_regulatory_reactions.values())

    def metabolic_regulatory_metabolites_gen(self):

        """
        Get metabolic regulatory genes
        :return: generator object of metabolic regulatory metabolites values
        """

        return (val for val in self.metabolic_regulatory_metabolites.values())

    def regulatory_conditions_gen(self):

        """
        Get regulatory conditions
        :return: generator object of regulatory conditions values
        """

        return (val for val in self.regulatory_conditions.values())

    def add_regulator(self, identifier,
                      name=None,
                      expression_coef=None,
                      aliases=None,
                      regulatory_interactions=None,
                      targets=None,
                      regulators=None,
                      symbol=None):

        """


        :param identifier:
        :param name:
        :param expression_coef:
        :param aliases:
        :param regulatory_interactions:
        :param targets:
        :param regulators:
        :param symbol:
        :return:
        """

        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        if identifier in self.regulators:

            regulator = self.regulators[identifier]

            if name is not None:
                regulator.name = name
            if expression_coef is not None:
                regulator.expression_coef = expression_coef
            if aliases is not None:
                regulator.aliases = aliases
            if regulatory_interactions is not None:
                regulator.regulatory_interactions.update(regulatory_interactions)
            if targets is not None:
                regulator.targets.update(targets)
            if regulators is not None:
                regulator.regulators.update(regulators)
            if symbol is not None:
                regulator.symbol = symbol

            return regulator

        if identifier in self.targets:

            regulator = self.targets[identifier]

            if name is not None:
                regulator.name = name
            if expression_coef is not None:
                regulator.expression_coef = expression_coef
            if aliases is not None:
                regulator.aliases = aliases
            if regulatory_interactions is not None:
                regulator.regulatory_interactions.update(regulatory_interactions)
            if targets is not None:
                regulator.targets.update(targets)
            if regulators is not None:
                regulator.regulators.update(regulators)
            if symbol is not None:
                regulator.symbol = symbol
            regulator.is_regulator = True

            return regulator

        regulator = IntegratedRegulatoryVariable(
            identifier=identifier,
            name=name,
            expression_coef=expression_coef,
            aliases=aliases,
            regulatory_interactions=regulatory_interactions,
            targets=targets,
            regulators=regulators,
            model=self,
            cbm_model=self.cbm_simulation_interface,
            symbol=symbol,
            is_regulator=True,
            is_target=False)

        self.regulatory_variables[regulator.id] = regulator

        return regulator

    def add_target(self, identifier,
                   name=None,
                   expression_coef=None,
                   aliases=None,
                   regulatory_interactions=None,
                   targets=None,
                   regulators=None,
                   symbol=None):

        """

        :param identifier:
        :param name:
        :param expression_coef:
        :param aliases:
        :param regulatory_interactions:
        :param targets:
        :param regulators:
        :param symbol:
        :return:
        """
        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        if identifier in self.targets:

            target = self.targets[identifier]

            if name is not None:
                target.name = name
            if expression_coef is not None:
                target.expression_coef = expression_coef
            if aliases is not None:
                target.aliases = aliases
            if regulatory_interactions is not None:
                target.regulatory_interactions.update(regulatory_interactions)
            if targets is not None:
                target.targets.update(targets)
            if regulators is not None:
                target.regulators.update(regulators)
            if symbol is not None:
                target.symbol = symbol

            return target

        if identifier in self.regulators:

            target = self.regulators[identifier]

            if name is not None:
                target.name = name
            if expression_coef is not None:
                target.expression_coef = expression_coef
            if aliases is not None:
                target.aliases = aliases
            if regulatory_interactions is not None:
                target.regulatory_interactions.update(regulatory_interactions)
            if targets is not None:
                target.targets.update(targets)
            if regulators is not None:
                target.regulators.update(regulators)

            target.is_target = True

            return target

        target = IntegratedRegulatoryVariable(
            identifier=identifier,
            name=name,
            expression_coef=expression_coef,
            aliases=aliases,
            regulatory_interactions=regulatory_interactions,
            targets=targets,
            regulators=regulators,
            model=self,
            cbm_model=self.cbm_simulation_interface,
            symbol=symbol,
            is_regulator=False,
            is_target=True)

        self.regulatory_variables[target.id] = target

        return target

    def add_regulatory_interaction(self, identifier,
                                   rule,
                                   target,
                                   target_aliases=None,
                                   regulators=None):

        """

        :param identifier:
        :param rule:
        :param target:
        :param target_aliases:
        :param regulators:
        :return:
        """

        if not isinstance(identifier, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        if not isinstance(rule, str):
            raise TypeError("The regulatory rule must be a string. Use str() to convert")

        if not isinstance(target, str) and not isinstance(target, RegulatoryVariable):
            raise TypeError("Target must be a string or a RegulatoryVariable object")

        if not target_aliases:
            target_aliases = []

        if identifier in self.regulatory_interactions:
            def add_reg_inter_warn(message):
                warn(message, UserWarning, stacklevel=2)

            add_reg_inter_warn("The ID already exists in the dictionary of regulatory interactions of the regulatory "
                               "model. Ignoring the alteration")
            return

        return IntegratedRegulatoryInteraction.from_str(identifier=identifier,
                                                        rule=rule,
                                                        target=target,
                                                        target_aliases=target_aliases,
                                                        regulators=regulators,
                                                        model=self,
                                                        cbm_model=self.cbm_simulation_interface)

    def __build_gprs__(self):

        self._gprs = {}

        # A GPR object can also be created and stored under the gprs variable. This object can be much likely
        # designed as the regulatory interaction one
        # So far, the gprs_evaluater contains a list of sympy symbols as arguments, a sympy expression and sympy
        # lambda function for a given gpr rule
        self._gprs_evaluator = {}

        for reac_id in self.reactions:

            gpr = self._cbm_simulation_interface.get_gpr(reac_id)

            if not gpr:
                gpr = ''

            self._gprs[reac_id] = gpr

            self._gprs_evaluator[reac_id] = boolean_rule_from_str(gpr)

    def __infer_metabolic_regulatory_variables__(self):

        # stores all gene-like regulatory variables available in the cbm and regulatory model
        self._metabolic_regulatory_genes = {}

        # stores all reaction-like  regulatory variables available in the cbm and regulatory model
        self._metabolic_regulatory_reactions = {}

        # stores all metabolite-like  regulatory variables available in the cbm and regulatory model
        self._metabolic_regulatory_metabolites = {}

        # stores all regulatory variables (regulators only) not available in the cbm model but available in the
        # regulatory one
        self._regulatory_conditions = {}

        # map between all regulatory and metabolic entities in both models, as they might be linked by alias
        self._aliases_map = {}

        for variable in self.regulatory_variables_gen():

            if variable.cbm_model_id in self.genes:
                self._metabolic_regulatory_genes[variable.id] = variable
                self._aliases_map[variable.cbm_model_id] = variable.id

            elif variable.cbm_model_id in self.reactions:
                self._metabolic_regulatory_reactions[variable.id] = variable
                self._aliases_map[variable.cbm_model_id] = variable.id

            elif variable.cbm_model_id in self.metabolites:
                self._metabolic_regulatory_metabolites[variable.id] = variable
                self._aliases_map[variable.cbm_model_id] = variable.id

            else:
                if variable.id not in self.targets:
                    self._regulatory_conditions[variable.id] = variable

    def decode_metabolic_state(self, state):

        """
        Method responsible for decoding the RFBA metabolic state, namely the state of all metabolic genes associated
        at least with one reaction in the GPRs rule.

        :param state: dict, key is the id of the regulatory variable (metabolic target) while value can be 0, 1 or float (reactions and metabolites predicates)
        :return: SimulationResult object.
        """

        # This function decodes a given metabolic state (state of the metabolic genes) and simulates the cbm model
        # with the resulting constraints

        constraints = {}

        # for every reaction that has a valid gpr rule
        for rxn, (_, _, tree, variables, _, conditionals) in self._gprs_evaluator.items():

            if tree.value is None:
                continue

            # creates a state map (dict) for that specific rule
            state_map = {}
            true_list = []

            for arg in variables:
                if arg in self._aliases_map:
                    if self._aliases_map[arg] in state:
                        state_map[arg] = state[self._aliases_map[arg]]
                        continue
                state_map[arg] = 1

            new_state_map = {}

            for key, value in state_map.items():

                if key not in conditionals:

                    if value:
                        true_list.append(key)

                else:
                    new_state_map[key] = value

            # solving the rule
            evaluator = BooleanEvaluator(true_list, new_state_map)
            res = tree.evaluate(evaluator.f_operand, evaluator.f_operator)
            res = solution_decode(res)

            # if the reaction is off
            if not bool(res):
                constraints[rxn] = (0.0, 0.0)

        return constraints

    def simulate_cbm_model(self, objective=None, maximize=True, method=None, constraints=None):

        """

        Standard CBM model simulation in mewpy. See simulate method of the Simulation object for further detail

        :param objective: dict
        :param maximize: bool
        :param method: SimulationMethod object
        :param constraints: dict, additional constraints focused on the reactions
        :return: SimulationResult object, it contains the result of the simulate method.
        """

        self.objective = objective
        self.maximize = maximize
        self.cbm_simulation_method = method

        # broader try, as sometimes an infeasible exception is thrown by mewpy simulate interface. Sepcial when the
        # initial state of the model is fulfilled by zeros
        try:
            return self.cbm_simulation_interface.simulate(objective=self.objective, method=self.cbm_simulation_method,
                                                          maximize=self.maximize, constraints=constraints)

        except:
            return SimulationResult(self.cbm_model,
                                    0.0,
                                    {},
                                    status=SStatus.INFEASIBLE)

    def simulate(self):

        raise NotImplementedError("IntegratedModel is a base classe without clear instructions over how regulatory "
                                  "constraints should be used for simulation. Instantiate RFBA or SRFBA instead")

    def __populate_solution__(self):

        for regvar in self.regulatory_variables_gen():
            regvar._solution = regvar.expression_coef
            regvar.expression_coef = self.initial_state[regvar.id]

    def update_aliases_from_tabular_format_file(self, aliases_file,
                                                sep=None,
                                                id_col=None,
                                                aliases_cols=None,
                                                header=None):
        """

        Updating the aliases of the regulatory variables in the integrated model allows to link regulatory variables
        to metabolic entities (e.g. genes, reactions or metabolites) present in the regulatory rules.

        Linking metabolic and regulatory entities is crucial for the integrated analysis to work properly.

        :param aliases_file: str, path to aliases file in the tabular format (txt, csv, ...)
        :param sep: str, separator
        :param id_col: int, index of the regulatory variables identifiers
        :param aliases_cols: int, index of the regulatory variables aliases
        :param header: int or None, index of the header row or None if there is none
        :return:
        """

        if not isinstance(aliases_file, str):
            raise TypeError("File or file path must be a string")

        df = read_tabular_aliases(aliases_file, sep=sep, id_col=id_col, aliases_cols=aliases_cols, header=header)

        aliases_cols = [col for col in df.columns if col != 'ids']

        for row in df.index:

            new_aliases = list(df.loc[row, aliases_cols])

            if row in self.regulatory_variables:
                self.regulatory_variables[row].aliases = new_aliases
                self.regulatory_variables[row].cbm_model = None
                self.regulatory_variables[row].cbm_model = self.cbm_simulation_interface

        self.__infer_metabolic_regulatory_variables__()

    # def essential_regulatory_variables(self):
    #
    #     # TODO: Implement a method to find all essential regulatory variables following the approach essential_genes
    #
    #     pass
