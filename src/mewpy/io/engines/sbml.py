import os
from functools import partial
from math import ceil
from typing import Union, TYPE_CHECKING

from mewpy.mew.algebra import Expression, Symbol, Or, And, NoneAtom, Float
from mewpy.model import RegulatoryModel, MetabolicModel
from mewpy.io.dto import DataTransferObject, VariableRecord, History, FunctionTerm, CompartmentRecord
from mewpy.util.constants import ModelConstants

from .engine import Engine
from .engines_utils import (build_symbolic,
                            ASTNODE_BOOLEAN_VALUES, ASTNODE_RELATIONAL_OPERATORS,
                            ASTNODE_NAME, ASTNODE_BOOLEAN_OPERATORS, ASTNODE_VALUES,
                            pattern_notes,
                            f_id, fs_id,
                            F_GENE, F_SPECIE, F_REACTION, F_SPECIE_REV, F_GENE_REV, F_REACTION_REV, F_TRANSITION,
                            F_TRANSITION_REV,
                            convert_fbc, get_sbml_doc_to_write, get_sbml_doc_to_read,
                            UNIT_ID, UNITS,
                            add_sbml_parameter, LOWER_BOUND_ID, UPPER_BOUND_ID, ZERO_BOUND_ID,
                            BOUND_MINUS_INF, BOUND_PLUS_INF,
                            SBO_DEFAULT_FLUX_BOUND, SBO_FLUX_BOUND,
                            get_sbml_lb_id, get_sbml_ub_id, write_sbml_doc,
                            set_math, set_gpr,
                            expression_warning, sbml_warning)

if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


class RegulatorySBML(Engine):

    def __init__(self, io, config, model=None):

        super().__init__(io, config, model)

    @property
    def model_type(self):
        return 'regulatory'

    @property
    def model(self):

        if self._model is None:

            identifier = self.dto.id

            if not identifier:
                identifier = self.get_identifier()

            return RegulatoryModel(identifier=identifier)

        return self._model

    @staticmethod
    def _parse_coefficients_note(notes):

        if not notes:
            return set()

        coefficients = set()

        for row in notes.split('\n'):

            row = row.strip()

            if row.startswith('<p>') and row.endswith('</p>'):
                st, level = row.split(':')

                st = st.strip('<p>')
                level = level.strip('</p>')

                if 'state' in st.lower():
                    st = int(st.lower().replace('state', '').strip())

                lb = ModelConstants.REACTION_LOWER_BOUND
                ub = ModelConstants.REACTION_UPPER_BOUND

                level = level.replace('+inf', f'{ub}').replace('-inf', f'{lb}')

                min_coef, max_coef = level.split(',')

                min_coef = float(min_coef[1:])
                coefficients.add(min_coef)

                max_coef = float(max_coef[0:-1])
                coefficients.add(max_coef)

        return coefficients

    @staticmethod
    def _parse_initial_level_note(notes):

        if not notes:
            return

        for row in notes.split('\n'):

            row = row.strip()

            if row.startswith('<p>') and row.endswith('</p>'):
                parameter, level = row.split(':')

                parameter = parameter.strip('<p>')
                level = level.strip('</p>')

                if 'initial' in parameter.lower():
                    return float(level.replace(' ', ''))

    def _update_qual_species_coefficients(self, qual_species, value):

        # all species are created with the following coefficients: {0:0, 1: maximum level}.
        # but intermediate states can be encoded into the transitions and thus added to the coefficients dict

        value = float(value)

        if qual_species not in self.dto.variables:
            raise SyntaxError

        coefficients = self.dto.variables[qual_species].coefficients

        if min(coefficients) > value > max(coefficients):
            raise SyntaxError

        coefficients.add(value)

    def _parse_values(self, ast_node, ast_node_qual_species, inputs_thresholds):

        node_type = ast_node.getType()

        if node_type in ASTNODE_BOOLEAN_VALUES:

            return ASTNODE_BOOLEAN_VALUES[node_type]

        elif ast_node.getType() == ASTNODE_NAME:

            # In this case, the right child is encoded in the list of inputs

            input_id = fs_id(ast_node.getName(), (F_SPECIE, F_GENE, F_REACTION, F_TRANSITION))
            value = inputs_thresholds.get(input_id, None)

            if value is None:
                raise SyntaxError

            regulator_id = fs_id(ast_node_qual_species.getName(), (F_SPECIE, F_GENE, F_REACTION))
            self._update_qual_species_coefficients(regulator_id, value)

            return Float(value)

        elif ast_node.getType() in ASTNODE_VALUES:

            # In this case, the right child is encoded in the regulators current coefficients

            value = ast_node.getValue()

            if value is None:
                raise SyntaxError

            regulator_id = fs_id(ast_node_qual_species.getName(), (F_SPECIE, F_GENE, F_REACTION))
            self._update_qual_species_coefficients(regulator_id, value)

            return Float(value)

        else:
            raise SyntaxError

    def _parse_symbols(self, ast_node):

        ast_node_id = fs_id(ast_node.getName(), (F_SPECIE, F_GENE, F_REACTION))

        if ast_node_id not in self.dto.variables:
            raise SyntaxError

        return Symbol(value=ast_node_id)

    def _parse_math_node(self, ast_node, inputs_thresholds):

        node_type = ast_node.getType()

        if node_type in ASTNODE_BOOLEAN_OPERATORS:

            Operator = ASTNODE_BOOLEAN_OPERATORS[node_type]
            variables = []

            for child in range(ast_node.getNumChildren()):
                variable = self._parse_math_node(ast_node.getChild(child), inputs_thresholds)

                variables.append(variable)

            return Operator(variables=variables)

        elif node_type in ASTNODE_RELATIONAL_OPERATORS:

            Operator = ASTNODE_RELATIONAL_OPERATORS[node_type]
            symbol = self._parse_symbols(ast_node.getChild(0))
            value = self._parse_values(ast_node.getChild(1), ast_node.getChild(0), inputs_thresholds)

            return Operator(variables=[symbol, value])

        elif node_type == ASTNODE_NAME:

            return self._parse_symbols(ast_node)

        else:

            raise SyntaxError

    def parse_math_node(self, ast_node, inputs_thresholds):

        """
        Reads and parses a node of type math ASTNode into an algebraic expression.

        :param ast_node: ASTNode, math ASTNode from functionTerm.getMath()
        :param inputs_thresholds:
        :return: Symbolic
        """

        if not inputs_thresholds:
            inputs_thresholds = {}

        try:

            symbolic = self._parse_math_node(ast_node, inputs_thresholds)

        except SyntaxError:

            self.warnings.append(partial(expression_warning,
                                         f'{ast_node} cannot be parsed. Assigning empty expression instead'))

            symbolic = NoneAtom()

        return symbolic

    def get_identifier(self):

        if os.path.exists(self.io):
            _, identifier = os.path.split(self.io)
            return os.path.splitext(identifier)[0]

        return 'model'

    def _open_to_read(self):

        self._dto = DataTransferObject()

        # -----------------------------------------------------------------------------
        # Doc, Model and Qual Plugin
        # -----------------------------------------------------------------------------

        self.dto.doc = get_sbml_doc_to_read(self.io)

        self.dto.model = self.dto.doc.getModel()

        if self.dto.model is None:
            raise OSError(f'{self.io} is not a valid input. Model SBML section is missing. '
                          f'Provide a correct path or file handler')

        self.dto.qual_plugin = self.dto.model.getPlugin('qual')

        if self.dto.qual_plugin is None:
            raise OSError(f'Although {self.io} is a valid SBML file, the qual plugin was not detected. Thus, '
                          f'regulatory interactions cannot be determined')

        identifier = self.dto.model.getIdAttribute()

        if not identifier:
            self.warnings.append(partial(sbml_warning, 'Model identifier is not encoded in the SBML file'))

            identifier = self.get_identifier()

        self.dto.id = identifier
        self.dto.name = self.dto.model.getName()

    def _open_to_write(self):

        self._dto = DataTransferObject()

        # -----------------------------------------------------------------------------
        # Doc, Model and FBC Plugin
        # -----------------------------------------------------------------------------

        self.dto.doc = get_sbml_doc_to_write(self.io,
                                             level=3,
                                             version=1,
                                             packages=('qual',),
                                             packages_version=(1,),
                                             packages_required=(True,),
                                             sbo_term=False)

        if self.dto.model is None:

            self.dto.model = self.dto.doc.createModel()

        else:

            self.dto.model = self.dto.doc.getModel()

        self.dto.qual_plugin = self.dto.model.getPlugin('qual')

        if self.model.id is not None:
            self.dto.model.setId(self.model.id)
            self.dto.model.setMetaId('meta_' + self.model.id)

        else:
            self.dto.model.setMetaId('meta_model')

        if self.model.name is not None:
            self.dto.model.setName(self.model.name)

    def open(self, mode='r'):

        if mode == 'r':

            return self._open_to_read()

        elif mode == 'w':

            return self._open_to_write()

        else:
            raise ValueError(f'{mode} mode is not recognized. Try one of the following: r, w')

    def parse(self):

        if self.dto is None:
            raise OSError('SBML file is not open')

        if self.dto.id is None:
            raise OSError('SBML file is not open')

        if self.dto.doc is None:
            raise OSError('SBML file is not open')

        if self.dto.model is None:
            raise OSError(f'SBML file is not open')

        if self.dto.qual_plugin is None:
            raise OSError(f'SBML file is not open')

        self.dto.level = self.dto.model.getLevel()
        self.dto.version = self.dto.model.getVersion()

        has_history = self.dto.model.isSetModelHistory()

        if has_history:

            history = self.dto.model.getModelHistory()

            created = None
            if history.isSetCreatedDate():
                created = history.getCreatedDate()

            creators = history.getListCreators()

            self.dto.history = History(data=created, creators=creators)

        # -----------------------------------------------------------------------------
        # Compartments
        # -----------------------------------------------------------------------------
        for compartment in self.dto.model.getListOfCompartments():
            self.dto.compartments[compartment.getIdAttribute()] = CompartmentRecord(id=compartment.getIdAttribute(),
                                                                                    name=compartment.getName())

        # -----------------------------------------------------------------------------
        # Qualitative Species
        # -----------------------------------------------------------------------------
        for qual_species in self.dto.qual_plugin.getListOfQualitativeSpecies():

            identifier = fs_id(qual_species.getIdAttribute(), (F_SPECIE, F_GENE, F_REACTION))
            name = qual_species.getName()

            if not name:
                name = identifier

            aliases = {identifier, name}
            compartment = qual_species.getCompartment()

            constant = qual_species.getConstant()

            # coefficients can be encoded into the notes section. It is not advised though
            notes = qual_species.getNotesString()

            coefficients = self._parse_coefficients_note(notes)

            if coefficients:

                self.warnings.append(partial(sbml_warning, f'Are the {identifier} coefficients encoded in the notes '
                                                           f'section? Coefficients must be hard coded during the '
                                                           f'math nodes of each transition.'))

                maximum_level = max(coefficients)
                minimum_level = min(coefficients)

            else:

                # coefficients might be updated later during transition parsing
                maximum_level = float(qual_species.getMaxLevel())
                minimum_level = 0.0

                coefficients = {minimum_level, maximum_level}

            # finding the active coefficient
            if qual_species.isSetInitialLevel():
                active_coefficient = float(qual_species.getInitialLevel())

            else:

                active_coefficient = self._parse_initial_level_note(notes)

                if active_coefficient is None:
                    self.warnings.append(partial(sbml_warning, f'{identifier} initial level was not found. '
                                                               f'Setting active/initial coefficient to the '
                                                               f'minimum value'))

                    active_coefficient = 0.0

            if minimum_level > active_coefficient > maximum_level:

                raise ValueError(f'Initial level is higher/lower than the minimum/maximum level for the {identifier} '
                                 f'qual species')

            else:

                coefficients.add(active_coefficient)

            variable = VariableRecord(id=identifier,
                                      name=name,
                                      aliases=aliases,
                                      compartment=compartment,
                                      constant=constant,
                                      notes=notes,
                                      coefficients=coefficients,
                                      active_coefficient=active_coefficient)

            self.dto.variables[identifier] = variable

        # -----------------------------------------------------------------------------
        # Interactions/List of Transitions
        # -----------------------------------------------------------------------------

        # iterate over all transitions/interactions
        for transition in self.dto.qual_plugin.getListOfTransitions():

            identifier = f_id(transition.getIdAttribute(), F_TRANSITION)
            name = transition.getName()

            if not name:
                name = identifier

            # -----------------------------------------------------------------------------
            # Regulators/List of Inputs
            # -----------------------------------------------------------------------------

            # regulators can be encoded into the list of inputs by the qualitative species or by an input identifier.
            # A regulator associated with an input identifier must also be associated with the qualitative species
            # and sometimes with the threshold level.
            # The threshold level stands for the coefficient of the regulator in a given function term.
            # In this case, the ASTNode math parsing will have to replace the input id the is enclosed by the tag <ci>
            # by the threshold level.

            # A regulator encoded into the list of inputs using only the qualitative species is a simpler approach.
            # In this case, the threshold level is regularly encoded in a <cn> tag of type integer or rational

            # Nevertheless, the second option is more powerful, as we can set multiple non-integer values

            # to be used during functional term parsing
            regulators = {}
            # to be used during functional term parsing
            inputs_thresholds = {}

            for regulator_input in transition.getListOfInputs():

                if not regulator_input.isSetQualitativeSpecies():
                    raise ValueError(f"Qualitative species for {identifier}'s transition input is not set")

                regulator_id = fs_id(regulator_input.getQualitativeSpecies(), (F_SPECIE, F_GENE, F_REACTION))

                regulator_record = self.dto.variables[regulator_id]

                regulators[regulator_id] = regulator_record

                self.dto.regulators[regulator_id] = regulator_record

                self.variables[regulator_id].add('regulator')

                # input identifier
                if regulator_input.isSetId():

                    input_id = fs_id(regulator_input.getIdAttribute(), (F_SPECIE, F_GENE, F_REACTION, F_TRANSITION))

                else:
                    input_id = regulator_id

                # input threshold level
                if regulator_input.isSetThresholdLevel():
                    input_threshold = regulator_input.getThresholdLevel()

                    inputs_thresholds[input_id] = input_threshold

            if inputs_thresholds:
                self.warnings.append(partial(sbml_warning, f'{identifier} threshold levels detected. It is '
                                                           f'recommended to encode input levels directly in the math '
                                                           f'node'))

            # -----------------------------------------------------------------------------
            # Targets/List of Outputs
            # -----------------------------------------------------------------------------

            # In the case of outputs/targets list, only the qualitative species really matters
            targets = {}
            for target_output in transition.getListOfOutputs():

                if not target_output.isSetQualitativeSpecies():
                    raise ValueError(f"Qualitative species for {identifier}'s transition output is not set")

                target_id = fs_id(target_output.getQualitativeSpecies(), (F_SPECIE, F_GENE, F_REACTION))

                target_record = self.dto.variables[target_id]

                targets[target_id] = target_record

                self.dto.targets[target_id] = target_record

                self.variables[target_id].add('target')

            # -----------------------------------------------------------------------------
            # Function Terms
            # -----------------------------------------------------------------------------

            # A function term specify how regulators determine the state of the targets
            function_terms = {}

            # -----------------------------------------------------------------------------
            # Default term
            # -----------------------------------------------------------------------------

            # The default term does not have an expression, as it just refers to the base coefficient value
            # for the target.
            default_term = transition.getDefaultTerm()

            if default_term is None:

                self.warnings.append(partial(sbml_warning, f'Default function term '
                                                           f'not set for transition {identifier}.'
                                                           f'Setting default function term of zero.'))

                result_level = 0.0

            else:

                result_level = float(default_term.getResultLevel())

                if result_level is None:
                    result_level = 0.0

                    self.warnings.append(partial(sbml_warning, 'Default function term result level not found. '
                                                               'Setting to zero'))

            function_terms[result_level] = FunctionTerm(id=f'{identifier}_{result_level}',
                                                        symbolic=NoneAtom(),
                                                        coefficient=result_level)

            # -----------------------------------------------------------------------------
            # Math based function terms
            # -----------------------------------------------------------------------------

            # The other function terms define how the outputs/targets can take multiple states.
            # These terms encode the conditions and logic for a given target being active or inactive
            for function_term in transition.getListOfFunctionTerms():

                # result level of the function term
                result_level = function_term.getResultLevel()

                if result_level is None:
                    result_level = 0

                    self.warnings.append(partial(sbml_warning, f'Function term result level not found '
                                                               f'for term in transition {identifier}. '
                                                               f'Setting to zero'))

                # a transition can only have one function term per result level
                if result_level in function_terms:
                    raise ValueError(f'Function term {identifier} with {result_level} result level has already '
                                     f'been used by other function terms')

                # a transition can only have function terms to which the result levels are less or equal than the
                # maximum levels of all its outputs. _update_qual_species_coefficients will check these issues
                for target in targets:
                    self._update_qual_species_coefficients(target, result_level)

                # math node of the function term
                if function_term.isSetMath():

                    math_node = function_term.getMath()

                    # It will parse the symbolic algebraic expression in the ASTMath node. The numeric atoms will be
                    # replaced by the state coefficient according to the inputs' threshold levels or qualitative
                    # species' coefficients
                    symbolic = self.parse_math_node(math_node, inputs_thresholds)

                else:

                    self.warnings.append(partial(sbml_warning, f'Function term math node was'
                                                               f'not set for transition {identifier} with '
                                                               f'{result_level}.'
                                                               f'Setting function term equal to the result level.'))

                    symbolic = NoneAtom()

                function_terms[result_level] = FunctionTerm(id=f'{identifier}_{result_level}',
                                                            symbolic=symbolic,
                                                            coefficient=result_level)

            # -----------------------------------------------------------------------------
            # Interaction record
            # -----------------------------------------------------------------------------

            # Setting an interaction per target/output

            for target in targets.values():
                interaction_id = f'{target.id}_interaction'

                interaction_record = VariableRecord(id=interaction_id,
                                                    name=identifier,
                                                    aliases={interaction_id, identifier, name, target.id},
                                                    target=target,
                                                    function_terms=function_terms,
                                                    regulators=regulators)

                self.dto.interactions[interaction_id] = interaction_record

                self.variables[interaction_id].add('interaction')

    def read(self,
             model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
             variables=None):

        if not model:
            model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = self.model

        if not variables:
            variables = self.variables

        if self.dto.id:
            model._id = self.dto.id

        if self.dto.name:
            model.name = self.dto.name

        model.compartments = {compartment.id: compartment.name
                              for compartment in self.dto.compartments.values()}

        processed_vars = set()

        for interaction_id, interaction_record in self.dto.interactions.items():

            target_record = interaction_record.target

            target, warning = target_record.to_variable(model=model,
                                                        types=variables.get(target_record.id, {'target'}),
                                                        name=target_record.name,
                                                        aliases=target_record.aliases,
                                                        coefficients=target_record.coefficients,
                                                        active_coefficient=target_record.active_coefficient)

            if warning:
                self.warnings.append(partial(sbml_warning, warning))

            processed_vars.add(target_record.id)

            regulators_records = interaction_record.regulators

            regulators = {}

            for regulator_id, regulator_record in regulators_records.items():

                regulator, warn = regulator_record.to_variable(model=model,
                                                               types=variables.get(regulator_id, {'regulator'}),
                                                               name=regulator_record.name,
                                                               aliases=regulator_record.aliases,
                                                               coefficients=regulator_record.coefficients,
                                                               active_coefficient=regulator_record.active_coefficient)

                if warn:
                    self.warnings.append(partial(sbml_warning, warn))

                regulators[regulator_id] = regulator

                processed_vars.add(regulator_id)

            regulatory_events = {}

            for func_term in interaction_record.function_terms.values():
                expression_regulators = {symbol.name: regulators[symbol.name]
                                         for symbol in func_term.symbolic.atoms(symbols_only=True)}

                regulatory_events[func_term.coefficient] = Expression(symbolic=func_term.symbolic,
                                                                      variables=expression_regulators)

            interaction, warn = interaction_record.to_variable(model=model,
                                                               types=variables.get(interaction_id, {'interaction'}),
                                                               name=interaction_record.name,
                                                               aliases=interaction_record.aliases,
                                                               target=target,
                                                               regulatory_events=regulatory_events)

            if warn:
                self.warnings.append(partial(sbml_warning, warn))

            model.add(interaction, 'interaction', comprehensive=True)

        if len(processed_vars) != len(self.dto.variables):

            for variable_id, variable_record in self.dto.variables.items():

                if variable_id not in processed_vars:

                    variable, warn = variable_record.to_variable(model=model,
                                                                 types=variables.get(variable_id, {'regulator'}),
                                                                 name=variable_record.name,
                                                                 aliases=variable_record.aliases,
                                                                 coefficients=variable_record.coefficients,
                                                                 active_coefficient=variable_record.active_coefficient)

                    if warn:
                        self.warnings.append(partial(sbml_warning, warn))

                    model.add(variable, 'regulator')

        return model

    def _reverse_f_id(self, variable):

        if variable.is_reaction():

            return f_id(variable.id, F_REACTION_REV)

        elif variable.is_metabolite():

            return f_id(variable.id, F_SPECIE_REV)

        elif variable.is_gene():

            return f_id(variable.id, F_GENE_REV)

        elif variable.is_interaction():

            if variable.target:

                target_id = self._reverse_f_id(variable.target)

                return f_id(target_id, F_TRANSITION_REV)

            else:

                return f_id(variable.id, F_TRANSITION_REV)

        else:

            return variable.id

    @staticmethod
    def _expression_replace(expression_string, replacements):

        for key, value in replacements.items():
            expression_string = expression_string.replace(key, value)

        return expression_string

    def write(self):

        if self.dto is None:
            raise OSError('SBML file is not open')

        if self.dto.doc is None:
            raise OSError('SBML file is not open')

        if self.dto.model is None:
            raise OSError(f'SBML file is not open')

        # -----------------------------------------------------------------------------
        # Compartments
        # -----------------------------------------------------------------------------
        default_compartment = None

        for compartment_id, compartment_name in self.model.compartments.items():
            sbml_compartment = self.dto.model.createCompartment()
            sbml_compartment.setId(compartment_id)
            sbml_compartment.setName(compartment_name)
            sbml_compartment.setConstant(True)

            if default_compartment is None:
                default_compartment = compartment_id

        if default_compartment is None:
            default_compartment = 'e'
            sbml_compartment = self.dto.model.createCompartment()
            sbml_compartment.setId('e')
            sbml_compartment.setName('extracellular')
            sbml_compartment.setConstant(True)

        # -----------------------------------------------------------------------------
        # Variables/Qualitative Species - Targets
        # -----------------------------------------------------------------------------
        processed_regulators = []

        for target in self.model.yield_targets():

            qual_species = self.dto.qual_plugin.createQualitativeSpecies()

            target_id = self._reverse_f_id(target)
            qual_species.setId(target_id)

            if hasattr(target, 'compartment'):

                if target.compartment is None:
                    compartment = default_compartment

                else:
                    compartment = target.compartment

            else:
                compartment = default_compartment

            qual_species.setCompartment(compartment)

            qual_species.setConstant(False)
            qual_species.setName(target.name)

            qual_species.setInitialLevel(int(ceil(target.coefficient.active_coefficient)))
            qual_species.setMaxLevel(int(ceil(target.coefficient.maximum_coefficient)))

            processed_regulators.append(target.id)

        # -----------------------------------------------------------------------------
        # Variables/Qualitative Species - Regulators
        # -----------------------------------------------------------------------------

        for regulator in self.model.yield_regulators():

            if regulator.id not in processed_regulators:

                qual_species = self.dto.qual_plugin.createQualitativeSpecies()

                regulator_id = self._reverse_f_id(regulator)
                qual_species.setId(regulator_id)

                if hasattr(regulator, 'compartment'):

                    if regulator.compartment is None:
                        compartment = default_compartment

                    else:
                        compartment = regulator.compartment

                else:
                    compartment = default_compartment

                qual_species.setCompartment(compartment)

                qual_species.setConstant(False)
                qual_species.setName(regulator.name)

                qual_species.setInitialLevel(int(ceil(regulator.coefficient.active_coefficient)))
                qual_species.setMaxLevel(int(ceil(regulator.coefficient.maximum_coefficient)))

        # -----------------------------------------------------------------------------
        # Transitions - Interactions
        # -----------------------------------------------------------------------------
        for interaction in self.model.yield_interactions():

            transition = self.dto.qual_plugin.createTransition()

            transition_id = self._reverse_f_id(interaction)

            transition.setId(transition_id)
            transition.setName(interaction.name)

            # -----------------------------------------------------------------------------
            # Output
            # -----------------------------------------------------------------------------

            if interaction.target:
                output = transition.createOutput()

                target_id = self._reverse_f_id(interaction.target)

                output.setId(f'{target_id}_out')
                output.setQualitativeSpecies(target_id)

            # -----------------------------------------------------------------------------
            # Inputs
            # -----------------------------------------------------------------------------

            for regulator in interaction.yield_regulators():
                regulator_id = self._reverse_f_id(regulator)

                input_id = f'{regulator_id}_input'
                reg_input = transition.createInput()
                reg_input.setId(input_id)
                reg_input.setQualitativeSpecies(regulator_id)

            # -----------------------------------------------------------------------------
            # Default Term
            # -----------------------------------------------------------------------------
            default_term = transition.createDefaultTerm()
            default_assigned = False

            for coefficient, expression in interaction.regulatory_events.items():

                # default term assignment
                if expression.is_none and not default_assigned:
                    default_term.setResultLevel(int(coefficient))

                    default_assigned = True

                    continue

                # -----------------------------------------------------------------------------
                # Function Terms
                # -----------------------------------------------------------------------------
                function_term = transition.createFunctionTerm()
                function_term.setResultLevel(int(coefficient))

                # operators to be replaced
                replacements = {'&': '&&',
                                '|': '||',
                                '~': '!',
                                '=': '==',
                                '<==': '<=',
                                '>==': '>='}

                # reverse ids for the variables that must be replaced
                symbols_reverse_ids = {variable.id: self._reverse_f_id(variable)
                                       for variable in expression.variables.values()}

                replacements.update(symbols_reverse_ids)

                expression_string = expression.to_string()

                expression_string = self._expression_replace(expression_string=expression_string,
                                                             replacements=replacements)

                if not expression_string:
                    self.warnings.append(partial(sbml_warning, f'Empty expression to be set as a function term in '
                                                               f'transition {transition_id}'))

                    continue

                warning = set_math(transition_id, expression_string, function_term)

                if warning:
                    self.warnings.append(partial(sbml_warning, warning))

        write_sbml_doc(self.io, self.dto.doc)

    def close(self):

        if hasattr(self.io, 'close'):
            self.io.close()

    def clean(self):
        self._dto = None


class MetabolicSBML(Engine):

    def __init__(self, io, config, model=None):
        super().__init__(io, config, model)

    @property
    def model_type(self):
        return 'metabolic'

    @property
    def model(self):

        if self._model is None:

            identifier = self.dto.id

            if not identifier:
                identifier = self.get_identifier()

            return MetabolicModel(identifier=identifier)

        return self._model

    def parse_notes(self, notes):

        notes_store = {}

        for match in pattern_notes.finditer(notes):

            try:
                key, value = match.group("content").split(":", 1)

            except ValueError:
                self.warnings.append(partial(sbml_warning,
                                             "Unexpected content format {}.".format(match.group("content"))))
                continue

            value = value.strip()

            if value:
                notes_store[key.strip()] = value

        return notes_store

    @staticmethod
    def parse_symbols(fbc_association):

        gene_id = f_id(fbc_association.getGeneProduct(), F_GENE)

        symbol = Symbol(gene_id)

        return symbol

    def parse_leaves(self, fbc_association):

        if fbc_association.isGeneProductRef():

            return self.parse_symbols(fbc_association)

        else:

            return

    def _parse_gpa(self, fbc_association):

        if fbc_association.isFbcOr():

            args = [self._parse_gpa(child)
                    for child in fbc_association.getListOfAssociations()]

            return Or(variables=args)

        elif fbc_association.isFbcAnd():

            args = [self._parse_gpa(child)
                    for child in fbc_association.getListOfAssociations()]

            return And(variables=args)

        else:

            return self.parse_leaves(fbc_association)

    def parse_gpa(self, fbc_association):

        """
        Reads and parses a node of type math ASTNode into a boolean algebraic expression.

        :param fbc_association: FBCAssociation
        :return: Expression
        """

        try:

            symbolic = self._parse_gpa(fbc_association)

        except SyntaxError:

            self.warnings.append(partial(expression_warning,
                                         f'{fbc_association} cannot be parsed. Assigning empty expression instead'))

            symbolic = NoneAtom()

        return symbolic

    def get_identifier(self):

        if os.path.exists(self.io):
            _, identifier = os.path.split(self.io)
            return os.path.splitext(identifier)[0]

        return 'model'

    def _open_to_read(self):

        self._dto = DataTransferObject()

        # -----------------------------------------------------------------------------
        # Doc, Model and FBC Plugin
        # -----------------------------------------------------------------------------

        self.dto.doc = get_sbml_doc_to_read(self.io)

        self.dto.model = self.dto.doc.getModel()

        if self.dto.model is None:
            raise OSError(f'{self.io} is not a valid input. Model SBML section is missing. '
                          f'Provide a correct path or file handler')

        identifier = self.dto.model.getIdAttribute()

        if not identifier:
            self.warnings.append(partial(sbml_warning, 'Model identifier is not encoded in the SBML file'))

            identifier = self.get_identifier()

        self.dto.id = identifier
        self.dto.name = self.dto.model.getName()

    def _open_to_write(self):

        self._dto = DataTransferObject()

        # -----------------------------------------------------------------------------
        # Doc, Model and FBC Plugin
        # -----------------------------------------------------------------------------

        self.dto.doc = get_sbml_doc_to_write(self.io,
                                             level=3,
                                             version=1,
                                             packages=('fbc',),
                                             packages_version=(2,),
                                             packages_required=(False,),
                                             sbo_term=True)

        if self.dto.model is None:

            self.dto.model = self.dto.doc.createModel()

        else:

            self.dto.model = self.dto.doc.getModel()

        # fbc plugin is added by get_sbml_doc_to_write
        self.dto.fbc_plugin = self.dto.model.getPlugin('fbc')
        self.dto.fbc_plugin.setStrict(True)

        if self.model.id is not None:
            self.dto.model.setId(self.model.id)
            self.dto.model.setMetaId('meta_' + self.model.id)

        else:
            self.dto.model.setMetaId('meta_model')

        if self.model.name is not None:
            self.dto.model.setName(self.model.name)

    def open(self, mode='r'):

        if mode == 'r':

            return self._open_to_read()

        elif mode == 'w':

            return self._open_to_write()

        else:
            raise ValueError(f'{mode} mode is not recognized. Try one of the following: r, w')

    def parse(self):

        if self.dto is None:
            raise OSError('SBML file is not open')

        if self.dto.id is None:
            raise OSError('SBML file is not open')

        if self.dto.doc is None:
            raise OSError('SBML file is not open')

        if self.dto.model is None:
            raise OSError(f'SBML file is not open')

        self.dto.fbc_plugin = self.dto.model.getPlugin("fbc")

        if not self.dto.fbc_plugin:
            self.warnings.append(partial(sbml_warning, "SBML model does not have fbc plugin"))

        else:
            if not self.dto.fbc_plugin.isSetStrict():
                self.warnings.append(partial(sbml_warning, 'SBML model fbc plugin is not set to strict. It must '
                                                           'fbc:strict="true"'))

            doc_fbc = self.dto.doc.getPlugin("fbc")
            fbc_version = doc_fbc.getPackageVersion()

            # fbc 1 to 2. If fails, an import error is launched
            if fbc_version == 1:
                self.warnings.append(partial(sbml_warning, 'Models should be encoded using fbc version 2. Converting '
                                                           'fbc v1 to fbc v2'))

                convert_fbc(self.dto.doc)

        # -----------------------------------------------------------------------------
        # Model id, name, etc
        # -----------------------------------------------------------------------------

        self.dto.level = self.dto.model.getLevel()
        self.dto.version = self.dto.model.getVersion()

        has_history = self.dto.model.isSetModelHistory()

        if has_history:

            history = self.dto.model.getModelHistory()

            created = None
            if history.isSetCreatedDate():
                created = history.getCreatedDate()

            creators = history.getListCreators()

            self.dto.history = History(created, creators)

        # -----------------------------------------------------------------------------
        # Compartments
        # -----------------------------------------------------------------------------

        for compartment in self.dto.model.getListOfCompartments():
            self.dto.compartments[compartment.getIdAttribute()] = CompartmentRecord(id=compartment.getIdAttribute(),
                                                                                    name=compartment.getName())

        # -----------------------------------------------------------------------------
        # Metabolites and Extracellular metabolites
        # -----------------------------------------------------------------------------

        metabolites = {}
        extracellular_metabolites = {}

        if self.dto.model.getNumSpecies() == 0:
            self.warnings.append(partial(sbml_warning, "SBML model does not have species/metabolites"))

        for met in self.dto.model.getListOfSpecies():

            met_id = f_id(met.getIdAttribute(), F_SPECIE)
            met_name = met.getName()
            met_aliases = {met_id, met_name}

            if not met_name:
                met_name = met_id

            met_notes = met.getNotesString()
            met_annotation = met.getAnnotationString()
            met_compartment = met.getCompartment()

            met_fbc = met.getPlugin("fbc")
            if met_fbc:
                met_charge = met_fbc.getCharge()
                met_formula = met_fbc.getChemicalFormula()
            else:
                if met.isSetCharge():
                    met_charge = met.getCharge()
                    met_formula = None
                else:
                    met_charge = None
                    met_formula = None

            met_record = VariableRecord(id=met_id,
                                        name=met_name,
                                        aliases=met_aliases,
                                        notes=met_notes,
                                        annotation=met_annotation,
                                        compartment=met_compartment,
                                        charge=met_charge,
                                        formula=met_formula)

            self.variables[met_id].add('metabolite')

            if met.getBoundaryCondition() is True:
                # extracellular metabolites
                extracellular_metabolites[met_id] = met_record

            metabolites[met_id] = met_record

            self.dto.variables[met_id] = met_record

        self.dto.metabolites = metabolites
        self.dto.extracellular_metabolites = extracellular_metabolites

        # -----------------------------------------------------------------------------
        # Genes
        # -----------------------------------------------------------------------------

        genes = {}

        # parsing genes encoded in the gene products section of the sbml
        if self.dto.fbc_plugin:

            if self.dto.fbc_plugin.getNumGeneProducts() == 0:
                self.warnings.append(partial(sbml_warning, "SBML model fbc plugin does not have gene products/genes"))

            for gene in self.dto.fbc_plugin.getListOfGeneProducts():

                gene_id = f_id(gene.getIdAttribute(), F_GENE)
                gene_name = gene.getName()

                if not gene_name:
                    gene_name = gene_id

                gene_aliases = {gene_id, gene_name}

                gene_notes = gene.getNotesString()
                gene_annotation = gene.getAnnotationString()

                gene_record = VariableRecord(id=gene_id,
                                             name=gene_name,
                                             aliases=gene_aliases,
                                             notes=gene_notes,
                                             annotation=gene_annotation)

                self.variables[gene_id].add('gene')

                genes[gene_id] = gene_record
                self.dto.variables[gene_id] = gene_record

        self.dto.genes = genes

        # -----------------------------------------------------------------------------
        # Reactions
        # -----------------------------------------------------------------------------

        reactions = {}

        if self.dto.model.getNumReactions() == 0:
            self.warnings.append(partial(sbml_warning, "SBML model does not have reactions"))

        for reaction in self.dto.model.getListOfReactions():

            rxn_id = f_id(reaction.getIdAttribute(), F_REACTION)
            rxn_name = reaction.getName()

            if not rxn_name:
                rxn_name = rxn_id

            rxn_aliases = {rxn_id, rxn_name}

            rxn_notes = self.parse_notes(reaction.getNotesString())
            rxn_annotation = reaction.getAnnotationString()

            # ------------------------------------------------
            # Reaction fbc plugin is of the upmost importance to parse reaction properties
            # ------------------------------------------------

            rxn_fbc = reaction.getPlugin("fbc")

            # ------------------------------------------------
            # Reaction Bounds
            # ------------------------------------------------

            if rxn_fbc:

                rxn_bounds = [ModelConstants.REACTION_LOWER_BOUND, ModelConstants.REACTION_UPPER_BOUND]

                # bounds in fbc parameters section
                bounds_ids = (rxn_fbc.getLowerFluxBound(), rxn_fbc.getUpperFluxBound())

                for i, bound_id in enumerate(bounds_ids):

                    if bound_id:

                        bound_parameter = self.dto.model.getParameter(bound_id)

                        if bound_parameter and bound_parameter.getConstant() and bound_parameter.getValue() is not None:

                            rxn_bounds[i] = bound_parameter.getValue()

                        else:

                            self.warnings.append(
                                partial(sbml_warning, f"Incorrect {bound_parameter} bound for {reaction} reaction. "
                                                      f"Set to default"))

                    else:

                        self.warnings.append(partial(sbml_warning, f"Bound for {reaction} reaction not found in the "
                                                                   f"SBML model fbc plugin. Set to "
                                                                   f"default. Try to set all bounds explicitly on all "
                                                                   f"reactions"))

            elif reaction.isSetKineticLaw():

                self.warnings.append(partial(sbml_warning, f"{reaction} reaction fbc plugin not found. This might "
                                                           f"hinder reaction parsing"))

                self.warnings.append(partial(
                    sbml_warning, f"Bounds have been detected in kinetic laws for {reaction} reaction. "
                                  f"Try to set all bounds explicitly on all reactions using the fbc plugin, as mewpy "
                                  f"can miss sometimes kinetic laws"))

                # bounds encoded in the kinetic law. Not advised
                kinetic_law = reaction.getKineticLaw()

                rxn_bounds = [ModelConstants.REACTION_LOWER_BOUND, ModelConstants.REACTION_UPPER_BOUND]

                # parameters of the kinetic law
                kinetic_parameters = ('LOWER_BOUND', 'UPPER_BOUND')

                for i, kinetic_parameter in enumerate(kinetic_parameters):

                    bound = kinetic_law.getParameter(kinetic_parameter)

                    if bound:
                        rxn_bounds[i] = bound.getValue()

                    else:
                        self.warnings.append(partial(sbml_warning,
                                                     f"{kinetic_parameter} has not been detected. "
                                                     f"{kinetic_parameter} has been set to default. "
                                                     f"Try to set all bounds explicitly on all reactions using "
                                                     f"the fbc plugin"))

            else:

                self.warnings.append(partial(sbml_warning, f"{reaction} reaction fbc plugin not found. This might "
                                                           f"hinder reaction parsing"))

                self.warnings.append(partial(sbml_warning,
                                             f"Bounds have not been detected. Bounds have been set to default. Try to "
                                             f"set all bounds explicitly on all reactions using the fbc plugin"))

                rxn_bounds = [ModelConstants.REACTION_LOWER_BOUND, ModelConstants.REACTION_UPPER_BOUND]

            # ------------------------------------------------
            # Reaction stoichiometry and metabolites
            # ------------------------------------------------

            rxn_stoichiometry = {}
            rxn_metabolites = {}
            rxn_products = {}
            rxn_reactants = {}

            # reactants
            for reactant in reaction.getListOfReactants():

                reactant_id = f_id(reactant.getSpecies(), F_SPECIE)

                if reactant_id not in self.dto.metabolites:
                    raise ValueError(f"{reactant_id} reactant of {reaction} reaction "
                                     f"is not listed as specie in the SBML Model.")

                reactant_record = self.dto.metabolites[reactant_id]

                rxn_metabolites[reactant_id] = reactant_record
                rxn_reactants[reactant_id] = reactant_record
                rxn_stoichiometry[reactant_id] = -reactant.getStoichiometry()

            # products
            for product in reaction.getListOfProducts():

                product_id = f_id(product.getSpecies(), F_SPECIE)

                if product_id not in self.dto.metabolites:
                    raise ValueError(f"{product_id} product of {reaction} reaction "
                                     f"is not listed as specie in the SBML Model.")

                product_record = self.dto.metabolites[product_id]

                rxn_metabolites[product_id] = product_record
                rxn_products[product_id] = product_record
                rxn_stoichiometry[product_id] = product.getStoichiometry()

            # ------------------------------------------------
            # Reaction Gene-Protein-Reaction rule
            # ------------------------------------------------

            symbolic = NoneAtom()

            if rxn_fbc:

                # if the gpr is encoded in the reaction fbc plugin.
                # It is advised.

                gpa = rxn_fbc.getGeneProductAssociation()

                if gpa is not None:
                    symbolic = self.parse_gpa(gpa.getAssociation())

            else:

                self.warnings.append(partial(sbml_warning,
                                             "Please use fbc plugin fbc:gpr to encode gprs in the future, "
                                             "as parsing gprs from notes might be troublesome"))

                # Else the gpr parsing tries to find the gpr rule (string) within the notes

                gpr_rule = rxn_notes.get('GENE ASSOCIATION',
                                         rxn_notes.get('GENE_ASSOCIATION', None))

                if gpr_rule is None:

                    self.warnings.append(partial(sbml_warning,
                                                 "GPR was not found within the reaction's notes section"))

                else:

                    gpr_rule = ' '.join(f_id(child, F_GENE) for child in gpr_rule.split(' '))

                    symbolic, warning = build_symbolic(expression=gpr_rule)

                    if warning:
                        self.warnings.append(partial(expression_warning, warning))

            rxn_genes = {}
            for symbol in symbolic.atoms(symbols_only=True):

                if rxn_fbc:

                    if symbol.name in self.dto.genes:

                        rxn_genes[symbol.name] = self.dto.genes[symbol.name]

                        continue

                    else:

                        self.warnings.append(partial(sbml_warning,
                                                     f'{symbol.name} is not listed in the SBML model fbc plugin'))

                gene_record = VariableRecord(id=symbol.name,
                                             name=symbol.name,
                                             aliases={symbol.name, symbol.value})

                self.variables[symbol.name].add('gene')

                self.dto.variables[symbol.name] = gene_record
                self.dto.genes[symbol.name] = gene_record

                rxn_genes[symbol.name] = gene_record

            # -----------------------------------------------------------------------------
            # GPR Function term
            # -----------------------------------------------------------------------------
            function_term = FunctionTerm(id='gpr_term', symbolic=symbolic, coefficient=1)

            # ------------------------------------------------
            # Building reaction record and assigning to containers
            # ------------------------------------------------

            reaction_record = VariableRecord(id=rxn_id,
                                             name=rxn_name,
                                             aliases=rxn_aliases,
                                             notes=rxn_notes,
                                             annotation=rxn_annotation,
                                             bounds=tuple(rxn_bounds),
                                             genes=rxn_genes,
                                             gpr=function_term,
                                             metabolites=rxn_metabolites,
                                             products=rxn_products,
                                             reactants=rxn_reactants,
                                             stoichiometry=rxn_stoichiometry)

            self.variables[rxn_id].add('reaction')

            reactions[rxn_id] = reaction_record

            self.dto.variables[rxn_id] = reaction_record

        self.dto.reactions = reactions

        # Some transport reactions have fictitious extracellular metabolites, encoded with the prefix _b,
        # such as M_a_b. In this case, a exchange reaction is created for each metabolite

        extracellular_reactions = {}

        for extracellular_met in self.dto.extracellular_metabolites.values():
            # ------------------------------------------------
            # Building reaction record and assigning to containers
            # ------------------------------------------------

            reaction_record = VariableRecord(id=f'EX_{extracellular_met.id}',
                                             name=f'EX_{extracellular_met.id}',
                                             bounds=(ModelConstants.REACTION_LOWER_BOUND,
                                                     ModelConstants.REACTION_UPPER_BOUND),
                                             metabolites={extracellular_met.id: extracellular_met},
                                             reactants={extracellular_met.id: extracellular_met},
                                             stoichiometry={extracellular_met.id: -1})

            self.variables[f'EX_{extracellular_met.id}'].add('reaction')

            extracellular_reactions[f'EX_{extracellular_met.id}'] = reaction_record

            self.dto.variables[f'EX_{extracellular_met.id}'] = reaction_record

            self.warnings.append(partial(sbml_warning,
                                         f'EX_{extracellular_met.id} reaction added '
                                         f'for metabolite {extracellular_met.id}'))

        self.dto.extracellular_reactions = extracellular_reactions
        self.dto.reactions.update(extracellular_reactions)

        # -----------------------------------------------------------------------------
        # Objective
        # -----------------------------------------------------------------------------

        model_objective = {}

        if self.dto.fbc_plugin:

            objectives = self.dto.fbc_plugin.getListOfObjectives()

            if objectives is None:
                self.warnings.append(partial(sbml_warning, "listOfObjectives element not found"))

            elif objectives.size() == 0:
                self.warnings.append(partial(sbml_warning, "No objective in listOfObjectives"))

            elif not objectives.getActiveObjective():
                self.warnings.append(partial(sbml_warning, "No active objective in listOfObjectives"))

            else:
                active_objectives = objectives.getActiveObjective()
                objective = self.dto.fbc_plugin.getObjective(active_objectives)

                direction = objective.getType()

                for flux_objective in objective.getListOfFluxObjectives():

                    flux_objective_id = f_id(flux_objective.getReaction(), F_REACTION)

                    objective_rxn = self.dto.reactions.get(flux_objective_id, None)

                    if objective_rxn is None:
                        self.warnings.append(partial(sbml_warning,
                                                     f"Objective {flux_objective_id} reaction "
                                                     f"not found in the SBML model"))

                        continue

                    coef = flux_objective.getCoefficient()

                    if direction == 'minimize':
                        coef = -coef

                    model_objective[flux_objective_id] = coef

        else:
            self.warnings.append(partial(sbml_warning,
                                         f"Objective might be encoded in kinetic laws of a given reaction. However, "
                                         f"mewpy does not handle kinetic laws. The objective has not been set. Try "
                                         f"to set the objective explicitly on the fbc plugin"))

        if len(model_objective) == 0:
            self.warnings.append(partial(sbml_warning, "No objective found for the model"))

        self.dto.objective = model_objective

    def read(self,
             model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
             variables=None):

        if not model:
            model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = self.model

        if not variables:
            variables = self.variables

        if self.dto.id:
            model._id = self.dto.id

        if self.dto.name:
            model.name = self.dto.name

        model.compartments = {compartment.id: compartment.name
                              for compartment in self.dto.compartments.values()}

        processed_metabolites = set()
        processed_genes = set()

        for rxn_id, rxn_record in self.dto.reactions.items():

            genes = {}

            for gene_id, gene_record in rxn_record.genes.items():

                gene, warning = gene_record.to_variable(model=model,
                                                        types=variables.get(gene_id, {'gene'}),
                                                        name=gene_record.name,
                                                        aliases=gene_record.aliases)

                if warning:
                    self.warnings.append(partial(sbml_warning, warning))

                genes[gene_id] = gene

                processed_genes.add(gene_id)

            stoichiometry = {}

            for met_id, met_record in rxn_record.metabolites.items():

                met, warning = met_record.to_variable(model=model,
                                                      types=variables.get(met_id, {'metabolite'}),
                                                      name=met_record.name,
                                                      aliases=met_record.aliases,
                                                      compartment=met_record.compartment,
                                                      charge=met_record.charge,
                                                      formula=met_record.formula)

                if warning:
                    self.warnings.append(partial(sbml_warning, warning))

                coef = rxn_record.stoichiometry[met_id]

                stoichiometry[met] = coef

                processed_metabolites.add(met_id)

            gpr = Expression(symbolic=rxn_record.gpr.symbolic, variables=genes)

            rxn, warning = rxn_record.to_variable(model=model,
                                                  types=variables.get(rxn_id, {'reaction'}),
                                                  bounds=rxn_record.bounds,
                                                  gpr=gpr,
                                                  stoichiometry=stoichiometry)

            if warning:
                self.warnings.append(partial(sbml_warning, warning))

            model.add(rxn, 'reaction')

        to_append = []

        for met_id, met_record in self.dto.metabolites.items():

            if met_id not in processed_metabolites:
                met, warning = met_record.to_variable(model=model,
                                                      types=variables.get(met_id, {'metabolite'}),
                                                      name=met_record.name,
                                                      aliases=met_record.aliases,
                                                      compartment=met_record.compartment,
                                                      charge=met_record.charge,
                                                      formula=met_record.formula)

                if warning:
                    self.warnings.append(partial(sbml_warning, warning))

                to_append.append(met)

        for gene_id, gene_record in self.dto.genes.items():

            if gene_id not in processed_genes:
                gene, warning = gene_record.to_variable(model=model,
                                                        types=variables.get(gene_id, {'gene'}),
                                                        name=gene_record.name,
                                                        aliases=gene_record.aliases)

                if warning:
                    self.warnings.append(partial(sbml_warning, warning))

                to_append.append(gene)

        model.add(to_append)

        model.objective = self.dto.objective

        return model

    def write(self):

        if self.dto is None:
            raise OSError('SBML file is not open')

        if self.dto.doc is None:
            raise OSError('SBML file is not open')

        if self.dto.model is None:
            raise OSError(f'SBML file is not open')

        # -----------------------------------------------------------------------------
        # Units
        # -----------------------------------------------------------------------------
        units = self.config.get('units', False)
        unit_definition = None
        if units:

            unit_definition = self.dto.model.createUnitDefinition()
            unit_definition.setId(UNIT_ID)

            for unit in UNITS:
                _unit = unit_definition.createUnit()
                _unit.setKind(unit.kind)
                _unit.setExponent(unit.exponent)
                _unit.setScale(unit.scale)
                _unit.setMultiplier(unit.multiplier)

        # -----------------------------------------------------------------------------
        # Constants and parameters
        # -----------------------------------------------------------------------------
        add_sbml_parameter(sbml_model=self.dto.model,
                           parameter_id=LOWER_BOUND_ID,
                           value=ModelConstants.REACTION_LOWER_BOUND,
                           constant=True,
                           sbo=SBO_DEFAULT_FLUX_BOUND)

        add_sbml_parameter(sbml_model=self.dto.model,
                           parameter_id=UPPER_BOUND_ID,
                           value=ModelConstants.REACTION_UPPER_BOUND,
                           constant=True,
                           sbo=SBO_DEFAULT_FLUX_BOUND)

        add_sbml_parameter(sbml_model=self.dto.model,
                           parameter_id=ZERO_BOUND_ID,
                           value=0,
                           constant=True,
                           sbo=SBO_DEFAULT_FLUX_BOUND)

        add_sbml_parameter(sbml_model=self.dto.model,
                           parameter_id=BOUND_MINUS_INF,
                           value=-float("Inf"),
                           constant=True,
                           sbo=SBO_DEFAULT_FLUX_BOUND)

        add_sbml_parameter(sbml_model=self.dto.model,
                           parameter_id=BOUND_PLUS_INF,
                           value=float("Inf"),
                           constant=True,
                           sbo=SBO_DEFAULT_FLUX_BOUND)

        # -----------------------------------------------------------------------------
        # Compartments
        # -----------------------------------------------------------------------------
        for compartment_id, compartment_name in self.model.compartments.items():
            sbml_compartment = self.dto.model.createCompartment()
            sbml_compartment.setId(compartment_id)
            sbml_compartment.setName(compartment_name)
            sbml_compartment.setConstant(True)

        # -----------------------------------------------------------------------------
        # Metabolites
        # -----------------------------------------------------------------------------
        for metabolite in self.model.yield_metabolites():
            sbml_specie = self.dto.model.createSpecies()

            sbml_specie.setId(f_id(metabolite.id, F_SPECIE_REV))
            sbml_specie.setConstant(False)
            sbml_specie.setBoundaryCondition(False)
            sbml_specie.setHasOnlySubstanceUnits(False)
            sbml_specie.setName(metabolite.name)
            sbml_specie.setCompartment(metabolite.compartment)
            specie_fbc = sbml_specie.getPlugin('fbc')
            specie_fbc.setCharge(metabolite.charge)
            specie_fbc.setChemicalFormula(metabolite.formula)

        # -----------------------------------------------------------------------------
        # Genes
        # -----------------------------------------------------------------------------
        for gene in self.model.yield_genes():
            gene_product = self.dto.fbc_plugin.createGeneProduct()

            gene_product.setId(f_id(gene.id, F_GENE_REV))
            gene_product.setName(gene.name)
            gene_product.setLabel(gene.id)

        # -----------------------------------------------------------------------------
        # Objective
        # -----------------------------------------------------------------------------
        objective = self.dto.fbc_plugin.createObjective()
        objective.setId('obj')
        objective.setType('maximize')
        self.dto.fbc_plugin.setActiveObjectiveId('obj')

        for reaction, coefficient in self.model.objective.items():
            flux_objective = objective.createFluxObjective()
            flux_objective.setReaction(f_id(reaction.id, F_REACTION_REV))
            flux_objective.setCoefficient(coefficient)

        # -----------------------------------------------------------------------------
        # Reactions
        # -----------------------------------------------------------------------------
        for reaction in self.model.yield_reactions():

            sbml_reaction = self.dto.model.createReaction()
            sbml_reaction.setId(f_id(reaction.id, F_REACTION_REV))
            sbml_reaction.setName(reaction.name)
            sbml_reaction.setFast(False)
            sbml_reaction.setReversible(reaction.reversibility)

            # -----------------------------------------------------------------------------
            # Stoichiometry
            # -----------------------------------------------------------------------------
            for metabolite, st in reaction.stoichiometry.items():

                if st < 0:
                    sbml_reactant = sbml_reaction.createReactant()
                    sbml_reactant.setSpecies(f_id(metabolite.id, F_SPECIE_REV))
                    sbml_reactant.setStoichiometry(-st)
                    sbml_reactant.setConstant(True)

                else:
                    sbml_product = sbml_reaction.createProduct()
                    sbml_product.setSpecies(f_id(metabolite.id, F_SPECIE_REV))
                    sbml_product.setStoichiometry(st)
                    sbml_product.setConstant(True)

            # -----------------------------------------------------------------------------
            # Bounds
            # -----------------------------------------------------------------------------
            sbml_rxn_fbc = sbml_reaction.getPlugin('fbc')

            lb_parameter_id = get_sbml_lb_id(sbml_model=self.dto.model,
                                             reaction=reaction,
                                             unit_definition=unit_definition)

            if lb_parameter_id is None:
                lb_parameter_id = f'{f_id(reaction.id, F_REACTION_REV)}_lb'

                add_sbml_parameter(sbml_model=self.dto.model,
                                   parameter_id=lb_parameter_id,
                                   value=reaction.lower_bound,
                                   sbo=SBO_FLUX_BOUND,
                                   constant=True,
                                   unit_definition=unit_definition)

            sbml_rxn_fbc.setLowerFluxBound(lb_parameter_id)

            ub_parameter_id = get_sbml_ub_id(sbml_model=self.dto.model,
                                             reaction=reaction,
                                             unit_definition=unit_definition)

            if ub_parameter_id is None:
                ub_parameter_id = f'{f_id(reaction.id, F_REACTION_REV)}_lb'

                add_sbml_parameter(sbml_model=self.dto.model,
                                   parameter_id=ub_parameter_id,
                                   value=reaction.upper_bound,
                                   sbo=SBO_FLUX_BOUND,
                                   constant=True,
                                   unit_definition=unit_definition)

            sbml_rxn_fbc.setUpperFluxBound(ub_parameter_id)

            # -----------------------------------------------------------------------------
            # Reaction Gene-Protein-Reaction rule
            # -----------------------------------------------------------------------------
            set_gpr(self, sbml_warning, reaction, sbml_rxn_fbc)

        write_sbml_doc(self.io, self.dto.doc)

    def close(self):

        if hasattr(self.io, 'close'):
            self.io.close()

    def clean(self):
        self._dto = None
