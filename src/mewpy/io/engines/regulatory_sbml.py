import os
from functools import partial
from typing import Union, TYPE_CHECKING

from mewpy.io.dto import DataTransferObject, VariableRecord, History, FunctionTerm, CompartmentRecord
from mewpy.germ.algebra import Expression, Symbol, NoneAtom, Float
from mewpy.germ.models import RegulatoryModel, MetabolicModel
from mewpy.util.constants import ModelConstants
from .engine import Engine
from .engines_utils import (ASTNODE_BOOLEAN_VALUES, ASTNODE_RELATIONAL_OPERATORS,
                            ASTNODE_NAME, ASTNODE_BOOLEAN_OPERATORS, ASTNODE_VALUES,
                            f_id, fs_id,
                            F_GENE, F_SPECIE, F_REACTION, F_SPECIE_REV, F_GENE_REV, F_REACTION_REV, F_TRANSITION,
                            F_TRANSITION_REV,
                            get_sbml_doc_to_write, get_sbml_doc_to_read,
                            write_sbml_doc,
                            set_math, expression_warning, sbml_warning)

if TYPE_CHECKING:
    from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel


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
                _, level = row.split(':')

                # st = st.strip('<p>')
                # level = level.strip('</p>')
                #
                # if 'state' in st.lower():
                #     st = int(st.lower().replace('state', '').strip())

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

            # finding the default coefficient
            if qual_species.isSetInitialLevel():
                active_coefficient = float(qual_species.getInitialLevel())

            else:

                active_coefficient = self._parse_initial_level_note(notes)

                if active_coefficient is None:
                    self.warnings.append(partial(sbml_warning, f'{identifier} initial level was not found. '
                                                               f'Setting default/initial coefficient to the '
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
                                                        coefficients=target_record.coefficients)

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
                                                               coefficients=regulator_record.coefficients)

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

            model.add(interaction, comprehensive=True)

        if len(processed_vars) != len(self.dto.variables):

            for variable_id, variable_record in self.dto.variables.items():

                if variable_id not in processed_vars:

                    variable, warn = variable_record.to_variable(model=model,
                                                                 types=variables.get(variable_id, {'regulator'}),
                                                                 name=variable_record.name,
                                                                 aliases=variable_record.aliases,
                                                                 coefficients=variable_record.coefficients)

                    if warn:
                        self.warnings.append(partial(sbml_warning, warn))

                    model.add(variable)

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

            qual_species.setInitialLevel(int(round(min(target.coefficients))))
            qual_species.setMaxLevel(int(round(max(target.coefficients))))

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

                qual_species.setInitialLevel(int(round(min(regulator.coefficients))))
                qual_species.setMaxLevel(int(round(max(regulator.coefficients))))

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
