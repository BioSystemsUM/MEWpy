import os
from functools import partial
from typing import Union, TYPE_CHECKING

from mewpy.io.dto import DataTransferObject, VariableRecord, History, FunctionTerm, CompartmentRecord
from mewpy.germ.algebra import Expression, Symbol, Or, And, NoneAtom
from mewpy.germ.models import RegulatoryModel, MetabolicModel
from mewpy.util.constants import ModelConstants
from .engine import Engine
from .engines_utils import (build_symbolic,
                            pattern_notes,
                            f_id, F_GENE, F_SPECIE, F_REACTION, F_SPECIE_REV, F_GENE_REV, F_REACTION_REV, convert_fbc,
                            get_sbml_doc_to_write, get_sbml_doc_to_read,
                            UNIT_ID, UNITS,
                            add_sbml_parameter, LOWER_BOUND_ID, UPPER_BOUND_ID, ZERO_BOUND_ID,
                            BOUND_MINUS_INF, BOUND_PLUS_INF,
                            SBO_DEFAULT_FLUX_BOUND, SBO_FLUX_BOUND,
                            get_sbml_lb_id, get_sbml_ub_id, write_sbml_doc,
                            set_gpr,
                            expression_warning, sbml_warning)

if TYPE_CHECKING:
    from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel


class MetabolicSBML(Engine):
    """
    SBML engine for metabolic models.
    Credits to the cobrapy team (https://github.com/opencobra/cobrapy), as this SBML engine is heavily inspired in their
    SBML parser.
    """

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

            model.add(rxn)

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

        model.add(*to_append)

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
