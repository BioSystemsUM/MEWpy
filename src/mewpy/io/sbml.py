import os
from collections import OrderedDict
from math import isinf, isnan

from libsbml import AssignmentRule, SBMLReader

from ..model.kinetic import ODEModel, Compartment, Metabolite, KineticReaction, Rule


def load_sbml(filename):
    """ Loads an SBML file.

    :param filename: SBML file path, str.
    :returns: SBMLModel

    """

    if not os.path.exists(filename):
        raise IOError("Model file was not found")

    reader = SBMLReader()
    document = reader.readSBML(str(filename))
    sbml_model = document.getModel()

    if sbml_model is None:
        document.printErrors()
        raise IOError(f'Failed to load model {filename}.')

    return sbml_model


def load_ODEModel(filename):
    sbml_model = load_sbml(filename)
    ode_model = ODEModel(sbml_model.getId())
    _load_compartments(sbml_model, ode_model)
    _load_metabolites(sbml_model, ode_model)
    # adds constants and boundaries to metabolites
    for species in sbml_model.getListOfSpecies():
        m = ode_model.metabolites[species.getId()]
        if not hasattr(m, "constant"):
            m.constant = species.getConstant()
        if not hasattr(m, "boundary"):
            m.boundary = species.getBoundaryCondition()
    # load_reactions(sbml_model, ode_model)
    _load_concentrations(sbml_model, ode_model)
    _load_global_parameters(sbml_model, ode_model)
    _load_ratelaws(sbml_model, ode_model)
    _load_assignment_rules(sbml_model, ode_model)
    return ode_model


def extract_metadata(sbml_elem, elem):

    sboterm = sbml_elem.getSBOTermID()
    if sboterm:
        elem.metadata['SBOTerm'] = sboterm

    notes = sbml_elem.getNotes()
    if notes:
        recursive_node_parser(notes, elem.metadata)

    annotation = sbml_elem.getAnnotationString()
    if annotation:
        elem.metadata['XMLAnnotation'] = annotation


def recursive_node_parser(node, cache):
    node_data = node.getCharacters()
    if ':' in node_data:
        key, value = node_data.split(':', 1)
        cache[key.strip()] = value.strip()

    for i in range(node.getNumChildren()):
        recursive_node_parser(node.getChild(i), cache)


def _load_compartments(sbml_model, model):
    for compartment in sbml_model.getListOfCompartments():
        model.add_compartment(_load_compartment(compartment))


def _load_compartment(compartment):
    size = compartment.getSize()
    if isnan(size) or isinf(size):
        size = 1.0

    comp = Compartment(compartment.getId(), compartment.getName(), False, size)
    extract_metadata(compartment, comp)
    return comp


def _load_metabolites(sbml_model, model):
    for species in sbml_model.getListOfSpecies():
        model.add_metabolite(_load_metabolite(species))


def _load_metabolite(species):
    metabolite = Metabolite(species.getId(), species.getName(), species.getCompartment())
    try:
        fbc_species = species.getPlugin('fbc')
        if fbc_species.isSetChemicalFormula():
            formula = fbc_species.getChemicalFormula()
            metabolite.metadata['FORMULA'] = formula

        if fbc_species.isSetCharge():
            charge = fbc_species.getCharge()
            metabolite.metadata['CHARGE'] = str(charge)
    except Exception:
        pass
    extract_metadata(species, metabolite)
    return metabolite


def _load_concentrations(sbml_model, odemodel):
    for species in sbml_model.getListOfSpecies():
        odemodel.set_concentration(species.getId(), species.getInitialConcentration())


def _load_global_parameters(sbml_model, odemodel):
    for parameter in sbml_model.getListOfParameters():
        odemodel.set_global_parameter(parameter.getId(), parameter.getValue(), parameter.getConstant())


def _load_ratelaws(sbml_model, odemodel):
    for reaction in sbml_model.getListOfReactions():
        formula = reaction.getKineticLaw().getFormula()
        parameters = OrderedDict()
        modifiers = []
        stoichiometry = OrderedDict()

        for reactant in reaction.getListOfReactants():
            m_id = reactant.getSpecies()
            coeff = -reactant.getStoichiometry()

            if m_id not in stoichiometry:
                stoichiometry[m_id] = coeff
            else:
                stoichiometry[m_id] += coeff

        for product in reaction.getListOfProducts():
            m_id = product.getSpecies()
            coeff = product.getStoichiometry()
            if m_id not in stoichiometry:
                stoichiometry[m_id] = coeff
            else:
                stoichiometry[m_id] += coeff
            if stoichiometry[m_id] == 0.0:
                del stoichiometry[m_id]

        for parameter in reaction.getKineticLaw().getListOfParameters():
            parameters[parameter.getId()] = parameter.getValue()
        for modifier in reaction.getListOfModifiers():
            m_id = modifier.getSpecies()
            modifiers.append(m_id)

        law = KineticReaction(reaction.getId(), formula, name=reaction.getName(), 
                              stoichiometry=stoichiometry,
                              parameters=parameters, modifiers=modifiers, 
                              reversible=reaction.getReversible())
        odemodel.set_ratelaw(reaction.getId(), law)


def _load_assignment_rules(sbml_model, odemodel):
    for rule in sbml_model.getListOfRules():
        if isinstance(rule, AssignmentRule):
            r_id = rule.getVariable()
            odemodel.set_assignment_rule(r_id, Rule(r_id, rule.getFormula()))
