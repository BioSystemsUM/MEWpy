from collections import OrderedDict
from libsbml import AssignmentRule, SBMLReader
import os
from ..model.kinetic import ODEModel, KineticReaction, Rule


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
    from reframed.io.sbml import load_compartments, load_metabolites, load_reactions
    ode_model = ODEModel(sbml_model.getId())
    load_compartments(sbml_model, ode_model)
    load_metabolites(sbml_model, ode_model)
    # adds constants and boundaries to metabolites
    for species in sbml_model.getListOfSpecies():
        m = ode_model.metabolites[species.getId()]
        if not hasattr(m, "constant"):
            m.constant = species.getConstant()
        if not hasattr(m, "boundary"):
            m.boundary = species.getBoundaryCondition()
    load_reactions(sbml_model, ode_model)
    _load_concentrations(sbml_model, ode_model)
    _load_global_parameters(sbml_model, ode_model)
    _load_ratelaws(sbml_model, ode_model)
    _load_assignment_rules(sbml_model, ode_model)
    return ode_model


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
        substrates = []
        products = []
        modifiers = []
        for parameter in reaction.getKineticLaw().getListOfParameters():
            parameters[parameter.getId()] = parameter.getValue()
        for reactant in reaction.getListOfReactants():
            m_id = reactant.getSpecies()
            substrates.append(m_id)
        for product in reaction.getListOfProducts():
            m_id = product.getSpecies()
            products.append(m_id)
        for modifier in reaction.getListOfModifiers():
            m_id = modifier.getSpecies()
            modifiers.append(m_id)
        law = KineticReaction(reaction.getId(), formula, parameters=parameters,
                              substrates=substrates, products=products, modifiers=modifiers)
        odemodel.set_ratelaw(reaction.getId(), law)


def _load_assignment_rules(sbml_model, odemodel):
    for rule in sbml_model.getListOfRules():
        if isinstance(rule, AssignmentRule):
            r_id = rule.getVariable()
            odemodel.set_assignment_rule(r_id, Rule(r_id, rule.getFormula()))
