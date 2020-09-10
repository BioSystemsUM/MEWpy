import libsbml
import warnings
from sympy.logic.boolalg import And, Or, Not
from sympy.core.relational import StrictGreaterThan, StrictLessThan, GreaterThan, LessThan, Eq
from sympy import sympify, symbols
import os

DIR = os.path.dirname(os.path.realpath(__file__))
PATH = os.path.join(DIR, '../../../examples/models/regulation/e_coli_core.xml')

ASTNODE_BOOLEAN_OPERATORS = {
    libsbml.AST_LOGICAL_AND: And,
    libsbml.AST_LOGICAL_OR: Or,
    libsbml.AST_LOGICAL_NOT: Not}

ASTNODE_RELATIONAL_OPERATORS = {
    libsbml.AST_RELATIONAL_EQ: Eq,
    libsbml.AST_RELATIONAL_LEQ: LessThan,
    libsbml.AST_RELATIONAL_GEQ: GreaterThan,
    libsbml.AST_RELATIONAL_LT: StrictLessThan,
    libsbml.AST_RELATIONAL_GT: StrictGreaterThan
}

ASTNODE_VALUES = {
    libsbml.AST_INTEGER: lambda x: int(x.getValue()),
    libsbml.AST_REAL: lambda x: x.getValue(),
    libsbml.AST_REAL_E: lambda x: x.getValue(),
    libsbml.AST_RATIONAL: lambda x: x.getValue(),
    libsbml.AST_CONSTANT_FALSE: 0,
    libsbml.AST_CONSTANT_TRUE: 1
}

ASTNODE_NAME = libsbml.AST_NAME


def readMathNode(ast_node):

    """
    Reads and parses a node of type math ASTNode into a boolean sympy expression.

    :param ast_node: ASTNode, math ASTNode from ft.getMath()
    :return: sympy expression object or None
    """

    def _process_rule(node):

        node_type = node.getType()

        if node_type in ASTNODE_BOOLEAN_OPERATORS:

            # TODO: If op == Not and node.getNumChildren() > 1, something is wrong in the sbml file

            op = ASTNODE_BOOLEAN_OPERATORS[node_type]
            args = [_process_rule(node.getChild(child)) for child in range(node.getNumChildren())]

            # noinspection PyArgumentList
            sp = op(*args, evaluate=False)

            return sp

        elif node_type in ASTNODE_RELATIONAL_OPERATORS:

            # TODO: If node.getNumChildren() > 2, something is wrong in the sbml file
            # TODO: Process equalities

            op = ASTNODE_RELATIONAL_OPERATORS[node_type]
            args = [_process_rule(node.getChild(child)) for child in range(node.getNumChildren())]

            # noinspection PyArgumentList
            sp = op(*args, evaluate=False)

            return sp

        elif node_type == ASTNODE_NAME:

            # TODO: Verify if the name appears in the list of inputs and species
            # TODO: see if it is possible to change the name of the node. If not, create a aliases_map to find the
            #  variable

            return symbols(node.getName())

        elif node_type in ASTNODE_VALUES:

            # TODO: It should be a state for the previous variable, zero (which means not variable), or one (which
            #  means that the variable directly influences the outcome of the outputs)

            return sympify(ASTNODE_VALUES[node_type](node))

        else:

            return

    try:
        return _process_rule(ast_node)
    except Exception as e:
        print(e)
        print('Moving on')
        return

doc = libsbml.readSBML(PATH)

model = doc.getModel()

model_qual = model.getPlugin("qual")

model_id = model.getIdAttribute()

valid_id = libsbml.SyntaxChecker.isValidSBMLSId(model_id)

model_name = model.getName()

level = model.getLevel()

version = model.getVersion()

has_history = model.isSetModelHistory()

if has_history:

    history = model.getModelHistory()

    if history.isSetCreatedDate():
        created = history.getCreatedDate()

    creators = history.getListCreators()

# TODO: parsing compartments to create metabolite aliases with the respective compartment annotation. Either way,
#  the ID reflect the compartment and be equal to the one in the CBM model
compartments = {}
for compartment in model.getListOfCompartments():
    compartments[compartment.getIdAttribute()] = compartment.getName()

# TODO: parsing qual species IDs
# TODO: parsing states in notes

qual_species = {}
for qs in model_qual.getListOfQualitativeSpecies():
    qs_id = qs.getIdAttribute()
    name = qs.getName()
    if not name:
        name = qs_id
    compartment = qs.getCompartment()
    constant = qs.getConstant()
    initialLevel = qs.getInitialLevel()
    maxlevel = qs.getMaxLevel()
    notes = qs.getNotesString()
    qual_species[qs.getIdAttribute()] = (qs_id, name, compartment, constant, initialLevel, maxlevel, notes)

interactions = {}

for qt in model_qual.getListOfTransitions():

    qt_id = qt.getIdAttribute()
    name = qt.getName()
    if not name:
        name = qt_id

    # regulators
    inputs = {}
    for reg in qt.getListOfInputs():
        reg_id = reg.getIdAttribute()
        reg_name = reg.getName()
        specie_id = qual_species[reg.getQualitativeSpecies()][0]
        inputs[specie_id] = (reg_id, reg_name, specie_id)

    # TODO: if a given variable is the output of more than one transition, these transitions must be combined into a
    #  single one. Alternatively, launch a warning to repair the sbml as we do not handle this

    # TODO: Caution, there may be more than one output for a given transition

    # targets
    outputs = {}
    for tg in qt.getListOfOutputs():
        tg_id = tg.getIdAttribute()
        tg_name = tg.getName()
        specie_id = qual_species[tg.getQualitativeSpecies()][0]
        outputs[specie_id] = (tg_id, tg_name, specie_id)

    function_terms = {'default': ('', qt.getDefaultTerm().getResultLevel())}

    for ft in qt.getListOfFunctionTerms():

        ft_result = ft.getResultLevel()

        ft_rule = readMathNode(ft.getMath())

        if ft_rule is None:
            ft_result = ''

        function_terms['interaction_' + str(ft_result)] = (ft_rule, ft_result)

    interactions[qt_id] = (qt_id, name, inputs, outputs, function_terms)

# TODO: workflow:
#  Create empty model ->
#  Add regulatory variables (parsing) ->
#  Add regulatory interactions from sympy object