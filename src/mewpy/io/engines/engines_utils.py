import os
import re
import warnings
from collections import namedtuple
from functools import partial
from typing import Tuple

import libsbml

from mewpy.util.constants import ModelConstants
from mewpy.germ.algebra import (And, Or, Not, Less, Greater, LessEqual, GreaterEqual, BoolFalse, BoolTrue,
                                NoneAtom, Symbolic, Equal, parse_expression)


def build_symbolic(expression) -> Tuple[Symbolic, str]:
    try:

        return parse_expression(expression), ''

    except SyntaxError:

        return NoneAtom(), f'{expression} cannot be parsed. Assigning empty expression instead'


def get_sbml_doc_to_read(io):
    if isinstance(io, str):

        if os.path.exists(io):
            doc = libsbml.readSBMLFromFile(io)

        else:
            raise OSError(f'{io} is not a valid input. Provide the path or file handler')

    elif hasattr(io, 'read'):

        doc = libsbml.readSBMLFromString(io.read())

    else:

        raise OSError(f'{io} is not a valid input. Provide the path or file handler')

    return doc


def get_sbml_doc_to_write(io,
                          level,
                          version,
                          packages,
                          packages_version,
                          packages_required,
                          sbo_term=False):
    doc = None

    if isinstance(io, str):
        # string

        try:
            doc = get_sbml_doc_to_read(io)

        except OSError:

            pass

    elif hasattr(io, 'read'):

        # filepath handler

        try:
            doc = get_sbml_doc_to_read(io)

        except OSError:

            pass

    else:

        raise OSError(f'{io} is not a valid input. Provide the path or file handler')

    if doc is None:

        # only supporting SBML level 3 version 1 with fbc version 2 or qual version 1
        sbml_namespace = libsbml.SBMLNamespaces(level, version)

        for package, package_version in zip(packages, packages_version):
            sbml_namespace.addPackageNamespace(package, package_version)

        doc = libsbml.SBMLDocument(sbml_namespace)

    for package, required in zip(packages, packages_required):
        doc.setPackageRequired(package, required)

    if sbo_term:
        doc.setSBOTerm(SBO_FBA_FRAMEWORK)

    return doc


# -----------------------------------------------------------------------------
# SBML AST NODES
# -----------------------------------------------------------------------------
ASTNODE_BOOLEAN_OPERATORS = {
    libsbml.AST_LOGICAL_AND: And,
    libsbml.AST_LOGICAL_OR: Or,
    libsbml.AST_LOGICAL_NOT: Not}

ASTNODE_RELATIONAL_OPERATORS = {
    libsbml.AST_RELATIONAL_LEQ: LessEqual,
    libsbml.AST_RELATIONAL_GEQ: GreaterEqual,
    libsbml.AST_RELATIONAL_LT: Less,
    libsbml.AST_RELATIONAL_GT: Greater,
    libsbml.AST_RELATIONAL_EQ: Equal
}

ASTNODE_VALUES = {
    libsbml.AST_INTEGER,
    libsbml.AST_REAL,
    libsbml.AST_REAL_E,
    libsbml.AST_RATIONAL,
}

ASTNODE_BOOLEAN_VALUES = {
    libsbml.AST_CONSTANT_FALSE: BoolFalse(),
    libsbml.AST_CONSTANT_TRUE: BoolTrue(),
}

ASTNODE_NAME = libsbml.AST_NAME

# -----------------------------------------------------------------------------
# CONSTANTS
# -----------------------------------------------------------------------------
LOWER_BOUND_ID = 'mewpy_lb'
UPPER_BOUND_ID = 'mewpy_ub'
ZERO_BOUND_ID = 'mewpy_zero_b'

BOUND_MINUS_INF = 'minus_inf'
BOUND_PLUS_INF = 'plus_inf'

# -----------------------------------------------------------------------------
# SBO
# -----------------------------------------------------------------------------
SBO_FBA_FRAMEWORK = "SBO:0000624"
SBO_DEFAULT_FLUX_BOUND = "SBO:0000626"
SBO_FLUX_BOUND = "SBO:0000625"
SBO_EXCHANGE_REACTION = "SBO:0000627"

# -----------------------------------------------------------------------------
# UNITS
# -----------------------------------------------------------------------------
UNIT_ID = 'mmol_per_gDW_per_hr'

Unit = namedtuple('Unit', ['kind', 'scale', 'multiplier', 'exponent'])
UNITS = (Unit(kind=libsbml.UNIT_KIND_MOLE, scale=-3, multiplier=1, exponent=1),
         Unit(kind=libsbml.UNIT_KIND_GRAM, scale=0, multiplier=1, exponent=-1),
         Unit(kind=libsbml.UNIT_KIND_SECOND, scale=0, multiplier=3600, exponent=-1))

# IMPORTANT NOTE: SOME FUNCTIONS FOR PARSING SBML ENTITIES (namely, metabolic entities) HAVE BEEN HEAVILY
# INSPIRED BY THE TALENTED PEOPLE DEVELOPING COBRAPY. CHECK THE SOURCE: https://github.com/opencobra/cobrapy

# -----------------------------------------------------------------------------
# Functions for id replacements (import/export)
# -----------------------------------------------------------------------------
SBML_DOT = "__SBML_DOT__"

# -----------------------------------------------------------------------------
# Note pattern
# -----------------------------------------------------------------------------
pattern_notes = re.compile(r'<(?P<prefix>(\w+:)?)p[^>]*>(?P<content>.*?)</(?P=prefix)p>', re.IGNORECASE | re.DOTALL)

pattern_to_sbml = re.compile(r'([^0-9_a-zA-Z])')

pattern_from_sbml = re.compile(r'__(\d+)__')


def _escape_non_alphanum(non_ascii):
    """converts a non alphanumeric character to a string representation of
    its ascii number """
    return "__" + str(ord(non_ascii.group())) + "__"


def _number_to_chr(digit_str):
    """converts an ascii number to a character """
    return chr(int(digit_str.group(1)))


def _clip(sid, prefix):
    """Clips a prefix from the beginning of a string if it exists."""
    return sid[len(prefix):] if sid.startswith(prefix) else sid


def _f_gene(sid, prefix='G_'):
    """Clips gene prefix from id."""
    sid = sid.replace(SBML_DOT, ".")
    sid = pattern_from_sbml.sub(_number_to_chr, sid)
    return _clip(sid, prefix)


def _f_gene_rev(sid, prefix='G_'):
    """Adds gene prefix to id."""
    sid = pattern_to_sbml.sub(_escape_non_alphanum, sid)
    return prefix + sid.replace(".", SBML_DOT)


def _f_specie(sid, prefix='M_'):
    """Clips specie/metabolite prefix from id."""
    sid = pattern_from_sbml.sub(_number_to_chr, sid)
    return _clip(sid, prefix)


def _f_specie_rev(sid, prefix='M_'):
    """Adds specie/metabolite prefix to id."""
    sid = pattern_to_sbml.sub(_escape_non_alphanum, sid)
    return prefix + sid


def _f_reaction(sid, prefix='R_'):
    """Clips reaction prefix from id."""
    sid = pattern_from_sbml.sub(_number_to_chr, sid)
    return _clip(sid, prefix)


def _f_reaction_rev(sid, prefix='R_'):
    """Adds reaction prefix to id."""
    sid = pattern_to_sbml.sub(_escape_non_alphanum, sid)
    return prefix + sid


def _f_transition(sid, prefix='TR_'):
    """Clips transition prefix from id."""
    sid = pattern_from_sbml.sub(_number_to_chr, sid)
    return _clip(sid, prefix)


def _f_transition_rev(sid, prefix='TR_'):
    """Adds transition prefix to id."""
    sid = pattern_to_sbml.sub(_escape_non_alphanum, sid)
    return prefix + sid


F_GENE = 'F_GENE'
F_GENE_REV = 'F_GENE_REV'
F_SPECIE = 'F_SPECIE'
F_SPECIE_REV = 'F_SPECIE_REV'
F_REACTION = 'F_REACTION'
F_REACTION_REV = 'F_REACTION_REV'
F_TRANSITION = 'F_TRANSITION'
F_TRANSITION_REV = 'F_TRANSITION_REV'

F_REPLACE = {
    F_GENE: _f_gene,
    F_GENE_REV: _f_gene_rev,
    F_SPECIE: _f_specie,
    F_SPECIE_REV: _f_specie_rev,
    F_REACTION: _f_reaction,
    F_REACTION_REV: _f_reaction_rev,
    F_TRANSITION: _f_transition,
    F_TRANSITION_REV: _f_transition_rev,

}


def f_id(identifier, f_type):
    return F_REPLACE.get(f_type, lambda x: x)(identifier)


def fs_id(identifier, f_types):
    for f_type in f_types:
        identifier = f_id(identifier, f_type)

    return identifier


def add_sbml_parameter(sbml_model, parameter_id, value, constant=True, sbo=None, unit_definition=None):
    parameter = sbml_model.createParameter()
    parameter.setId(parameter_id)
    parameter.setValue(value)
    parameter.setConstant(constant)

    if sbo:
        parameter.setSBOTerm(sbo)

    if unit_definition is not None:
        parameter.setUnits(unit_definition.getId())


def get_sbml_lb_id(sbml_model, reaction, unit_definition=None):
    value = reaction.lower_bound

    if value == ModelConstants.REACTION_LOWER_BOUND:
        return LOWER_BOUND_ID

    elif value == 0:
        return ZERO_BOUND_ID

    elif value == -float('Inf'):
        return BOUND_MINUS_INF

    else:

        return None


def get_sbml_ub_id(sbml_model, reaction, unit_definition=None):
    value = reaction.upper_bound

    if value == ModelConstants.REACTION_UPPER_BOUND:
        return UPPER_BOUND_ID

    elif value == 0:
        return ZERO_BOUND_ID

    elif value == float('Inf'):
        return BOUND_PLUS_INF

    else:

        return None


def set_gpr(engine, warning, reaction, sbml_rxn_fbc):
    if not reaction.gpr.is_none:

        gpr = str(reaction.gpr)

        if gpr:

            gpr = gpr.replace('&', 'and').replace('|', 'or')

            gpa = sbml_rxn_fbc.createGeneProductAssociation()
            op = gpa.setAssociation(gpr, True, False)

            if op is None:
                engine.warnings.append(partial(warning, f"Could not set {gpr} GPR for {reaction.id}"))

            elif isinstance(op, int) and op != libsbml.LIBSBML_OPERATION_SUCCESS:
                engine.warnings.append(partial(warning, f"Could not set {gpr} GPR for {reaction.id}"))


def set_math(identifier, expression, function_term):
    math = libsbml.parseL3Formula(expression)

    op = function_term.setMath(math)

    if op is None:
        return f'Could not set {expression} expression for {identifier}'

    elif isinstance(op, int) and op != libsbml.LIBSBML_OPERATION_SUCCESS:
        return f'Could not set {expression} expression for {identifier}'

    return ''


def convert_fbc(doc):
    conversion_properties = libsbml.ConversionProperties()
    conversion_properties.addOption("convert fbc v1 to fbc v2", True, "Convert FBC-v1 model to FBC-v2")

    result = doc.convert(conversion_properties)

    if result != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise ImportError("Conversion of SBML fbc plugin version 1 to 2 failed")


def write_sbml_doc(io, doc):
    if isinstance(io, str):
        libsbml.writeSBMLToFile(doc, io)

    elif hasattr(io, 'write'):
        sbml = libsbml.writeSBMLToString(doc)
        io.write(sbml)


# -----------------------------------------------------------------------------
# Warning stuff
# -----------------------------------------------------------------------------

def prom_warning(message):
    return warnings.warn(message, Warning, stacklevel=2)


def coregflux_warning(message):
    return warnings.warn(message, Warning, stacklevel=2)


def sbml_warning(message):
    return warnings.warn(message, Warning, stacklevel=2)


def cobra_warning(message):
    return warnings.warn(message, Warning, stacklevel=2)


def csv_warning(message):
    return warnings.warn(message, Warning, stacklevel=2)


def expression_warning(message):
    return warnings.warn(message, Warning, stacklevel=2)
