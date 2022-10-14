"""
Algebra Constants

This module contains all algebra constants, namely operators and operand variables.
In addition, one can also find in this module the corresponding string representations of each operator and operand.
String representations will be used during expression parsing and evaluation, using python built-in eval in safe mode.

This module also contains common escape characters which will be replaced by the descriptive naming during expression
parsing.
"""

# Algebra operators and operands objects
from .symbolic import (BoolTrue,
                       BoolFalse,
                       And,
                       Or,
                       Not,
                       Equal,
                       NotEqual,
                       Inequality,
                       Greater,
                       GreaterEqual,
                       Less,
                       LessEqual,
                       Integer,
                       Float,
                       Zero,
                       One,
                       Symbol,
                       NoneAtom)

# Boolean algebra escape chars to be replaced
BOOLEAN_ESCAPE_CHARS = {'-': '_dash_',
                        ',': '_comma_',
                        ';': '_semicolon_',
                        '[': '_lbracket_',
                        ']': '_rbracket_',
                        ':': '_colon_',
                        '+': '_plus_',
                        '/': '_slash_',
                        '\\': '_backslash_',
                        '$': '_dollar_',
                        '%': '_percentage_',
                        '"': '_quotes_',
                        '(e)': '_e_',
                        '(c)': '_c_',
                        '(p)': '_p_'}

# Relational algebra escape chars to be replaced
RELATIONAL_ESCAPE_CHARS = BOOLEAN_ESCAPE_CHARS

# Main string representations for Boolean algebra operators
AND = '&'
OR = '|'
NOT = '~'

# Other string representations for Boolean algebra operators
BOOLEAN_OPERATORS = {'not': NOT,
                     '!': NOT,
                     'and': AND,
                     'or': OR,
                     AND: AND,
                     OR: OR,
                     NOT: NOT}

# Main string representations for Boolean algebra operators (set)
GLOBAL_BOOLEAN_OPERATORS = {AND, OR, NOT}

# Main string representations for Boolean algebra operands (except symbolic variables)
TRUE = '1'
FALSE = '0'

# Other string representations for Boolean algebra operands (except symbolic variables)
BOOLEAN_STATES = {'on': TRUE,
                  'off': FALSE,
                  'true': TRUE,
                  'false': FALSE,
                  TRUE: TRUE,
                  FALSE: FALSE}

# Main string representations for Boolean algebra operands (set)
GLOBAL_BOOLEAN_STATES = {TRUE, FALSE}

# Main string representations for Strict Relational algebra operators
GREATER = '>'
LESS = '<'
EQUAL = '='
NOT_EQUAL = '!='

# Other string representations for Strict Relational algebra operators
RELATIONAL_OPERATORS = {'greater than': GREATER,
                        'less than': LESS,
                        'greater': GREATER,
                        'less': LESS,
                        '==': EQUAL,
                        'equal': EQUAL,
                        'not equal': NOT_EQUAL,
                        GREATER: GREATER,
                        LESS: LESS,
                        EQUAL: EQUAL,
                        NOT_EQUAL: NOT_EQUAL}

# Main string representations for Strict Relational algebra operators (set)
GLOBAL_RELATIONAL_OPERATORS = {GREATER, LESS, EQUAL, NOT_EQUAL}

# Main string representations for Non-Strict Relational algebra operators
GREATER_EQUAL = '>='
LESS_EQUAL = '=>'

# Other string representations for Non-Strict Relational algebra operators
RELATIONAL_EQUAL_OPERATORS = {'greater than or equal': GREATER_EQUAL,
                              'less than or equal': LESS_EQUAL,
                              'greater or equal': GREATER_EQUAL,
                              'less or equal': LESS_EQUAL,
                              GREATER_EQUAL: GREATER_EQUAL,
                              LESS_EQUAL: LESS_EQUAL}

# Main string representations for Non-Strict Relational algebra operators (set)
GLOBAL_RELATIONAL_EQUAL_OPERATORS = {GREATER_EQUAL, LESS_EQUAL}

# Main string representations for Non-Strict and Strict algebra operands (except symbolic variables)
RELATIONAL_STATES = BOOLEAN_STATES

# Main string representations for all algebra operands (except symbolic variables)
GLOBAL_RELATIONAL_STATES = GLOBAL_BOOLEAN_STATES

# Main string representations for Non-Strict and Strict algebra operators and operands including empty, numeric and
# symbolic variables
GLOBAL_MEWPY_OPERATORS = {'BoolFalse': BoolFalse,
                          'BoolTrue': BoolTrue,
                          'And': And,
                          'Or': Or,
                          'Not': Not,
                          'Equal': Equal,
                          'NotEqual': NotEqual,
                          'Inequality': Inequality,
                          'Greater': Greater,
                          'GreaterEqual': GreaterEqual,
                          'Less': Less,
                          'LessEqual': LessEqual,
                          'Integer': Integer,
                          'Float': Float,
                          'Zero': Zero,
                          'One': One,
                          'Symbol': Symbol,
                          'NoneAtom': NoneAtom}
