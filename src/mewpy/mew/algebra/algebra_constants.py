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

# Boolean escape chars to be replaced
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

RELATIONAL_ESCAPE_CHARS = BOOLEAN_ESCAPE_CHARS

# Boolean string operator symbols
AND = '&'
OR = '|'
NOT = '~'

BOOLEAN_OPERATORS = {'not': NOT,
                     '!': NOT,
                     'and': AND,
                     'or': OR,
                     AND: AND,
                     OR: OR,
                     NOT: NOT}

GLOBAL_BOOLEAN_OPERATORS = {AND, OR, NOT}

# Boolean states
TRUE = '1'
FALSE = '0'

BOOLEAN_STATES = {'on': TRUE,
                  'off': FALSE,
                  'true': TRUE,
                  'false': FALSE,
                  TRUE: TRUE,
                  FALSE: FALSE}

GLOBAL_BOOLEAN_STATES = {TRUE, FALSE}

# Strict relational operators
GREATER = '>'
LESS = '<'
EQUAL = '='
NOT_EQUAL = '!='

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

GLOBAL_RELATIONAL_OPERATORS = {GREATER, LESS, EQUAL, NOT_EQUAL}

# Relational operators
GREATER_EQUAL = '>='
LESS_EQUAL = '=>'

RELATIONAL_EQUAL_OPERATORS = {'greater than or equal': GREATER_EQUAL,
                              'less than or equal': LESS_EQUAL,
                              'greater or equal': GREATER_EQUAL,
                              'less or equal': LESS_EQUAL,
                              GREATER_EQUAL: GREATER_EQUAL,
                              LESS_EQUAL: LESS_EQUAL}

GLOBAL_RELATIONAL_EQUAL_OPERATORS = {GREATER_EQUAL, LESS_EQUAL}

RELATIONAL_STATES = BOOLEAN_STATES

GLOBAL_RELATIONAL_STATES = GLOBAL_BOOLEAN_STATES

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
