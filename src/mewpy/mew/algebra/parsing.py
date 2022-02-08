import re
from io import StringIO
from token import NAME, OP, NUMBER
from tokenize import generate_tokens, untokenize
from typing import List

from .symbolic import NoneAtom, Symbolic
from .algebra_constants import (BOOLEAN_STATES,
                                BOOLEAN_OPERATORS,
                                TRUE,
                                FALSE,
                                RELATIONAL_STATES,
                                RELATIONAL_OPERATORS,
                                RELATIONAL_EQUAL_OPERATORS,
                                GLOBAL_RELATIONAL_EQUAL_OPERATORS,
                                GLOBAL_RELATIONAL_OPERATORS,
                                BOOLEAN_ESCAPE_CHARS,
                                RELATIONAL_ESCAPE_CHARS,
                                GLOBAL_MEWPY_OPERATORS)

_relational_ops = set()
_relational_ops.update(GLOBAL_RELATIONAL_OPERATORS)
_relational_ops.update(GLOBAL_RELATIONAL_EQUAL_OPERATORS)


class ExpressionParser:

    def __init__(self,
                 expression=None,

                 stringify_expression='',
                 parsed_expression='',
                 symbolic_expression='',
                 tokenized_expression=None,

                 tokens=None,
                 filters=None,
                 transformations=None,
                 escape_chars=None,
                 replaces=None,

                 symbols=None,
                 aliases=None,

                 is_boolean=False,
                 is_true=False,
                 is_false=False,
                 is_and=False,
                 is_or=False,
                 is_not=False,

                 is_relational=False,
                 is_equal=False,
                 is_not_equal=False,
                 is_inequality=False,
                 is_greater=False,
                 is_greater_equal=False,
                 is_less=False,
                 is_less_equal=False,

                 is_numeric=False,
                 is_integer=False,
                 is_float=False,
                 is_one=False,
                 is_zero=False,

                 is_symbol=False,
                 is_atom=False,

                 is_none=False):

        """
        Internal use only!

        The ExpressionParser object is to be used only by the next parsing, filtering, transformation
        and processing methods.
        It infers the algebra expression nature, and holds corresponding transforming, filtering and replacing methods
        and dictionaries to be applied to the raw expression.

        """

        if aliases is None:
            aliases = {}

        if symbols is None:
            symbols = {}

        if replaces is None:
            replaces = {}

        if escape_chars is None:
            escape_chars = {}

        if tokens is None:
            tokens = []

        if filters is None:
            filters = []

        if transformations is None:
            transformations = []

        if tokenized_expression is None:
            tokenized_expression = []

        self.expression = expression

        self.stringify_expression = stringify_expression
        self.parsed_expression = parsed_expression
        self.symbolic_expression = symbolic_expression
        self.tokenized_expression = tokenized_expression

        self.tokens = tokens
        self.filters = filters
        self.transformations = transformations
        self.escape_chars = escape_chars
        self.replaces = replaces

        self.symbols = symbols
        self.aliases = aliases

        self.is_boolean = is_boolean
        self.is_true = is_true
        self.is_false = is_false
        self.is_and = is_and
        self.is_or = is_or
        self.is_not = is_not

        self.is_relational = is_relational
        self.is_equal = is_equal
        self.is_not_equal = is_not_equal
        self.is_inequality = is_inequality
        self.is_greater = is_greater
        self.is_greater_equal = is_greater_equal
        self.is_less = is_less
        self.is_less_equal = is_less_equal

        self.is_numeric = is_numeric
        self.is_integer = is_integer
        self.is_float = is_float
        self.is_one = is_one
        self.is_zero = is_zero

        self.is_symbol = is_symbol
        self.is_atom = is_atom

        self.is_none = is_none

    def build(self):

        for token in self.tokenized_expression:

            if token.lower() in BOOLEAN_OPERATORS:

                if self.is_boolean:
                    continue

                self.is_boolean = True
                diff = set(boolean_filters) - set(self.filters)
                self.filters = self.filters + list(diff)
                diff = set(boolean_transformations) - set(self.transformations)
                self.transformations = self.transformations + list(diff)
                self.escape_chars.update(BOOLEAN_ESCAPE_CHARS)

            elif token.lower() in RELATIONAL_OPERATORS or token.lower() in RELATIONAL_EQUAL_OPERATORS:

                if self.is_relational:
                    continue

                self.is_relational = True
                diff = set(relational_filters) - set(self.filters)
                self.filters = self.filters + list(diff)
                diff = set(relational_transformations) - set(self.transformations)
                self.transformations = self.transformations + list(diff)
                self.escape_chars.update(RELATIONAL_ESCAPE_CHARS)

            elif token.lower() in BOOLEAN_STATES or token.lower() in RELATIONAL_STATES:

                if self.is_one or self.is_true or self.is_false or self.is_zero:
                    continue

                bool_state = BOOLEAN_STATES.get(token.lower(), token.lower())

                if bool_state == TRUE:
                    self.is_one = True
                    self.is_true = True

                if bool_state == FALSE:
                    self.is_zero = True
                    self.is_false = True

                diff = set(state_filters) - set(self.filters)
                self.filters = self.filters + list(diff)
                diff = set(state_transformations) - set(self.transformations)
                self.transformations = self.transformations + list(diff)

            else:

                if self.is_symbol:
                    continue

                self.is_symbol = True
                diff = set(symbol_filters) - set(self.filters)
                self.filters = self.filters + list(diff)
                diff = set(symbol_transformations) - set(self.transformations)
                self.transformations = self.transformations + list(diff)
                self.escape_chars.update({**BOOLEAN_ESCAPE_CHARS, **RELATIONAL_ESCAPE_CHARS})

        if not self.is_boolean or self.is_relational or self.is_one or self.is_true or \
                self.is_zero or self.is_false or self.is_symbol:
            self.filters = all_filters
            self.transformations = all_transformations
            self.escape_chars = {**BOOLEAN_ESCAPE_CHARS, **RELATIONAL_ESCAPE_CHARS}


def tokenize(rule: str) -> List[str]:
    """
    Tokenizes a stringify expression.
    :param rule: stringify expression as string
    :return: it returns all tokens of the expression
    """
    return list(filter(lambda x: x != '', rule.replace('(', ' ( ').replace(')', ' ) ').split(' ')))


def escape_chars_filter(expression: ExpressionParser):
    """
    For a given expression it parses out the special chars by replacing them with the corresponding values.
    The global dictionary available in the algebra constants module is used according to the expression type.
    This global dictionary can be altered for adding or removing special chars.
    :param expression: the ExpressionParser object
    """

    for escape_char, replace in expression.escape_chars.items():

        if escape_char in expression.parsed_expression:
            expression.replaces[replace] = escape_char
            expression.parsed_expression = expression.parsed_expression.replace(escape_char, replace)


def digit_filter(expression):
    """
    For a given regulatory self it checks if all regulatory variables start with a digit (str.is_digit()) and
    parses out the regulatory variable by adding the prefix '_dg_'

    """

    expression.parsed_expression = '(' + expression.parsed_expression + ')'
    regexp = re.compile(r'[^a-zA-Z|^_][0-9]+[a-zA-Z]|[^a-zA-Z|^_][0-9]+_')
    res = list(regexp.finditer(expression.parsed_expression))

    new_rule = ''

    last_nd = 0
    for i in range(len(res)):
        st = res[i].start() + 1
        nd = res[i].end()
        new_rule = new_rule + expression.parsed_expression[last_nd:st] + '_dg_' + expression.parsed_expression[
                                                                                  st:nd]
        last_nd = res[i].end()
        expression.replaces['_dg_'] = ''

    new_rule = new_rule + expression.parsed_expression[last_nd:]

    expression.parsed_expression = new_rule[1:-1]


def bitwise_filter(expression):
    """
    Transforms python boolean operators and, or and not, or other propositional logic syntax (e.g. greater
    than) to corresponding bitwise ones
    """

    bit_replace = {}

    if expression.is_boolean:
        bit_replace.update(BOOLEAN_OPERATORS)

    if expression.is_relational:
        bit_replace.update(RELATIONAL_EQUAL_OPERATORS)

    if expression.is_one or expression.is_true or expression.is_zero or expression.is_false:
        bit_replace.update(RELATIONAL_STATES)
        bit_replace.update(BOOLEAN_STATES)

    bit_ = {}

    for key, val in bit_replace.items():

        if key == val:
            continue

        if val in bit_:
            bit_[val].update({key, key.upper(), key.lower(), key.title()})
        else:
            bit_[val] = {key, key.upper(), key.lower(), key.title()}

    for bit_val, bools in bit_.items():

        for python_bool in bools:
            expression.parsed_expression = re.sub(r'\b{}\b'.format(python_bool), bit_val,
                                                  expression.parsed_expression)


def relational_filter(expression):
    """
    For a given boolean self (string) finds all relational expressions and encloses all of them with
    parenthesis and replaces the ambiguous operators

    """

    # regex to match a relational of any type: x>0; x<0; x>=0; x=>0;reverse order 0>x ...; with or without white
    # spaces ([a-zA-Z0-9_]+>[0-9]+)|([a-zA-Z0-9_]+\str>\str[0-9]+)

    regex_str = r'|'.join([r'([a-zA-Z0-9_]+' + regex + r'[0-9.]+)' + r'|' +
                           r'([a-zA-Z0-9_]+\s' + regex + r'\s[0-9.]+)' + r'|' +
                           r'([0-9.]+' + regex + r'[a-zA-Z0-9_]+)' + r'|' +
                           r'([0-9.]+\s' + regex + r'\s[a-zA-Z0-9_]+)'
                           for regex in _relational_ops])

    regexp = re.compile(regex_str)
    res = list(regexp.finditer(expression.parsed_expression))

    for match in set(map(lambda x: x.group(), res)):
        expression.parsed_expression = re.sub(r'\b{}\b'.format(match), '(' + match + ')',
                                              expression.parsed_expression)


def symbolic_transform(expression):
    """
    Heavily inspired by the talented people contributing for the impressive sympy package

    :param expression:
    :return:
    """

    result = []

    expression.tokens.append((None, None))  # so zip traverses all tokens

    for tok, next_token in zip(expression.tokens, expression.tokens[1:]):

        token_number, token_val = tok

        if token_number == NAME:

            name = token_val

            if name == 'False':
                result.extend([(NAME, 'BoolFalse'),
                               (OP, '('),
                               (NAME, repr(str(name))),
                               (OP, ')')])

            elif name == 'True':
                result.extend([(NAME, 'BoolTrue'),
                               (OP, '('),
                               (NAME, repr(str(name))),
                               (OP, ')')])

            elif name == 'None':
                result.extend([(NAME, 'NoneAtom'),
                               (OP, '('),
                               (NAME, repr(str(name))),
                               (OP, ')')])

            else:
                result.extend([
                    (NAME, 'Symbol'),
                    (OP, '('),
                    (NAME, repr(str(token_val))),
                    (OP, ')'),
                ])

        else:

            result.append((token_number, token_val))

    expression.tokens = result


# FIXME: parsing does not handle minus so far. BOOLEAN_ESCAPE_CHARS = {'-': '_dash_',
# FIXME: parsing does not handle scientific notation so far. BOOLEAN_ESCAPE_CHARS = {'-': '_dash_',
def numeric_transform(expression):
    """

    Heavily inspired by the talented people contributing for the impressive sympy package

    :param expression:
    :return:
    """

    result = []

    for token_number, token_val in expression.tokens:

        if token_number == NUMBER:

            number = token_val

            if number in ('1.0', '1'):

                seq = [(NAME, 'One'),
                       (OP, '('),
                       (NUMBER, repr(str(number))),
                       (OP, ')')]

            elif number in ('0.0', '0'):

                seq = [(NAME, 'Zero'),
                       (OP, '('),
                       (NUMBER, repr(str(number))),
                       (OP, ')')]

            elif '.' in number or (('e' in number or 'E' in number)
                                   and not (number.startswith('0x') or number.startswith('0X'))):

                seq = [(NAME, 'Float'),
                       (OP, '('),
                       (NUMBER, repr(str(number))),
                       (OP, ')')]

            else:

                seq = [(NAME, 'Integer'),
                       (OP, '('),
                       (NUMBER, number),
                       (OP, ')')]

            result.extend(seq)

        else:

            result.append((token_number, token_val))

    expression.tokens = result


all_filters = [escape_chars_filter, digit_filter, bitwise_filter, relational_filter]
all_transformations = [symbolic_transform, numeric_transform]

boolean_filters = [escape_chars_filter, digit_filter, bitwise_filter]
boolean_transformations = [symbolic_transform, numeric_transform]

relational_filters = [escape_chars_filter, digit_filter, bitwise_filter, relational_filter]
relational_transformations = [symbolic_transform, numeric_transform]

symbol_filters = [escape_chars_filter, digit_filter]
symbol_transformations = [symbolic_transform]

state_filters = [bitwise_filter]
state_transformations = [symbolic_transform, numeric_transform]


def apply_filters(expression, filters=None):
    if not filters:
        filters = expression.filters

    sorted_filters = []

    for _filter in sorted_filters:
        if callable(_filter):
            if _filter is escape_chars_filter:
                sorted_filters.insert(0, _filter)
            else:
                sorted_filters.append(_filter)

    for _filter in filters:
        _filter(expression)

    return expression


def apply_transformations(expression, transformations=None):
    """
    Heavily inspired by the talented people contributing for the impressive sympy package

    :param expression:
    :param transformations:
    :return:
    """

    if not transformations:
        transformations = expression.transformations

    sorted_transformations = []

    for transform in transformations:
        if callable(transform):
            if transform is symbolic_transform:
                sorted_transformations.insert(0, transform)
            else:
                sorted_transformations.append(transform)

    for transform in sorted_transformations:
        transform(expression)

    expression.symbolic_expression = untokenize(expression.tokens)

    return expression


def evaluate_expression(symbolic_expression, local_dict=None, global_dict=None):
    if not local_dict:
        local_dict = {}

    if not global_dict:
        global_dict = GLOBAL_MEWPY_OPERATORS

    return eval(symbolic_expression, global_dict, local_dict)


def parse_expression(expression: str) -> Symbolic:
    """
    Parsing an expression encoded into a string object to a Symbolic object that provides symbolic algebra evaluation
    with Symbolic-based symbols using regular python binary operators or custom functions

    For more details, consult the Expression object.

    :param expression: an expression as string object containing symbols and python binary (or not) operators
    :return: it returns a Symbolic object that allows for symbolic algebra evaluation
    """
    if not expression:

        expr_parser = ExpressionParser()
        expr_parser.expression = NoneAtom()
        expr_parser.is_none = True

    elif isinstance(expression, str):

        expr_parser = ExpressionParser(expression=None,
                                       stringify_expression=expression,
                                       parsed_expression=expression,
                                       tokenized_expression=tokenize(expression))

        expr_parser.build()

        expr_parser = apply_filters(expr_parser)

        expr_parser.tokenized_expression = tokenize(expr_parser.parsed_expression)

        input_code = StringIO(expr_parser.parsed_expression.strip())

        for token_number, token_val, _, _, _ in generate_tokens(input_code.readline):
            expr_parser.tokens.append((token_number, token_val))

        expr_parser = apply_transformations(expr_parser)

        if len(expr_parser.tokenized_expression) < 1:
            expr_parser.expression = NoneAtom()
            expr_parser.is_none = True

        else:

            expr_parser.expression = evaluate_expression(expr_parser.symbolic_expression)

    else:

        raise ValueError('Expression could not be parsed')

    for element in expr_parser.expression:

        if element.is_symbol:

            old_reg_name = ''.join(element.value)

            for replace, special_char in expr_parser.replaces.items():
                old_reg_name = old_reg_name.replace(replace, special_char)

            element.name = old_reg_name

    return expr_parser.expression
