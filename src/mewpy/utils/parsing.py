from abc import abstractmethod
from operator import add, sub, mul, truediv, pow
import sys
import re
from warnings import warn
from sympy.parsing.sympy_parser import parse_expr
from sympy.logic.boolalg import BooleanTrue, BooleanFalse
from sympy.core.numbers import Zero, One
from sympy import lambdify, Symbol


# Boolean operator symbols
S_AND = '&'
S_OR = '|'
S_NOT = '~'
S_ON = '1'
S_OFF = '0'
S_GREATER = '>'
S_LESS = '<'
S_EQUAL = '='
S_GREATER_THAN_EQUAL = '>='
S_LESS_THAN_EQUAL = '=>'
# Empty leaf symbol
EMPTY_LEAF = '@'

BOOLEAN_OPERATORS = [S_AND, S_OR, S_NOT]

BOOLEAN_STATES = [S_ON, S_OFF]

BOOLEAN_RELATIONALS = [S_GREATER, S_LESS, S_EQUAL]
BOOLEAN_EQUAL_RELATIONALS = [S_GREATER_THAN_EQUAL, S_LESS_THAN_EQUAL]

# TODO: This might be problematic, but it seems to be the better way. A disclaimer/instructions of how genes and
#  special vars must be formatted can also work. If the user provides a list of all genes/regulatory variables in the
#  TRN or search them in the ids column, this can also be some work arounds.

# Special chars to be replaced
BOOLEAN_SPECIAL_CHARS = {'-': '_dash_',
                         # '.' : '_dot_',
                         ',': '_comma_',
                         ';': '_semicolon_',
                         '[': '_lbracket_',
                         ']': '_rbracket_',
                         ':': '_colon_',
                         '+': '_plus_',
                         '/': '_slash_',
                         '\\"': '_backslash_',
                         '=': '_equal_',
                         '$': '_dollar_',
                         '%': '_percentage_',
                         '"': '_quotes_',
                         '(e)': '_e_',
                         '(c)': '_c_',
                         '(p)': '_p_',
                         '(E)': '_E_',
                         '(C)': '_C_',
                         '(P)': '_P_',
                         'lambda': '_lambda_function_',
                         'Add': '_Add_',
                         'add': '_add_',
                         'Mul': '_Mul_',
                         'mul': '_mul_',
                         'Pow': '_Pow_',
                         'pow': '_pow_',
                         }

# Increases the system recursion limit for long expressions
sys.setrecursionlimit(100000)


def evaluate_expression(expression, variables):
    """ Evaluates a logical expression (containing variables, 'and','or','(' and ')') 
        against the presence (True) or absence (False) of propositions within a list.
    """
    expression = expression.replace('(', '( ').replace(')', ' )')
    # Symbol conversion not mandatory
    expression = expression.replace('and', '&').replace('or', '|')
    tokens = expression.split()
    sentence = []
    for token in tokens:
        if token in "&|()":
            sentence.append(token)
        elif token in variables:
            sentence.append("True")
        else:
            sentence.append("False")
    proposition = ' '.join(sentence)
    return eval(proposition)


def evaluate_expression_tree(expression, variables):
    """ Evaluates a logical expression (containing variables, 'and','or', 'not','(' , ')') 
        against the presence (True) or absence (False) of propositions within a list.
        Assumes the correctness of the expression.
    """
    t = build_tree(expression, Boolean)
    evaluator = BooleanEvaluator(variables)
    res = t.evaluate(evaluator.f_operand, evaluator.f_operator)
    return res


class Node(object):
    """
    Binary syntax tree node
    """

    def __init__(self, value, left=None, right=None):
        if left is not None and not isinstance(left, Node):
            raise ValueError("Invalid left element")
        if right is not None and not isinstance(right, Node):
            raise ValueError("Invalid right element")
        self.value = value
        self.left = left
        self.right = right

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if self.is_leaf():
            return str(self.value)
        else:
            return f"{str(self.value)} ( {str(self.left)} , {str(self.right)} )"

    def is_leaf(self):
        return not self.left and not self.right

    def is_empty_leaf(self):
        return self.value == EMPTY_LEAF

    def is_unary(self):
        return (self.left.is_empty_leaf() and not self.right.is_empty_leaf()) or (
                not self.left.is_empty_leaf() and self.right.is_empty_leaf())

    def is_binary(self):
        return not self.left.is_empty_leaf() and not self.right.is_empty_leaf()

    def get_operands(self):
        if self.is_leaf():
            if self.value == EMPTY_LEAF:
                return set()
            else:
                return {self.value}
        else:
            return self.left.get_operands().union(self.right.get_operands())

    def print_node(self, level=0):
        tabs = ""
        for _ in range(level):
            tabs += "\t"
        if self.is_leaf():
            if self.value != EMPTY_LEAF:
                print(tabs, f"|____ \033[1;31m{self.value}\033[0;0m")
            else:
                pass
        else:
            print(tabs, f"|____[\033[;1m\033[94m{self.value}\033[0;0m]")
        if self.left is not None:
            self.left.print_node(level + 1)
        if self.right is not None:
            self.right.print_node(level + 1)

    def evaluate(self, f_operand=None, f_operator=None):
        """
        Evaluates the expression using the f_operand and f_operator mapping functions
        """
        if f_operand is None or f_operator is None:
            return eval(str(self))
        elif self.is_leaf():
            return f_operand(self.value)
        else:
            return f_operator(self.value)(self.left.evaluate(f_operand, f_operator),
                                          self.right.evaluate(f_operand, f_operator))

    def get_conditions(self):
        """
        Retrieves the propositional conditions

        return: The list of conditions
        """
        ops = self.get_operands()
        return set([i for i in ops if is_condition(i)])



class Syntax:
    """Defines an interface for the tree syntax parsing 
       with operators, their precedence and associativity.
    """

    @abstractmethod
    def is_operator(self, op):
        raise NotImplementedError

    @abstractmethod
    def is_greater_precedence(self, op1, op2):
        raise NotImplementedError

    @abstractmethod
    def associativity(self, op):
        raise NotImplementedError

    @staticmethod
    def arity(op):
        return 2

    @staticmethod
    def replace():
        return {}


class Arithmetic(Syntax):

    @staticmethod
    def is_operator(op):
        return op in ['+', '-', '*', '/', '^']

    @staticmethod
    def is_greater_precedence(op1, op2):
        pre = {'+': 0, '-': 0, '*': 1, '/': 1, '^': 2}
        return pre[op1] >= pre[op2]

    @staticmethod
    def associativity(op):
        ass = {'+': 0, '-': 0, '*': 0, '/': 0, '^': 1}
        return ass[op]

    @staticmethod
    def arity(op):
        ar = {'+': 2, '-': 2, '*': 2, '/': 2, '^': 2}
        return ar[op]


class ArithmeticEvaluator:

    @staticmethod
    def f_operator(op):
        operators = {'+': add, '-': sub, '*': mul, '/': truediv, '^': pow}
        if op in operators.keys():
            return operators[op]
        else:
            raise ValueError(f"Operator {op} not defined")

    @staticmethod
    def f_operand(op):
        if is_number(op):
            return float(op)
        else:
            raise ValueError(f"Operator {op} not defined")


class Boolean(Syntax):
    """
    A boolean syntax parser where (NOT) is considered to
    be defined as a binary operator (NOT a) = (EMPTY NOT a)
    """

    @staticmethod
    def is_operator(op):
        return op in [S_AND, S_OR, S_NOT]

    @staticmethod
    def is_greater_precedence(op1, op2):
        pre = {S_AND: 0, S_OR: 0, S_NOT: 1}
        return pre[op1] >= pre[op2]

    @staticmethod
    def associativity(op):
        ass = {S_AND: 0, S_OR: 0, S_NOT: 1}
        return ass[op]

    @staticmethod
    def arity(op):
        ar = {S_AND: 2, S_OR: 2, S_NOT: 1}
        return ar[op]

    @staticmethod
    def replace():
        r = {'not': [EMPTY_LEAF, S_NOT], 'and': [S_AND], 'or': [S_OR], '&': [S_AND], '|': [S_OR], '~': [EMPTY_LEAF, S_NOT]}
        return r

    SYMPY_REPLACE = {'not': S_NOT,
                      '!': S_NOT,
                      'and': S_AND,
                      'or': S_OR,
                      'on': S_ON,
                      'off': S_OFF,
                      'true': S_ON,
                      'false': S_OFF}

    _SYMPY_RELATIONAL_REPLACE = {'greater than or equal': S_GREATER_THAN_EQUAL,
                                 'less than or equal': S_LESS_THAN_EQUAL,
                                 'greater or equal': S_GREATER_THAN_EQUAL,
                                 'less or equal': S_LESS_THAN_EQUAL,
                                 'greater than': S_GREATER,
                                 'less than': S_LESS,
                                 'greater': S_GREATER,
                                 'less': S_LESS,
                                 '==': S_EQUAL,
                                 'equal': S_EQUAL}


class BooleanEvaluator:
    """A boolean evaluator

        arguments:
            true_list (list): operands evaluated as True. Operands not in the list are evaluated as False
            variables (dict): a dictionary mapping symbols to values. Used to evaluate conditions.
        """

    def __init__(self, true_list=[], variables={}):
        self.true_list = true_list
        self.vars = variables


    def f_operator(self, op):
        operators = {S_AND: lambda x, y: x and y, S_OR: lambda x, y: x or y, S_NOT: lambda x, y: not y}
        if op in operators.keys():
            return operators[op]
        else:
            raise ValueError(f"Operator {op} not defined")

    def f_operand(self, op):
        if op.upper() == 'TRUE' or op == '1' or op in self.true_list:
            return True
        elif is_condition(op):
            return eval(op, None, self.vars)
        else:
            return False

    def set_true_list(self, true_list):
        self.true_list = true_list


class GeneEvaluator:
    """An evaluator for genes

        arguments:
            genes_value (dic): maps genes to values.
                               If a gene is not in the dictionary, a default value of 1 is considered.
            and_operator : function to be applied instead of (and', '&') (not case sensitive)
            or_operator : function to be applied instead of ('or','|') (not case sensitive)
    """

    def __init__(self, genes_value, and_operator, or_operator, prefix=""):

        self.genes_value = genes_value
        self.and_operator = and_operator
        self.or_operator = or_operator
        self.prefix = prefix

    def f_operator(self, op):
        operators = {S_AND: self.and_operator, S_OR: self.or_operator}
        if op in operators.keys():
            return operators[op]
        else:
            raise ValueError(f"Operator {op} not defined")

    def f_operand(self, op):
        if op[len(self.prefix):] in self.genes_value:
            return self.genes_value[op]
        else:
            return 1

def build_tree(exp, rules):
    """ Builds a parsing syntax tree for basic mathematical expressions

        parameters:
            exp (str) : the expression to be parsed
            rules (Syntax): defines the syntax rules
    """
    replace_dic = rules.replace()
    exp_ = tokenize_infix_expression(exp)
    exp_list = []
    for i in exp_:
        if i.lower() in replace_dic:
            exp_list.extend(replace_dic[i.lower()])
        else:
            exp_list.append(i)
    stack = []
    tree_stack = []
    predecessor = None
    for i in exp_list:
        if not (rules.is_operator(i) or i in ['(', ')']):
            if predecessor and not (rules.is_operator(predecessor) or predecessor in ['(', ')']):
                s = tree_stack[-1].value
                tree_stack[-1].value = s + " " + i
            else:
                t = Node(i)
                tree_stack.append(t)
        elif rules.is_operator(i):
            if not stack or stack[-1] == '(':
                stack.append(i)

            elif rules.is_greater_precedence(i, stack[-1]) and rules.associativity(i) == 1:
                stack.append(i)

            else:
                while stack and stack[-1] != '(' and rules.is_greater_precedence(stack[-1], i) and rules.associativity(
                        i) == 0:
                    popped_item = stack.pop()
                    t = Node(popped_item)
                    t1 = tree_stack.pop()
                    t2 = tree_stack.pop()
                    t.right = t1
                    t.left = t2
                    tree_stack.append(t)
                stack.append(i)

        elif i == '(':
            stack.append('(')

        elif i == ')':
            while stack[-1] != '(':
                popped_item = stack.pop()
                t = Node(popped_item)
                t1 = tree_stack.pop()
                t2 = tree_stack.pop()
                t.right = t1
                t.left = t2
                tree_stack.append(t)
            stack.pop()
        predecessor = i

    while stack:
        popped_item = stack.pop()
        t = Node(popped_item)
        # t1 = None if len(tree_stack) < 1 else tree_stack.pop()
        # t2 = None if len(tree_stack) < 1 else tree_stack.pop()
        t1 = tree_stack.pop()
        t2 = tree_stack.pop()
        t.right = t1
        t.left = t2
        tree_stack.append(t)

    t = tree_stack.pop()

    return t


def tokenize_infix_expression(exp):
    return list(filter(lambda x: x != '', exp.replace('(', ' ( ').replace(')', ' ) ').split(' ')))


def is_number(token):
    """ Returns True if the token is a number
    """
    return token.replace('.', '', 1).replace('-', '', 1).isnumeric()


def is_condition(token):
    """ Returns True if the token is a condition
    """
    regexp = re.compile(r'>|<|=')
    return bool(regexp.search(token))


# generic regulatory rule parser
def generic_regulatory_rule_parser(rule):
    exp = tokenize_infix_expression(rule)
    parsed_str = ''.join(exp)

    return parsed_str, exp


# special chars filter
def special_chars_filter(rule, special_chars=None):
    """
    For a given rule it parses out the special chars by replacing them with the corresponding values in the
    special_chars argument. By default, the global dictionary BOOL_SPECIAL_CHARS is used. This global dictionary can be
    altered for adding or removing special chars.


    :param rule: str, the regulatory rule
    :param special_chars: dict or None, special (or not) chars to be replaced. BOOL_SPECIAL_CHARS by default.
    special_chars argument must have unique keys and unique values!!
    :returns: str, the new parsed rule
    :returns: dict, the occurred replacements {replace: special_char}
    """

    if special_chars is None:
        special_chars = BOOLEAN_SPECIAL_CHARS

    replaces = {}

    for special_char, replace in special_chars.items():

        if special_char in rule:
            replaces[replace] = special_char
            rule = rule.replace(special_char, replace)

    return rule, replaces


# initial digit filter
def starts_with_digit_filter(rule):
    """
    For a given regulatory rule it checks if all regulatory variables start with a digit (str.is_digit()) and parses
    out the regulatory variable by adding the prefix '_st_dg_'

    :param rule: str, the regulatory rule
    :return: tuple, the new parsed regulatory variable as a string object and the occurred replacements as dict object
    """

    replaces = {}
    rule = '(' + rule + ')'
    regexp = re.compile(r'[^a-zA-Z|^_][0-9]+[a-zA-Z]|[^a-zA-Z|^_][0-9]+_')
    res = list(regexp.finditer(rule))

    new_rule = ''

    last_nd = 0
    for i in range(len(res)):
        st = res[i].start() + 1
        nd = res[i].end()
        new_rule = new_rule + rule[last_nd:st] + '_st_dg_' + rule[st:nd]
        last_nd = res[i].end()
        replaces['_st_dg_'] = ''

    new_rule = new_rule + rule[last_nd:]

    return new_rule[1:-1], replaces


# relational filter
def single_relational_filter(expression):
    """
    For a given regulatory expression checks if it is a relational expression and identifies the regulatory
    variable. Otherwise, returns false

    :param expression: str, the regulatory expression
    :return: tuple, the regulatory variable as a string object and the condition enclosed by parenthesis, or False
    """

    replace_dict = Boolean._SYMPY_RELATIONAL_REPLACE

    for condition in replace_dict:

        if condition in expression:

            expression = expression.replace(condition, replace_dict[condition])

        elif condition.upper() in expression:

            expression = expression.replace(condition.upper(), replace_dict[condition])

        elif condition.title() in expression:

            expression = expression.replace(condition.title(), replace_dict[condition])

        else:

            continue

    # first let's handle GreaterEqual LessEqual
    regex_str = '|'.join([regex for regex in BOOLEAN_EQUAL_RELATIONALS])
    regexp = re.compile(regex_str)
    res = list(regexp.finditer(expression))

    if len(res) == 1:

        return expression[0:res[0].start()], '(' + expression + ')'

    elif len(res) == 2:

        return expression[res[0].end():res[1].start()], '(' + expression + ')'

    elif len(res) > 2:

        raise ValueError("The condition {} cannot be parsed".format(expression))

    else:

        # Then let's handle Greater Less and Equal
        regex_str = '|'.join([regex for regex in BOOLEAN_RELATIONALS])
        regexp = re.compile(regex_str)
        res = list(regexp.finditer(expression))

        if len(res) == 1:
            return expression[0:res[0].start()], '(' + expression + ')'

        if len(res) == 2:
            return expression[res[0].end():res[1].start()], '(' + expression + ')'

        if len(res) > 2:
            raise ValueError("The condition {} cannot be parsed".format(expression))

    return False


# symbol filter
def symbol_filter(expression, local_dict):
    """
    For a given regulatory expression it identifies the regulatory variable and parses it out as a sympy symbol.

    :param expression: str, the regulatory expression
    : param local_dict: dict, mapping between the variables in the expression and a Symbol Sympy's object
    :return: Symbol, the regulatory variable as a sympy Symbol object
    """

    sympy_expr = parse_expr(expression, evaluate=False, local_dict=local_dict)

    def get_symbol(atoms):

        for atom in atoms:

            if not atom.is_Atom:
                return get_symbol(atom.atoms())

            if atom.is_Atom and atom.is_Symbol:
                return atom

        return

    return get_symbol(sympy_expr.atoms())


# variable filter
def variable_from_str(expression, filter=True, special_chars=None, replaces=None):
    """
    For a given expression (string) it identifies the variable and parses it out as a sympy symbol.
    If more than one variable is present in the regulatory expression, a BaseException is raised

    It also applies the starts_with_digit_filter and special_chars_filter if filter is True. If the replaces
    dictionary is passed or reconstructed during filtering, the function can return the regulatory variable alias.
    Otherwise, the alias is equal to the regulatory variable ID depending on the filter condition.

    :param expression: str, expression
    :param filter: bool, if filter, special_chars_filter and starts_with_digit_filter are applied
    :param special_chars: dict, dict or None, special (or not) chars to be replaced. BOOL_SPECIAL_CHARS is used by
    default. special_chars argument must have unique keys and unique values!!
    :param replaces: dict, replaces occurred for aliases construction. None is the default
    :returns regulatory variable: Sympy Symbol, regulatory variable
    :returns alias: dict, the new variable name map to the old one
    """

    if filter:

        expression, new_replaces = special_chars_filter(expression, special_chars)
        expression, new_replaces_d = starts_with_digit_filter(expression)
        new_replaces.update(new_replaces_d)

    else:

        new_replaces = {}

    if replaces is None:
        replaces = {}

    replaces.update(new_replaces)

    sympify = parse_expr(expression, evaluate=False)
    aliases = {}

    if len(sympify.atoms()) > 1:
        raise BaseException("Something went wrong processing the variable {}".format(expression))

    if not isinstance(sympify, Symbol):
        raise BaseException("Something went wrong processing the variable {}".format(expression))

    old_reg_name = ''.join(sympify.name)

    for replace, special_char in replaces.items():
        old_reg_name = old_reg_name.replace(replace, special_char)

    aliases[sympify.name] = old_reg_name

    return sympify, aliases


# # get tree
# def get_tree(rule):
#
#     try:
#
#         return build_tree(rule, Boolean)
#
#     except:
#
#         return Node(rule)


# relational filter
def relational_filter(rule):
    """
    For a given boolean expression (string) finds all relational expressions and encloses all of them with
    parenthesis and replaces the ambiguous operators

    :param rule: str, the boolean expression
    :return: str, the new boolean expression having all relational enclosed by parenthesis
    """

    for relational in Boolean._SYMPY_RELATIONAL_REPLACE:

        if relational in rule:

            rule = rule.replace(relational, Boolean._SYMPY_RELATIONAL_REPLACE[relational])

        elif relational.upper() in rule:

            rule = rule.replace(relational.upper(), Boolean._SYMPY_RELATIONAL_REPLACE[relational])

        elif relational.title() in rule:

            rule = rule.replace(relational.title(), Boolean._SYMPY_RELATIONAL_REPLACE[relational])

        else:

            continue

    # regex to match a relational of any type: x>0; x<0; x>=0; x=>0;reverse order 0>x ...; with or without white spaces
    # ([a-zA-Z0-9_]+>[0-9]+)|([a-zA-Z0-9_]+\s>\s[0-9]+)

    regex_str = '|'.join(['([a-zA-Z0-9_]+' + regex + '[0-9.]+)' + '|' +
                          '([a-zA-Z0-9_]+\s' + regex + '\s[0-9.]+)' + '|' +
                          '([0-9.]+' + regex + '[a-zA-Z0-9_]+)' + '|' +
                          '([0-9.]+\s' + regex + '\s[a-zA-Z0-9_]+)'
                          for regex in BOOLEAN_EQUAL_RELATIONALS + BOOLEAN_RELATIONALS])

    regexp = re.compile(regex_str)
    res = list(regexp.finditer(rule))

    for match in set(map(lambda x: x.group(), res)):
        rule = re.sub(r'\b{}\b'.format(match), '(' + match + ')', rule)

    return rule


def bitwise_rule_filter(rule):
    """
    Transforms python boolean operators and, or and not to corresponding bitwise ones
    :param rule: str, boolean rule
    :return: str, new parsed boolean rule
    """

    for bool, bit_val in Boolean.SYMPY_REPLACE.items():

        rule = re.sub(r'\b{}\b'.format(bool), bit_val, rule)
        rule = re.sub(r'\b{}\b'.format(bool.upper()), bit_val, rule)
        rule = re.sub(r'\b{}\b'.format(bool.title()), bit_val, rule)

    return rule


# Boolean regulatory rule from str
def boolean_rule_from_str(rule):
    """
    It creates a new parsed regulatory rule, list all its elements (operators and operands) and creates a sympy
    expression.
    It also identifies all regulatory variables and returns them.

    The regulatory variables might be filtered and transformed with the special characters and starts with digit
    filters. The old regulatory variable ID is kept as an alias

    :param rule: str, boolean rule
    :returns rule: str, parsed boolean rule
    :return elements: list, all elements in the rule (operators and operands)
    :return sympify: Sympy Boolean object Or, And or Not, Sympy Symbol object, or Sympy Integral object
    :return symbols: dict, all regulatory variables in the regulatory rule. Keys stand for their new IDs and
    values stand for their sympy symbols
    :return symbols: dict, all regulatory variables in the regulatory rule. Keys stand for their new IDs and
    values stand for their aliases
    """

    regex_str = '|'.join(['(' + regex + '[0-9.]+)' + '|' +
                          '(' + regex + '\s[0-9.]+)' + '|'
                          for regex in BOOLEAN_EQUAL_RELATIONALS + BOOLEAN_RELATIONALS])
    regexp = re.compile(regex_str)

    rule, replaces = special_chars_filter(rule)
    rule, new_replaces = starts_with_digit_filter(rule)
    replaces.update(new_replaces)

    rule = bitwise_rule_filter(rule)
    rule = relational_filter(rule)

    elements = tokenize_infix_expression(rule)

    aliases = {}

    if len(elements) < 1:
        return rule, elements, Node(None), [], aliases, []

    if len(elements) == 1:

        if elements[0] == '1' or elements[0] == '0':
            return rule, elements, build_tree(exp=rule, rules=Boolean), [], [], []

    bool_tree = build_tree(exp=rule, rules=Boolean)

    operands = bool_tree.get_operands()
    conditionals = bool_tree.get_conditions()
    new_operands = []
    new_conditionals = []

    for op in operands:

        res = list(regexp.finditer(op))

        new_op = op

        for match in set(map(lambda x: x.group(), res)):
            new_op = re.sub(r'\b{}\b'.format(match), '', op)

        new_operands.append(new_op)
        if op in conditionals:
            new_conditionals.append(new_op)

        old_reg_name = ''.join(new_op)

        for replace, special_char in replaces.items():
            old_reg_name = old_reg_name.replace(replace, special_char)

        aliases[new_op] = [old_reg_name]

    return rule, elements, bool_tree, new_operands, aliases, new_conditionals


# # Boolean regulatory rule from str
# def boolean_rule_from_str(rule):
#     """
#     It creates a new parsed regulatory rule, list all its elements (operators and operands) and creates a sympy
#     expression.
#     It also identifies all regulatory variables and returns them.
#
#     The regulatory variables might be filtered and transformed with the special characters and starts with digit
#     filters. The old regulatory variable ID is kept as an alias
#
#     :param rule: str, boolean rule
#     :returns rule: str, parsed boolean rule
#     :return elements: list, all elements in the rule (operators and operands)
#     :return sympify: Sympy Boolean object Or, And or Not, Sympy Symbol object, or Sympy Integral object
#     :return symbols: dict, all regulatory variables in the regulatory rule. Keys stand for their new IDs and
#     values stand for their sympy symbols
#     :return symbols: dict, all regulatory variables in the regulatory rule. Keys stand for their new IDs and
#     values stand for their aliases
#     """
#
#     rule, replaces = special_chars_filter(rule)
#     rule, new_replaces = starts_with_digit_filter(rule)
#     replaces.update(new_replaces)
#
#     rule = bitwise_rule_filter(rule)
#     rule = relational_filter(rule)
#
#     elements = tokenize_infix_expression(rule)
#
#     symbols = {}
#     aliases = {}
#
#     if len(elements) < 1:
#         return rule, elements, None, lambda x: None, symbols, aliases
#
#     else:
#
#         sympify = parse_expr(rule, evaluate=False)
#
#         if isinstance(sympify, BooleanFalse) or isinstance(sympify, Zero):
#
#             return rule, elements, sympify, lambda x: 0, symbols, aliases
#
#         elif isinstance(sympify, BooleanTrue) or isinstance(sympify, One):
#
#             return rule, elements, sympify, lambda x: 1, symbols, aliases
#
#         else:
#
#             for sb in sympify.atoms():
#
#                 if isinstance(sb, Symbol):
#
#                     symbols[sb.name] = sb
#
#                     old_reg_name = ''.join(sb.name)
#
#                     for replace, special_char in replaces.items():
#                         old_reg_name = old_reg_name.replace(replace, special_char)
#
#                     aliases[sb.name] = [old_reg_name]
#
#     return rule, elements, sympify, sympify.subs, symbols, aliases


# # Boolean regulatory rule parser
# def bool_regulatory_rule_parser(rule):
#     """
#     It creates a new parsed regulatory rule, list all its elements (operators and operands). It also identifies all
#     regulatory variables and returns them in a dict object
#
#     :param rule: str, regulatory rule
#     :returns parsed regulatory rule: str, parsed regulatory rule
#     :return elements: list, all elements in the rule (operators and operands)
#     :returns regulatory variables: dict, all regulatory variables in the regulatory rule (keys are the new id and
#     values are tuples (Sympy symbol, alias))
#     """
#
#     replace_dic = Boolean.SYMPY_REPLACE
#
#     new_rule, replaces = special_chars_filter(rule)
#     new_rule, new_replaces = starts_with_digit_filter(new_rule)
#     replaces.update(new_replaces)
#
#     exp = tokenize_infix_expression(new_rule)
#
#     parsed_str = ''
#     exp_elements = []
#     reg_vars = {}
#
#     for op in exp:
#
#         if op.lower() in replace_dic:
#
#             operator = replace_dic[op.lower()]
#             exp_elements.append(operator)
#             parsed_str += operator
#
#         elif op in BOOLEAN_OPERATORS + BOOLEAN_STATES:
#
#             exp_elements.append(op)
#             parsed_str += op
#
#         elif op in [')', '(']:
#
#             exp_elements.append(op)
#             parsed_str += op
#
#         else:
#
#             condition = single_relational_filter(op)
#
#             if condition:
#
#                 exp_elements.append(condition[1])
#                 parsed_str += condition[1]
#
#                 op = condition[0]
#
#             else:
#
#                 exp_elements.append(op)
#                 parsed_str += op
#
#             reg_var = variable_from_str(op, filter=False, replaces=replaces)
#
#             if not reg_var:
#                 def op_message(message):
#                     warn(message, UserWarning, stacklevel=2)
#
#                 op_message("The {} could not be parsed into a Sympy Symbol object. This may raise problems in "
#                            "the future".format(op))
#
#             else:
#                 reg_vars[reg_var[0].name] = reg_var
#
#     return parsed_str, exp_elements, reg_vars
#
#
# # Boolean sympy expression
# def bool_sympy_expression(rule, local_dict):
#     """
#     Wrapper for a sympy expression obtained from a parsed string. Input must be parsed with bool_regulatory_rule_parser
#
#     :param rule: str, parsed regulatory rule
#     :param local_dict: dict, mapping between the variables in the expression and teh respecitve Symbol Sympy's objects
#     :return: Sympy Boolean object Or, And or Not, Sympy Symbol object, or Sympy Integral object
#     """
#
#     if not rule:
#         return None
#     return parse_expr(rule, evaluate=False, local_dict=local_dict)
#
#
# # Boolean sympy lambda function
# def bool_lambdify(expression, args):
#     """
#         Wrapper for lambdify obtained from a bool_sympy_expression. Input must be parsed with
#         bool_regulatory_rule_parser and then a sympy expression must be created with bool_sympy_expression_from_str
#
#         Unfornatelly, lambdify does not handle Not operator in the latter versions
#
#         :param expression: sympy core object, parsed regulatory sympy expression
#         :param args: tuple, sympy Symbols of the regulatory sympy expression
#         :return: Sympy Boolean object Or, And or Not, Sympy Symbol object, or Sympy Integral object
#         """
#
#     return lambdify(args, expression, 'sympy', dummify=False)


def test_1():
    # aritmetic example

    t = build_tree(" 1 + 2 + 3 + ( 2 * 2 )", Arithmetic)
    res = t.evaluate(ArithmeticEvaluator.f_operand,
                     ArithmeticEvaluator.f_operator)
    print(t, ' = ', res)

    # boolean example

    t = build_tree(" A or not (B)", Boolean)
    print(t)
    print(type(t))
    t.print_node()
    print(t.get_operands())
    # propositions in the list are evaluated as True and the remaining are False
    evaluator = BooleanEvaluator(['a', 'c'])
    res = t.evaluate(evaluator.f_operand, evaluator.f_operator)
    print('value = ', res)

    # gene values example
    t = build_tree("( a or not b ) and (c)", Boolean)
    t.print_node()
    value_map = {'a': 1, 'b': 2, 'c': 3, 'd': 4}

    def op_and(x, y): return min(x, y)

    def op_or(x, y): return max(x, y)

    evaluator = GeneEvaluator(value_map, op_and, op_or)
    res = t.evaluate(evaluator.f_operand, evaluator.f_operator)
    print('Evaluation ', res)


def test_2():
    from urllib.request import urlretrieve
    from cobra.io import read_sbml_model
    import random

    path, _ = urlretrieve('http://bigg.ucsd.edu/static/models/RECON1.xml')
    model = read_sbml_model(path)
    ogpr = model.reactions.ATPS4m.gene_name_reaction_rule
    gpr = ogpr
    print(gpr)
    t = build_tree(gpr, Boolean)
    print(t)
    genes = list(t.get_operands())
    print("GENES:\n", genes)
    print("Evaluations:")
    evaluator = BooleanEvaluator(genes)
    res = t.evaluate(evaluator.f_operand, evaluator.f_operator)
    print(evaluator.true_list, " ==> ", res)
    for _ in range(20):
        l = []
        n = random.randint(0, len(genes))
        for _ in range(n):
            i = random.randint(0, len(genes) - 1)
            l.append(genes[i])
        evaluator.set_true_list(l)
        res = t.evaluate(evaluator.f_operand, evaluator.f_operator)
        print(evaluator.true_list, " ==> ", res)


if __name__ == '__main__':
    # test_1()
    # test_2()

    # rule = '(fuc-L(e)>0)'
    # rule = '(not fuc-L(e)>0)'
    # rule = '(not fuc-L(e)>0.1)'
    # rule = '(not fuc-L(e) > 0.1)'
    rule = '((NOT (o2(e)>0)) AND (NOT (no3(e)>0)) AND (NOT (no2(e)>0)) AND (NOT (tmao(e)>0)) AND (NOT (dmso(e)>0)) ' \
           'AND (for(e)>0))'
    # rule = '(RpoS OR ' \
    #        'Fnr OR ' \
    #        '((NOT IclR) AND ' \
    #        '((NOT( (glcn(e)>0) OR ' \
    #        '(glc-D(e)>0) OR ' \
    #        '(arab-L(e)>0) OR ' \
    #        '(xyl-D(e)>0) OR ' \
    #        '(rib-D(e)>0) OR ' \
    #        '(mal-L(e)>0) OR ' \
    #        '(glyc(e)>0) OR ' \
    #        '(sbt-D(e)>0) OR ' \
    #        '(lac-D(e)>0) OR ' \
    #        '(man(e)>0) OR ' \
    #        '(succ(e)>0))))))'
    # rule = '(NOT IclR) AND ((NOT ArcA) OR FruR > 0)'
    rule = 'S AND (NOT g22)'
    # rule = rule.replace('(e)', '_e_')
    # t = build_tree(rule, Boolean)
    # print(rule)
    # t.print_node()
    # print(t.get_operands())

    # res = bool_regulatory_rule_parser(rule)
    # local_dict = {key: val[0] for key, val in res[2].items()}
    # res_sympy = bool_sympy_expression(res[0], local_dict)
    # arguments = list(map(lambda x: x[0], res[2].values()))
    # res_sympy_function = bool_lambdify(res_sympy, arguments)
    # # state_map = {symbol.name: 1 for symbol in arguments}
    # state_map = {'IclR': 0, 'ArcA': 1, 'FruR': 1}
    # print(res_sympy_function(**state_map))
    #
    # rule = 'symbols'
    # res = bool_regulatory_rule_parser(rule)
    # # local_dict = {key: val[0] for key, val in res[2].items()}
    # local_dict = {}
    # res_sympy = bool_sympy_expression(res[0], local_dict)
    # arguments = list(map(lambda x: x[0], res[2].values()))
    # res_sympy_function = bool_lambdify(res_sympy, arguments)
    # # state_map = {symbol.name: 1 for symbol in arguments}
    # # state_map = {'IclR': 0, 'ArcA': 1, 'FruR': 1}
    # # print(res_sympy_function(**state_map))

    rule, elements, tree, variables, aliases = boolean_rule_from_str(rule)
