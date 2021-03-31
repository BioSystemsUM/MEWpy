import re
import sys
from abc import abstractmethod
from operator import add, sub, mul, truediv, pow

from sympy import Symbol
from sympy.parsing.sympy_parser import parse_expr

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
    The evaluation is achieved usind python native eval function.

    :param str expression: The expression to be evaluated.
    :param list variables: List of variables to be evaluated as True.
    :returns: A boolean evaluation of the expression.

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
    """ Evaluates a logical expression (containing variables, 'and','or', 'not','(' , ')') against the presence (True)
    or absence (False) of propositions within a list.
    Assumes the correctness of the expression. The evaluation is achieved by means of a parsing tree.

    :param str expression: The expression to be evaluated.
    :param list variables: List of variables to be evaluated as True.
    :returns: A boolean evaluation of the expression.

    """
    t = build_tree(expression, Boolean)
    evaluator = BooleanEvaluator(variables)
    res = t.evaluate(evaluator.f_operand, evaluator.f_operator)
    return res


class Node(object):
    """
    Binary syntax tree node.

    :param value: The node value.
    :param left: The left node or None.
    :param right: The right node or None.
    """

    def __init__(self, value, left=None, right=None, tp=0):
        """Binary tree.  Empty leafs with None or 'EMPTY_LEAF' are always left sided.

        Args:
            value : Node value
            left (Node, optional): [description]. Defaults to None.
            right (Node, optional): [description]. Defaults to None.
            tp (int, optional): 1-Unary function, 2-binary function. Defaults to 0.

        Raises:
            ValueError
        """
        if left is not None and not isinstance(left, Node):
            raise ValueError("Invalid right element")
        if right is not None and not isinstance(right, Node):
            raise ValueError("Invalid right element")
        self.value = value
        self.left = left
        self.right = right
        self.tp = tp

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if self.is_leaf():
            return str(self.value)
        else:
            return f"{str(self.value)} ( {str(self.left)} , {str(self.right)} )"

    def is_leaf(self):
        """
        :returns: True if the node is a leaf False otherwise. Both left and right are None.
        """
        return not self.left and not self.right

    def is_empty_leaf(self):
        """
        :returns: True if the node is an empty leaf False otherwise.
        """
        return self.value == EMPTY_LEAF

    def is_unary(self):
        """
        :returns: True if the node is a unary operation False otherwise.
        """
        return (self.left.is_empty_leaf() and not self.right.is_empty_leaf()) or (
            not self.left.is_empty_leaf() and self.right.is_empty_leaf())

    def is_binary(self):
        """
        :returns: True if the node is a binary operation False otherwise.
        """
        return not self.left.is_empty_leaf() and not self.right.is_empty_leaf()

    def get_operands(self):
        if self.is_leaf():
            if self.value == EMPTY_LEAF:
                return set()
            else:
                return {self.value}
        else:
            return self.left.get_operands().union(self.right.get_operands())

    def get_parameters(self):
        if self.is_leaf():
            if self.value != EMPTY_LEAF and not is_number(self.value):
                return {self.value}
            else:
                return set()
        else:
            return self.left.get_parameters().union(self.right.get_parameters())

    def print_node(self, level=0):
        tabs = ""
        for _ in range(level):
            tabs += "\t"
        if self.is_leaf():
            if self.value != EMPTY_LEAF:
                print(tabs, f"|____{self.value}")
            else:
                pass
        else:
            print(tabs, f"|____{self.value}")
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

    def replace(self, r_map: dict):
        """Apply the mapping replacing to the tree

        Args:
            map (dict): replacement mapping

        Returns:
            Node: new Tree with replace entries
        """
        if self.is_leaf():
            v = r_map[self.value] if self.value in r_map.keys() else self.value
            return Node(v, None, None)
        else:
            return Node(self.value, self.left.replace(r_map), self.right.replace(r_map), self.tp)

    def to_infix(self) -> str:
        if self.is_leaf():
            if self.value == EMPTY_LEAF:
                return ''
            else:
                return str(self.value)
        elif self.tp == 2:
            return ''.join([self.value, '( ', self.left.to_infix(), ' , ', self.right.to_infix(), ' )'])
        elif self.tp == 1:
            return ''.join([self.value, '( ', self.right.to_infix(), ' )'])
        else:
            return ''.join(['( ', self.left.to_infix(), ' ', self.value, ' ', self.right.to_infix(), ' )'])

    def copy(self):
        if self.is_leaf():
            return Node(self.value.copy(), None, None)
        else:
            return Node(self.value.copy(), self.left.copy(), self.right.copy(), self.tp)


class Syntax:
    """Defines an interface for the tree syntax parsing
       with operators, their precedence and associativity.
    """

    @ abstractmethod
    def is_operator(self, op):
        raise NotImplementedError

    @ abstractmethod
    def is_greater_precedence(self, op1, op2):
        raise NotImplementedError

    @ abstractmethod
    def associativity(self, op):
        raise NotImplementedError

    @ staticmethod
    def arity(op):
        return 2

    @ staticmethod
    def replace():
        return {}


class Arithmetic(Syntax):
    """Defines a basic arithmetic sintax.
    """

    @ staticmethod
    def is_operator(op):
        return op in ['+', '-', '*', '/', '^']

    @ staticmethod
    def is_greater_precedence(op1, op2):
        pre = {'+': 0, '-': 0, '*': 1, '/': 1, '^': 2}
        return pre[op1] >= pre[op2]

    @ staticmethod
    def associativity(op):
        ass = {'+': 0, '-': 0, '*': 0, '/': 0, '^': 1}
        return ass[op]

    @ staticmethod
    def arity(op):
        ar = {'+': 2, '-': 2, '*': 2, '/': 2, '^': 2}
        return ar[op]


class ArithmeticEvaluator:

    @ staticmethod
    def f_operator(op):
        operators = {'+': add, '-': sub, '*': mul, '/': truediv, '^': pow, 'pow': pow}
        if op in operators.keys():
            return operators[op]
        else:
            raise ValueError(f"Operator {op} not defined")

    @ staticmethod
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

    @ staticmethod
    def is_operator(op):
        return op in [S_AND, S_OR, S_NOT]

    @ staticmethod
    def is_greater_precedence(op1, op2):
        pre = {S_AND: 0, S_OR: 0, S_NOT: 1}
        return pre[op1] >= pre[op2]

    @ staticmethod
    def associativity(op):
        ass = {S_AND: 0, S_OR: 0, S_NOT: 1}
        return ass[op]

    @ staticmethod
    def arity(op):
        ar = {S_AND: 2, S_OR: 2, S_NOT: 1}
        return ar[op]

    @ staticmethod
    def replace():
        r = {'not': [EMPTY_LEAF, S_NOT], 'and': [S_AND],
             'or': [S_OR], '&': [S_AND], '|': [S_OR], '~': [EMPTY_LEAF, S_NOT]}
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
    """A boolean evaluator.

    :param list true_list: Operands evaluated as True. Operands not in the list are evaluated as False
    :param dict variables: A dictionary mapping symbols to values. Used to evaluate conditions.

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
    """An evaluator for genes expression.

    :param genes_value: A dictionary mapping genes to values. Gene not in the dictionary have a value of 1.
    :param and_operator: function to be applied instead of (and', '&') (not case sensitive)
    :param or_operator: function to be applied instead of ('or','|') (not case sensitive)
    """

    def __init__(self, genes_value, and_operator=min, or_operator=max, prefix=""):

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


def tokenize_function(exp: str):
    p = 0
    s = -1
    tokens = []
    i = 0
    while i < len(exp):
        if exp[i] == '(':
            s += 1
            if s == 0:
                p = i+1
            if not tokens:
                tokens = [exp[:i]]
        elif exp[i] == ')':
            s -= 1
            if s == -1:
                tokens.append(exp[p:i])
        elif exp[i] == ',' and s == 0:
            tokens.append(exp[p:i])
            p = i+1
        i += 1
    if not tokens:
        return [exp]
    else:
        return tokens


# Tree
def build_tree(exp, rules):
    """ Builds a parsing syntax tree for basic mathematical expressions

    :param exp: the expression to be parsed
    :param rules: Sintax definition rules
    """
    replace_dic = rules.replace()
    exp_ = tokenize_infix_expression(exp)
    exp_list = []
    for token in exp_:
        if token.lower() in replace_dic:
            exp_list.extend(replace_dic[token.lower()])
        else:
            exp_list.append(token)
    stack = []
    tree_stack = []
    predecessor = None
    i = 0
    while i < len(exp_list):
        token = exp_list[i]
        if not (rules.is_operator(token) or token in ['(', ')']):
            if i < len(exp_list)-2 and exp_list[i+1] == '(':
                s = 1
                p = i+2
                while p < len(exp_list) and s > 0:
                    if exp_list[p] == '(':
                        s += 1
                    elif exp_list[p] == ')':
                        s -= 1
                    p += 1
                token = ' '.join(exp_list[i:p])
                i = p - 1

            if predecessor and not (rules.is_operator(predecessor) or predecessor in ['(', ')']):
                s = tree_stack[-1].value
                tree_stack[-1].value = s + " " + token
            else:
                if '(' in token:
                    f = tokenize_function(token)
                    if len(f) == 2:
                        t = Node(f[0], Node(EMPTY_LEAF), build_tree(f[1], rules), 1)
                    elif len(f) == 3:
                        t = Node(f[0], build_tree(f[1], rules), build_tree(f[2], rules), 2)
                    else:
                        t = Node(token)
                else:
                    t = Node(token)
                tree_stack.append(t)
        elif rules.is_operator(token):
            if not stack or stack[-1] == '(':
                stack.append(token)

            elif rules.is_greater_precedence(token, stack[-1]) and rules.associativity(token) == 1:
                stack.append(token)

            else:
                while stack and stack[-1] != '(' and rules.is_greater_precedence(stack[-1], token) \
                        and rules.associativity(token) == 0:
                    popped_item = stack.pop()
                    t = Node(popped_item)
                    t1 = tree_stack.pop()
                    t2 = tree_stack.pop()
                    t.right = t1
                    t.left = t2
                    tree_stack.append(t)
                stack.append(token)

        elif token == '(':
            stack.append('(')

        elif token == ')':
            while stack[-1] != '(':
                popped_item = stack.pop()
                t = Node(popped_item)
                t1 = tree_stack.pop()
                t2 = tree_stack.pop()
                t.right = t1
                t.left = t2
                tree_stack.append(t)
            stack.pop()
        predecessor = token
        i += 1

    while stack:
        popped_item = stack.pop()
        t = Node(popped_item)
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
    :param special_chars: dict or None, special (or not) chars to be replaced. BOOL_SPECIAL_CHARS by default \
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
    :param local_dict: dict, mapping between the variables in the expression and a Symbol Sympy's object
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
    :param special_chars: dict, dict or None, special (or not) chars to be replaced. BOOL_SPECIAL_CHARS is used by \
        default. special_chars argument must have unique keys and unique values!!
    :param replaces: dict, replaces occurred for aliases construction. None is the default.
    :return: regulatory variable Sympy Symbol, regulatory variable
    :return: alias: dict, the new variable name map to the old one

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
                          r'([a-zA-Z0-9_]+\s' + regex + r'\s[0-9.]+)' + '|' +
                          '([0-9.]+' + regex + '[a-zA-Z0-9_]+)' + '|' +
                          r'([0-9.]+\s' + regex + r'\s[a-zA-Z0-9_]+)'
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
    :return symbols: dict, all regulatory variables in the regulatory rule. Keys stand for their new IDs \
        and values stand for their sympy symbols
    :return symbols: dict, all regulatory variables in the regulatory rule. Keys stand for their new IDs \
        and values stand for their aliases
    """

    regex_str = '|'.join(['(' + regex + '[0-9.]+)' + '|' +
                          '(' + regex + r'\s[0-9.]+)' + '|'
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
