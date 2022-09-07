"""
##############################################################################
File containing string exoression parsing utility functions.
Expressions are decomposed into parsing trees for infix evaluation.

Authors: Vitor Pereira
##############################################################################
"""
import re
import sys
from abc import abstractmethod
from operator import add, sub, mul, truediv, pow

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


def maybe_fn(f, v1, v2):
    """Maybe evaluator: if one of the arguments is None, it retuns the value of the other\
    argument. If both arguments are None, it returns None. If both arguments are not None\
    it returns the evaluation f(v1,v2).

    :param f: a function
    :param v1: the first argument
    :param v2: the second argument  
    """
    if v1 is None:
        return v2
    elif v2 is None:
        return v1
    else:
        return f(v1, v2)


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
        """Parameters are all non numeric symbols in an expression"""
        if self.is_leaf():
            if self.value != EMPTY_LEAF and not is_number(self.value):
                return {self.value}
            else:
                return set()
        else:
            return self.left.get_parameters().union(self.right.get_parameters())

    def print_node(self, level=0):
        """Prints a parsing tree of the expression."""
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
            return maybe_fn(f_operator(self.value),
                            self.left.evaluate(f_operand, f_operator),
                            self.right.evaluate(f_operand, f_operator)
                            )

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

    def __init__(self, genes_value, and_operator=min, or_operator=max, prefix="", unexpressed_value=1):

        self.genes_value = genes_value
        self.and_operator = and_operator
        self.or_operator = or_operator
        self.prefix = prefix
        self.unexpressed_value = unexpressed_value

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
            return self.unexpressed_value


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
