from abc import abstractmethod
from operator import add, sub, mul, truediv, pow
import sys
import re

# Boolean operator symbols
S_AND = '&'
S_OR  = '|'
S_NOT = '~' 
# Empty leaf symbol
EMPTY_LEAF = '@'

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
    sentence=[]
    for token in tokens:
        if token in "&|()":
            sentence.append(token)
        elif token in variables:
            sentence.append("True")
        else:
            sentence.append("False") 
    proposition=' '.join(sentence)    
    return eval(proposition)




def evaluate_expression_tree(expression, variables):
    """ Evaluates a logical expression (containing variables, 'and','or', 'not','(' , ')') 
        against the presence (True) or absence (False) of propositions within a list.
        Assumes the correctness of the expression.
    """
    t = build_tree(expression, Boolean)
    evaluator = BooleanEvaluator(variables)
    res = t.evaluate(evaluator.f_operand,evaluator.f_operator)
    return res




class Node(object):
    """
    Binary sintax tree node
    """
    def __init__(self,value,left = None,right = None):
        if left is not None and not isinstance(left,Node):
            raise ValueError("Invalid left element")
        if right is not None and not isinstance(right,Node):
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
        return (not self.left and not self.right )
    
    
    def is_empty_leaf(self):
        return self.value == EMPTY_LEAF


    def is_unary(self):
        return (self.left.is_empty_leaf() and not self.right.is_empty_leaf()) or  (not self.left.is_empty_leaf() and self.right.is_empty_leaf())


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
        
    
    def print_node(self, level = 0):
        tabs = ""
        for _ in range(level): tabs += "\t"
        if self.is_leaf():
            if self.value != EMPTY_LEAF:
                print(tabs, f"|____ \033[1;31m{self.value}\033[0;0m")    
            else: 
                pass
        else:
            print(tabs, f"|____[\033[;1m\033[94m{self.value}\033[0;0m]")
        if (self.left != None): 
            self.left.print_node(level+1)
        if (self.right != None): 
            self.right.print_node(level+1)

    

    
    
    def evaluate(self, f_operand = None, f_operator = None):
        """
        Evaluates the expression using the f_operand and f_operator mapping functions
        """
        if f_operand is None or f_operator is None:
            return eval(str(self))
        elif self.is_leaf():
            return f_operand(self.value)
        else:
            return f_operator(self.value)(self.left.evaluate(f_operand,f_operator),
			                              self.right.evaluate(f_operand,f_operator))




class Syntax:
    """Defines an interface for the tree syntax parsing 
       with operatores, their precedence and associativity.
    """
    @abstractmethod
    def is_operator(self,op):
        raise NotImplementedError

    @abstractmethod
    def is_greater_precedence(self,op1, op2):
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



class Aritmetic(Syntax):

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


class AritmeticEvaluator:

    @staticmethod
    def f_operator(op):
        operators = {'+':add, '-':sub, '*':mul, '/':truediv, '^':pow}
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
        pre = {S_AND: 0, S_OR: 0 , S_NOT : 1 }
        return pre[op1] >= pre[op2]

    @staticmethod
    def associativity(op):
        ass = {S_AND: 0, S_OR: 0 , S_NOT : 1}
        return ass[op]
    
    @staticmethod
    def arity(op):
        ar = {S_AND: 2, S_OR: 2, S_NOT: 1}
        return ar[op]

    @staticmethod
    def replace():
        r = {}
        r['not'] = [EMPTY_LEAF, S_NOT]
        r['and'] = [S_AND]
        r['or'] =  [S_OR]
        return r


class BooleanEvaluator:
    """A boolean evaluator
    
        arguments:
            true_list (list): operands evaluated as True. Operands not in the list are evaluated as False
        """
    def __init__(self, true_list = []):
        self.true_list = true_list

    def f_operator(self,op):
        operators = {S_AND: lambda x,y: x and y , S_OR: lambda x,y: x or y, S_NOT:lambda x , y : not y}
        if op in operators.keys():
            return operators[op]
        else:
            raise ValueError(f"Operator {op} not defined")

    def f_operand(self,op):
        if op.upper() == 'TRUE' or op =='1' or op in self.true_list:
            return True
        else:
            return False
    
    def set_true_list(self, l):
        self.true_list = l




class GeneEvaluator:
    """An evaluator for genes
    
        arguments:
            genes_value (dic): maps genes to values. 
                               If a gene is not in the dictionary, a default value of 1 is considered.
            and_operator : function to be applied instead of (and', '&') (not case sentitive) 
            or_operator : function to be applied instead of ('or','|') (not case sentitive)
    """
    def __init__(self, genes_value, and_operator, or_operator, prefix = ""):
        
        self.genes_value = genes_value
        self.and_operator = and_operator
        self.or_operator = or_operator
        self.prefix =prefix

    def f_operator(self,op):
        operators = {S_AND: self.and_operator, S_OR: self.or_operator}
        if op in operators.keys():
            return operators[op]
        else:
            raise ValueError(f"Operator {op} not defined")

    def f_operand(self,op):
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
    exp_list =[]
    for i in exp_:
        if i.lower() in replace_dic:
            exp_list.extend(replace_dic[i.lower()])
        else:
            exp_list.append(i)
    stack = []
    tree_stack = []
    for i in exp_list:
        if not (rules.is_operator(i) or i in ['(', ')']):
            t = Node(i)
            tree_stack.append(t)

        elif rules.is_operator(i):
            if not stack or stack[-1] == '(':
                stack.append(i)

            elif rules.is_greater_precedence(i, stack[-1]) and rules.associativity(i) == 1:
                stack.append(i)

            else:
                while stack and stack[-1] != '(' and rules.is_greater_precedence(stack[-1], i) and rules.associativity(i) == 0:
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



def test_1():
    # aritmetic example
    
    t = build_tree(" 1 + 2 + 3 + ( 2 * 2 )", Aritmetic)
    res = t.evaluate(AritmeticEvaluator.f_operand,AritmeticEvaluator.f_operator)
    print(t,' = ', res)
    
    # boolean example
    
    t = build_tree(" A or not (B)", Boolean)
    print(t)
    t.print_node()
    print(t.get_operands())
    # propositions in the list are evaluated as True and the remaining are False
    evaluator = BooleanEvaluator(['a','c'])
    res = t.evaluate(evaluator.f_operand,evaluator.f_operator)
    print('value = ',res)
    

    # gene values example
    t = build_tree("( a or not b ) and (c)", Boolean)    
    t.print_node()
    value_map = {'a':1,'b':2,'c':3,'d':4}
    op_and = lambda x,y: min(x,y)
    op_or = lambda x,y: max(x,y)
    evaluator = GeneEvaluator( value_map, op_and , op_or)
    res = t.evaluate( evaluator.f_operand, evaluator.f_operator)
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
    print("GENES:\n",genes)
    print("Evaluations:")
    evaluator = BooleanEvaluator(genes)
    res = t.evaluate(evaluator.f_operand,evaluator.f_operator)
    print(evaluator.true_list," ==> ",res)
    for _ in range(20):
        l = []
        n = random.randint(0,len(genes))
        for _ in range(n):
            i = random.randint(0,len(genes)-1)
            l.append(genes[i])
        evaluator.set_true_list(l)
        res = t.evaluate(evaluator.f_operand,evaluator.f_operator)
        print(evaluator.true_list," ==> ",res)
        

        
			


if __name__ == '__main__':
    
    test_1()
    #test_2()
    




