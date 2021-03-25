from .expression import Expression
from .symbolic import *
from .parsing import parse_expression
from .algebra_utils import solution_decode

# TODO: I have implemented an algebra sub-package that allows smooth modeling operations with variables
#  and algebra expressions of any kind. It is similar to sympy expression handling. It is also similar to
#  the boolean syntax parser and probably slower. However, it allows me to handle variables in this expressions
#  and so dynamically retrieve regulators/genes objects for reactions/interactions/targets in my models.
#  I don't know if there is an easy way to combine all expression parsing. Thus, I should create a method to create a
#  Node of a boolean syntax parser from my symbolic mathematics
