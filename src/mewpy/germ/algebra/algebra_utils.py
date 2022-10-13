from typing import Union, Dict, Any


def solution_decode(solution: Union[int, bool], decoder: Dict[Any, int] = None) -> int:
    """
    Decode a solution from expression or symbolic evaluation.
    Solution decoding is sometimes not straightforward due to python bitwise evaluations

    :param solution: a boolean or int result
    :param decoder: a custom dictionary for decoding the solution (key) into a given output (value)
    :return: it returns a solution decoded into the integer 0 or 1
    """

    if not decoder:
        decoder = {True: 1,
                   False: 0,
                   1: 1,
                   0: 0,
                   -1: 1,
                   -2: 0
                   }

    return decoder.get(solution, solution)


def _walk(symbolic, reverse=False):
    """
    Internal use!
    Helping function to walk through an algebra expression implemented by the symbolic objects.
    This function is heavily inspired by the sympy method preorder_traversal and postorder_traversal.
    Credits to sympy contributors: https://github.com/sympy/sympy

    :param symbolic: symbolic-like object
    :param reverse: order of the traversal
    :return:
    """
    if reverse:

        for child in symbolic.variables:
            for subtree in _walk(child, reverse):
                yield subtree

        yield symbolic

    else:

        yield symbolic

        for child in symbolic.variables:
            yield from _walk(child, reverse)
