from typing import Union

from mewpy.solvers import get_default_solver
from mewpy.solvers.sglobal import __MEWPY_solvers__ as solvers
from mewpy.solvers.solver import Solver


integer_coefficients = ((0, 0), (1, 1), (0.0, 0.0), (1.0, 1.0), (0, 1), (0.0, 1.0))


def get_solver_instance(solver: Union[str, Solver] = None) -> Solver:
    """
    It returns a new empty mewpy solver instance. However, if a solver instance is provided,
    it only checks if it is a mewpy solver.
    :param solver: Solver, CplexSolver, GurobiSolver or OptLangSolver instance or name of the solver
    :return: a mewpy solver instance
    """
    if solver is None:
        solver_name = get_default_solver()

        SolverType = solvers[solver_name]

        solver = SolverType()

    elif isinstance(solver, str):

        SolverType = solvers.get(solver, None)

        if SolverType is None:
            raise ValueError(f'{solver} is not listed as valid solver. Check the valid solvers: {solvers}')

        solver = SolverType()

    elif isinstance(solver, Solver):

        pass

    else:
        raise ValueError(f'Invalid solver {solver}. Check the valid solvers: {solvers}')

    return solver


class Node:

    def __init__(self, value, length=None, idxes=None):

        if not length:
            length = 0

        if not idxes:
            idxes = None

        self._next = None
        self._previous = None

        self.value = value
        self.length = length
        self.idxes: slice = idxes

    def __str__(self):
        return self.value

    @property
    def next(self):
        return self._next

    @property
    def previous(self):
        return self._previous

    def unlink(self):
        self._next = None
        self._previous = None


class LinkedList:

    def __init__(self, *args):

        if args:
            head = args[0]
            tail = args[-1]
            nodes = list(args)
            nodes.append(None)
            nodes = list(zip(nodes[:-1], nodes[1:]))

        else:
            head = None
            tail = None
            nodes = []

        self._data = {}

        self._head = head
        self._tail = tail

        for node, next_node in nodes:
            node._next = next_node
            if next_node:
                next_node._previous = node

        self.build_data()

    @property
    def data(self):
        return self._data

    def __len__(self):

        res = 0

        if self._tail:

            res = self._data.get(self._tail.value).idxes.stop

        elif self._head:

            res = self._data.get(self._tail.value).idxes.stop

        if res > 0:
            return res

        return 0

    def __hash__(self):
        return self._data.__hash__()

    def __eq__(self, other):
        return self._data.__eq__(other)

    def __contains__(self, item):
        return self._data.__contains__(item)

    def __getitem__(self, item):

        return self._data.__getitem__(item).idxes

    def __setitem__(self, key, value):

        raise NotImplementedError('Linked lists do not support item setting. Try pop or add')

    def get(self, value, default=None):

        node = self._data.get(value, None)

        if node:
            return node.idxes

        return default

    def keys(self, unique=True):

        if unique:
            yield from self._data.keys()

            return

        for key, node in self._data.items():

            if node.idxes.stop - node.idxes.start > 1:

                for i in range(node.idxes.start, node.idxes.stop):
                    yield f'{key}_{i}'

            else:
                yield key

    def values(self):

        return (node.idxes for node in self._data.values())

    def items(self):

        return ((key, node.idxes) for key, node in self._data.items())

    def traverse(self):

        node = self._head

        while node is not None:
            yield node
            node = node.next

    def map(self, function):

        node = self._head

        while node is not None:
            function(node)
            node = node.next

    def get_node(self, value, default=None):

        return self._data.get(value, default)

    def build_data(self):

        self._data = {}

        node = self._head

        start = 0
        while node is not None:
            stop = start + node.length

            node.idxes = slice(start, stop)

            self._data[node.value] = node

            start = stop

            node = node.next

    def extend(self, nodes):

        for node in nodes:
            self.add(node)

    def add(self, node):

        if isinstance(node, (tuple, list)):
            node = Node(node[0], node[1])

        elif isinstance(node, dict):
            node = Node(node['value'], node['length'])

        elif isinstance(node, Node):
            pass

        else:
            raise TypeError('Node must be a tuple, list, dict(value=val, length=len) or Node instance')

        if node.value in self.data:
            raise ValueError('Node value is already in linked list')

        if not self._head:

            node._previous = None
            node._next = None

            self._head = node
            self._tail = node

            if not node.idxes:
                node.idxes = slice(0, node.length)

            self.data[node.value] = node

        else:

            if not node.idxes:
                # noinspection PyProtectedMember
                node.idxes = slice(self._tail.idxes.stop, self._tail.idxes.stop + node.length)

            self.data[node.value] = node

            node._previous = self._tail
            node._next = None

            self._tail._next = node
            self._tail = node

    def pop(self, value):

        if isinstance(value, Node):
            # noinspection PyUnresolvedReferences
            value = Node.value

        node = self._data.pop(value)
        previous_node = node.previous
        next_node = node.next

        if previous_node and next_node:

            previous_node._next = next_node
            next_node._previous = previous_node

            node._next = None
            node._previous = None

            start = previous_node.idxes.stop

        elif previous_node and not next_node:

            previous_node._next = None

            node._next = None
            node._previous = None

            self._tail = previous_node

            return node

        elif not previous_node and next_node:

            next_node._previous = None

            node._next = None
            node._previous = None

            self._head = next_node

            start = 0

        else:

            node._next = None
            node._previous = None

            self._head = None
            self._tail = None

            self._data = {}

            return node

        _node = next_node

        while _node is not None:
            stop = start + _node.length

            _node.idxes = slice(start, stop)

            self._data[_node.value] = _node

            start = stop

            _node = _node.next

        return node

    def clear(self):

        self._data = {}

        self.map(lambda n: n.unlink())

        self._tail = None
        self._head = None
