import functools
import re
import types
import time
from collections import Iterable


class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""


class Timer:
    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError("Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError("Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None
        print(f"Elapsed time: {elapsed_time:0.6f} seconds")

    def __enter__(self):
        """Start a new timer as a context manager"""
        self.start()
        return self

    def __exit__(self, *exc_info):
        """Stop the context manager timer"""
        self.stop()


class Singleton(object):
    _instance = None

    def __new__(class_, *args, **kwargs):
        if not isinstance(class_._instance, class_):
            class_._instance = object.__new__(class_, *args, **kwargs)
        return class_._instance


def copy_func(f):
    """Based on http://stackoverflow.com/a/6528148/190597 (Glenn Maynard)"""
    g = types.FunctionType(f.__code__, f.__globals__, name=f.__name__,
                           argdefs=f.__defaults__,
                           closure=f.__closure__)
    g = functools.update_wrapper(g, f)
    g.__kwdefaults__ = f.__kwdefaults__
    return g


class Dispatcher:

    def __init__(self):
        """
        Dispatcher for the simulate method of the Simulation interface
        It allows a simplification of the if else chain of methods provided as input to the simulate method

        based on https://stackoverflow.com/questions/36836161/singledispatch-based-on-value-instead-of-type
        (Ilja Everilä)

        """

        # weak ref key for the garbage collector
        self.registry = {}

    def __get__(self, instance, owner):

        if instance is None:
            return self

        return self.dispatch(instance, owner)

    def dispatch(self, instance, owner):

        def wrapper(state, *args, **kwargs):

            method = self.registry.get(state).__get__(instance, owner)

            return method(*args, **kwargs)

        return wrapper

    def register(self, state):

        def wrapper(method):

            self.registry[state] = method

            return method

        return wrapper


def iterable(obj, is_string=False):
    if isinstance(obj, Iterable):

        if is_string and isinstance(obj, str):
            return (obj,)

        return obj

    return (obj,)


def generator(container):
    return (value for value in container.values())


# Taken from the talented team responsible for developing cobrapy!!!!
chemical_formula_re = re.compile('([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)')


def elements(formula):
    all_elements = re.findall(chemical_formula_re, formula)
    atoms = {}
    for atom, count in all_elements:
        if not count:
            count = '1'
        atoms[atom] = atoms.get(atom, 0) + int(count)
    return atoms
