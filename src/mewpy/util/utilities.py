import joblib
import contextlib
import functools
import re
import types
import time
from collections.abc import Iterable
from .constants import atomic_weights
from warnings import warn

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

    def find(self, pattern=None, sort=False):
        """A user friendly method to find metabolites, reactions or genes in the model.

        :param pattern: The pattern which can be a regular expression, defaults to None in which case all entries are listed.
        :type pattern: str, optional
        :param sort: if the search results should be sorted, defaults to False
        :type sort: bool, optional
        :return: the search results
        :rtype: pandas dataframe
        """
        values = list(self.keys())
        if pattern:
            import re
            if isinstance(pattern, list):
                patt = '|'.join(pattern)
                re_expr = re.compile(patt)
            else:
                re_expr = re.compile(pattern)
            values = [x for x in values if re_expr.search(x) is not None]
        if sort:
            values.sort()

        import pandas as pd
        data = [{'attribute':x,'value':self.get(x)} for x in values]
        
        if data:
            df = pd.DataFrame(data)
            df = df.set_index(df.columns[0])
        else: 
            df = pd.DataFrame()
        return df

    def __repr__(self) -> str:
        return str(self.find())


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

def molecular_weight(formula, element=None):
    elems = elements(formula)
    if element:
        mw = elems.get(element, 0) * atomic_weights.get(element, 0)
    else:
        missing = set(elems) - set(atomic_weights)
        if missing:
            warn(f"Atomic weight not listed for elements: {missing}")

        mw = sum(atomic_weights.get(elem, 0) * n for elem, n in elems.items())

    return mw

@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()


def get_all_subclasses(cls):
    '''Returns all subclasses of a class'''
    all_subclasses = []

    for subclass in cls.__subclasses__():
        all_subclasses.append(subclass)
        all_subclasses.extend(get_all_subclasses(subclass))

    return all_subclasses

def is_notebook() -> bool:
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter
