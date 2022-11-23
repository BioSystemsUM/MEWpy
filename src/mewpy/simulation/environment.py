""" 
##############################################################################
Environment module
Addapted from REFRAMED

Author: Vitor Pereira
##############################################################################
"""
from collections import OrderedDict
from math import inf
from types import FunctionType
from warnings import warn
from . import get_simulator
from .simulation import Simulator


class Environment(OrderedDict):
    """ This class represents the exchange of compounds between an organism and the environment. """

    def __init__(self):
        OrderedDict.__init__(self)

    def __str__(self):
        entries = (f"{r_id}\t{lb}\t{ub}" for r_id, (lb, ub) in self.items())
        return '\n'.join(entries)

    def __repr__(self):
        return str(self)

    def get_compounds(self, fmt_func=None):
        """
        Return the list of compounds in the growth medium for this environment.
        Args:
            fmt_func (str or function): python format string (see Notes)
        Returns:
            list: compounds in the medium
        Notes:
            The format function parameter is used to convert exchange reaction ids to metabolite ids.
            By default assumes BiGG notation ('R_EX_h2o_e' is transformed to 'h2o').
        """

        if fmt_func is None:
            def fmt_func(x):
                if x[:2] == 'R_':
                    return x[5:-2]
                else:
                    return x[3:-2]
        elif not isinstance(fmt_func, FunctionType):
            raise RuntimeError("fmt_func argument must be a string or function.")

        compounds = []

        for r_id, (lb, _) in self.items():
            if lb < 0:
                compounds.append(fmt_func(r_id))

        return compounds

    def apply(self, model, exclusive=True, inplace=True, warning=True, prefix=None):
        """
        Apply environmental conditions to a given model
        Args:
            model (REFRAME or COBRAPY model or Simulator): model
            exclusive (bool): block uptake of any model compounds not specified in this environment (default: True)
            warning (bool): print warning for exchange reactions not found in the model (default: True)
            inplace (bool): apply to model, otherwise return a constraints dict (default: True)
        """
        if not isinstance(model, Simulator):
            sim = get_simulator(model)
        else:
            sim = model

        if prefix:
            e = {f"{prefix}{rxn_id}": (lb, ub) for rxn_id, (lb, ub) in self.items()}
        else:
            e = self

        if exclusive:
            env = Environment.empty(model)
            env.update(e)
        else:
            env = e

        if not inplace:
            constraints = {}

        for r_id, (lb, ub) in env.items():
            if r_id in sim.reactions:
                if inplace:
                    sim.set_reaction_bounds(r_id, lb, ub)
                else:
                    constraints[r_id] = (lb, ub)
                    
        if not inplace:
            return constraints

    def simplify(self, inplace=False):
        """ Keep only uptake reactions for the respective medium. """

        if inplace:
            env = self
        else:
            env = self.copy()

        to_remove = []

        for r_id, (lb, _) in env.items():
            if lb >= 0:
                to_remove.append(r_id)

        for r_id in to_remove:
            del env[r_id]

        if not inplace:
            return env

    @staticmethod
    def from_reactions(reactions, max_uptake=10.0):
        """
        Create an environment from list of uptake reactions
        Arguments:
            reactions (list): exchange reactions
            max_uptake (float): maximum uptake rate for given compounds (default: 10.0)
        Returns:
            Environment
        """

        env = Environment()
        for r_id in reactions:
            env[r_id] = (-max_uptake, inf)

        return env

    @staticmethod
    def from_compounds(compounds, fmt_func=None, max_uptake=10.0, prefix=''):
        """
        Initialize environment from list of medium compounds
        Arguments:
            compounds (list): List of compounds present in the medium
            fmt_func (str or function): python format string (see Notes)
            max_uptake (float): maximum uptake rate for given compounds (default: 10.0)
        Returns:
            Environment
        Notes:
            Please provide a formatting string or function to convert compound identifiers to exchange reactions.
            Default format uses BiGG notation: "R_EX_{}_e" (for example "h2o" becomes "R_EX_h2o_e").
        """

        if fmt_func is None:
            def fmt_func(x):
                return f"{prefix}EX_{x}_e"
        elif isinstance(fmt_func, str):
            fmt_str = fmt_func
            def fmt_func(x): return fmt_str.format(x)
        elif not isinstance(fmt_func, FunctionType):
            raise RuntimeError("fmt_func argument must be a string or function.")

        reactions = map(fmt_func, compounds)

        return Environment.from_reactions(reactions, max_uptake=max_uptake)

    @staticmethod
    def from_model(model):
        """
        Extract environmental conditions from a given model
        Arguments:
            model (REFRAME or COBRAPY model or Simulator): model
        Returns:
            Environment: environment from provided model
        """

        if not isinstance(model, Simulator):
            sim = get_simulator(model)
        else:
            sim = model

        env = Environment()

        for r_id in sim.get_exchange_reactions():
            env[r_id] = sim.get_reaction_bounds(r_id)

        return env

    @staticmethod
    def from_defaults(model, max_uptake=10.0, max_secretion=inf, inplace=False):
        """
        Generate default environmental conditions for a given model
        Arguments:
            model (REFRAME or COBRAPY model or Simulator): model
            max_uptake (float): maximum uptake rate (default: 10.0)
            max_secretion (float): maximum secretion rate (default: 1000.0)
            inplace (bool): apply to model (default: False)
        Returns:
            Environment: Default environment for provided model
        """
        if not isinstance(model, Simulator):
            sim = get_simulator(model)
        else:
            sim = model

        env = Environment()

        for r_id in sim.get_exchange_reactions():
            env[r_id] = (-max_uptake, max_secretion)

        if inplace:
            env.apply(model, exclusive=False, inplace=True)
        else:
            return env

    @staticmethod
    def complete(model, max_uptake=10.0, inplace=False):
        """
        Generate a complete growth medium for a given model
        Arguments:
            model (REFRAME or COBRAPY model or Simulator): model
            max_uptake (float): maximum uptake rate (default: 1000.0)
            inplace (bool): apply to model (default: False)
        Returns:
            Environment: complete medium for provided model
        """

        return Environment.from_defaults(model, max_uptake=max_uptake, inplace=inplace)

    @staticmethod
    def empty(model, inplace=False):
        """
        Generate an empty growth medium for a given model
        Arguments:
            model (REFRAME or COBRAPY model or Simulator): model
            inplace (bool): apply to model (default: False)
        Returns:
            Environment: empty medium for provided model
        """

        return Environment.from_defaults(model, max_uptake=0, inplace=inplace)
