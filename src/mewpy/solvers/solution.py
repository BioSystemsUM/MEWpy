""" Adapted from Daniel Machado's REFRAMED"""
from enum import Enum
import re
from ..simulation import get_simulator
from ..simulation.simulation import Simulator


class Status(Enum):
    """ Enumeration of possible solution status. """
    OPTIMAL = 'Optimal'
    UNKNOWN = 'Unknown'
    SUBOPTIMAL = 'Suboptimal'
    UNBOUNDED = 'Unbounded'
    INFEASIBLE = 'Infeasible'
    INF_OR_UNB = 'Infeasible or Unbounded'


class Solution(object):
    """ Stores the results of an optimization.

    Instantiate without arguments to create an empty Solution representing a failed optimization.
    """

    def __init__(self, status=Status.UNKNOWN, message=None, fobj=None, values=None,
                 shadow_prices=None, reduced_costs=None):
        self.status = status
        self.message = message
        self.fobj = fobj
        self.values = values
        self.shadow_prices = shadow_prices
        self.reduced_costs = reduced_costs

    def __str__(self):
        return f"Objective: {self.fobj}\nStatus: {self.status.value}\n"

    def __repr__(self):
        return str(self)

    def to_dataframe(self):
        """ Convert reaction fluxes to *pandas.DataFrame*

        Returns:
            pandas.DataFrame: flux values
        """
        try:
            import pandas as pd
        except ImportError:
            raise RuntimeError("Pandas is not installed.")

        return pd.DataFrame(self.values.values(), columns=["value"], index=self.values.keys())
