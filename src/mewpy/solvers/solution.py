# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
##############################################################################
LP Solution interface

Adapted by Vitor Pereira from Daniel Machado's REFRAMED
https://github.com/cdanielmachado/reframed
##############################################################################
"""
from mewpy.simulation import SStatus, get_simulator, Simulator, SimulationResult
from enum import Enum
import re


class Status(Enum):
    """ Enumeration of possible solution status. """
    OPTIMAL = 'Optimal'
    UNKNOWN = 'Unknown'
    SUBOPTIMAL = 'Suboptimal'
    UNBOUNDED = 'Unbounded'
    INFEASIBLE = 'Infeasible'
    INF_OR_UNB = 'Infeasible or Unbounded'

status_mapping={
    Status.OPTIMAL : SStatus.OPTIMAL,
    Status.UNKNOWN : SStatus.UNKNOWN,
    Status.SUBOPTIMAL : SStatus.SUBOPTIMAL,
    Status.UNBOUNDED : SStatus.UNBOUNDED,
    Status.INFEASIBLE : SStatus.INFEASIBLE,
    Status.INF_OR_UNB : SStatus.INF_OR_UNB
}

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

def to_simulation_result(model, objective_value, constraints, sim, solution, method=None):
    res = SimulationResult(model.model if isinstance(model, Simulator) else model,
                           objective_value,
                           status= status_mapping[solution.status],
                           fluxes=solution.values,
                           envcond=sim.environmental_conditions,
                           model_constraints=sim._constraints.copy(),
                           simul_constraints=constraints,
                           method=method
                           )                           
    return res

def print_values(value_dict, pattern=None, sort=False, abstol=1e-9):

    values = [(key, value) for key, value in value_dict.items() if abs(value) > abstol]

    if pattern:
        re_expr = re.compile(pattern)
        values = [x for x in values if re_expr.search(x[0]) is not None]

    if sort:
        values.sort(key=lambda x: x[1])

    entries = (f'{r_id:<12} {val: .6g}'for (r_id, val) in values)

    print('\n'.join(entries))


def print_balance(values, m_id, model, sort=False, percentage=False, abstol=1e-9):

    sim = get_simulator(model)
    inputs = sim.get_metabolite_producers(m_id)
    outputs = sim.get_metabolite_consumers(m_id)

    fwd_in = [(r_id, sim.get_reaction(r_id).stoichiometry[m_id] * values[r_id], '--> o')
              for r_id in inputs if values[r_id] > 0]
    rev_in = [(r_id, sim.get_reaction(r_id).stoichiometry[m_id] * values[r_id], 'o <--')
              for r_id in outputs if values[r_id] < 0]
    fwd_out = [(r_id, sim.get_reaction(r_id).stoichiometry[m_id] * values[r_id], 'o -->')
               for r_id in outputs if values[r_id] > 0]
    rev_out = [(r_id, sim.get_reaction(r_id).stoichiometry[m_id] * values[r_id], '<-- o')
               for r_id in inputs if values[r_id] < 0]

    flux_in = [x for x in fwd_in + rev_in if x[1] > abstol]
    flux_out = [x for x in fwd_out + rev_out if -x[1] > abstol]

    if sort:
        flux_in.sort(key=lambda x: x[1], reverse=True)
        flux_out.sort(key=lambda x: x[1], reverse=False)

    if percentage:
        turnover = sum([x[1] for x in flux_in])
        flux_in = [(x[0], x[1] / turnover, x[2]) for x in flux_in]
        flux_out = [(x[0], x[1] / turnover, x[2]) for x in flux_out]
        print_format = '[ {} ] {:<12} {:< 10.2%}'
    else:
        print_format = '[ {} ] {:<12} {:< 10.6g}'

    lines = (print_format.format(x[2], x[0], x[1]) for x in flux_in + flux_out)

    print('\n'.join(lines))
