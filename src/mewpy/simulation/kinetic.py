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
Kinetic simulation module

Author: Vitor Pereira
##############################################################################
"""
from multiprocessing import Process, Manager
from collections import OrderedDict
from mewpy.simulation.simulation import SimulationResult, SimulationInterface
from mewpy.model.kinetic import ODEModel
from mewpy.solvers import (KineticConfigurations,
                           SolverConfigurations,
                           ODEStatus,
                           ode_solver_instance)
import warnings
import numpy as np
from typing import List, Dict, Tuple, Union, TYPE_CHECKING

if TYPE_CHECKING:
    import pandas


def kinetic_solve(model: ODEModel,
                  y0: List[float],
                  time_steps: List[float],
                  parameters: Dict[str, float] = None,
                  factors: Dict[str, float] = None
                  ) -> Tuple[ODEStatus,
                             Dict['str', float],
                             Dict['str', float],
                             List[float],
                             List[float]]:
    """Kinetic solve method that invokes an available ODE solver. 

    :param model: The kinetic model
    :type model: ODEModel
    :param y0: vector of initial concentrations
    :type y0: List[float]
    :param time_steps: integration time steps
    :type time_steps: List[float]
    :param parameters: Parameters to be modified, defaults to None
    :type parameters: Dict[str,float], optional
    :param factors: factors to be applied to parameters, defaults to None
    :type factors: Dict[str, float], optional
    :return: _description_
    :rtype: _type_
    """

    rates = OrderedDict()
    f = model.get_ode(r_dict=rates, params=parameters, factors=factors)
    solver = ode_solver_instance(f, KineticConfigurations.SOLVER_METHOD)

    C, t, y = solver.solve(y0, time_steps)

    for c in C:
        if c < -1 * SolverConfigurations.RELATIVE_TOL:
            return ODEStatus.ERROR, {}, {}

    # values bellow solver precision will be set to 0
    rates.update({k: 0 for k, v in rates.items() if (
                  v < SolverConfigurations.ABSOLUTE_TOL
                  and v > - SolverConfigurations.ABSOLUTE_TOL)})
    conc = OrderedDict(zip(model.metabolites.keys(), C))
    return ODEStatus.OPTIMAL, rates, conc, t, y


class KineticThread(Process):
    """
    Solves the ODE inside a thread enabling to impose a timeout limit
    with thread.join(timeout)
    """

    def __init__(self,
                 model: ODEModel,
                 initial_concentrations: List[float] = None,
                 time_steps: List[float] = None,
                 parameters: Dict[str, float] = None,
                 factors: Dict[str, float] = None) -> None:
        """
        TSolves the ODE inside a thread enabling to impose a timeout limit
        with thread.join(timeout)

        :param model: The kinetic model
        :type model: ODEModel
        :param initial_concentrations: A list of initial concentrations, defaults to None
        :type initial_concentrations: List[float], optional
        :param time_steps: List of integration time steps, defaults to None
        :type time_steps: List[float], optional
        :param parameters: Kinetic parameters to be modified, defaults to None
        :type parameters: Dict[str, float], optional
        :param factors: Factors to be applied to kinetic parameters, defaults to None
        :type factors: Dict[str, float], optional
        """

        Process.__init__(self, daemon=False)
        self.model = model
        self.parameters = parameters
        self.factors = factors
        self.initial_concentrations = initial_concentrations
        self.time_steps = time_steps

        self.result = Manager().dict()
        self.result['status'] = None
        self.result["rates"] = None
        self.result["concentrations"] = None
        self.result["t"] = None
        self.result["y"] = None

    def run(self):
        try:
            status, rates, concentrations, t, y = kinetic_solve(self.model,
                                                                self.initial_concentrations,
                                                                self.time_steps,
                                                                self.parameters,
                                                                self.factors)
            self.result['status'] = status
            self.result["rates"] = rates
            self.result["concentrations"] = concentrations
            self.result["t"] = t
            self.result["y"] = y
        except Exception:
            warnings.warn('Timeout')
        return


class KineticSimulationResult(SimulationResult):

    def __init__(self,
                 model: ODEModel,
                 status: ODEStatus,
                 factors: Dict[str, float] = None,
                 rates: Dict[str, float] = None,
                 concentrations: List[float] = None,
                 t: List[float] = None,
                 y: List[float] = None) -> None:
        """Result class of a kinetic simulation

        :param model: The kinetic model
        :type model: ODEModel
        :param status: The solve status
        :type status: ODEStatus
        :param factors: factors used in the simulation, defaults to None
        :type factors: Dict[str, float], optional
        :param rates: _description_, defaults to None
        :type rates: Dict[str, float], optional
        :param concentrations: _description_, defaults to None
        :type concentrations: List[float], optional
        :param t: integration time points, defaults to None
        :type t: List[float], optional
        :param y: _description_, defaults to None
        :type y: List[float], optional
        """
        super(KineticSimulationResult, self).__init__(model, None, fluxes=rates, status=status)
        self.factors = factors
        self.concentrations = concentrations
        self.t = t
        self.y = y
        self.m_indexes = {k: v for v, k in enumerate(concentrations.keys())}

    def get_y(self, m_id):
        if m_id in self.m_indexes:
            return np.array(self.y).T[:, self.m_indexes[m_id]]
        else:
            raise ValueError(f"Unknown metabolite {m_id}")

    def get_concentrations(self, format: str = None) -> Union["pandas.DataFrame", Dict[str, float]]:
        """_summary_

        :param format:The output format ("df" or None), defaults to None
        :type format: str, optional
        :return: the steady-state metabolite concentrations
        :rtype: _type_
        """
        if format and format == 'df':
            import pandas as pd
            return pd.DataFrame(self.concentrations)
        else:
            return self.concentrations

    def plot(self, met: List[str] = None, size: Tuple[int, int] = None):
        import matplotlib.pyplot as plt
        if size:
            plt.rcParams["figure.figsize"] = size
        if not met:
            _mets = list(self.concentrations.keys())
        elif isinstance(met, str):
            _mets = [met]
        elif isinstance(met, list) and len(met) <= 4:
            _mets = met
        else:
            raise ValueError('fluxes should be a reaction identifier,' 
                             'a list of reaction identifiers or None.')
        ax = plt.subplot()
        if len(_mets) != 2:
            for k in _mets:
                ax.plot(self.t, self.get_y(k), label=k)
            if len(_mets) == 1:
                ax.set_ylabel(self.model.get_metabolite(_mets[0]).name)
            else:
                ax.set_ylabel('Concentrations')
                plt.legend()
        else:
            ax.plot(self.t, self.get_y(_mets[0]), label=_mets[0])
            ax2 = plt.twinx(ax)
            ax2.plot(self.t, self.get_y(_mets[1]), label=_mets[1], color='r')
            ax.set_ylabel(self.model.get_metabolite(_mets[0]).name, color='b')
            ax2.set_ylabel(self.model.get_metabolite(_mets[1]).name, color='r')

        ax.set_xlabel('Time')
        return ax


class KineticSimulation(SimulationInterface):

    def __init__(self,
                 model: ODEModel,
                 parameters: Dict[str, float] = None,
                 t_points: List[float] = [0, 1e9],
                 timeout: int = KineticConfigurations.SOLVER_TIMEOUT) -> None:
        """Class that runs kinetic simulations

        :param model: The kinetic model
        :type model: ODEModel
        :param parameters: Dictionary of modified kinetic parameter, defaults to None
             in which case the parameter values in the model are used.
        :type parameters: Dict[str, float], optional
        :param t_points: the integration time points or span, defaults to [0, 1e9]
        :type t_points: List[float], optional
        :param timeout: The integration timeout, defaults to KineticConfigurations.SOLVER_TIMEOUT
        :type timeout: int, optional
        """
        if not isinstance(model, ODEModel):
            raise ValueError('model is not an instance of ODEModel.')
        self.model = model
        self.t_points = t_points
        self.timeout = timeout
        self.parameters = parameters if parameters else dict()

    def get_initial_concentrations(self, initcon: Dict[str, float] = None):
        values = []
        _initcon = initcon if initcon else dict()
        for i, m in enumerate(self.model.metabolites):
            try:
                values.append(_initcon.get(m, self.model.concentrations[m]))
            except:
                values.append(None)
        return values

    def set_time(self, start: int, end: int, steps: int):
        """
        This function sets the time parameters for the model.

        :param int start: the start time - usually 0
        :param int end: the end time (default is 100)
        :param int steps: the number of timepoints for the output
        """
        self.t_points = np.linspace(start, end, steps)

    def get_time_points(self):
        """Returns the time point or span."""
        return self.t_points

    def simulate(self,
                 parameters: Dict[str, float] = None,
                 initcon: List[float] = None,
                 factors: Dict[str, float] = None,
                 t_points: List[float] = None) -> KineticSimulationResult:
        """
        Solve an initial value problem for a system of ODEs.

        :param dict parameters: Parameters to be modified. Default None
        :para dict initcon: initial conditions, metabolite concentrations. Default None
        :param dict factors: Modification over the kinetic model.
        :param list t_points: Times at which to store the computed solution,\
            must be sorted and lie within t_span. Default None, in which case the number of
            time steps is defined by SolverConfigurations.N_STEPS.
        :returns: Returns a kineticSimulationResult with the steady-state flux distribution and concentrations.
        """

        _factors = factors if factors is not None else {}
        initConcentrations = self.get_initial_concentrations(initcon)

        status = None
        sstateRates = None
        sstateConc = None
        t = None
        y = None
        params = self.parameters
        if parameters:
            params.update(parameters)

        time_steps = t_points if t_points else self.get_time_points()

        if len(time_steps) == 2:
            time_steps = np.linspace(time_steps[0],
                                     time_steps[1],
                                     num=SolverConfigurations.N_STEPS,
                                     endpoint=True)

        if self.timeout:
            try:
                th = KineticThread(self.model,
                                   initial_concentrations=initConcentrations,
                                   time_steps=time_steps,
                                   parameters=params,
                                   factors=_factors)

                th.start()
                th.join(self.timeout)
                status = th.result['status']
                sstateRates = th.result['rates']
                sstateConc = th.result['concentrations']
                t = th.result['t']
                y = th.result['y']
            except AssertionError as e:
                raise AssertionError(f"{str(e)}. Installing ray for multiprocessing will solve this issue.")
            except Exception as e:
                warnings.warn(str(e))
        else:
            status, sstateRates, sstateConc, t, y = kinetic_solve(self.model,
                                                                  initConcentrations,
                                                                  time_steps,
                                                                  params,
                                                                  _factors)

        return KineticSimulationResult(self.model, status, factors=_factors, rates=sstateRates,
                                       concentrations=sstateConc, t=t, y=y)
