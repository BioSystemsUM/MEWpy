from multiprocessing import Process, Manager
from collections import OrderedDict
from mewpy.simulation.simulation import SimulationResult, SimulationInterface
from mewpy.model.kinetic import ODEModel
from mewpy.solvers import KineticConfigurations, SolverConfigurations, ODEStatus, ode_solver_instance
import warnings
import numpy as np

def kinetic_solve(model, y0, time_steps, parameters=None, factors=None):

    rates = OrderedDict()
    f = model.get_ode(r_dict=rates, params=parameters, factors=factors)
    solver = ode_solver_instance(f, KineticConfigurations.SOLVER_METHOD)

    C, t, y = solver.solve(y0, time_steps)

    for c in C:
        if c < -1 * SolverConfigurations.RELATIVE_TOL:
            return ODEStatus.ERROR, {}, {}

    # values bellow solver precision will be set to 0
    rates.update({k: 0 for k, v in rates.items() if
                  v < SolverConfigurations.ABSOLUTE_TOL and v > - SolverConfigurations.ABSOLUTE_TOL})
    conc = OrderedDict(zip(model.metabolites.keys(), C))
    return ODEStatus.OPTIMAL, rates, conc, t, y


class KineticThread(Process):
    """
    Solves the ODE inside a thread enabling to impose a timeout limit with thread.join(timeout)
    """

    def __init__(self, model, initial_concentrations=None, time_steps=None, parameters=None, factors=None):
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

    def __init__(self, model, status, factors=None, rates=None, concentrations=None, t=None, y=None):
        super(KineticSimulationResult, self).__init__(model, None, fluxes=rates, status=status)
        self.factors = factors
        self.concentrations = concentrations
        self.t = t
        self.y = y
        self.m_indexes = {k: v for v, k in enumerate(concentrations.keys())} 

    def get_y(self, m_id):
        if m_id in self.m_indexes:
            return np.array(self.y).T[:,self.m_indexes[m_id]]
        else:
            raise ValueError(f"Unknown metabolite {m_id}")

    def get_ss_concentrations(self, format=None):
        if format and format=='df':
            import pandas as pd
            return pd.DataFrame(self.concentrations)
        else:
            return self.concentrations

    def plot(self, met=None, size:tuple=None):    
        import matplotlib.pyplot as plt
        if size:
            plt.rcParams["figure.figsize"] = size
        if not met:
            _mets = list(self.concentrations.keys()) 
        elif isinstance(met,str):
            _mets =[met]
        elif isinstance(met,list) and len(met)<=4:
            _mets = met
        else:
            raise ValueError('fluxes should be a reaction identifier, a list of reaction identifiers or None.')
        ax = plt.subplot()
        if len(_mets)!=2: 
            for k in _mets:
                ax.plot(self.t, self.get_y(k), label=k)
            if len(_mets)==1:
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

    def __init__(self, model,  parameters=None, t_points=[0, 1e9], 
                 timeout=KineticConfigurations.SOLVER_TIMEOUT):
        if not isinstance(model, ODEModel):
            raise ValueError('model is not an instance of ODEModel.')
        self.model = model
        self.t_points = t_points
        self.timeout = timeout
        self.parameters = parameters if parameters else dict()

    def get_initial_concentrations(self, initcon=None):
        values = []
        _initcon = initcon if initcon else dict()
        for i, m in enumerate(self.model.metabolites):
            try:
                values.append(_initcon.get(m, self.model.concentrations[m]))
            except :
                values.append(None)
        return values

    def set_time(self, start, end, steps):
        """
        This function sets the time parameters for the model.  This is how long the model will simulate

        Args:
            start (int): the start time - usually 0
            end (int): the end time (default is 100)
            steps (int): the number of timepoints for the output
        """
        self.t_points = np.linspace(start, end, steps)


    def get_time_points(self):
        """Returns the time point or span."""
        return self.t_points

    def simulate(self, parameters=None, initcon=None, factors=None, t_points=None):
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
        
        if len(time_steps)==2:
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
