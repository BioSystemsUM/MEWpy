import threading
from collections import OrderedDict
import warnings
from ..simulation.simulation import SimulationResult, SimulationInterface
from ..model.kinetic import ODEModel
from ..solvers import KineticConfigurations, SolverConfigurations, ODEStatus, ode_solver_instance


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


class KineticThread(threading.Thread):
    """
    Solves the ODE inside a thread enabling to impose a timeout limit with thread.join(timeout)
    """

    def __init__(self, model, initial_concentrations=None, time_steps=None, parameters=None, factors=None):
        super(KineticThread, self).__init__()
        self.model = model
        self.parameters = parameters
        self.factors = factors
        self.initial_concentrations = initial_concentrations
        self.time_steps = time_steps
        self.status = None
        self.rates = None
        self.concentrations = None
        self.t= None
        self.y = None

    def run(self):
        try:
            self.status, self.rates, self.concentrations, self.t, self.y = kinetic_solve(self.model,
                                                                         self.initial_concentrations,
                                                                         self.time_steps,
                                                                         self.parameters,
                                                                         self.factors)
        except Exception as e:
            warnings.warn(f"{e}")
        return


class kineticSimulationResult(SimulationResult):

    def __init__(self, model, status, factors=None, rates=None, concentrations=None, t=None, y=None):
        super(kineticSimulationResult, self).__init__(model, None, fluxes=rates, status=status)
        self.factors = factors
        self.concentations = concentrations
        self.t = t
        self.y = y


class KineticSimulation(SimulationInterface):

    def __init__(self, model,  parameters=None, tSteps=[0, 1e9], timeout=KineticConfigurations.SOLVER_TIMEOUT):
        if not isinstance(model, ODEModel):
            raise ValueError('model is not an instance of ODEModel.')
        self.model = model
        self.tSteps = tSteps
        self.timeout = timeout
        self.parameters = parameters if parameters else dict()

    def get_initial_concentrations(self, initcon=None):
        values = []
        _initcon = initcon if initcon else dict()
        for i, m in enumerate(self.model.metabolites):
            try:
                values.append(_initcon.get(m, self.model.concentrations[m]))
            except:
                values.append(None)
        return values

    def get_time_steps(self):
        return self.tSteps

    def simulate(self, parameters=None, initcon=None, factors=None):
        """
        This method preform the phenotype simulation of the kinetic model, using the solverId method and applying
        the modifications present in the instance of overrideSimulProblem.

        :param dict factors: Modification over the kinetic model.
        :returns: Returns a kineticSimulationResult with the steady-state flux distribution and concentrations.
        """

        _factors = factors if factors is not None else {}
        initConcentrations = self.get_initial_concentrations(initcon)

        status = None
        sstateRates = None
        sstateConc = None
        params = self.parameters
        if parameters:
            params.update(parameters)

        if self.timeout:
            th = KineticThread(self.model,
                               initial_concentrations=initConcentrations,
                               time_steps=self.get_time_steps(),
                               parameters=params,
                               factors=_factors
                               )

            th.start()
            th.join(self.timeout)
            if th.is_alive():
                th.join()
            status = th.status
            sstateRates = th.rates
            sstateConc = th.concentrations
            t = th.t
            y = th.y

        else:
            status, sstateRates, sstateConc, t, y = kinetic_solve(self.model,
                                                            initConcentrations,
                                                            self.get_time_steps(),
                                                            params,
                                                            _factors)

        return kineticSimulationResult(self.model, status, factors=_factors, rates=sstateRates,
                                       concentrations=sstateConc, t=t, y=y )
