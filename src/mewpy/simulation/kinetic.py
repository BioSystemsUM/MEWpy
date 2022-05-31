import threading
from collections import OrderedDict
import warnings
from ..simulation.simulation import SimulationResult
from ..model.kinetic import ODEModel
from ..solvers import KineticConfigurations, SolverConfigurations, ODEStatus, ode_solver_instance


def kinetic_solve(model, y0, time_steps, parameters=None, factors=None):

    rates = OrderedDict()
    f = model.get_ode(r_dict=rates, params=parameters, factors=factors)
    solver = ode_solver_instance(f, KineticConfigurations.SOLVER_METHOD)

    C, t, X = solver.solve(y0, time_steps)

    for c in C:
        if c < -1 * SolverConfigurations.RELATIVE_TOL:
            return ODEStatus.ERROR, {}, {}

    # values bellow solver precision will be set to 0
    rates.update({k: 0 for k, v in rates.items() if
                  v < SolverConfigurations.ABSOLUTE_TOL and v > - SolverConfigurations.ABSOLUTE_TOL})
    conc = OrderedDict(zip(model.metabolites.keys(), C))
    return ODEStatus.OPTIMAL, rates, conc


class KineticThread(threading.Thread):
    """
    Solves the ODE inside a thread enabling to impose a timeout limit with thread.join(timeout)
    """

    def __init__(self, model, initial_concentrations=None, time_steps=None, parameters=None, final_factors=None):
        super(KineticThread, self).__init__()
        self.model = model
        self.parameters = parameters
        self.final_factors = final_factors
        self.initial_concentrations = initial_concentrations
        self.time_steps = time_steps
        self.status = None
        self.rates = None
        self.concentrations = None

    def run(self):
        try:
            self.status, self.rates, self.concentrations = kinetic_solve(self.model,
                                                                         self.initial_concentrations,
                                                                         self.time_steps,
                                                                         self.parameters,
                                                                         self.final_factors)
        except Exception as e:
            warnings.warn(f"{e}")
        return


class kineticSimulationResult(SimulationResult):

    def __init__(self, model, status, factors=None, rates=None, concentrations=None):
        super(kineticSimulationResult, self).__init__(model, None, fluxes=rates, status=status)
        self.factors = factors
        self.concentations = concentrations


class KineticSimulation:

    def __init__(self, model, parameters=None, tSteps=[0, 1e9], timeout=KineticConfigurations.SOLVER_TIMEOUT):
        if not isinstance(model, ODEModel):
            raise ValueError('model is not an instance of ODEModel.')
        self.model = model
        self.parameters = parameters
        self.tSteps = tSteps
        self.timeout = timeout

    def __getstate__(self):
        state = OrderedDict(self.__dict__.copy())
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def get_initial_concentrations(self):
        values = []
        for i, m in enumerate(self.model.metabolites):
            try:
                values.append(self.model.concentrations[m])
            except:
                values.append(None)
        return values

    def get_time_steps(self):
        return self.tSteps

    def simulate(self, factors=None):
        """
        This method preform the phenotype simulation of the kinetic model, using the solverId method and applying
        the modifications present in the instance of overrideSimulProblem.

        :param dict factors: Modification over the kinetic model.
        :returns: Returns a kineticSimulationResult with the steady-state flux distribution and concentrations.
        """

        final_factors = factors if factors is not None else {}
        # update initial concentrations when a [enz] is changed: == 0, up or down regulated
        initConcentrations = self.get_initial_concentrations()

        status = None
        sstateRates = None
        sstateConc = None

        if self.timeout:

            th = KineticThread(self.model, parameters=self.parameters,
                               final_factors=final_factors,
                               initial_concentrations=initConcentrations,
                               time_steps=self.get_time_steps())

            th.start()
            th.join(self.timeout)
            if th.is_alive():
                th.join()
            status = th.status
            sstateRates = th.rates
            sstateConc = th.concentrations

        else:
            status, sstateRates, sstateConc = kinetic_solve(self.model,
                                                            initConcentrations,
                                                            self.get_time_steps(),
                                                            self.parameters,
                                                            final_factors)

        return kineticSimulationResult(self.model, status, factors=final_factors, rates=sstateRates,
                                       concentrations=sstateConc)
