
from mewpy.utils.ode import KineticConfigurations, SolverConfigurations, ODEStatus, ODESpySolver
from mewpy.simulation.simulation import SimulationResult
from mewpy.model.kinetic import ODEModel
from collections import OrderedDict
import threading
import warnings


def kinetic_solve(model, finalParameters, finalFactors, initialConc, timePoints):
    """
    Private function: auxiliary function required to avoid the pickling the solver.solve function

    """
    finalRates = OrderedDict()
    f = model.get_ode(r_dict=finalRates, params=finalParameters, factors=finalFactors)

    def func(x, t):
        return f(t, x)

    solver = ODESpySolver(KineticConfigurations.SOLVER_METHOD).get_solver(func)
    solver.set_initial_condition(list(initialConc.values()))
    try:
        X, _ = solver.solve(timePoints)
        # if solver returns a solution where any concentration is negative
        for c in X[1]:
            if c < -1*SolverConfigurations.RELATIVE_TOL:
                return ODEStatus.ERROR, {}, {}

    except Exception as e:
        warnings.warn(f"ODE Solver error: {e}.")
        return ODEStatus.ERROR, {}, {}

    # values bellow solver precision will be set to 0
    finalRates.update({k: 0 for k, v in finalRates.items() if
                       v < SolverConfigurations.ABSOLUTE_TOL and v > - SolverConfigurations.ABSOLUTE_TOL})
    conc = OrderedDict(zip(model.metabolites.keys(), X[1]))
    return ODEStatus.OPTIMAL, finalRates, conc


class KineticThread(threading.Thread):
    """
    Solves the ODE inside a thread enabling to impose a timeout limit with thread.join(timeout)

    """
    def __init__(self, model, parameters=None, final_factors=None, initial_concentrations=None, time_steps=None):
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
            self.status, self.rates, self.concentrations = kinetic_solve(self.model, self.parameters,
                                                                         self.final_factors,
                                                                         self.initial_concentrations,
                                                                         self.time_steps)
        except Exception as e:
            warnings.warn(f"{e}")
        return


class kineticSimulationResult(SimulationResult):
    def __init__(self, model, status, factors=None, rates=None, concentrations=None):
        super(kineticSimulationResult, self).__init__(model, None, fluxes=rates, status=status)
        self.factors = factors
        self.concentations = concentrations


class KineticSimulation:

    def __init__(self, model, parameters=None, tSteps=[0, 1e9], timeout=KineticConfigurations.SOLVER_TIMEOUT,
                 solver=KineticConfigurations.SOLVER, method=KineticConfigurations.SOLVER_METHOD):
        if not isinstance(model, ODEModel):
            raise ValueError('model is not an instance of ODEModel.')
        self.model = model
        self.parameters = parameters
        self.tSteps = tSteps
        self.timeout = timeout
        self.solver = solver
        self.method = method

    def __getstate__(self):
        state = OrderedDict(self.__dict__.copy())
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def get_initial_concentrations(self):
        return self.model.concentrations

    def get_time_steps(self):
        return self.tSteps

    def simulate(self, factors=None):
        """
        This method preform the phenotype simulation of the kinetic model, using the solverId method and applying
        the modifications present in the instance of overrideSimulProblem.

        :param dict factores: Modification over the kinetic model.
        :returns: Returns a kineticSimulationResult with the steady-state flux distribution and concentrations.

        """

        final_factors = factors if factors is not None else {}
        # update initial concentrations when a [enz] is changed: == 0, up or down regulated
        initConcentrations = self.get_initial_concentrations().copy()

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
            status = th.status
            sstateRates = th.rates
            sstateConc = th.concentrations

        else:
            status, sstateRates, sstateConc = kinetic_solve(self.model, self.parameters,
                                                            final_factors,
                                                            initConcentrations,
                                                            self.get_time_steps())

        return kineticSimulationResult(self.model, status, factors=final_factors, rates=sstateRates,
                                       concentrations=sstateConc)
