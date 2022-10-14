import warnings
from mewpy.model.kinetic import ODEModel
from mewpy.solvers import KineticConfigurations
from mewpy.simulation import get_simulator
from mewpy.simulation.simulation import Simulator
from mewpy.simulation.kinetic import KineticSimulation
from collections import OrderedDict
from warnings import warn
import pandas as pd
import numpy as np
from numpy.random import normal
from re import search
from tqdm import tqdm


warnings.filterwarnings('ignore', 'Timeout')


def sample(vmaxs, sigma=0.1):
    k = vmaxs.keys()
    f = np.exp(normal(0, sigma, len(vmaxs)))
    v = np.array(list(vmaxs.values()))
    r = list(v*f)
    return dict(zip(k, r))


class HybridSimulation:

    def __init__(self, kmodel, cbmodel, gDW=564.0, envcond=dict(),
                 mapping=dict(), t_points=[0, 1e9],
                 timeout=KineticConfigurations.SOLVER_TIMEOUT):

        if not isinstance(kmodel, ODEModel):
            raise ValueError('model is not an instance of ODEModel.')
        
        if not isinstance(cbmodel, Simulator):
            self.sim = get_simulator(cbmodel, envcond=envcond)
        else:
            self.sim = cbmodel

        self.kmodel = kmodel
        self.mapping = mapping
        self.t_points = t_points
        self.timeout = timeout
        self.models_verification()
        self.gDW = gDW

    def __getstate__(self):
        state = OrderedDict(self.__dict__.copy())
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def get_kinetic_model(self):
        return self.kmodel

    def get_mapping(self):
        return self.mapping

    def models_verification(self):
        """
        Function that verifies if it's possible to perform the Hibrid Simulation.
        """
        kmodel = self.get_kinetic_model()
        cbmodel = self.sim.model
        mapping = self.get_mapping()

        for k,v in mapping.items():
            if k not in kmodel.ratelaws.keys():
                raise ValueError(f"could not find reaction {k} in the kinetic model ")
            if v[0] not in cbmodel.reactions:
                raise ValueError(f"could not find reaction {v[0]} in the steady-state model ")
        
        return True
        
    def mapping_conversion(self, fluxes):
        """
        Function that converts the kinetic fluxes into constraint-based fluxes.

        :param fluxes: kinetic simulation fluxes
        :type fluxes: dict
        :return: kinetic fluxes compatible with the constraint-based model
        """
        mapping = self.get_mapping()
        flxs = dict()
        for k, value in fluxes.items():
            if k in mapping.keys():
                v = mapping[k]
                flxs[v[0]] = v[1]*value * 3600/self.gDW
        if len(flxs) != 0:
            return flxs
        else:
            raise warn('Mapping not done properly, please redo mapping')

    def mapping_bounds(self, lbs, ubs):
        """
        Function that converts the kinetic bounds into constraint-based flux bounds.

        :param lbs: kinetic lower bounds
        :type fluxes: dict
        :param ubs: kinetic upper bounds
        :type fluxes: dict
        :return: constraints
        """
        mapping = self.get_mapping()
        flxs = dict()
        for k, value in lbs.items():
            if k in mapping.keys():
                v = mapping[k]
                a = v[1]*value * 3600/self.gDW
                b = v[1]*ubs[k] * 3600/self.gDW
                flxs[v[0]] = (a, b) if a < b else (b, a)
        if len(flxs) != 0:
            return flxs
        else:
            raise warn('Mapping not done properly, please redo mapping')

    def nsamples(self, vmaxs, n=1, sigma=0.1):
        """
        Generates n fluxes samples varying vmax values on a log-norm distribtution
        with mean 0 and std sigma.

        """
        kmodel = self.get_kinetic_model()
        ksample = []
        ksim = KineticSimulation(model=kmodel, t_points=self.t_points, timeout=self.timeout)
        for _ in tqdm(range(n)):
            v = sample(vmaxs, sigma=sigma)
            try:
                res = ksim.simulate(parameters=v)
                if res.fluxes:
                    ksample.append(res.fluxes)
            except Exception as e:
                warn.warning(str(e))
        df = pd.DataFrame(ksample)
        # drop any NaN if exist
        df.dropna()
        return df
    
    
    def simulate(self, objective=None, 
                       initcond=None, 
                       parameters=None, 
                       constraints=None, 
                       method='pFBA'):
        """
        This method performs a phenotype simulation hibridizing a kinetic and a constraint-based model.

        :param objective: steady-state optimization objective.
        :type objective: dict, optional
        :param initcond:List of numbers in which the kinetic simulation fluxes will be scaled.
        :param initcond: list, optional
        :param parameters: Kinetic simulation parameters.
        :type parameters: dict, optional
        :param constraints: Constraint-based model simulation constraints.
        :type constraints: dict, optional
        :param method: the phenotype simulation method
        :type method: str. Default 'pFBA'
        :returns: Returns the solution of the hibridization.
        """
        mapp = self.models_verification()
        kmodel = self.get_kinetic_model()

        ksim = KineticSimulation(model=kmodel, t_points=self.t_points, timeout=self.timeout)
        result = ksim.simulate(parameters=parameters, initcon=initcond)
        fluxes = result.fluxes
        
        if constraints is None:
            constraints = dict()

        if not mapp:
            fluxes = self.mapping_conversion(fluxes)
        constraints.update(fluxes)
        if objective:
            solution = self.sim.simulate(objective= objective, method=method, constraints=constraints)
        else:
            solution = self.sim.simulate(method=method, constraints=constraints)
        return solution


    def simulate_distribution(self, df, q1=0.1, q2=0.9, objective=None, method='pFBA', constraints=None):
        """Runs a pFBA on the steady-state model with fluxes constrained to ranges between the q1-th and q2-th percentile
        of fluxes distributions sampled from the kinetic model.
        The kinetic flux distributions are provided as panda dataframes.     
        """
        const = dict()
        lbs = df.quantile(q1, axis=0).to_dict()
        ubs = df.quantile(q2, axis=0).to_dict()
        if constraints:
            const.update(constraints)
        k_const = self.mapping_bounds(lbs, ubs)
        const.update(k_const)
        if objective:
            solution = self.sim.simulate(method=method, constraints=const, objective=objective)
        else:
            solution = self.sim.simulate(method=method, constraints=const)
        solution.kinetic_constraints = k_const
        return solution


