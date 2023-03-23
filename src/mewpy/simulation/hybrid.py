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
Hybrid kinetic/constraint-based simulation module

Author: Vitor Pereira
Contributors: Mariana Pereira
##############################################################################
"""
from mewpy.model.kinetic import ODEModel
from mewpy.solvers import KineticConfigurations
from mewpy.simulation import get_simulator, Simulator
from mewpy.simulation.kinetic import KineticSimulation
from mewpy.solvers import solver_instance
from mewpy.util.utilities import AttrDict
from math import inf
from collections import OrderedDict
import warnings
from warnings import warn
import pandas as pd
import numpy as np
from numpy.random import normal
import itertools
from tqdm import tqdm

from typing import Tuple, Dict, Union, List, TYPE_CHECKING

if TYPE_CHECKING:
    from cobra import Model
    from reframed import CBModel


warnings.filterwarnings('ignore', 'Timeout')


def _partial_lMOMA(model, reactions: dict, biomass: str, constraints=None):
    """
    Run a (linear version of) Minimization Of Metabolic Adjustment (lMOMA) 
    simulation using fluxes from the Kinetic Simulation:

    :param model: a COBRAPY or REFRAMED model, or an instance of Simulator
    :param reactions: dictionary of reactions whose sum of fluxes is to be minimized
    :type reactions: dict
    :param biomass: name of the biomass reaction
    :type biomass: str
    :param constraints: constraints to be imposed, defaults to None
    :type constraints: dict, optional
    """

    simul = get_simulator(model)

    solver = solver_instance(simul)

    _constraints = simul._environmental_conditions.copy()

    if constraints:
        _constraints.update(constraints)

    _reactions = {k: v for k, v in reactions.items()}

    bio_ref = simul.simulate({biomass: 1}, constraints=constraints, slim=True)

    for r_id in _reactions.keys():
        d_pos, d_neg = r_id + '_d+', r_id + '_d-'
        solver.add_variable(d_pos, 0, inf, update=False)
        solver.add_variable(d_neg, 0, inf, update=False)
        solver.update()

    bio_plus = biomass + '_d+'
    bio_minus = biomass + '_d-'
    solver.add_variable(bio_plus, 0, inf, update=False)
    solver.add_variable(bio_minus, 0, inf, update=False)
    solver.update()

    for r_id in _reactions.keys():
        d_pos, d_neg = r_id + '_d+', r_id + '_d-'
        solver.add_constraint('c' + d_pos, {r_id: -1, d_pos: 1}, '>', -_reactions[r_id], update=False)
        solver.add_constraint('c' + d_neg, {r_id: 1, d_neg: 1}, '>', _reactions[r_id], update=False)
        solver.update()

    solver.add_constraint('c' + bio_plus, {biomass: -1, bio_plus: 1}, '>', -bio_ref, update=False)
    solver.add_constraint('c' + bio_minus, {biomass: 1, bio_minus: 1}, '>', bio_ref, update=False)
    solver.update()

    objective = dict()
    for r_id in _reactions.keys():
        d_pos, d_neg = r_id + '_d+', r_id + '_d-'
        objective[d_pos] = 1
        objective[d_neg] = 1

    objective[bio_plus] = 1
    objective[bio_minus] = 1

    solution = solver.solve(objective, minimize=True, constraints=constraints)

    return solution


def sample(vmaxs: Dict[str, float], sigma: float = 0.1):
    k = vmaxs.keys()
    f = np.exp(normal(0, sigma, len(vmaxs)))
    v = np.array(list(vmaxs.values()))
    r = list(v*f)
    return dict(zip(k, r))


class HybridSimulation:

    def __init__(self,
                 kmodel: ODEModel,
                 cbmodel: Union[Simulator, "Model", "CBModel"],
                 gDW: float = 564.0,
                 D: float = 0.278e-4,
                 envcond: Dict[str, Union[float, Tuple[float, float]]] = dict(),
                 mapping: Dict[str, Tuple[str, int]] = dict(),
                 t_points: List[Union[float, int]] = [0, 1e9],
                 timeout: int = KineticConfigurations.SOLVER_TIMEOUT):
        """_summary_

        :param kmodel: The kinetic model
        :type kmodel: ODEModel
        :param cbmodel: The constraint-based model
        :type cbmodel: A COBRApy or REFRAMED model
        :param gDW: The volume of the organims, defaults to 564.0 to E. coli
        :type gDW: float, optional
        :param envcond: The medium, defaults to dict()
        :type envcond: dict, optional
        :param mapping: a mapping from kinetic to CB reactions, defaults to dict()
        :type mapping: dict, optional
        :param t_points: Time point or span for kinetic integration, defaults to [0, 1e9]
        :type t_points: list, optional
        :param timeout: The ODE solver timeout, defaults to KineticConfigurations.SOLVER_TIMEOUT. If 0 no timeout is set.
        :type timeout: int, optional
        :raises ValueError: The models are not of the required types.
        """

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
        self.D = D

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

        for k, v in mapping.items():
            if k not in kmodel.ratelaws.keys():
                raise ValueError(f"could not find reaction {k} in the kinetic model ")
            if v[0] not in cbmodel.reactions:
                raise ValueError(f"could not find reaction {v[0]} in the steady-state model ")

        return True

    def unit_conv(self, value):
        return value*self.D*3600/self.gDW

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
                flxs[v[0]] = self.unit_conv(v[1]*value)
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
                a = self.unit_conv(v[1]*value)
                b = self.unit_conv(v[1]*ubs[k])
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
                 amplitude=None,
                 method='pFBA',
                 **kwargs):
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
        :param amplitude: the amplitude range centered in the flux value. Default None, in which case
            partial lMOMA is applied
        :type amplitude: float 
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

        if mapp:
            _fluxes = self.mapping_conversion(fluxes)

            # add a tolerance (v-A/2,v+A/2) to the flux rate
            if amplitude:
                c = dict()
                for k, v in _fluxes.items():
                    a, _ = self.sim.get_reaction_bounds(k)
                    lb = v-amplitude/2
                    ub = v+amplitude/2
                    if a < 0:
                        c[k] = (lb, ub)
                    else:
                        c[k] = (max(0, lb), ub)

            # uses a partial lMOMA to find the closest solution
            # in the feasable space.
            else:
                # assumes growth as model objective
                biomass = [*self.sim.objective][0]
                s = _partial_lMOMA(self.sim, _fluxes, biomass)
                c = {k: s.values[k] for k in _fluxes.keys()}
            constraints.update(c)
        else:
            raise ValueError('Could not mapp reactions.')

        if objective:
            solution = self.sim.simulate(objective=objective, method=method, constraints=constraints, **kwargs)
        else:
            solution = self.sim.simulate(method=method, constraints=constraints, **kwargs)
        return solution

    def simulate_distribution(self,
                              df,
                              q1=0.1,
                              q2=0.9,
                              objective=None,
                              method='pFBA',
                              constraints=None,
                              **kwargs):
        """
        Runs a pFBA on the steady-state model with fluxes constrained to ranges 
        between the q1-th and q2-th percentile of fluxes distributions sampled 
        from the kinetic model.
        The kinetic flux distributions are provided as panda dataframes.     

        :param df: _description_
        :type df: _type_
        :param q1: _description_, defaults to 0.1
        :type q1: float, optional
        :param q2: _description_, defaults to 0.9
        :type q2: float, optional
        :param objective: _description_, defaults to None
        :type objective: _type_, optional
        :param method: _description_, defaults to 'pFBA'
        :type method: str, optional
        :param constraints: _description_, defaults to None
        :type constraints: _type_, optional
        :return: _description_
        :rtype: _type_
        """
        const = dict()
        lbs = df.quantile(q1, axis=0).to_dict()
        ubs = df.quantile(q2, axis=0).to_dict()
        if constraints:
            const.update(constraints)
        k_const = self.mapping_bounds(lbs, ubs)
        const.update(k_const)
        if objective:
            solution = self.sim.simulate(method=method, constraints=const, objective=objective, **kwargs)
        else:
            solution = self.sim.simulate(method=method, constraints=const, **kwargs)
        solution.kinetic_constraints = k_const
        return solution


class Mapper:
    def __init__(self, vmax, forward, backward, sense=1) -> None:
        """
        A reaction mapper.

        :param (str) vmax: the vmax parameter identifier.
        :param (list) forward: a list of tupples (proteins,kcat) for forward reactions.
        :param (list) backward: a list of tupples (proteins,kcat) for backward reactions.
        :param (int) sense: if the kinetic and GECKO reactions have the same sense (1) or reverse (-1).
        """
        self.vmax_id = vmax
        self.sense = sense
        self.forward = AttrDict(forward)
        self.backward = AttrDict(backward)

    @property
    def proteins(self):
        proteins = []
        proteins.extend(list(self.forward.keys()))
        proteins.extend(list(self.backward.keys()))
        return list(set(proteins))


class Map(AttrDict):

    def intersection(self, rid1, rid2):
        """
        Identifies common proteins
        """
        p1 = self.get(rid1).proteins
        p2 = self.get(rid2).proteins
        return list(set(p1).intersection(set(p2)))

    def intersections(self):
        """ 
        Identifies kinetic reactions that use
        a same enzyme.
        """
        rxns = list(self.keys())
        comb = list(itertools.combinations(rxns, 2))
        combinations = []
        for a, b in comb:
            if self.intersection(a, b):
                combinations.append((a, b))
        return combinations

    @property
    def proteins(self):
        prots = []
        for v in self.values():
            prots.extend(v.proteins)
        return prots


def read_map(jsonfile: str):
    """
    Reads kinetic to GECKO mapping json files.

    :param (str) jsonfile: the json file name.
    :returns: an instance of Map

    The json file is expected to have the structure:
        { kinetic_reaction1:
            [ vmax,
              [(protein,kcat)], # forward
              [(protein,kcat)], # forward
              sense # 1 or -1
            ]
        }
    """
    import json
    with open(jsonfile) as f:
        mapp = json.load(f)
    m = dict()
    for k, v in mapp.items():
        m[k] = Mapper(v[0], v[1], v[2], v[3])
    return Map(m)


def hasNaN(values):
    import math
    for x in values:
        if math.isnan(x):
            return True
    return False


class HybridGeckoSimulation:

    def __init__(self,
                 kmodel: ODEModel,
                 cbmodel: Union[Simulator, "Model", "CBModel"],
                 gDW: float = 564.0,
                 D: float = 0.278e-4,
                 envcond: Dict[str, Union[float, Tuple[float, float]]] = dict(),
                 enzyme_mapping: Map = None,
                 protein_prefix: str = 'R_draw_prot_',
                 t_points: List[Union[float, int]] = [0, 1e9],
                 timeout: int = KineticConfigurations.SOLVER_TIMEOUT):
        """
        Hybrid Gecko Simulation.

        :param (ODEModel) kmodel: the kinetic model.
        :param (GeckoModel) cbmodel: the GECKO model.
        :param (float) gDW: the cell volume. Default E. coli from [1].
        :param (dict) envcond: the medium definition.
        :param (Map) enzyme_mapping: An instance of Map.
        :param (str) protein_prefix: the draw protein pseudo reaction prefix, 
               e.g. for protein XXXXX 'R_draw_prot_XXXXX'.
        :param (array like) t_points: the integrative time points or span. Default [0, 1e9].
        :param (int) timeout: The integration timeout. If timeout=0, no timeout is applied.     


        [1] Chassagnole et. al, Dynamic Modeling of Central Carbon Metabolism of Escherichia 
        coli,(2002). DOI: 10.1002/bit.10288  
        """

        if not isinstance(kmodel, ODEModel):
            raise ValueError('model is not an instance of ODEModel.')

        if not isinstance(cbmodel, Simulator):
            self.sim = get_simulator(cbmodel, envcond=envcond)
        else:
            self.sim = cbmodel

        self.kmodel = kmodel
        self.t_points = t_points
        self.gDW = gDW
        self.D = D
        self.enzyme_mapping = enzyme_mapping
        self.protein_prefix = protein_prefix
        self.timeout = timeout

    def unit_conv(self, value):
        return value*self.D*3600/self.gDW

    def simulate(self, objective=None,
                 initcond=None,
                 parameters=None,
                 constraints=None,
                 method='pFBA',
                 apply_lb=True,
                 lb_tolerance=0.05,
                 **kwargs):
        """
        Runs a hybrid simulation on GECKO models by defining enzymatic 
        constraints that limit enzyme usage in function of vmax, fluxes and kcat values.

        :param objective: the optimization objective.
        :type obective: dict, optional
        :param parameters: Kinetic simulation parameters.
        :type parameters: dict, optional
        :param constraints: Constraint-based model simulation constraints.
        :type constraints: dict, optional
        :param method: the phenotype simulation method
        :type method: str. Default 'pFBA'
        :param maximize: The optimization direction (True: maximize, False:minimize).
        :type maximize: bool. Default True.
        :param apply_lb: If the lb of pseudo draw reactions are to be constrained using 
            the kinetic flux rate values.
        :type apply_lb: bool. Default True.
        :param lb_tolerance: A tolerance for the lb, ie, the lb is set to the value obtained 
            from the kinetic flux rate less the tolerance (or 0 if negative).
        :type lb_tolerance: float. Default 0.05.         
        :returns: Returns the solution of the hibridization.
        """

        # Solve the kinetic model
        ksim = KineticSimulation(model=self.kmodel, t_points=self.t_points, timeout=self.timeout)
        result = ksim.simulate(parameters=parameters, initcon=initcond)
        fluxes = result.fluxes

        params = self.kmodel.merge_constants()
        if parameters:
            params.update(parameters)

        enzymatic_constraints = dict()
        for krxn, mapper in self.enzyme_mapping.items():
            # identify the sense of the reaction
            if fluxes[krxn] > 0:
                sense = mapper.sense
            else:
                sense = -1*mapper.sense

            # A same enzyme may have different kcats
            # for each sense
            if sense > 0:
                proteins = mapper.forward
            else:
                proteins = mapper.backward

            vmax_value = params.get(mapper.vmax_id)
            flux = fluxes[krxn]

            if vmax_value:
                for protein, kcat in proteins.items():
                    # Units:
                    #    vmax:  mM/s
                    #    kcat:  1/h
                    #    gDW:   gDW/L
                    max_enzyme_usage = vmax_value * 3600 * self.D / (kcat * self.gDW)
                    if apply_lb:
                        min_enzyme_usage = max(0, abs(flux) * 3600 * self.D / (kcat * self.gDW)-lb_tolerance)
                    else:
                        min_enzyme_usage = 0
                    draw_p = f"{self.protein_prefix}{protein}"

                    # For promiscuous enzymes, the ub of enzyme usage is
                    # the sum of usages for each reaction, and the lb is
                    # the minimum usage of all reactions.
                    if draw_p in enzymatic_constraints:
                        lb, ub = enzymatic_constraints[draw_p]
                        enzymatic_constraints[draw_p] = (min(min_enzyme_usage, lb),
                                                         ub+max_enzyme_usage)

                    else:
                        enzymatic_constraints[draw_p] = (min_enzyme_usage,
                                                         max_enzyme_usage)
        if constraints is None:
            constraints = dict()
        constraints.update(enzymatic_constraints)
        if objective:
            solution = self.sim.simulate(objective=objective,
                                         method=method,
                                         constraints=constraints,
                                         **kwargs)
        else:
            solution = self.sim.simulate(method=method,
                                         constraints=constraints,
                                         **kwargs)

        return solution
