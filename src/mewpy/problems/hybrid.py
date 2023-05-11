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
Hybrid Kinetic/Constraint-Based Optimization Problems

Author: Vitor Pereira
##############################################################################   
"""
from mewpy.problems import GeckoOUProblem, GOUProblem
from mewpy.solvers import KineticConfigurations
from mewpy.simulation.kinetic import KineticSimulation
from mewpy.simulation.hybrid import Map
from re import search
from typing import Union, TYPE_CHECKING, Dict, Tuple, List

if TYPE_CHECKING:
    from cobra.core import Model
    from reframed.core.cbmodel import CBModel
    from mewpy.optimization.evaluation import EvaluationFunction
    from mewpy.simulation.simulation import Simulator
    from mewpy.model.kinetic import ODEModel

class HybridGOUProblem(GOUProblem):

    def __init__(self,
                 model: Union["Model", "CBModel"],
                 hconstraints: Dict[str, Tuple[float, float]],
                 fevaluation: List["EvaluationFunction"] = None,
                 **kwargs):
        """ Overrides GOUProblem by applying constraints resulting from
        sampling a kinetic model.

        :param model: The constraint metabolic model.
        :param dict hconstraints: The hybrid constraints definind kinetic model solution space. 
        :param list fevaluation: A list of callable EvaluationFunctions.
        Optional:

        :param OrderedDict envcond: Environmental conditions.
        :param OrderedDict constraints: Additional constraints to be applied to the model.
        :param int candidate_min_size: The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
        :param int candidate_max_size: The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
        :param list target: List of modification target genes.
        :param list non_target: List of non target genes. Not considered if a target list is provided.
        :param float scalefactor: A scaling factor to be used in the LP formulation.
        :param dic reference: Dictionary of flux values to be used in the over/under expression values computation.
        :param tuple operators: (and, or) operations. Default (MIN, MAX).
        :param list levels: Over/under expression levels (Default EAConstants.LEVELS).
        :param boolean twostep: If deletions should be applied before identifiying reference flux values.

        """
        super().__init__(model, fevaluation, **kwargs)
        self.hconstraints=hconstraints

    def solution_to_constraints(self, candidate):
        constraints = super().solution_to_constraints(candidate)
        # Apply the hybrid contraints:
        # Finds the intersection ranges of the kinetic and steady-state 
        # constraints genetic modifications space.
        # If not overlapping, the genetic modifications are not applied. 
        for r,v in  self.hconstraints.items():
            if r in constraints.keys():
                x = constraints[r]
                if isinstance(x,tuple):
                    lb,ub = x
                else:
                    lb=x
                    ub=x  
                hlb,hub = v
                l = max(lb, hlb)
                u = min(ub, hub)
                if l<=u:
                    constraints[r] = (l, u)
                else:
                    constraints[r] = (hlb, hub)
            else:
                constraints[r]=v
        return constraints    


class HybridGeckoOUProblem(GeckoOUProblem):

    def __init__(self,
                 kmodel: "ODEModel",
                 cbmodel: Union["Simulator", "Model", "CBModel"],
                 enzyme_mapping: Map,
                 fevaluation=None,
                 gDW: float = 564.0,
                 t_points: List[Union[float, int]] = [0, 1e9],
                 timeout: int = KineticConfigurations.SOLVER_TIMEOUT,
                 **kwargs):
        """ Overrides GeckoOUProblem by applying constraints resulting from kinetic simulations.
        
        :param kmodel: The kinetic model.
        :param cmodel: A GECKO model.
        :param Map enzyme_mapping: The mapping between kinetic and GECKO model 
            (instance of mewpy.simulation.hybrid.Map).
        :param list fevaluation: A list of callable EvaluationFunctions.
        
        Optional:
        :param float gDW: the organims MW in gDW. Default 564.0
        :param list t_points: Time point or span. Default [0, 1e9].
        :param int timeout: ODE solver timeout. Default KineticConfigurations.SOLVER_TIMEOUT. 
         
        Optional from GECKO Problem:

        :param OrderedDict envcond: Environmental conditions.
        :param OrderedDict constraints: Additional constraints to be applied to the model.
        :param int candidate_min_size: The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
        :param int candidate_max_size: The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
        :param list target: List of modification target genes.
        :param list non_target: List of non target genes. Not considered if a target list is provided.
        :param float scalefactor: A scaling factor to be used in the LP formulation.
        :param dic reference: Dictionary of flux values to be used in the over/under expression values computation.
        :param tuple operators: (and, or) operations. Default (MIN, MAX).
        :param list levels: Over/under expression levels (Default EAConstants.LEVELS).
        :
        """

        super().__init__(cbmodel, fevaluation, **kwargs)

        self.kmodel = kmodel
        self.gDW = gDW
        self.enzyme_mapping = enzyme_mapping
        self.t_points = t_points
        self.timeout = timeout

        # additional optional parameters

        # initial concentrations
        self.initcond = kwargs.get('initcond', None)
        # constraint the lower bound
        self.apply_lb = kwargs.get('apply_lb', True)
        # the lb tolerance
        self.lb_tolerance = kwargs.get('lb_tolerance', 0.05)

        self.ksim = KineticSimulation(model=self.kmodel, t_points=self.t_points, timeout=self.timeout)

        self.vmaxs = None

    def decode(self, candidate):
        """
        Decodes a candidate, an integer set, into a dictionary of constraints
        """
        decoded_candidate = dict()
        for idx, lv_idx in candidate:
            try:
                decoded_candidate[self.target_list[idx]] = self.levels[lv_idx]
            except IndexError:
                raise IndexError(
                    f"Index out of range: {idx} from {len(self.target_list[idx])}")
        return decoded_candidate


    def _build_target_list(self):
        """ Generates a target list, set of Vmax variables and proteins.
            It expects Vmax variables to be defined using "rmax"/"Rmax" or "vmar"/"Vmax"
            substrings, e.g., 'rmaxPGM' or vmax_PK'. For other notations, a user should
            provide a target list that includes maximum velocities variables identifiers
            and genes associated to reactions not modeled by kinetic laws.
        """
        p = list(self.kmodel.get_parameters(exclude_compartments=True))
        vmaxs = []
        for k in p:
            if search(r'(?i)[rv]max', k):
                vmaxs.append(k)
        self.vmaxs = vmaxs[:]
        vmaxs = set(vmaxs)
        enz = self.enzyme_mapping.proteins
        proteins = [x for x in self.simulator.proteins if x not in enz]

        # remove IDs defined as non modification targets
        if self.non_target:
            vmaxs = vmaxs - set(self.non_target)
            proteins = proteins - set(self.non_target)

        # remove IDs included in the partial solution
        if self._partial_solution:
            vmaxs = vmaxs - set(self._partial_solution.keys())
            proteins = proteins - set(self._partial_solution.keys())

        # the target list
        self._trg_list = list(vmaxs)+list(proteins)

    def solution_to_constraints(self, 
                                candidate: Dict[str, float]
                                ) -> Dict[str, Union[float, Tuple[float, float]]]:
        """Converts a dictionary of modifications to metabolic constraints.

        :param candidate: the genetic modifications 
        :type candidate: Dict[str,float]
        :return: a dictionary of metabolic constraints
        """
        # cb constraints
        cb_candidate = {"{}{}".format(self.prot_prefix, k): v 
                          for k, v in candidate.items() 
                          if k not in self.vmaxs}
        constraints = super().solution_to_constraints(cb_candidate)

        # kinetic modifications
        factors = {k: v for k, v in candidate.items() if k in self.vmaxs}
        result = self.ksim.simulate(factors=factors, initcon=self.initcond)
        fluxes = result.fluxes
        params = self.kmodel.merge_constants()

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
                    
                    # TODO: include dil rate
                    max_enzyme_usage = vmax_value * 3600 / (kcat * self.gDW)
                    if self.apply_lb:
                        min_enzyme_usage = max(0, abs(flux) * 3600 / (kcat * self.gDW)-self.lb_tolerance)
                    else:
                        min_enzyme_usage = 0
                    draw_p = f"{self.prot_prefix}{protein}"

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

            constraints.update(enzymatic_constraints)
            return constraints
