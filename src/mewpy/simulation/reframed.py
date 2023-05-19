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
Simulation for REFRAMED models, wraps a REFRAMED CBModel to share a common
interface to be used by MEWpy.

Author: Vítor Pereira    
##############################################################################
"""

import logging
from collections import OrderedDict
import numpy as np

from reframed.cobra.simulation import FBA, pFBA, MOMA, lMOMA, ROOM
from reframed.core.cbmodel import CBModel
from reframed.solvers import set_default_solver
from reframed.solvers import solver_instance
from reframed.solvers.solution import Solution
from reframed.solvers.solution import Status as s_status
from reframed.core.model import ReactionType

from . import SimulationMethod, SStatus, get_default_solver
from .simulation import Simulator, SimulationResult, ModelContainer
from mewpy.model.gecko import GeckoModel
from mewpy.util.constants import ModelConstants
from mewpy.util.utilities import elements, AttrDict
from tqdm import tqdm

LOGGER = logging.getLogger(__name__)

solver_map = {'gurobi': 'gurobi', 'cplex': 'cplex', 'glpk': 'optlang'}

reaction_type_map = {
    "ENZ": ReactionType.ENZYMATIC,
    "TRP": ReactionType.TRANSPORT,
    "EX": ReactionType.EXCHANGE,
    "SINK": ReactionType.SINK,
    "OTHER": ReactionType.OTHER,
}

# TODO: missing proteins and set objective implementations
class CBModelContainer(ModelContainer):
    """ A basic container for REFRAMED models.

    :param model: A metabolic model.

    """
    _r_prefix = 'R_'
    _m_prefix = 'M_'
    _g_prefix = 'G_'

    def __init__(self, model: CBModel = None):
        self.model = model

    @property
    def id(self) -> str:
        return self.model.id

    @id.setter
    def id(self,sid:str):
        self.model.id=sid

    @property
    def reactions(self):
        return list(self.model.reactions.keys())

    def get_reaction(self, r_id:str):
        if r_id not in self.reactions:
            raise ValueError(f"Reactions {r_id} does not exist")
        rxn = self.model.reactions[r_id]
        res = {'id': r_id, 'name': rxn.name, 'lb': rxn.lb, 'ub': rxn.ub, 'stoichiometry': rxn.stoichiometry}
        res['gpr'] = str(rxn.gpr) if rxn.gpr is not None else None
        res['annotations'] = rxn.metadata
        return AttrDict(res)

    @property
    def genes(self):
        return list(self.model.genes.keys())

    def get_gene(self, g_id:str):
        g = self.model.genes[g_id]
        gr = self.get_gene_reactions()
        r = gr.get(g_id,[])
        res = {'id': g_id, 'name': g.name, 'reactions': r}
        return AttrDict(res)

    @property
    def metabolites(self):
        return list(self.model.metabolites.keys())

    def get_metabolite(self, m_id:str):
        met = self.model.metabolites[m_id]
        res = {'id': m_id, 
               'name': met.name, 
               'compartment': met.compartment, 
               'formula': met.metadata.get('FORMULA', '')}
        return AttrDict(res)

    @property
    def compartments(self):
        return self.model.compartments

    def get_compartment(self, c_id:str):
        c = self.model.compartments[c_id]
        res = {'id': c_id, 'name': c.name, 'external': c.external}
        return AttrDict(res)

    def get_exchange_reactions(self):
        return self.model.get_exchange_reactions()

    def get_gene_reactions(self):
        """
        :returns: a map of genes to reactions.
        """
        if not self._gene_to_reaction:
            gr = OrderedDict()
            for rxn_id in self.reactions:
                rxn = self.model.reactions[rxn_id]
                if rxn.gpr:
                    genes = rxn.gpr.get_genes()
                    for g in genes:
                        if g in gr.keys():
                            gr[g].append(rxn_id)
                        else:
                            gr[g] = [rxn_id]
            self._gene_to_reaction = gr
        return self._gene_to_reaction

    @property
    def medium(self):

        def is_active(rxn):
            """Determine if a boundary reaction permits flux towards creating
            metabolites
            """
            reaction = self.model.reactions[rxn]
            return ((bool(reaction.get_products()) and (reaction.ub > 0)) or
                    (bool(reaction.get_substrates) and (reaction.lb < 0)))

        def get_active_bound(rxn):
            """For an active boundary reaction, return the relevant bound"""
            reaction = self.model.reactions[rxn]
            if reaction.get_substrates():
                return -reaction.lb
            elif reaction.get_products():
                return reaction.ub

        return {rxn: get_active_bound(rxn) for rxn in self.get_exchange_reactions()
                if is_active(rxn)}


class Simulation(CBModelContainer, Simulator):
    """Generic Simulation class for cobra Model.
       Defines the simulation conditions, and makes available a set of methods.

    :param model: An metabolic model instance.

    Optional:
    :param dic envcond: Dictionary of environmental conditions.
    :param dic constraints: A dictionary of reaction constraints.
    :param solver: An instance of the LP solver.
    :param dic reference: A dictionary of the wild type flux values.

    """

    # TODO: the parent init call is missing ... super() can resolve the mro of the simulation diamond inheritance
    def __init__(self, model: CBModel,
                 envcond=None,
                 constraints=None,
                 solver=None,
                 reference=None,
                 reset_solver=ModelConstants.RESET_SOLVER):

        if not isinstance(model, CBModel):
            raise ValueError(
                "Model is None or is not an instance of REFRAMED CBModel")

        self.model = model
        set_default_solver(solver_map[get_default_solver()])
        # keep track on reaction bounds changes
        self._environmental_conditions = OrderedDict() if envcond is None else envcond
        self._constraints = dict() if constraints is None else {
            k: v for k, v in constraints.items() if k not in list(self._environmental_conditions.keys())}

        self.solver = solver
        self._reference = reference
        self._gene_to_reaction = None
        self.solver = solver
        self._reset_solver = reset_solver
        self.reverse_sintax = [('_b', '_f')]
        self._m_r_lookup = None

        self.__status_mapping = {
            s_status.OPTIMAL: SStatus.OPTIMAL,
            s_status.UNBOUNDED: SStatus.UNBOUNDED,
            s_status.INFEASIBLE: SStatus.INFEASIBLE,
            s_status.INF_OR_UNB: SStatus.INF_OR_UNB,
            s_status.UNKNOWN: SStatus.UNKNOWN,
            s_status.SUBOPTIMAL: SStatus.SUBOPTIMAL
        }
        # apply the env. cond. and additional constraints to the model
        for r_id, bounds in self._environmental_conditions.items():
            self._set_model_reaction_bounds(r_id, bounds)
        for r_id, bounds in self._constraints.items():
            self._set_model_reaction_bounds(r_id, bounds)

        # if modifications on the envirenment are permited
        # during simulations
        self._allow_env_changes = False
        self.biomass_reaction = None 
        try: 
            self.biomass_reaction = model.biomass_reaction
        except:
            pass
    
    def _set_model_reaction_bounds(self, r_id, bounds):
        if isinstance(bounds, tuple):
            lb = bounds[0]
            ub = bounds[1]
        elif isinstance(bounds, (int, float)):
            lb = bounds
            ub = bounds
        else:
            raise ValueError(f"Invalid bounds definition {bounds}")
        self.model.set_flux_bounds(r_id, lb, ub)

    @property
    def environmental_conditions(self):
        return self._environmental_conditions.copy()

    @environmental_conditions.setter
    def environmental_conditions(self, a):
        raise ValueError("Can not change environmental conditions. Use set_reaction_bounds instead")

    @property
    def objective(self):
        return self.model.get_objective()

    @objective.setter
    def objective(self, objective):
        a = self.model.get_objective()
        d = {k: 0 for k in a}
        if isinstance(objective, str):
            d[objective] = 1
        elif isinstance(objective, dict):
            d.update(objective)
        else:
            raise ValueError(
                'The objective must be a reaction identifier or a dictionary of \
                reaction identifier with respective coeficients.')

        self.model.set_objective(d)

    def set_objective(self, reaction_id: str):
        self.model.set_objective({reaction_id: 1})

    def update(self):
        """Updates the model
        """
        self.model.update()

    def add_compartment(self, comp_id, name=None, external=False):
        """ Adds a compartment

            :param str comp_id: Compartment ID
            :param str name: Compartment name, default None
            :param bool external: If the compartment is external, default False.
        """
        from reframed.core.model import Compartment
        comp = Compartment(comp_id,name=name,external=external)
        self.model.add_compartment(comp)

    def add_metabolite(self, id, formula=None, name=None, compartment=None):
        """Adds a metabolite

        :param id: [description]
        :type id: [type]
        :param formula: [description], defaults to None
        :type formula: [type], optional
        :param name: [description], defaults to None
        :type name: [type], optional
        :param compartment: [description], defaults to None
        :type compartment: [type], optional
        """
        from reframed.core.model import Metabolite
        meta = Metabolite(id, name=name, compartment=compartment)
        meta.metadata['FORMULA'] = formula
        self.model.add_metabolite(meta)

    def add_gene(self,id,name):
        from reframed.core.cbmodel import Gene
        g = Gene(id,name)
        self.model.add_gene(g)

    def add_reaction(self, 
                     rxn_id,
                     name=None,
                     stoichiometry=None,
                     reversible=True,
                     lb=ModelConstants.REACTION_LOWER_BOUND,
                     ub=ModelConstants.REACTION_UPPER_BOUND, 
                     gpr= None,
                     objective=0,
                     replace=True,
                     annotations={},
                     reaction_type=None):
        """Adds a reaction to the model

        :param rxn_id: The reaction identifier
        :type rxn_id: str
        :param name: The name of the reaction, defaults to None
        :type name: str, optional
        :param stoichiometry: The stoichiometry of the reaction, defaults to None
        :type stoichiometry: dict[str,float], optional
        :param reversible: If the reaction is reversible or not, defaults to True
        :type reversible: bool, optional
        :param lb: The reaction lower bound, defaults to ModelConstants.REACTION_LOWER_BOUND
        :type lb: float, optional
        :param ub: The reaction upper bound, defaults to ModelConstants.REACTION_UPPER_BOUND
        :type ub: float, optional
        :param gpr: The gene-protein-reactio rule, defaults to None
        :type gpr: str, optional
        :param objective: the objective coeficient, defaults to 0
        :type objective: float, optional
        :param replace: If replaces the reaction with a same identifier, defaults to True
        :type replace: bool, optional
        """
        from reframed.core.cbmodel import CBReaction
        from reframed.io.sbml import parse_gpr_rule
        if gpr and isinstance(gpr,str):
            gpr = parse_gpr_rule(gpr)
        _reaction_type = reaction_type_map.get(reaction_type, None)
        reaction = CBReaction(rxn_id,
                              stoichiometry=stoichiometry,
                              name=name,
                              lb=lb,
                              ub=ub, 
                              gpr_association=gpr,
                              reversible=reversible,
                              objective=objective,
                              reaction_type=_reaction_type)

        reaction.metadata = annotations
        self.model.add_reaction(reaction, replace=replace)

    def remove_reaction(self, r_id:str):
        """Removes a reaction from the model.

        Args:
            r_id (str): The reaction identifier.
        """
        self.model.remove_reaction(r_id)
    
    def remove_reactions(self, rxn_ids):
        """_summary_

        Args:
            rxn_ids (List[str]): _description_
        """
        for r_id in rxn_ids:
            self.model.remove_reactions(r_id)
    

    def update_stoichiometry(self, rxn_id:str, stoichiometry:dict):
        rxn = self.model.reactions[rxn_id]
        rxn.stoichiometry = OrderedDict(stoichiometry)
        self.model._needs_update = True
        self.solver = None

    def get_uptake_reactions(self):
        """
        :returns: The list of uptake reactions.

        """
        drains = self.get_exchange_reactions()
        reacs = [r for r in drains if self.model.reactions[r].reversible or
                 ((self.model.reactions[r].lb is None or self.model.reactions[r].lb < 0)
                  and len(self.model.reactions[r].get_substrates()) > 0) or
                 ((self.model.reactions[r].ub is None or self.model.reactions[r].ub > 0)
                  and len(self.model.reactions[r].get_products())) > 0]
        return reacs

    def get_transport_reactions(self):
        """
        :returns: The list of transport reactions.
        """
        transport_reactions = []
        for rx in self.reactions:
            s_set = set()
            p_set = set()
            s = self.model.reactions[rx].get_substrates()
            for x in s:
                c = self.model.metabolites[x].compartment
                s_set.add(c)
            p = self.model.reactions[rx].get_products()
            for x in p:
                c = self.model.metabolites[x].compartment
                p_set.add(c)
            if len(s) == 1 and len(p) == 1 and len(p_set.intersection(s_set)) == 0:
                transport_reactions.append(rx)
        return transport_reactions

    def get_transport_genes(self):
        """Returns the list of genes that only catalyze transport reactions.
        """
        trp_rxs = self.get_transport_reactions()
        r_g = self.get_gene_reactions()
        genes = []
        for g, rxs in r_g.items():
            if set(rxs).issubset(set(trp_rxs)):
                genes.append(g)
        return genes

    def reverse_reaction(self, reaction_id:str):
        """
        Identify if a reaction is reversible and returns the
        reverse reaction if it is the case.

        :param reaction_id: A reaction identifier.
        :return: A reverse reaction identifier or None

        """

        # TODO: ... use regex instead.

        rxn = self.model.reactions[reaction_id]
        reactions = self.model.reactions
        if rxn.lb < 0:
            rxn.reversible = True
            return reaction_id
        # The model might have been converted to irreversible by REFRAMED in which case reversible reactions
        # are decoupled into forward (reaction_id+'_f') and backward (reaction_id+'_b') reactions
        # or migth be using some other identifier which must be included in self.reverse_sufix
        else:
            for a, b in self.reverse_sintax:
                n = len(reaction_id) - len(a)
                m = len(reaction_id) - len(b)
                if reaction_id[n:] == a and reactions[reaction_id[:n] + b]:
                    return reaction_id[:n] + b
                elif reaction_id[m:] == b and reactions[reaction_id[:m] + a]:
                    return reaction_id[:m] + a
                else:
                    continue
            return None

    def metabolite_reaction_lookup(self, force_recalculate=False):
        """ Return the network topology as a nested map from metabolite to reaction to coefficient.
        :return: a dictionary lookup table
        """
        if force_recalculate:
            self.model._m_r_lookup=None
        return self.model.metabolite_reaction_lookup()

    def metabolite_elements(self, metabolite_id):
        formula = self.model.metabolites[metabolite_id].metadata['FORMULA']
        return elements(formula)

    def get_reaction_bounds(self, reaction):
        """
        Returns the bounds for a given reaction.
        :param reaction: str, reaction ID
        :return: lb(s), ub(s), tuple
        """
        lb, ub = self.model.reactions[reaction].lb, self.model.reactions[reaction].ub
        #return lb if lb > -np.inf else ModelConstants.REACTION_LOWER_BOUND,\
        #    ub if ub < np.inf else ModelConstants.REACTION_UPPER_BOUND
        return lb,ub

    def set_reaction_bounds(self, reaction, lb, ub, track=True):
        """
        Sets the bounds for a given reaction.
        :param reaction: str, reaction ID
        :param float lb: lower bound 
        :param float ub: upper bound
        :param bool track: if the changes are to be logged. Default True
        """
        if track:
            if reaction in self.get_uptake_reactions():
                self._environmental_conditions[reaction] = (lb, ub)
            else:
                self._constraints[reaction] = (lb, ub)
        self.model.set_flux_bounds(reaction, lb, ub)

    def find_bounds(self):
        """
        Return the median upper and lower bound of the metabolic model.
        Bounds can vary from model to model. Cobrapy defaults to (-1000, 1000).
        """
        lower_bounds = np.asarray([rxn.lb for rxn in self.model.reactions.values()], dtype=float)
        upper_bounds = np.asarray([rxn.ub for rxn in self.model.reactions.values()], dtype=float)
        lower_bound = np.nanmedian(lower_bounds[lower_bounds != 0.0])
        upper_bound = np.nanmedian(upper_bounds[upper_bounds != 0.0])
        if np.isnan(lower_bound):
            LOGGER.warning("Could not identify a median lower bound.")
            lower_bound = ModelConstants.REACTION_LOWER_BOUND
        if np.isnan(upper_bound):
            LOGGER.warning("Could not identify a median upper bound.")
            upper_bound = ModelConstants.REACTION_UPPER_BOUND
        return lower_bound, upper_bound

    def find_unconstrained_reactions(self):
        """Return list of reactions that are not constrained at all."""
        lower_bound, upper_bound = self.find_bounds()
        return [
            rxn
            for rxn in self.model.reactions
            if rxn.lb <= lower_bound and rxn.ub >= upper_bound
        ]

    # Simulate
    def simulate(self, objective=None, method=SimulationMethod.FBA,
                 maximize=True, constraints=None, reference=None,
                 scalefactor=None, solver=None, slim=False, shadow_prices=False):
        '''
        Simulates a phenotype when applying a set constraints using the specified method.

        :param dic objective: The simulation objective. If none, the model objective is considered.
        :param method: The SimulationMethod (FBA, pFBA, lMOMA, etc ...)
        :param boolean maximize: The optimization direction
        :param dic constraints: A dictionary of constraints to be applied to the model.
        :param dic reference: A dictionary of reaction flux values.
        :param float scalefactor: A positive scaling factor for the solver. Default None.
        :param solver: An instance of the solver.
        '''

        if not objective:
            objective = self.model.get_objective()

        simul_constraints = OrderedDict()
        if constraints:
            if not self._allow_env_changes:
                simul_constraints.update({k: v for k, v in constraints.items()
                                        if k not in list(self._environmental_conditions.keys())})
            else:
                simul_constraints.update(constraints)

        if solver is not None:
            a_solver = solver
        elif self._reset_solver:
            a_solver = solver_instance(self.model)
        else:
            if self.solver is None:
                self.solver = solver_instance(self.model)
            a_solver = self.solver

        # scales the model if a scalling factor is defined.
        # ... scalling should be implemented at the solver level.
        if scalefactor:
            for _, rxn in self.model.reactions.items():
                rxn.lb = rxn.lb * scalefactor
                rxn.ub = rxn.ub * scalefactor
            if simul_constraints:
                for idx, constraint in simul_constraints.items():
                    if isinstance(constraint, (int, float)):
                        simul_constraints[idx] = constraint * scalefactor
                    elif isinstance(constraint, tuple):
                        simul_constraints[idx] = tuple(
                            x * scalefactor for x in constraint)
                    else:
                        raise ValueError("Could not scale the model")

        # TODO: simplifly ...using python >=3.10 cases
        if method in [SimulationMethod.lMOMA, SimulationMethod.MOMA, SimulationMethod.ROOM] and reference is None:
            reference = self.reference

        if method == SimulationMethod.FBA:
            get_values = not slim
            solution = FBA(self.model, objective=objective, minimize=not maximize,
                           constraints=simul_constraints, solver=a_solver, get_values=get_values)
        elif method == SimulationMethod.pFBA:
            solution = pFBA(self.model, objective=objective, minimize=not maximize,
                            constraints=simul_constraints, solver=a_solver, obj_frac=0.999)
        elif method == SimulationMethod.lMOMA:
            solution = lMOMA(self.model, constraints=simul_constraints,
                             reference=reference, solver=a_solver)
        elif method == SimulationMethod.MOMA:
            solution = MOMA(self.model, constraints=simul_constraints,
                            reference=reference, solver=a_solver)
        elif method == SimulationMethod.ROOM:
            solution = ROOM(self.model, constraints=simul_constraints,
                            reference=reference, solver=a_solver)
        # Special case in which only the simulation context is required without any simulatin result
        elif method == SimulationMethod.NONE:
            solution = Solution(status=s_status.UNKNOWN,
                                message=None, fobj=None, values=None)
        else:
            raise Exception(
                "Unknown method to perform the simulation.")

        # undoes the model scaling
        if scalefactor:
            for _, rxn in self.model.reactions.items():
                rxn.lb = rxn.lb / scalefactor
                rxn.ub = rxn.ub / scalefactor
            if solution.status in (s_status.OPTIMAL, s_status.SUBOPTIMAL):
                solution.fobj = solution.fobj / scalefactor
                for x, y in solution.values.items():
                    solution.values[x] = y / scalefactor

        if slim:
            return solution.fobj

        else:
            status = self.__status_mapping[solution.status]
            result = SimulationResult(self.model, solution.fobj, fluxes=solution.values, status=status,
                                      envcond=self.environmental_conditions, model_constraints=self._constraints.copy(),
                                      simul_constraints=constraints, maximize=maximize, method=method)
            return result

    def FVA(self, reactions=None, obj_frac=0.9, constraints=None, loopless=False, internal=None, solver=None,
            format='dict'):
        """ Flux Variability Analysis (FVA).

        :param model: An instance of a constraint-based model.
        :param float obj_frac: The minimum fraction of the maximum growth rate (default 0.9).\
            Requires that the objective value is at least the fraction times maximum objective value.\
            A value of 0.85 for instance means that the objective has to be at least at 85% percent of its maximum.
        :param list reactions: List of reactions to analyze (default: all).
        :param dic constraints: Additional constraints (optional).
        :param boolean loopless: Run looplessFBA internally (very slow) (default: false).
        :param list internal: List of internal reactions for looplessFBA (optional).
        :param solver: A pre-instantiated solver instance (optional)
        :param format: The return format: 'dict', returns a dictionary,'df' returns a data frame.
        :returns: A dictionary of flux variation ranges.

        """
        simul_constraints = {}
        if constraints:
            simul_constraints.update({k: v for k, v in constraints.items()
                                      if k not in list(self._environmental_conditions.keys())})
        if reactions is None:
            _reactions = self.reactions
        elif isinstance(reactions, str):
            _reactions = [reactions]
        elif isinstance(reactions, list):
            _reactions = reactions
        else:
            raise ValueError('Invalid reactions.')


        from reframed.cobra.variability import FVA
        res = FVA(self.model, obj_frac=obj_frac, reactions=_reactions,
                  constraints=simul_constraints, loopless=loopless, internal=internal, solver=solver)
        if format == 'df':
            import pandas as pd
            e = res.items()
            f = [[a, b, c] for a, [b, c] in e]
            df = pd.DataFrame(f, columns=['Reaction ID', 'Minimum', 'Maximum'])
            return df
        return res

    def create_empty_model(self,model_id:str):
        return Simulation(CBModel(model_id))


class GeckoSimulation(Simulation):

    def __init__(self, model: GeckoModel, envcond=None, constraints=None, solver=None, reference=None,
                 reset_solver=ModelConstants.RESET_SOLVER, protein_prefix=None):
        super(GeckoSimulation, self).__init__(
            model, envcond, constraints, solver, reference, reset_solver)
        self.protein_prefix = protein_prefix if protein_prefix else 'R_draw_prot_'
        self._essential_proteins = None

    @property
    def proteins(self):
        return self.model.proteins

    @property
    def protein_pool_exchange(self):
        return self.model.protein_pool_exchange

    @property
    def protein_rev_reactions(self):
        return self.model.protein_rev_reactions

    def essential_proteins(self, min_growth=0.01):
        if self._essential_proteins is not None:
            return self._essential_proteins
        wt_solution = self.simulate()
        wt_growth = wt_solution.objective_value
        self._essential_proteins = []
        proteins = self.model.proteins
        for p in tqdm(proteins):
            rxn = "{}{}".format(self.protein_prefix, p)
            res = self.simulate(constraints={rxn: 0})
            if res:
                if (res.status == SStatus.OPTIMAL and res.objective_value < wt_growth * min_growth) \
                        or res.status == SStatus.INFEASIBLE:
                    self._essential_proteins.append(rxn)
        return self._essential_proteins

    def protein_reactions(self, protein:str):
        """
        Returns the list of reactions associated with a protein.

        :params (str) protein: The protein identifier.
        :returns: A list of reaction identifiers
        """
        reactions = []
        for r_id, rxn in self.model.reactions.items():
            lsub = rxn.get_substrates()
            for m in lsub:
                if protein in m:
                    reactions.append(r_id)
        return reactions

    def reverse_reaction(self, reaction_id:str):
        """
        Identify if a reaction is reversible and returns the
        reverse reaction if it is the case.

        :param (str) reaction_id: The reaction identifier.
        :returns: A reaction identifier or None.
        """
        f, d = zip(*self.model.protein_rev_reactions.values())
        if reaction_id in f:
            return d[f.index(reaction_id)]
        elif reaction_id in d:
            return f[d.index(reaction_id)]
        else:
            return None

    def get_Kcats(self, protein:str):
        """ 
        Returns a dictionary of reactions and respective Kcat for a 
        specific protein/enzyme·

        :params (str) protein: the protein identifier.

        :returns: A dictionary of reactions and respective Kcat values.
        """
        import re
        re_expr = re.compile(f"{protein}_")
        values = [x for x in self.metabolites if re_expr.search(x) is not None]
        if len(values)==1:
            m_r = self.metabolite_reaction_lookup()
            r_d = m_r[values[0]]
            return {k: -1/v for k, v in r_d.items() if self.protein_prefix not in k}
        elif len(values)>1: 
            raise ValueError(f"More than one protein match {values}")
        else:
            raise ValueError(f"Protein {protein} not founded.")

    def set_Kcat(self, protein:str, reaction:str, kcat:float):
        """Alters an enzyme kcat for a given reaction.

        :param protein: The protein identifier
        :type protein: str
        :param reaction: The reaction identifier
        :type reaction: str
        :param kcat: The kcat value, a real positive value.
        :type kcat: float
        """
        if kcat <= 0:
            raise ValueError('kcat value needs to be positive.')
        rx = self.model.reactions[reaction]
        st = rx.stoichiometry
        mets =[x for x in list(st.keys()) if protein in x]
        if len(mets)==1:
            met=mets[0]
            st[met] = -1/kcat
            rx.stoichiometry = st
            self.model._needs_update=True
            self.solver=None
        else:
            LOGGER.warn(f'Could not identify {protein} ' 
                        f'protein specie in reaction {reaction}')
 