"""
Simulation for REFRAMED models
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

from . import SimulationMethod, SStatus, get_default_solver
from .simulation import Simulator, SimulationResult, ModelContainer
from ..model.gecko import GeckoModel
from ..util.constants import ModelConstants
from ..util.parsing import evaluate_expression_tree

LOGGER = logging.getLogger(__name__)

solver_map = {'gurobi': 'gurobi', 'cplex': 'cplex', 'glpk': 'optlang'}


class CBModelContainer(ModelContainer):
    """ A basic container for REFRAMED models.

    :param model: A metabolic model.

    """

    def __init__(self, model: CBModel):
        if not isinstance(model, CBModel):
            raise ValueError(
                "The model is not an instance of ReFramed CBModel")
        self.model = model

    @property
    def reactions(self):
        return list(self.model.reactions.keys())

    @property
    def genes(self):
        return list(self.model.genes.keys())

    @property
    def metabolites(self):
        return list(self.model.metabolites.keys())

    @property
    def compartments(self):
        return self.model.compartments

    def get_gpr(self, reaction_id):
        """Returns the gpr rule (str) for a given reaction ID.

        :param str reaction_id: The reaction identifier.
        :returns: A string representation of the GPR rule.

        """
        if reaction_id not in self.reactions:
            raise ValueError(f"Reactions {reaction_id} does not exist")
        reaction = self.model.reactions[reaction_id]
        if reaction.gpr:
            return str(reaction.gpr)
        else:
            return None

    def get_drains(self):
        return self.model.get_exchange_reactions()

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

        return {rxn: get_active_bound(rxn) for rxn in self.get_drains()
                if is_active(rxn)}


class Simulation(CBModelContainer, Simulator):
    """Generic Simulation class for cobra Model.
       Defines the simulation conditions, and makes available a set of methods.

    :param model: An metabolic model instance.

    Optional:

    :param objective: The model objective.
    :param dic envcond: Dictionary of environmental conditions.
    :param dic constraints: A dictionary of reaction constraints.
    :param solver: An instance of the LP solver.
    :param dic reference: A dictionary of the wild type flux values.

    """

    def __init__(self, model: CBModel, objective=None, envcond=None, constraints=None, solver=None, reference=None,
                 reset_solver=ModelConstants.RESET_SOLVER):

        if not isinstance(model, CBModel):
            raise ValueError(
                "Model is None or is not an instance of REFRAMED CBModel")
        self.model = model
        set_default_solver(solver_map[get_default_solver()])
        self.environmental_conditions = OrderedDict() if envcond is None else envcond
        self.constraints = OrderedDict() if constraints is None else constraints
        self.solver = solver
        self._essential_reactions = None
        self._essential_genes = None
        self._reference = reference
        self._gene_to_reaction = None
        self.solver = solver
        self._reset_solver = reset_solver
        self.reverse_sintax = [('_b', '_f')]
        self._index_metabolites_reactions = None
        self._m_r_lookup = None

        self.__status_mapping = {
            s_status.OPTIMAL: SStatus.OPTIMAL,
            s_status.UNBOUNDED: SStatus.UNBOUNDED,
            s_status.INFEASIBLE: SStatus.INFEASIBLE,
            s_status.INF_OR_UNB: SStatus.INF_OR_UNB,
            s_status.UNKNOWN: SStatus.UNKNOWN,
            s_status.SUBOPTIMAL: SStatus.SUBOPTIMAL
        }

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

    @property
    def reference(self):
        """The reference wild type reaction flux values.

        :returns: A dictionary of wild type reaction flux values.

        """
        if self._reference is None:
            self._reference = self.simulate(
                method=SimulationMethod.pFBA).fluxes
        return self._reference

    @property
    def essential_reactions(self, min_growth=0.01):
        """Essential reactions are those when knocked out enable a biomass flux value above a minimal growth defined as
        a percentage of the wild type growth.

        :param float min_growth: Minimal percentage of the wild type growth value. Default 0.01 (1%).
        :returns: A list of essential reactions.

        """
        if self._essential_reactions is not None:
            return self._essential_reactions
        wt_solution = self.simulate()
        wt_growth = wt_solution.objective_value
        reactions = self.model.reactions.keys()
        self._essential_reactions = []
        for rxn in reactions:
            res = self.simulate(constraints={rxn: 0})
            if res:
                if (res.status == SStatus.OPTIMAL and res.objective_value < wt_growth * min_growth) \
                        or res.status == SStatus.INFEASIBLE:
                    self._essential_reactions.append(rxn)
        return self._essential_reactions

    @property
    def essential_genes(self, min_growth=0.01):
        """Essential genes are those when deleted enable a biomass flux value above a minimal growth defined as
        a percentage of the wild type growth.

        :param float min_growth: Minimal percentage of the wild type growth value. Default 0.01 (1%).
        :returns: A list of essential genes.

        """
        if self._essential_genes is not None:
            return self._essential_genes
        self._essential_genes = []
        wt_solution = self.simulate()
        wt_growth = wt_solution.objective_value
        genes = self.model.genes
        for gene in genes:
            active_genes = set(self.model.genes) - set([gene])
            active_reactions = self.evaluate_gprs(active_genes)
            inactive_reactions = set(
                self.model.reactions) - set(active_reactions)
            gr_constraints = {rxn: 0 for rxn in inactive_reactions}
            res = self.simulate(constraints=gr_constraints)
            if res:
                if (res.status == SStatus.OPTIMAL and res.objective_value < wt_growth * min_growth) \
                        or res.status == SStatus.INFEASIBLE:
                    self._essential_genes.append(gene)
        return self._essential_genes

    def evaluate_gprs(self, active_genes):
        """Returns the list of active reactions for a given list of active genes.

        :param list active_genes: List of genes identifiers.
        :returns: A list of active reaction identifiers.

        """
        active_reactions = []
        reactions = self.model.reactions
        for r_id, reaction in reactions.items():
            if reaction.gpr:
                if evaluate_expression_tree(str(reaction.gpr), active_genes):
                    active_reactions.append(r_id)
            else:
                active_reactions.append(r_id)
        return active_reactions

    def update(self):
        """Updates the model
        """
        self.model.update()

    def add_reaction(self, reaction, replace=True):
        """Adds a reaction to the model

        Args:
            reaction: The reaction, a Reframed reaction, to be added.
            replace (bool, optional): If the reaction should be replaced in case it is already defined.\
            Defaults to True.
        """
        self.model.add_reaction(reaction, replace=replace)

    def remove_reaction(self, r_id):
        """Removes a reaction from the model.

        Args:
            r_id (str): The reaction identifier.
        """
        self.model.remove_reaction(r_id)

    def get_metabolite_reactions(self, metabolite):
        """List all reactions that have the metabolite as reactant or product.

        Args:
            metabolite (str): The metabolite identifier.

        Returns:
            list: List of reactions
        """
        if not self._index_metabolites_reactions:
            self.__index_metabolites_reactions__()
        return self._index_metabolites_reactions[metabolite]

    def get_metabolite_compartment(self, metabolite):
        """Returns the compartment of a metabolite.

        Args:
            metabolite (str): The metabolite identifier.

        Returns:
            str: The compartment identifier.
        """
        return self.model.metabolites[metabolite].compartment

    def get_uptake_reactions(self):
        """
        :returns: The list of uptake reactions.

        """
        drains = self.get_drains()
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
        r_g = self.gene_reactions()
        genes = []
        for g, rxs in r_g.items():
            if set(rxs).issubset(set(trp_rxs)):
                genes.append(g)
        return genes

    def reverse_reaction(self, reaction_id):
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

    def gene_reactions(self):
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

    def get_reactions_for_genes(self, genes):
        """
        Returns the list of reactions catalysed by a list of genes

        :param list genes: A list of gene IDs.
        :returns: A list of reaction identifiers.

        """
        if not self._gene_to_reaction:
            self.gene_reactions()
        reactions = []
        for gene in genes:
            reactions.extend(self._gene_to_reaction[gene])
        return reactions

    def get_reaction_metabolites(self, reaction):
        '''
        Returns all metabolites of a given reaction

        :param reaction: reaction (str)
        :return: metabolites (dict)

        '''
        return self.model.reactions[reaction].stoichiometry

    def is_reactant(self, reaction, metabolite):
        '''
        Returns if a metabolite is reactant into a given reaction

        :param reaction: reaction (str)
        :param metabolite: metabolite (str)
        :return: bool
        '''
        if metabolite not in self.model.reactions[reaction].stoichiometry:
            raise KeyError("{} not in {}".format(metabolite, reaction))
        return self.model.reactions[reaction].stoichiometry[metabolite] < 0.0

    def is_product(self, reaction, metabolite):
        '''
        Returns if a metabolite is product into a given reaction.

        :param reaction: reaction (str)
        :param metabolite: metabolite (str)
        :return: bool
        '''
        if metabolite not in self.model.reactions[reaction].stoichiometry:
            raise KeyError("{} not in {}".format(metabolite, reaction))
        return self.model.reactions[reaction].stoichiometry[metabolite] > 0.0

    def __index_metabolites_reactions__(self):
        self._index_metabolites_reactions = {}
        for reaction in self.reactions:
            metabolites = self.get_reaction_metabolites(reaction)
            for metabolite in metabolites:
                if metabolite in self._index_metabolites_reactions:
                    self._index_metabolites_reactions[metabolite].append(
                        reaction)
                else:
                    self._index_metabolites_reactions[metabolite] = [reaction]

    def metabolite_reaction_lookup(self, force_recalculate=False):
        """ Return the network topology as a nested map from metabolite to reaction to coefficient.
        :return: a dictionary lookup table
        """

        if not self._m_r_lookup or force_recalculate:
            self._m_r_lookup = OrderedDict([(m_id, OrderedDict()) for m_id in self.metabolites])

            for r_id, reaction in self.model.reactions.items():
                for m_id, coeff in reaction.stoichiometry.items():
                    self._m_r_lookup[m_id][r_id] = coeff

        return self._m_r_lookup

    def set_objective(self, reaction):
        self.model.set_objective({reaction: 1})

    def get_objective(self):
        return list(self.objective.keys())

    def get_S(self):
        """
        Returns the S matrix as a numpy array
        :return: S matrix, np.array
        """

        return np.array(self.model.stoichiometric_matrix())

    def get_reaction_bounds(self, reaction):
        """
        Returns the bounds for a given reaction.
        :param reaction: str, reaction ID
        :return: lb(s), ub(s), tuple
        """

        if reaction in self.constraints:
            lb, ub = self.constraints[reaction]
        elif reaction in self.environmental_conditions:
            lb, ub = self.environmental_conditions[reaction]
        else:
            lb, ub = self.model.reactions[reaction].lb, self.model.reactions[reaction].ub

        return lb if lb > -np.inf else -999999, ub if ub < np.inf else 999999

    def get_bounds(self):
        """
        Returns the whole set of lower and upper bounds as numpy arrays

        :return: lb(s), ub(s), tuple of lists

        """
        lbs, ubs = list(zip(*[self.get_reaction_bounds(reaction) for reaction in self.model.reactions]))
        return list(lbs), list(ubs)

    def find_bounds(self):
        """
        Return the median upper and lower bound of the metabolic model.
        Bounds can vary from model to model. Cobrapy defaults to (-1000, 1000).
        """
        lower_bounds = np.asarray([rxn.lb for rxn in self.model.reactions], dtype=float)
        upper_bounds = np.asarray([rxn.ub for rxn in self.model.reactions], dtype=float)
        lower_bound = np.nanmedian(lower_bounds[lower_bounds != 0.0])
        upper_bound = np.nanmedian(upper_bounds[upper_bounds != 0.0])
        if np.isnan(lower_bound):
            LOGGER.warning("Could not identify a median lower bound.")
            lower_bound = -1000.0
        if np.isnan(upper_bound):
            LOGGER.warning("Could not identify a median upper bound.")
            upper_bound = 1000.0
        return lower_bound, upper_bound

    def find_unconstrained_reactions(self):
        """Return list of reactions that are not constrained at all."""
        lower_bound, upper_bound = self.find_bounds()
        return [
            rxn
            for rxn in self.model.reactions
            if rxn.lb <= lower_bound and rxn.ub >= upper_bound
        ]

    def get_boundary_reaction(self, metabolite):
        """
        Finds the boundary reaction associated with an extracellular metabolite.
        If none is found, None is returned

        :param metabolite: str, metabolite ID

        :returns: reaction, str or  None
        """

        # reaction.reaction_type.value == 'exchange'
        # len(reaction.stoichiometry) == 1

        for reaction in self.get_metabolite_reactions(metabolite):

            if reaction in self.medium or self.model.reactions[reaction].reaction_type.value == 'exchange' \
                    or len(self.model.reactions[reaction].stoichiometry) == 1:
                return reaction

        return None

    # Simulate
    def simulate(self, objective=None, method=SimulationMethod.FBA,
                 maximize=True, constraints=None, reference=None,
                 scalefactor=None, solver=None):
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
            simul_constraints.update(constraints)
        if self.constraints:
            simul_constraints.update(self.constraints)
        if self.environmental_conditions:
            simul_constraints.update(self.environmental_conditions)

        a_solver = solver
        if not self._reset_solver and not a_solver:
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

        # TODO: simplifly ...
        if method in [SimulationMethod.lMOMA, SimulationMethod.MOMA, SimulationMethod.ROOM] and reference is None:
            reference = self.reference

        if method == SimulationMethod.FBA:
            solution = FBA(self.model, objective=objective, minimize=not maximize,
                           constraints=simul_constraints, solver=a_solver)
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

        status = self.__status_mapping[solution.status]

        result = SimulationResult(self.model, solution.fobj, fluxes=solution.values, status=status,
                                  envcond=self.environmental_conditions, model_constraints=self.constraints,
                                  simul_constraints=constraints, maximize=maximize, method=method)
        return result

    def FVA(self, obj_frac=0.9, reactions=None, constraints=None, loopless=False, internal=None, solver=None,
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
        _constraints = {}
        if self.environmental_conditions:
            _constraints.update(self.environmental_conditions)
        if constraints:
            _constraints.update(constraints)
        from reframed.cobra.variability import FVA
        res = FVA(self.model, obj_frac=obj_frac, reactions=reactions,
                  constraints=_constraints, loopless=loopless, internal=internal, solver=solver)
        if format == 'df':
            import pandas as pd
            e = res.items()
            f = [[a, b, c] for a, [b, c] in e]
            df = pd.DataFrame(f, columns=['Reaction ID', 'Minimum', 'Maximum'])
            return df
        return res


class GeckoSimulation(Simulation):

    def __init__(self, model: GeckoModel, objective=None, envcond=None, constraints=None, solver=None, reference=None,
                 reset_solver=ModelConstants.RESET_SOLVER, protein_prefix=None):
        super(GeckoSimulation, self).__init__(
            model, objective, envcond, constraints, solver, reference, reset_solver)
        self.protein_prefix = protein_prefix if protein_prefix else 'draw_prot_'
        self._essential_proteins = None

    @property
    def proteins(self):
        return self.model.proteins

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
        for p in proteins:
            rxn = "{}{}".format(self.protein_prefix, p)
            res = self.simulate(constraints={rxn: 0})
            if res:
                if (res.status == SStatus.OPTIMAL and res.objective_value < wt_growth * min_growth) \
                        or res.status == SStatus.INFEASIBLE:
                    self._essential_proteins.append(rxn)
        return self._essential_proteins

    def protein_reactions(self, protein):
        """
        Returns the list of reactions associated to a protein
        """
        reactions = []
        for r_id, rxn in self.model.reactions.items():
            lsub = rxn.get_substrates()
            for m in lsub:
                if protein in m:
                    reactions.append(r_id)
        return reactions

    def reverse_reaction(self, reaction_id):
        """
        Identify if a reaction is reversible and returns the
        reverse reaction if it is the case.

        :returns: A reaction identifier or None.
        """
        f, d = zip(*self.model.protein_rev_reactions.values())
        if reaction_id in f:
            return d[f.index(reaction_id)]
        elif reaction_id in d:
            return f[d.index(reaction_id)]
        else:
            return None
