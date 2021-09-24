"""
Simulation for COBRApy models
"""
import logging
from collections import OrderedDict

import numpy as np
from cobra.core.model import Model
from cobra.core.solution import Solution
from cobra.flux_analysis import pfba, moma, room

from . import get_default_solver, SimulationMethod, SStatus
from .simulation import Simulator, SimulationResult, ModelContainer
from ..util.constants import ModelConstants
from ..util.parsing import evaluate_expression_tree
from tqdm import tqdm


LOGGER = logging.getLogger(__name__)


class CobraModelContainer(ModelContainer):
    """ A basic container for COBRApy models.

    :param model: A metabolic model.

    """

    def __init__(self, model: Model):
        if not isinstance(model, Model):
            raise ValueError("The model is not an instance of cobrapy Model")
        self.model = model

    @property
    def id(self):
        return self.model.id

    @property
    def reactions(self):
        return [rxn.id for rxn in self.model.reactions]

    def get_reaction(self, r_id):
        rxn = self.model.reactions.get_by_id(r_id)
        stoichiometry = {met.id: val for met, val in rxn.metabolites.items()}
        res = {'id': r_id, 'name': rxn.name, 'lb': rxn.lower_bound,
               'ub': rxn.upper_bound, 'stoichiometry': stoichiometry}
        res['gpr'] = rxn.gene_reaction_rule if rxn.gene_reaction_rule is not None else None
        return res

    @property
    def genes(self):
        return [gene.id for gene in self.model.genes]

    def get_gene(self, g_id):
        g = self.model.genes.get_by_id(g_id)
        res = {'id': g_id, 'name': g.name}
        return res

    @property
    def metabolites(self):
        return [met.id for met in self.model.metabolites]

    def get_metabolite(self, m_id):
        met = self.model.metabolites.get_by_id(m_id)
        res = {'id': m_id, 'name': met.name, 'compartment': met.compartment, 'formula': met.formula}
        return res

    @property
    def medium(self):
        return self.model.medium

    @property
    def compartments(self):
        return self.model.compartments

    def get_compartment(self, c_id):
        c = self.model.compartments[c_id]
        from cobra.medium import find_external_compartment
        e = find_external_compartment(self.model)
        res = {'id': c_id, 'name': c, 'external': (e == c_id)}
        return res

    def get_gpr(self, reaction_id):
        """Returns the gpr rule (str) for a given reaction ID.

        :param str reaction_id: The reaction identifier.
        :returns: A string representation of the GPR rule.

        """
        if reaction_id not in self.reactions:
            raise ValueError(f"Reactions {reaction_id} does not exist")
        reaction = self.model.reactions.get_by_id(reaction_id)
        if reaction.gene_reaction_rule:
            return str(reaction.gene_reaction_rule)
        else:
            return None

    def get_substrates(self, rxn_id):
        reaction = self.model.reactions.get_by_id(rxn_id)
        return {k.id: v for k, v in reaction.metabolites.items() if v < 0}

    def get_products(self, rxn_id):
        reaction = self.model.reactions.get_by_id(rxn_id)
        return {k.id: v for k, v in reaction.metabolites.items() if v > 0}

    def get_exchange_reactions(self):
        rxns = [r.id for r in self.model.exchanges]
        return rxns


class Simulation(CobraModelContainer, Simulator):
    """Generic Simulation class for cobra Model.
       Defines the simulation conditions, and makes available a set of methods.

    :param model: An metabolic model instance.

    Optional:

    :param dic envcond: Dictionary of environmental conditions.
    :param dic constraints: A dictionary of reaction constraints.
    :param solver: An instance of the LP solver.
    :param dic reference: A dictionary of the wild type flux values.

    """

    def __init__(self, model: Model, envcond=None, constraints=None, solver=None, reference=None,
                 reset_solver=ModelConstants.RESET_SOLVER):

        if not isinstance(model, Model):
            raise ValueError("Model is incompatible or inexistent")

        self.model = model
        self.model.solver = get_default_solver()
        self.environmental_conditions = OrderedDict() if envcond is None else envcond
        self.constraints = OrderedDict() if constraints is None else constraints
        self.solver = solver
        self._essential_reactions = None
        self._essential_genes = None
        self._gene_to_reaction = None
        self._reference = reference
        self.__status_mapping = {
            'optimal': SStatus.OPTIMAL,
            'unbounded': SStatus.UNBOUNDED,
            'infeasible': SStatus.INFEASIBLE,
            'infeasible_or_unbounded': SStatus.INF_OR_UNB,
            'suboptimal': SStatus.SUBOPTIMAL,
            'unknown': SStatus.UNKNOWN
        }
        self.solver = solver
        self._reset_solver = reset_solver
        self.reverse_sintax = []
        self._m_r_lookup = None

        self._MAX_STR = 'maximize'
        self._MIN_STR = 'minimize'

    @property
    def objective(self):
        from cobra.util.solver import linear_reaction_coefficients
        d = dict(linear_reaction_coefficients(self.model))
        return {k.id: v for k, v in d.items()}

    @objective.setter
    def objective(self, objective):
        if isinstance(objective, str):
            self.model.objective = objective
        elif isinstance(objective, dict):
            from cobra.util.solver import set_objective
            linear_coef = {self.model.reactions.get_by_id(r_id): v for r_id, v in objective.items()}
            set_objective(self.model, linear_coef)
        else:
            raise ValueError(
                'The objective must be a reaction identifier or a dictionary of \
                reaction identifier with respective coeficients.')

    @property
    def reference(self):
        """The reference wild type reaction flux values.

        :returns: A dictionary of wild type reaction flux values.

        """
        if self._reference is None:
            self._reference = self.simulate(
                method=SimulationMethod.pFBA).fluxes
        return self._reference

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
        reactions = self.reactions
        self._essential_reactions = []
        for rxn in tqdm(reactions):
            res = self.simulate(constraints={rxn: 0})
            if res:
                if (res.status == SStatus.OPTIMAL and res.objective_value < wt_growth * min_growth) \
                        or res.status == SStatus.INFEASIBLE:
                    self._essential_reactions.append(rxn)
        return self._essential_reactions

    def essential_genes(self, min_growth=0.01):
        """Essential genes are those when deleted enable a biomass flux value above a minimal growth defined as a
        percentage of the wild type growth.

        :param float min_growth: Minimal percentage of the wild type growth value. Default 0.01 (1%).
        :returns: A list of essential genes.

        """
        if self._essential_genes is not None:
            return self._essential_genes
        self._essential_genes = []
        wt_solution = self.simulate()
        wt_growth = wt_solution.objective_value
        genes = self.genes
        for gene in tqdm(genes):
            active_genes = set(self.genes) - {gene}
            active_reactions = self.evaluate_gprs(active_genes)
            inactive_reactions = set(self.reactions) - set(active_reactions)
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
        for r_id in self.reactions:
            reaction = self.model.reactions.get_by_id(r_id)
            if reaction.gene_reaction_rule:
                if evaluate_expression_tree(str(reaction.gene_reaction_rule), active_genes):
                    active_reactions.append(r_id)
            else:
                active_reactions.append(r_id)
        return active_reactions

    def add_metabolite(self, id, formula=None, name=None, compartment=None):
        from cobra import Metabolite
        meta = Metabolite(id, formula=formula, name=name, compartment=compartment)
        self.model.add_metabolites([meta])

    def add_reaction(self, rxn_id, stoichiometry, lb=ModelConstants.REACTION_LOWER_BOUND,
                     ub=ModelConstants.REACTION_UPPER_BOUND, replace=True, *kwargs):
        """Adds a reaction to the model

        Args:
            rxn_id: The reaction identifier
            stoichiometry: The reaction stoichiometry, a dictionary of species: coefficient
            lb: Reaction flux lower bound, defaults to ModelConstants.REACTION_LOWER_BOUND
            ub: Reaction flux upper bound, defaults to ModelConstants.REACTION_UPPER_BOUND
            replace(bool, optional): If the reaction should be replaced in case it is already defined.\
                Defaults to True.
        """
        from cobra import Reaction
        reaction = Reaction(rxn_id)
        reaction.add_metabolites(stoichiometry)
        reaction.lower_bound = lb
        reaction.upper_bound = ub
        self.model.add_reaction(reaction)

    def remove_reaction(self, r_id):
        """Removes a reaction from the model.

        Args:
            r_id (str): The reaction identifier.
        """
        self.model.remove_reactions(r_id)

    def get_uptake_reactions(self):
        """
        :returns: The list of uptake reactions.
        """
        drains = self.get_exchange_reactions()
        rxns = [r for r in drains if self.model.reactions.get_by_id(r).reversibility
                or ((self.model.reactions.get_by_id(r).lower_bound is None
                     or self.model.reactions.get_by_id(r).lower_bound < 0)
                    and len(self.model.reactions.get_by_id(r).reactants) > 0)
                or ((self.model.reactions.get_by_id(r).upper_bound is None
                     or self.model.reactions.get_by_id(r).upper_bound > 0)
                    and len(self.model.reactions.get_by_id(r).products) > 0)
                ]
        return rxns

    def get_transport_reactions(self):
        """
        :returns: The list of transport reactions.
        """
        transport_reactions = []
        for rx in self.reactions:
            s_set = set()
            p_set = set()
            s = self.model.reactions.get_by_id(rx).reactants
            for x in s:
                s_set.add(x.compartment)
            p = self.model.reactions.get_by_id(rx).products
            for x in p:
                p_set.add(x.compartment)
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
        Identify if a reaction is reversible and returns the reverse reaction if it is the case.

        :param reaction_id: A reaction identifier.
        :returns: A reaction identifier or None.

        """
        rxn = self.model.reactions.get_by_id(reaction_id)
        reactions = self.reactions
        if rxn.lower_bound < 0:
            return reaction_id
        else:
            rev = rxn.reverse_id
            if rev in reactions:
                return rev
        return None

    def gene_reactions(self):
        """
        :returns: A map of genes to reactions.
        """
        if not self._gene_to_reaction:
            gr = OrderedDict()
            for rxn_id in self.reactions:
                rxn = self.model.reactions.get_by_id(rxn_id)
                genes = rxn.genes
                for g in genes:
                    if g in gr.keys():
                        gr[g.id].append(rxn_id)
                    else:
                        gr[g.id] = [rxn_id]
            self._gene_to_reaction = gr
        return self._gene_to_reaction

    def metabolite_reaction_lookup(self, force_recalculate=False):
        """ Return the network topology as a nested map from metabolite to reaction to coefficient.
        :return: a dictionary lookup table
        """

        if not self._m_r_lookup or force_recalculate:
            self._m_r_lookup = OrderedDict([(m_id, OrderedDict()) for m_id in self.metabolites])
            for _, reaction in enumerate(self.model.reactions):
                for m, coeff in reaction.metabolites.items():
                    self._m_r_lookup[m.id][reaction.id] = coeff

        return self._m_r_lookup

    def metabolite_elements(self, metabolite_id):
        return self.model.metabolites.get_by_id(metabolite_id).elements()

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
            lb, ub = self.model.reactions.get_by_id(reaction).bounds
        return lb if lb > -np.inf else ModelConstants.REACTION_LOWER_BOUND,\
            ub if ub < np.inf else ModelConstants.REACTION_UPPER_BOUND

    def set_reaction_bounds(self, reaction, lb=None, ub=None):
        """
        Sets the bounds for a given reaction.
        :param reaction: str, reaction ID
        :param float lb: lower bound 
        :param float ub: upper bound
        """
        rxn = self.model.reactions.get_by_id(r)
        rxn.bounds = (lb, ub)

    def find_bounds(self):
        """
        Return the median upper and lower bound of the metabolic model.
        Bounds can vary from model to model. Cobrapy defaults to (-1000, 1000).
        """
        lower_bounds = np.asarray([rxn.lower_bound for rxn in self.model.reactions], dtype=float)
        upper_bounds = np.asarray([rxn.upper_bound for rxn in self.model.reactions], dtype=float)
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
            rxn.id
            for rxn in self.model.reactions
            if rxn.lower_bound <= lower_bound and rxn.upper_bound >= upper_bound
        ]

    # The simulator
    def simulate(self, objective=None, method=SimulationMethod.FBA, maximize=True,
                 constraints=None, reference=None, scalefactor=None, solver=None):
        '''
        Simulates a phenotype when applying a set constraints using the specified method.

        :param dic objective: The simulation objective. If none, the model objective is considered.
        :param method: The SimulationMethod (FBA, pFBA, lMOMA, etc ...)
        :param boolean maximize: The optimization direction.
        :param dic constraints: A dictionary of constraints to be applied to the model.
        :param dic reference: A dictionary of reaction flux values.
        :param float scalefactor: A positive scaling factor for the solver. Default None.

        '''

        if not objective:
            objective = self.model.objective
        elif isinstance(objective, dict) and len(objective) > 0:
            objective = next(iter(objective.keys()))

        simul_constraints = OrderedDict()
        if constraints:
            simul_constraints.update(constraints)
        if self.constraints:
            simul_constraints.update(self.constraints)
        if self.environmental_conditions:
            simul_constraints.update(self.environmental_conditions)

        with self.model as model:
            model.objective = objective
            for rxn in list(simul_constraints.keys()):
                reac = model.reactions.get_by_id(rxn)
                # constraints defined as a tuple (lower_bound, upper_bound) or as a single float value
                if isinstance(simul_constraints.get(rxn), tuple):
                    reac.bounds = (simul_constraints.get(
                        rxn)[0], simul_constraints.get(rxn)[1])
                else:
                    reac.bounds = (simul_constraints.get(
                        rxn), simul_constraints.get(rxn))
            # NOTE: If working directly over optlang use 'max' and 'min'
            # such is the case with pytfa.core.Model... need to find some workaround
            objective_sense = self._MAX_STR if maximize else self._MIN_STR
            if method == SimulationMethod.FBA:
                solution = model.optimize(objective_sense=objective_sense)
            elif method == SimulationMethod.pFBA:
                # make fraction_of_optimum a configurable parameter?
                solution = pfba(model)
            elif method == SimulationMethod.lMOMA or method == SimulationMethod.MOMA:
                s = None
                if not maximize:
                    s = model.optimize(objective_sense=objective_sense)
                linear = True if method == SimulationMethod.lMOMA else False
                solution = moma(model, solution=s, linear=linear)
            elif method == SimulationMethod.ROOM:
                solution = room(model)
            # Special case in which only the simulation context is required without any simulation result
            elif method == SimulationMethod.NONE:
                solution = Solution(None, 'unknown', None)
                pass
            else:
                raise Exception(
                    "Unknown method to perform the simulation.")

        status = self.__status_mapping[solution.status]

        result = SimulationResult(model, solution.objective_value, fluxes=solution.fluxes.to_dict(OrderedDict),
                                  status=status, envcond=self.environmental_conditions,
                                  model_constraints=self.constraints,
                                  simul_constraints=constraints,
                                  maximize=maximize,
                                  method=method)
        return result

    def FVA(self, obj_frac=0.9, reactions=None, constraints=None, loopless=False, internal=None, solver=None,
            format='dict'):
        """ Flux Variability Analysis (FVA).

        :param model: An instance of a constraint-based model.
        :param float obj_frac: The minimum fraction of the maximum growth rate (default 0.9). Requires that the \
            objective value is at least the fraction times maximum objective value. A value of 0.85 for instance \
            means that the objective has to be at least at 85% percent of its maximum.
        :param list reactions: List of reactions to analyze (default: all).
        :param dic constraints: Additional constraints (optional).
        :param boolean loopless: Run looplessFBA internally (very slow) (default: false).
        :param list internal: List of internal reactions for looplessFBA (optional).
        :param solver: A pre-instantiated solver instance (optional).
        :param format: The return format: 'dict', returns a dictionary,'df' returns a data frame.
        :returns: A dictionary of flux variation ranges.

        """
        from cobra.flux_analysis.variability import flux_variability_analysis

        simul_constraints = {}
        if self.environmental_conditions:
            simul_constraints.update(self.environmental_conditions)
        if constraints:
            simul_constraints.update(constraints)

        with self.model as model:

            if simul_constraints:
                for rxn in list(simul_constraints.keys()):
                    reac = model.reactions.get_by_id(rxn)
                    if isinstance(simul_constraints.get(rxn), tuple):
                        reac.bounds = (simul_constraints.get(
                            rxn)[0], simul_constraints.get(rxn)[1])
                    else:
                        reac.bounds = (simul_constraints.get(
                            rxn), simul_constraints.get(rxn))

            df = flux_variability_analysis(
                model, reaction_list=reactions, loopless=loopless, fraction_of_optimum=obj_frac)

        variability = {}
        for r_id in reactions:
            variability[r_id] = [
                float(df.loc[r_id][0]), float(df.loc[r_id][1])]

        if format == 'df':
            import pandas as pd
            e = variability.items()
            f = [[a, b, c] for a, [b, c] in e]
            df = pd.DataFrame(f, columns=['Reaction ID', 'Minimum', 'Maximum'])
            return df
        else:
            return variability

    def set_objective(self, reaction):
        self.model.objective = reaction


class GeckoSimulation(Simulation):
    """
    Simulator for geckopy.gecko.GeckoModel
    """

    def __init__(self, model, envcond=None, constraints=None, solver=None, reference=None,
                 reset_solver=ModelConstants.RESET_SOLVER, protein_prefix=None):
        try:
            from geckopy.gecko import GeckoModel
            if not isinstance(model, GeckoModel):
                raise ValueError("The model is not an instance of geckopy.gecko.GeckoModel")
        except ImportError:
            raise RuntimeError("The geckopy package is not installed.")

        super(GeckoSimulation, self).__init__(
            model, envcond, constraints, solver, reference, reset_solver)
        self.protein_prefix = protein_prefix if protein_prefix else 'draw_prot_'
        self._essential_proteins = None
        self._protein_rev_reactions = None

    @property
    def proteins(self):
        return list(self.model.proteins)

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
                if (res.status == SStatus.OPTIMAL and res.objective_value < wt_growth * min_growth) or \
                        res.status == SStatus.INFEASIBLE:
                    self._essential_proteins.append(rxn)
        return self._essential_proteins

    def protein_reactions(self, protein):
        """
        Returns the list of reactions associated to a protein.
        """
        reactions = []
        for r_id, rxn in self.model.reactions.items():
            lsub = rxn.reactants
            for m in lsub:
                if protein in m:
                    reactions.append(r_id)
        return reactions

    def reverse_reaction(self, reaction_id):
        """
        Identify if a reaction is reversible and returns the reverse reaction if it is the case.

        :param reaction_id: A reaction identifier.
        :returns: The reverse reaction identifier if exists or None.
        """
        f, d = zip(*self.protein_rev_reactions.values())
        if reaction_id in f:
            return d[f.index(reaction_id)]
        elif reaction_id in d:
            return f[d.index(reaction_id)]
        else:
            return None

    @property
    def protein_rev_reactions(self):
        """
        Pairs of reverse reactions associated with a protein

        :returns: A dictionary which identifies for each protein (key) the list of reversible reactions pairs.

        """
        if not self._protein_rev_reactions:
            proteins = self.model.proteins
            reactions = self.reactions
            in_sub = {}
            for p in proteins:
                sub = []
                for r_id in reactions:
                    rxn = self.model.reactions.get_by_id(r_id)
                    lsub = rxn.reactants
                    for m in lsub:
                        if p in m.id:
                            sub.append(r_id)
                in_sub[p] = sub
            pairs = {}
            for k, s in in_sub.items():
                revs = [r for r in s if '_REV' in r]
                if len(revs) > 0:
                    for r in revs:
                        la = [a for a in s if r.replace('_REV', '') == a]
                        la.append(r)
                        if len(la) == 2:
                            if k in pairs.keys():
                                pairs[k].append((la[0], la[1]))
                            else:
                                pairs[k] = [(la[0], la[1])]
            self._protein_rev_reactions = pairs
        return self._protein_rev_reactions

    def getKcat(self, protein):
        """ Returns a dictionary of reactions and respective Kcat for a specific enzyme·
        """
        m_r = self.metabolite_reaction_lookup()
        r_d = m_r[protein]
        return {k: v for k, v in r_d.items() if v < 0}
