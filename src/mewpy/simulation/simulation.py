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
Defines an interface for the simulators, objects that wrap metabolic phenotipe
simulation toolboxes.

Author: Vitor Pereira
##############################################################################
"""
from abc import abstractmethod, ABC
from collections import OrderedDict
from enum import Enum
from joblib import Parallel, delayed
from tqdm import tqdm
from copy import deepcopy
import math

from ..util.parsing import evaluate_expression_tree
from ..util.process import cpu_count

from typing import List

class SimulationMethod(Enum):
    FBA = 'FBA'
    pFBA = 'pFBA'
    MOMA = 'MOMA'
    lMOMA = 'lMOMA'
    ROOM = 'ROOM'
    NONE = 'NONE'

    def __eq__(self, other):
        """Overrides equal to enable string name comparison.
        Allows to seamlessly use:
            SimulationMethod.FBA = SimulationMethod.FBA
            SimulationMethod.FBA = 'FBA'
        without requiring an additional level of comparison (SimulationMethod.FBA.name = 'FBA')
        """
        if isinstance(other, SimulationMethod):
            return super().__eq__(other)
        elif isinstance(other, str):
            return self.name == other
        else:
            return False

    def __hash__(self):
        return hash(self.name)


class SStatus(Enum):
    """ Enumeration of possible solution status. """
    OPTIMAL = 'Optimal'
    UNKNOWN = 'Unknown'
    SUBOPTIMAL = 'Suboptimal'
    UNBOUNDED = 'Unbounded'
    INFEASIBLE = 'Infeasible'
    INF_OR_UNB = 'Infeasible or Unbounded'

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name


class SimulationInterface(ABC):

    @abstractmethod
    def simulate(**kwargs):
        """Abstract method to run a simulation.

        :returns: A SimulationResult.

        """
        raise NotImplementedError


class ModelContainer(ABC):
    """Interface for Model container.
       Provides an abstraction from models implementations.
    """

    _r_prefix = ''
    _m_prefix = ''
    _g_prefix = ''


    @property
    def reactions(self):
        """
        :returns: A list of reaction identifiers.

        """
        raise NotImplementedError

    @property
    def genes(self):
        """
        :returns: A list of gene identifiers.

        """
        raise NotImplementedError

    @property
    def proteins(self):
        """
        :returns: A list of proteins identifiers.

        """
        raise NotImplementedError

    @property
    def metabolites(self):
        """
        :returns: A list of metabolite identifiers.

        """
        raise NotImplementedError

    @property
    def compartments(self):
        """
        :returns: A list of compartments identifiers.

        """
        raise NotImplementedError

    @property
    def medium(self):
        raise NotImplementedError

    def get_reaction(self, r_id):
        raise NotImplementedError

    def get_gene(self, g_id):
        raise NotImplementedError

    def get_compartment(self, c_id):
        raise NotImplementedError

    def get_reaction_metabolites(self, r_id):
        rxn = self.get_reaction(r_id)
        return rxn['stoichiometry']

    def get_substrates(self, rxn_id):
        met = self.get_reaction_metabolites(rxn_id)
        return {m_id: coeff for m_id, coeff in met.items() if coeff < 0}

    def get_products(self, rxn_id):
        met = self.get_reaction_metabolites(rxn_id)
        return {m_id: coeff for m_id, coeff in met.items() if coeff > 0}

    def get_exchange_reactions(self):
        return NotImplementedError

    def get_gene_reactions(self):
        raise NotImplementedError

    def get_gpr(self, rxn_id):
        """
        :returns: A string representation of the reaction GPR if exists None otherwise.
        """
        rxn = self.get_reaction(rxn_id)
        return rxn['gpr']

    def summary(self):
        print(f"Metabolites: {len(self.metabolites)}")
        print(f"Reactions: {len(self.reactions)}")
        print(f"Genes: {len(self.genes)}")

    def set_objective(self, reaction):
        raise NotImplementedError


class Simulator(ModelContainer, SimulationInterface):
    """
    Interface for simulators
    """

    @abstractmethod
    def simulate(self, objective=None, method=SimulationMethod.FBA, maximize=True, constraints=None, reference=None,
                 solver=None, **kwargs):
        """Abstract method to run a phenotype simulation.

        :returns: A SimulationResult.

        """
        raise NotImplementedError

    @abstractmethod
    def FVA(self, reactions=None, obj_frac=0, constraints=None, loopless=False, internal=None, solver=None):
        """ Abstract method to run Flux Variability Analysis (FVA).

        :returns: A dictionary of flux range values.

        """
        raise NotImplementedError

    def simulate_mp(self, objective=None, method=SimulationMethod.FBA, maximize=True, constraints_list=None,
                    reference=None, solver=None, jobs=None, desc="Parallel Simulation", **kwargs):
        """Parallel phenotype simulations.

        :param (dict) objective: The simulations objective. If none, the model objective is considered.
        :param (SimulationMethod) method: The SimulationMethod (FBA, pFBA, lMOMA, etc ...)
        :param (boolean) maximize: The optimization direction
        :param (dict) constraints: A dictionary of constraints to be applied to the model.
        :param (dict) reference: A dictionary of reaction flux values.
        :param (Solver) solver: An instance of the solver.
        :param (int) jobs: The number of parallel jobs.
        :param (str) desc: Description to present in tdqm. 
        """
        constraints_list = [None] if not constraints_list else constraints_list
        jobs = jobs if jobs else cpu_count()
        print(f"Using {jobs} jobs")
        from ..util.utilities import tqdm_joblib
        with tqdm_joblib(tqdm(desc=desc, total=len(constraints_list))) as progress_bar:
            res = Parallel(n_jobs=jobs)(delayed(simulate)(self.model, self.environmental_conditions, objective, method,
                                                          maximize, constraints, reference, solver, **kwargs) for constraints in constraints_list)
        return res

    @abstractmethod
    def get_reaction_bounds(self, r_id):
        raise NotImplementedError

    def set_environmental_conditions(self, medium):
        for k, v in medium.items():
            if isinstance(v, tuple):
                lb, ub = v
            elif isinstance(v, (float, int)):
                lb, ub = v, v
            else:
                raise ValueError(f"{v} is an inappropiate bound.")
            self.set_reaction_bounds(k, lb, ub, False)
        self._environmental_conditions = medium

    @abstractmethod
    def metabolite_reaction_lookup(self, force_recalculate=False):
        raise NotImplementedError

    def find(self, pattern=None, sort=False, find_in='r'):
        """A user friendly method to find metabolites, reactions or genes in the model.

        :param pattern: The pattern which can be a regular expression, defaults to None in which case all entries are listed.
        :type pattern: str, optional
        :param sort: if the search results should be sorted, defaults to False
        :type sort: bool, optional
        :param find_in: The search set: 'r' for reactions, 'm' for metabolites, 'g' for genes, defaults to 'r'
        :type find_in: str, optional
        :return: the search results
        :rtype: pandas dataframe
        """
        if find_in == 'm':
            values = self.metabolites
        elif find_in == 'g':
            values = self.genes
        else:
            values = self.reactions
        if pattern:
            import re
            if isinstance(pattern, list):
                patt = '|'.join(pattern)
                re_expr = re.compile(patt)
            else:
                re_expr = re.compile(pattern)
            values = [x for x in values if re_expr.search(x) is not None]
        if sort:
            values.sort()

        import pandas as pd
        if find_in == 'm':
            data = [self.get_metabolite(x) for x in values]
        elif find_in == 'g':
            data = [self.get_gene(x) for x in values]
        else:
            data = [self.get_reaction(x) for x in values]

        if data:
            df = pd.DataFrame(data)
            df = df.set_index(df.columns[0])
        else:
            df = pd.DataFrame()
        return df

    def find_genes(self, pattern=None, sort=False):
        return self.find(pattern=pattern, sort=sort, find_in='g')

    def find_metabolites(self, pattern=None, sort=False):
        return self.find(pattern=pattern, sort=sort, find_in='m')

    def find_reactions(self, pattern=None, sort=False):
        return self.find(pattern=pattern, sort=sort, find_in='r')

    def is_essential_reaction(self, rxn, min_growth=0.01):
        res = self.simulate(constraints={rxn: 0}, slim=True)
        return res is None or math.isnan(res) or res < min_growth

    def essential_reactions(self, min_growth=0.01):
        """Essential reactions are those when knocked out enable a biomass flux value above a minimal growth defined as
        a percentage of the wild type growth.

        :param float min_growth: Minimal percentage of the wild type growth value. Default 0.01 (1%).
        :returns: A list of essential reactions.
        """
        essential = getattr(self, '_essential_reactions', None)
        if essential is not None:
            return essential
        essential = []
        for rxn in tqdm(self.reactions):
            if self.is_essential_reaction(rxn, min_growth=min_growth):
                essential.append(rxn)
        self._essential_reactions = essential
        return self._essential_reactions

    def evaluate_gprs(self, active_genes):
        """Returns the list of active reactions for a given list of active genes.

        :param list active_genes: List of genes identifiers.
        :returns: A list of active reaction identifiers.

        """
        active_reactions = []
        reactions = self.reactions
        for r_id in reactions:
            gpr = self.get_gpr(r_id)
            if gpr:
                if evaluate_expression_tree(str(gpr), active_genes):
                    active_reactions.append(r_id)
            else:
                active_reactions.append(r_id)
        return active_reactions

    def is_essential_gene(self, gene, min_growth=0.01):
        gr = self.get_gene_reactions()
        active_genes = set(self.genes) - {gene}
        rxns = gr.get(gene, [])
        if not rxns:
            return False
        inactive_reactions = []
        for r_id in rxns:
            gpr = self.get_gpr(r_id)
            if gpr and not evaluate_expression_tree(str(gpr), active_genes):
                inactive_reactions.append(r_id)
        constraints = {rxn: 0 for rxn in inactive_reactions}
        res = self.simulate(constraints=constraints, slim=True)
        return res is None or math.isnan(res) or res < min_growth

    def essential_genes(self, min_growth=0.01):
        """Essential genes are those when deleted enable a biomass flux value above a minimal growth defined as
        a percentage of the wild type growth.

        :param float min_growth: Minimal percentage of the wild type growth value. Default 0.01 (1%).
        :returns: A list of essential genes.

        """
        essential = getattr(self, '_essential_genes', None)
        if essential is not None:
            return essential
        essential = []
        for g in tqdm(self.genes):
            if self.is_essential_gene(g, min_growth=min_growth):
                essential.append(g)
        self._essential_genes = essential
        return self._essential_genes

    @property
    def reference(self):
        """The reference wild type reaction flux values.

        :returns: A dictionary of wild type reaction flux values.

        """
        ref = getattr(self, '_reference', None)
        if ref is not None:
            return ref
        self._reference = self.simulate(method="pFBA").fluxes
        return self._reference

    def create_empty_model(self, model_id: str):
        return NotImplementedError

    def get_external_metabolites(self):
        external =[]
        ext_com = [c_id for c_id in self.compartments 
                   if self.get_compartment(c_id).external==True]
        for m_id in self.metabolites:
            if self.get_metabolite(m_id).compartment in ext_com:
                external.append(m_id)
        return external

    def blocked_reactions(self, constraints=None, reactions=None, abstol=1e-9):
        """ Find all blocked reactions in a model
        
        :param (dict) constraints: additional constraints (optional)
        :param (list) reactions: List of reactions which will be tested (default: None, test all reactions)
        :param (float) abstol: absolute tolerance (default: 1e-9)

        :returns: A list blocked reactions
        """

        variability = self.FVA(obj_frac=0, reactions=reactions, constraints=constraints)

        return [r_id for r_id, (lb, ub) in variability.items() if (abs(lb) + abs(ub)) < abstol]

    def get_metabolite_reactions(self,m_id: str) -> List[str]:
        """Returns the list or reactions that produce or consume a metabolite.

        :param m_id: the metabolite identifier
        :return: A list of reaction identifiers
        """
        m_r = self.metabolite_reaction_lookup()
        return list(m_r[m_id].keys())
    
    def get_metabolite_producers(self, m_id: str) -> List[str]:
        """Returns the list or reactions that produce a metabolite.

        :param m_id: the metabolite identifier
        :return: A list of reaction identifiers
        """
        m_r = self.metabolite_reaction_lookup()
        return [k for k, v in m_r[m_id].items() if v > 0]

    def get_metabolite_consumers(self, m_id):
        """Returns the list or reactions that consume a metabolite.

        :param m_id: the metabolite identifier
        :return: A list of reaction identifiers
        """
        m_r = self.metabolite_reaction_lookup()
        return [k for k, v in m_r[m_id].items() if v < 0]
    
    def copy(self):
        """Retuns a copy of the Simulator instance."""
        return deepcopy(self)


class SimulationResult(object):
    """Class that represents simulation results and performs operations over them."""

    def __init__(self, model, objective_value, fluxes=None, status=None, envcond=None, model_constraints=None,
                 simul_constraints=None, maximize=True, method=None, shadow_prices=None):
        """
        :param model: A model instance.
        :param objective_value: The phenotype simulation objective value.
        :param dict fluxes: A dictionary of reaction fluxes values.
        :param status: The LP status.
        :param dict envcond: The environmental conditions of the phenotype simulation.
        :param dict model_constraints: Possible persistent additional constraints.
        :param dict simul_constraints: The simulation constraints.
        :param boolean maximize: Optimization direction.
        :param SimulationMethod method: The phenotypic methos
        :param dict shadow_prices: shadow prices
        """
        self.model = model
        self.objective_value = objective_value
        self.status = status
        self.fluxes = fluxes
        self.maximize = maximize
        # Environmental conditions
        self.envcond = envcond if envcond else OrderedDict()
        # Persistent additional constraints not included in the environmental conditions
        self.model_constraints = model_constraints if model_constraints else OrderedDict()
        # Constraints specific to the simulation
        self.simulation_constraints = simul_constraints if simul_constraints else OrderedDict()
        self.method = method
        self.shadow_prices = shadow_prices

    def get_constraints(self):
        """
        :returns: All constraints applied during the simulation both persistent and simulation specific.
        """
        constraints = OrderedDict()
        constraints.update(self.envcond)
        constraints.update(self.model_constraints)
        constraints.update(self.simulation_constraints)
        return constraints

    def __repr__(self):
        return (f"objective: {self.objective_value}\nStatus: "
                f"{self.status}\nMethod:{self.method}")

    def __str__(self):
        return self.__repr__()
    
    def find(self, pattern=None, sort=False, shadow_prices=False, show_nulls=False):
        """Returns a dataframe of reactions and their fluxes matching a pattern or a list of patterns.

        :param pattern: a string or a list of strings. defaults to None
        :param sort: If the dataframe is to be sorted by flux rates. Defaults to False
        :type sort: bool, optional
        :return: returns a dataframe.
        :rtype: pandas.DataFrame
        """
        if shadow_prices:
            try:
                values = [(key, value) for key, value in self.shadow_prices.items()]
                columns = ['Metabolite', 'Shadow Price']
            except Exception:
                raise ValueError('No shadow prices')
        else:
            try:
                if show_nulls:
                    values = [(key, value) for key, value in self.fluxes.items()]
                else:
                    values = [(key, value) for key, value in self.fluxes.items() if value!=0.0]
                columns = ['Reaction ID', 'Flux rate']
            except Exception:
                raise ValueError('No fluxes')
        if pattern:
            import re
            if isinstance(pattern, list):
                patt = '|'.join(pattern)
                re_expr = re.compile(patt)
            else:
                re_expr = re.compile(pattern)
            values = [x for x in values if re_expr.search(x[0]) is not None]
        if sort:
            values.sort(key=lambda x: x[1])
        import pandas as pd
        df = pd.DataFrame(values, columns=columns)
        df = df.set_index(columns[0])
        return df

    @property
    def dataframe(self):
        return self.find()

    def get_net_conversion(self, biomassId=None):
        """Returns a string representation of the net conversion.

        :params str biosmassId: Biomass identifier (optional)
        :returns: A string representation of the net conversion.

        """

        left = ""
        right = ""
        firstLeft, firstRight = True, True
        from . import get_simulator
        sim = get_simulator(self.model)
        ssFluxes = self.fluxes
        for r_id in sim.reactions:
            fluxValue = ssFluxes[r_id]
            sub = list(sim.get_substrates(r_id).keys())
            prod = list(sim.get_products(r_id).keys())
            # if rId is a drain reaction
            if fluxValue != 0.0 and (len(sub) == 0 or len(prod) == 0):
                m = sub + prod
                if fluxValue < 0:
                    if firstLeft:
                        firstLeft = False
                    else:
                        left = left + " + "
                    left = left + str(-1 * fluxValue)
                    left = left + " " + m[0]
                else:
                    if firstRight:
                        firstRight = False
                    else:
                        right = right + " + "
                    right = right + str(fluxValue)
                    right = right + " " + m[0]

        if biomassId and biomassId in ssFluxes.keys():
            biomassFlux = ssFluxes[biomassId]
            if biomassFlux > 0:
                if firstRight:
                    firstRight = False
                else:
                    right = right + " + "
                right = right + str(biomassFlux)
                right = right + " " + biomassId

        return left + " --> " + right

    def get_metabolites_turnover(self, pattern=None, format='df'):
        """ Calculate metabolite turnovers.

         :param str format: the display format (pandas.DataFrame or dict).
             Default 'df' pandas.DataFrame

        Returns:
             dict or pandas.DataFrame: metabolite turnover rates
        """
        from . import get_simulator
        sim = get_simulator(self.model)

        if not self.fluxes:
            return None

        mets = None
        if pattern:
            import re
            if isinstance(pattern, list):
                patt = '|'.join(pattern)
                re_expr = re.compile(patt)
            else:
                re_expr = re.compile(pattern)
            mets = [x for x in sim.metabolites if re_expr.search(x) is not None]

        m_r_table = sim.metabolite_reaction_lookup()

        data = {m_id: 0.5*sum([abs(coeff * self.fluxes[r_id]) for r_id, coeff in neighbours.items()])
                for m_id, neighbours in m_r_table.items()}
        if mets is not None:
            data = {k: v for k, v in data.items() if k in mets}

        if format == 'df':
            import pandas as pd
            df = pd.DataFrame(list(data.values()), index=list(data.keys()), columns=['Turnover'])
            df.index.name = 'Metabolite'
            return df
        return data

    def get_metabolite(self, met_id, format='df'):
        """Displays the consumption/production of a metabolite in the reactions
        it participates.

        :param str met_id: the metabolite identifier.
        :param str format: the display format (pandas.DataFrame or dict).
             Default 'df' pandas.DataFrame

        Returns:
            dict or pandas.DataFrame: metabolite turnover rates
        """
        from . import get_simulator
        sim = get_simulator(self.model)

        if not self.fluxes:
            return None

        m_r = sim.metabolite_reaction_lookup()[met_id]
        data = {r_id: coeff * self.fluxes[r_id] for r_id, coeff in m_r.items()}
        if format == 'df':
            import pandas as pd
            df = pd.DataFrame(list(data.values()), index=list(data.keys()), columns=['Value'])
            df.index.name = 'Reaction'
            return df
        return data

    @classmethod
    def from_linear_solver(cls, solution):
        """Converts a solver Solution object to a SolutionResult object.

        :param solution: solution to be converted
        :type solution: mewpy.solver.solution.Solution
        :raises ValueError: if solution is not an instance of mewpy.solver.solution.Solution
        :return: an instance of SolutionResult
        :rtype: SolutionResult
        """
        from mewpy.solvers import Solution, Status
        smap = {Status.OPTIMAL: SStatus.OPTIMAL,
                Status.UNKNOWN: SStatus.UNKNOWN,
                Status.SUBOPTIMAL: SStatus.SUBOPTIMAL,
                Status.UNBOUNDED: SStatus.UNBOUNDED,
                Status.INFEASIBLE: SStatus.INFEASIBLE,
                Status.INF_OR_UNB: SStatus.INF_OR_UNB
                }
        if not isinstance(solution, Solution):
            raise ValueError('solution should be and instance of mewpy.solvers.solution.Solution')
        return cls(None, solution.fobj, fluxes=solution.values, status=smap[solution.status])


def simulate(model, envcond=None, objective=None, method=SimulationMethod.FBA, maximize=True, constraints=None, reference=None,
             solver=None, **kwargs):
    """Runs an FBA phenotype simulation.

    :param model: cobrapy, reframed, GERM constraint-base model
    :param (dict ) envcond : Environmental conditions, defaults to None
    :param (dict) objective: the FBA objective , defaults to None
    :param method: The FBA method, defaults to SimulationMethod.FBA
    :param maximize: optimization sense, defaults to True
    :param (dict) constraints: additional constraints , defaults to None
    :param (dict) reference: reference fluxes, defaults to None
    :param solver: a solver instance, defaults to None

    :returns: SimultationResult
    """
    from . import get_simulator
    sim = get_simulator(model, envcond=envcond)
    res = sim.simulate(objective=objective, method=method, maximize=maximize,
                       constraints=constraints, reference=reference,
                       solver=solver, **kwargs)
    return res
