from abc import abstractmethod, ABC
from collections import OrderedDict
from ..util.parsing import evaluate_expression_tree
from ..util.process import cpu_count
from . import SimulationMethod, SStatus
from tqdm import tqdm


class ModelContainer(ABC):
    """Interface for Model container.
       Provides an abstraction from models implementations.
    """

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
    def medium(self):
        raise NotImplementedError

    def get_exchange_reactions(self):
        return NotImplementedError

    def get_gpr(self, reaction_id):
        """
        :returns: A string representation of the reaction GPR if exists None otherwise.
        """
        raise NotImplementedError

    def summary(self):
        print(f"Metabolites: {len(self.metabolites)}")
        print(f"Reactions: {len(self.reactions)}")
        print(f"Genes: {len(self.genes)}")

    def set_objective(self, reaction):
        raise NotImplementedError


class Simulator(ModelContainer):
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
    def FVA(self, obj_frac=0, reactions=None, constraints=None, loopless=False, internal=None, solver=None):
        """ Abstract method to run Flux Variability Analysis (FVA).

        :returns: A dictionary of flux range values.

        """
        raise NotImplementedError

    def simulate_mp(self, objective=None, method=SimulationMethod.FBA, maximize=True, constraints_list=None,
                    reference=None, solver=None, n_mp=None, **kwargs):

        res = Parallel(n_jobs=cpu_count())(delayed(self.simulate)(objective, method, maximize, constraints, reference,
                                                                  solver, **kwargs) for constraints in constraints_list)
        return res

    @abstractmethod
    def get_reaction_bounds(self, r_id):
        raise NotImplementedError

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
        if find_in == 'g':
            data = [{'Gene': x} for x in values]
        else:
            data = [self.get_reaction(x) for x in values]
        df = pd.DataFrame(data)
        return df

    def essential_reactions(self, min_growth=0.01):
        """Essential reactions are those when knocked out enable a biomass flux value above a minimal growth defined as
        a percentage of the wild type growth.

        :param float min_growth: Minimal percentage of the wild type growth value. Default 0.01 (1%).
        :returns: A list of essential reactions.

        """
        essential = getattr(self, '_essential_reactions', None)
        if essential is not None:
            return essential
        wt_solution = self.simulate()
        wt_growth = wt_solution.objective_value
        self._essential_reactions = []
        print("Computing essential reactions:")
        for rxn in tqdm(self.reactions):
            res = self.simulate(constraints={rxn: 0},   slim=True)
            if not res or res < min_growth:
                self._essential_reactions.append(rxn)
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

    def essential_genes(self, min_growth=0.01):
        """Essential genes are those when deleted enable a biomass flux value above a minimal growth defined as
        a percentage of the wild type growth.

        :param float min_growth: Minimal percentage of the wild type growth value. Default 0.01 (1%).
        :returns: A list of essential genes.

        """
        essential = getattr(self, '_essential_genes', None)
        if essential is not None:
            return essential
        self._essential_genes = []
        wt_solution = self.simulate()
        wt_growth = wt_solution.objective_value
        genes = self.genes
        print("Computing essential genes:")
        for gene in tqdm(genes):
            active_genes = set(self.genes) - {gene}
            active_reactions = self.evaluate_gprs(active_genes)
            inactive_reactions = set(self.reactions) - set(active_reactions)
            gr_constraints = {rxn: 0 for rxn in inactive_reactions}
            res = self.simulate(constraints=gr_constraints, slim=True)
            if not res or res < min_growth:
                self._essential_genes.append(gene)
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


class SimulationResult(object):
    """Class that represents simulation results and performs operations over them."""

    def __init__(self, model, objective_value, fluxes=None, status=None, envcond=None, model_constraints=None,
                 simul_constraints=None, maximize=True, method=None):
        """
        :param model: A model instance.
        :param objective_value: The phenotype simulation objective value.
        :param dict fluxes: A dictionary of reaction fluxes values.
        :param status: The LP status.
        :param dict envcond: The environmental conditions of the phenotype simulation.
        :param dict model_constraints: Possible persistent additional constraints.
        :param dict simul_constraints: The simulation constraints.
        :param boolean maximize: Optimization direction.

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
                f"{self.status}\nConstraints: {self.get_constraints()}\nMethod:{self.method}")

    def __str__(self):
        return (f"objective: {self.objective_value}\nStatus: "
                f"{self.status}\nConstraints: {self.get_constraints()}\nMethod:{self.method}")

    def find(self, pattern=None, sort=False):
        """Returns a dataframe of reactions and their fluxes matching a pattern or a list of patterns.

        :param pattern: a string or a list of strings. defaults to None
        :param sort: If the dataframe is to be sorted by flux rates. Defaults to False
        :type sort: bool, optional
        :return: returns a dataframe.
        :rtype: pandas.DataFrame
        """
        values = [(key, value) for key, value in self.fluxes.items()]
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
        df = pd.DataFrame(values, columns=['Reaction ID', 'Flux rate'])
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

    def get_metabolites_turnover(self):
        """ Calculate metabolite turnover.

        Arguments:
            model: REFRAMED/Cobrapy model or Simulator that generated the solution

        Returns:
            dict: metabolite turnover rates
        """
        from . import get_simulator
        sim = get_simulator(self.model)

        if not self.values:
            return None

        m_r_table = sim.metabolite_reaction_lookup()
        t = {m_id: 0.5*sum([abs(coeff * self.fluxes[r_id]) for r_id, coeff in neighbours.items()])
             for m_id, neighbours in m_r_table.items()}
        return t
