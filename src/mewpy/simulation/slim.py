import logging
from typing import Union, Dict, Tuple, List, TYPE_CHECKING

import numpy as np
from slimpy.core import Model, MetabolicModel, RegulatoryModel
from slimpy.variables import Reaction

from scr.slimpy.analysis import FBA
from scr.slimpy.solution import ModelSolution
from . import SimulationMethod, SStatus
from .simulation import Simulator, SimulationResult, ModelContainer
from ..util.constants import ModelConstants
from ..util.utilities import Dispatcher

if TYPE_CHECKING:
    from ..solvers.solver import Solver
    from ..solvers import CplexSolver, GurobiSolver, OptLangSolver

LOGGER = logging.getLogger(__name__)


class SlimModelContainer(ModelContainer):

    def __init__(self, model: Union[Model, MetabolicModel, RegulatoryModel]):

        """
        A slim model container. It provides an interface to access slim model containers

        :type model: Union[Model, MetabolicModel, RegulatoryModel]
        :param model: A slim Model, MetabolicModel, RegulatoryModel or all together
        """

        if not isinstance(model, (Model, MetabolicModel, RegulatoryModel)):
            raise ValueError("The model is not an instance of slim.core.Model")

        self.model = model

    # -----------------------------------------------------------------------------
    # Metabolic static attributes
    # -----------------------------------------------------------------------------

    @property
    def reactions(self):

        if self.model.is_metabolic():
            return list(self.model.reactions.keys())

        return []

    @property
    def metabolites(self):

        if self.model.is_metabolic():
            return list(self.model.metabolites.keys())

        return []

    @property
    def genes(self):

        if self.model.is_metabolic():
            return list(self.model.genes.keys())

        return []

    @property
    def proteins(self):
        return []

    # -----------------------------------------------------------------------------
    # Metabolic dynamic attributes
    # -----------------------------------------------------------------------------

    @property
    def compartments(self):

        if self.model.is_metabolic():
            return self.model.compartments

        return {}

    def get_drains(self):

        if self.model.is_metabolic():
            return list(self.model.exchanges.keys())

        return []

    def get_gpr(self, reaction_id):

        reaction = self.model.get(reaction_id)

        if reaction:

            if not reaction.gpr.is_none:
                return reaction.gpr.to_string()

        return

    # -----------------------------------------------------------------------------
    # Metabolic static attributes setters
    # -----------------------------------------------------------------------------

    def set_objective(self, reaction):

        self.model.objective = reaction

    # -----------------------------------------------------------------------------
    # Regulatory static attributes
    # -----------------------------------------------------------------------------
    @property
    def interactions(self):

        if self.model.is_regulatory():
            return list(self.model.interactions.keys())

        return []

    @property
    def regulators(self):

        if self.model.is_regulatory():
            return list(self.model.regulators.keys())

        return []

    @property
    def targets(self):

        if self.model.is_regulatory():
            return list(self.model.targets.keys())

        return []

    @property
    def environmental_stimuli(self):

        if self.model.is_regulatory():
            return list(self.model.environmental_stimuli.keys())

        return []

    # -----------------------------------------------------------------------------
    # Model summary
    # -----------------------------------------------------------------------------

    def summary(self):

        if self.model.is_metabolic():
            print(f"Metabolites: {len(self.metabolites)}")
            print(f"Reactions: {len(self.reactions)}")
            print(f"Genes: {len(self.genes)}")

        if self.model.is_regulatory():
            print(f"Interactions: {len(self.interactions)}")
            print(f"Regulators: {len(self.regulators)}")
            print(f"Targets: {len(self.targets)}")


class Simulation(SlimModelContainer, Simulator):
    dispatcher = Dispatcher()

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 objective: str = None,
                 envcond: Dict[str, Tuple[Union[int, float], Union[int, float]]] = None,
                 constraints: Dict[str, Tuple[Union[int, float], Union[int, float]]] = None,
                 reference: Dict[str, Union[int, float]] = None,
                 reset_solver=ModelConstants.RESET_SOLVER):

        """
        Simulation supports simulation of a slim Model, MetabolicModel, RegulatoryModel or all.
        Additional environmental conditions and constraints can be set using this interface.
        The objective function can be altered too.
        A reference wild type reaction flux values can be provided for convenience.

        Simulation offers several methods to retrieve additional information from the model.

        :type model: Union[Model, MetabolicModel, RegulatoryModel]
        :param model: A slim Model, MetabolicModel, RegulatoryModel or all together

        Optional:

        :param objective: the model objective
        :param envcond: a dictionary of additional environmental conditions
        :param constraints: a dictionary of additional constraints
        :param reference: A dictionary of the wild type flux values
        :param reset_solver: Whether it is required to reset the solver or not after simulation.
        By default: ModelConstants.RESET_SOLVER
        """

        super().__init__(model)

        if not isinstance(model, (Model, MetabolicModel, RegulatoryModel)):
            raise ValueError("The model is not an instance of slim.core.Model")

        if not envcond:
            envcond = {}

        if not constraints:
            constraints = {}

        self.model = model

        self.environmental_conditions = envcond
        self.constraints = constraints

        self._essential_reactions = None
        self._essential_genes = None
        self._reference = reference
        self._reset_solver = reset_solver

        self._m_r_lookup = None

        self.__status_mapping = {
            'optimal': SStatus.OPTIMAL,
            'unbounded': SStatus.UNBOUNDED,
            'infeasible': SStatus.INFEASIBLE,
            'infeasible_or_unbounded': SStatus.INF_OR_UNB,
            'suboptimal': SStatus.SUBOPTIMAL,
            'unknown': SStatus.UNKNOWN
        }

    # -----------------------------------------------------------------------------
    # Metabolic interface
    # -----------------------------------------------------------------------------

    @property
    def objective(self) -> Dict[str, Union[float, int]]:
        return {var.id: coef for var, coef in self.model.objective.items()}

    @objective.setter
    def objective(self, value: Dict[str, Union[float, int]]):
        self.model.objective = value

    @property
    def reference(self) -> Dict[str, Union[float, int]]:
        """
        The reference wild type reaction flux values.

        :return: A dictionary of wild type reaction flux values.

        """
        if self._reference is None:
            self._reference = self.simulate(method=SimulationMethod.pFBA).fluxes
        return self._reference

    @property
    def essential_reactions(self) -> List[str]:

        """
        The set of knocked out reactions that impair a minimal growth rate predicted by the model with fba.
        The defined minimal percentage of the wild type growth rate is usually set to 0.01 (1%).

        :return: A list of essential reactions.

        """

        if self._essential_reactions is not None:
            return self._essential_reactions

        if not self.model.is_metabolic():
            return []

        wt_solution = self.simulate(method=SimulationMethod.FBA)
        wt_growth = wt_solution.objective_value

        self._essential_reactions = []

        reactions = self.reactions.copy()

        for rxn in reactions:

            res = self.simulate(constraints={rxn: 0})

            if res:

                if (res.status == SStatus.OPTIMAL and res.objective_value < wt_growth * 0.01) \
                        or res.status == SStatus.INFEASIBLE:
                    self._essential_reactions.append(rxn)

        return self._essential_reactions

    @property
    def essential_genes(self) -> List[str]:

        """
        The set of knocked out genes that impair a minimal growth rate predicted by the model with fba.
        The defined minimal percentage of the wild type growth rate is usually set to 0.01 (1%).

        :return: A list of essential genes.

        """

        if self._essential_genes is not None:
            return self._essential_genes

        if not self.model.is_metabolic():
            return []

        self._essential_genes = []

        wt_solution = self.simulate()
        wt_growth = wt_solution.objective_value

        values = {gene.id: 1.0 for gene in self.model.yield_genes()}

        for gene in self.model.yield_genes():

            values[gene.id] = 0.0

            constraints = {}

            for rxn in gene.yield_reactions():

                if not rxn.gpr.is_none:

                    res = rxn.gpr.evaluate(values=values)

                    if not res:
                        constraints[rxn.id] = 0

            res = self.simulate(constraints=constraints)

            values[gene.id] = 1.0

            if res:

                if (res.status == SStatus.OPTIMAL and res.objective_value < wt_growth * 0.01) \
                        or res.status == SStatus.INFEASIBLE:
                    self._essential_genes.append(gene)

        return self._essential_genes

    def evaluate_gprs(self, active_genes) -> List[str]:

        """
        Returns the list of active reactions for a given list of active genes.

        :param list active_genes: list of genes identifiers.
        :return: a list of active reaction identifiers.

        """

        if not self.model.is_metabolic():
            return []

        values = {gene.id: 1.0 if gene.id in active_genes else 0.0
                  for gene in self.model.yield_genes()}

        return [rxn.id for rxn in self.model.yield_reactions()
                if rxn.gpr.is_none or rxn.gpr.evaluate(values=values)]

    def add_reaction(self, reaction: Reaction, replace: bool = True, comprehensive=True):

        """
        Adds a reaction to the slim model

        :param reaction: Reaction object to be added
        :param replace: Whether to replace the reaction in the model if already exists
        :param comprehensive: Whether to add metabolites and genes associated with the reaction to the model

        :return:
        """

        if not self.model.is_metabolic():
            return

        if replace:

            self.model.add(reaction, 'reaction', comprehensive=comprehensive, history=False)

        else:

            if reaction.id not in self.model.reactions:
                self.model.add(reaction, 'reaction', comprehensive=comprehensive, history=False)

    def remove_reaction(self, reaction: Reaction, remove_orphans=True):

        """
        Removes a reaction from the slim model

        :param reaction: Reaction object to be removed
        :param remove_orphans: Whether to remove metabolites and genes orphans associated with the reaction in the model

        :return:
        """

        if not self.model.is_metabolic():
            return

        self.model.remove(reaction, 'reaction', remove_orphans=remove_orphans, history=False)

    def get_uptake_reactions(self) -> List[str]:
        """

        Get uptake reactions

        :return: The list of uptake reactions.

        """

        if not self.model.is_metabolic():
            return []

        return [rxn.id for rxn in self.model.yield_exchanges()
                if rxn.reversibility or rxn.lower_bound < 0]

    def get_transport_reactions(self) -> List[str]:
        """

        Get the transport reactions

        :return: The list of transport reactions.

        """

        ext = self.model.external_compartment

        return [rxn.id for rxn in self.model.yield_reactions()
                if ext in rxn.compartments and len(rxn.compartments) > 1]

    def get_transport_genes(self) -> List[str]:
        """

        Get genes associated with transport reactions

        :return: the list genes associated with transport reactions.

        """
        ext = self.model.external_compartment

        genes = set()

        for rxn in self.model.yield_reactions():

            if ext in rxn.compartments and len(rxn.compartments) > 1:
                genes.update(rxn.genes.keys())

        return list(genes)

    # TODO: The models provided to mewpy should always be reversible?
    def reverse_reaction(self, reaction_id: str) -> Union[str, None]:

        """
        Identify if a reaction is reversible and returns the reverse reaction if it is the case.

        :param reaction_id: a reaction identifier
        :return: a reverse reaction identifier or None

        """

        reaction = self.model.get(reaction_id)

        if reaction:

            if reaction.reversibility:
                return reaction_id

        return

    def metabolite_reaction_lookup(self, force_recalculate: bool = False) -> Dict[str, Dict[str, Union[int, float]]]:
        """
        A view of the S matrix for the model.
        That is, it provides the network topology as a nested map from metabolite to reaction to coefficient.

        :return: a dictionary lookup table
        """

        if not self._m_r_lookup or force_recalculate:

            metabolites = self.metabolites.copy()

            self._m_r_lookup = {met: {} for met in metabolites}

            for reaction in self.model.yield_reactions():

                for met, coef in reaction.stoichiometry.items():
                    self._m_r_lookup[met.id][reaction.id] = coef

        return self._m_r_lookup

    def get_reaction_bounds(self, reaction: str) -> Tuple[Union[int, float], Union[int, float]]:

        """
        Get the bounds for a given reaction.

        :param reaction: reaction identifier

        :return: a tuple with the reaction bounds
        """

        if reaction in self.constraints:
            lb, ub = self.constraints[reaction]

        elif reaction in self.environmental_conditions:
            lb, ub = self.environmental_conditions[reaction]

        else:
            lb, ub = self.model.get(reaction).bounds

        return lb if lb > -np.inf else -999999, ub if ub < np.inf else 999999

    def find_bounds(self) -> Tuple[Union[int, float], Union[int, float]]:
        """

        It finds the median upper and lower bound of the metabolic model.
        Bounds can vary from model to model. slim set reactions bounds to (-999999.0, 999999.0) by default.

        :return: a tuple of median upper and lower bounds
        """

        lbs, ubs = [], []

        for rxn in self.model.yield_reactions():
            lbs.append(rxn.lower_bound)
            ubs.append(rxn.upper_bound)

        lbs = np.asarray(lbs, dtype=float)
        ubs = np.asarray(ubs, dtype=float)

        lb = np.nanmedian(lbs[lbs != 0.0])
        ub = np.nanmedian(ubs[ubs != 0.0])

        if np.isnan(lb):
            LOGGER.warning("Could not identify a median lower bound. Setting the default")
            lb = -1000.0

        if np.isnan(ub):
            LOGGER.warning("Could not identify a median upper bound. Setting the default")
            ub = 1000.0

        return lb, ub

    def find_unconstrained_reactions(self) -> List[str]:

        """
        Finds a list of reactions that are not constrained in the model

        :return: list of unconstrained reactions
        """

        lb, ub = self.find_bounds()

        return [rxn for rxn in self.model.yield_reactions()
                if rxn.lower_bound <= lb and rxn.upper_bound >= ub]

    # -----------------------------------------------------------------------------
    # Simulation
    # -----------------------------------------------------------------------------

    @dispatcher.register(SimulationMethod.FBA)
    def _fba(self, model):

        return FBA(model=model).optimize()

    @dispatcher.register(SimulationMethod.pFBA)
    def _pfba(self, model):

        return FBA(model=model).optimize()

    @dispatcher.register(SimulationMethod.NONE)
    def _none(self, model):

        return ModelSolution(method=None, x=None, objective_value=None, status='unknown')

    def simulate(self,
                 objective: Dict[str, Union[int, float]] = None,
                 method: SimulationMethod = SimulationMethod.FBA,
                 maximize: bool = True,
                 constraints: Dict[str, Tuple[Union[int, float], Union[int, float]]] = None,
                 reference: Dict[str, Union[int, float]] = None,
                 scalefactor: float = None,
                 solver: Union['Solver', 'CplexSolver', 'GurobiSolver', 'OptLangSolver'] = None) -> SimulationResult:

        """

        Simulates a phenotype for a given objective and set of constraints using the specified method.
        Reference wild-type conditions are also accepted

        :param objective: The simulation objective. If none, the model objective is considered.
        :param method: The SimulationMethod (FBA, pFBA, lMOMA, etc ...).
        See available methods at mewpy.simulation.SimulationMethod
        :param maximize: The optimization direction
        :param constraints: A dictionary of constraints to be applied to the model
        :param reference: A dictionary of reaction flux values for lmoma, moma and room
        :param scalefactor: A positive scaling factor for the solver. Default None.
        :param solver: A mewpy solver instance

        :return: SimulationResult object that holds the model, objective, solution values, status,
        environmental conditions and further constraints
        """

        if not objective:
            objective = self.objective

        if not constraints:
            constraints = {}

        simulation_constraints = {**self.constraints, **self.environmental_conditions, **constraints}

        with self.model:

            if maximize:

                self.model.objective = objective

            else:
                self.model.objective = {key: -coef for key, coef in objective.items()}

            for rxn, cnt in simulation_constraints.items():

                if isinstance(cnt, tuple):

                    self.model.get(rxn).bounds = cnt

                else:

                    self.model.get(rxn).bounds = (cnt, cnt)

            solution: ModelSolution = self.dispatcher(method, self.model)

        status = self.__status_mapping[solution.status]

        return SimulationResult(model=self.model,
                                objective_value=solution.objective_value,
                                fluxes=solution.x,
                                status=status,
                                envcond=self.environmental_conditions,
                                model_constraints=self.constraints,
                                simul_constraints=constraints,
                                maximize=maximize,
                                method=method)

    def FVA(self,
            obj_frac: float = 0.9,
            reactions=None,
            constraints=None,
            loopless=False,
            internal=None,
            solver=None,
            format='dict'):

        """

        It performs a Flux Variability Analysis (FVA).

        :param obj_frac: The minimum fraction of the maximum growth rate (default 0.9).
        Requires that the objective value is at least the fraction times maximum objective value.
        A value of 0.85 for instance means that the objective has to be at least at 85% percent of its maximum.
        :param reactions: List of reactions to analyze (default: all).
        :param constraints: Additional constraints (optional).
        :param loopless: Run looplessFBA internally (very slow) (default: false).
        :param internal: List of internal reactions for looplessFBA (optional).
        :param solver: A pre-instantiated solver instance (optional)
        :param format: The return format: 'dict' to return a dictionary; 'df' to return a data frame.

        :return: A dictionary or data frame of flux variation ranges.
        """

        if not constraints:
            constraints = {}

        simulation_constraints = {**self.constraints, **self.environmental_conditions, **constraints}

        with self.model:
            for rxn, cnt in simulation_constraints.items():

                if isinstance(cnt, tuple):

                    self.model.get(rxn).bounds = cnt

                else:

                    self.model.get(rxn).bounds = (cnt, cnt)

            from slimpy.analysis import FVA
            fva = FVA(self.model)

            solution = fva.run(initial_objective=1, minimize_initial=False)

        if format == 'df':

            import pandas as pd
            e = solution.items()
            f = [[a, b, c] for a, [b, c] in e]
            df = pd.DataFrame(f, columns=['Reaction ID', 'Minimum', 'Maximum'])

            return df

        return solution
