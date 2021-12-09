from abc import abstractmethod, ABC
from collections import OrderedDict

from mewpy.simulation import SimulationMethod


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

    def __evaluator__(self, kwargs, candidate):
        res = self.simulate(constraints=candidate, **kwargs)
        return res

    def simulate_mp(self, objective=None, method=SimulationMethod.FBA, maximize=True, constraints_list=None,
                    reference=None,
                    solver=None, n_mp=None, **kwargs):
        from mewpy.util.process import get_fevaluator
        args = {}
        args['objective'] = objective
        args['method'] = method
        args['maximize'] = maximize
        args['reference'] = reference
        args.update(kwargs)
        from functools import partial
        func = partial(self.__evaluator__, args)

        mp_evaluator = get_fevaluator(func)
        res = mp_evaluator.evaluate(constraints_list, None)
        return res

    @abstractmethod
    def get_reaction_bounds(self, r_id):
        raise NotImplementedError

    @abstractmethod
    def metabolite_reaction_lookup(self, force_recalculate=False):
        raise NotImplementedError


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
