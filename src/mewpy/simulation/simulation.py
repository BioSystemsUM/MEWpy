from abc import abstractclassmethod
from collections import OrderedDict
from mewpy.simulation import SimulationMethod


class ModelContainer:
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

    def get_drains(self):
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

    @abstractclassmethod
    def simulate(self, objective=None, method=SimulationMethod.FBA, maximize=True, constraints=None, reference=None,
                 solver=None, **kwargs):
        """Abstract method to run a phenotype simulation.

        :returns: A SimulationResult.

        """
        raise NotImplementedError

    @abstractclassmethod
    def FVA(self, obj_frac=0, reactions=None, constraints=None, loopless=False, internal=None, solver=None):
        """ Abstract method to run Flux Variability Analysis (FVA).

        :returns: A dictionary of flux range values.

        """
        raise NotImplementedError


class SimulationResult(object):
    """Class that represents simulation results and performs operations over them."""

    def __init__(self, model, objective_value, fluxes=None, status=None, envcond=None, model_constraints=None,
                 simul_constraints=None, maximize=True):
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
        return "objective: {}\nStatus: {}".format(self.objective_value, self.status)

    def __str__(self):
        return "objective: {}\nStatus: {}".format(self.objective_value, self.status)

    def get_net_conversion(self, biomassId=None):
        """Returns a string representation of the net conversion.

        :params str biosmassId: Biomass identifier (optional)
        :returns: A string representation of the net conversion.

        """

        left = ""
        right = ""
        firstLeft, firstRight = True, True

        ssFluxes = self.fluxes
        reactions = self.model.reactions
        for r_id in reactions.keys():
            fluxValue = ssFluxes[r_id]
            sub = reactions[r_id].get_substrates()
            prod = reactions[r_id].get_products()
            # if rId is a drain reaction
            if fluxValue != 0.0 and (len(sub) == 0 or len(prod) == 0):
                m = sub + prod
                if fluxValue < 0:
                    if firstLeft:
                        firstLeft = False
                    else:
                        left = left+" + "
                    left = left+str(-1*fluxValue)
                    left = left+" "+m[0]
                else:
                    if firstRight:
                        firstRight = False
                    else:
                        right = right+" + "
                    right = right+str(fluxValue)
                    right = right+" "+m[0]

        if biomassId and biomassId in ssFluxes.keys():
            biomassFlux = ssFluxes[biomassId]
            if biomassFlux > 0:
                if firstRight:
                    firstRight = False
                else:
                    right = right+" + "
                right = right + str(biomassFlux)
                right = right + " " + biomassId

        return left + " --> " + right
