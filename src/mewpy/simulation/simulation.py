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
        Returns a list of reaction identifiers
        """
        raise NotImplementedError

    @property
    def genes(self):
        """
        Returns a list of gene identifiers
        """
        raise NotImplementedError

    @property
    def proteins(self):
        """
        Returns a list of proteins identifiers
        """
        raise NotImplementedError

    @property
    def metabolites(self):
        """
        Returns a list of metabolite identifiers
        """
        raise NotImplementedError

    @property
    def medium(self):
        raise NotImplementedError

    def get_drains(self):
        return NotImplementedError

    def get_gpr(self, reaction_id):
        """
        Returns a reaction gpr if exists None otherwise.
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
    def simulate(self, objective=None, method=SimulationMethod.FBA, maximize=True, constraints=None, reference=None, solver=None, **kwargs):
        """
        Returns a SimulationResult
        """
        raise NotImplementedError

    @abstractclassmethod
    def FVA(self, obj_frac=0, reactions=None, constraints=None, loopless=False, internal=None, solver=None):
        """ Run Flux Variability Analysis (FVA).
        """
        raise NotImplementedError


class SimulationResult(object):

    def __init__(self, model, objective_value, fluxes=None, status=None, envcond=None, model_constraints=None, simul_constraints=None, maximize=True):
        """
            class that represents simulation results and performs
            operations over them.
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
        Return all constraints applyed during the simulation both persistent and simulation specific
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
        '''
           Returs a string representation of the net conversion

           args:
                biosmassId (str) : optional  
        '''

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
