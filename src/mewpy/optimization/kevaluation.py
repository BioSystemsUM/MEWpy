
"""
Kinetic specific evaluation functions
"""

from .evaluation import KineticEvaluationFunction
from ..util.ode import ODEStatus



class KineticTargetFlux(KineticEvaluationFunction):
    """ Target Flux evaluation function. 
    The fitness value is the flux value of the identified reaction.

    
    :param reaction: (str) The reaction identifier whose flux value is to be used as fitness. Default None in which case the model objective is considered.
    
    """

    def __init__(self, reaction, maximize = True):
        super(KineticTargetFlux, self).__init__(maximize=maximize, worst_fitness=0.0)
        self.reaction = reaction
       
    def get_fitness(self, simul_results, candidate, **kwargs):
        if simul_results.status == ODEStatus.ERROR:
            return self.worst_fitness
        else:
            try :
                return simul_results.fluxes[self.reaction]
            except Exception as e:
                print(e)
                return self.worst_fitness

    def required_simulations(self):
        """
        If the evaluation function requires a pre-simulation to compute fitness values
        """
        return []

    def short_str(self):
        return "KTargetFlux"

    def method_str(self):
        return "Kinetic Target Flux"

