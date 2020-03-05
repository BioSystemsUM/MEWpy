from abc import ABC,abstractmethod
from reframed.core.cbmodel import CBModel


#### Depricated
class Model(ABC):
    """ Defines the interface required for a model
        to be compatible with mewpy
    """

    @abstractmethod
    def get_reactions(self):
        raise NotImplementedError

    @abstractmethod
    def get_compartemens(self):
        raise NotImplementedError

    @abstractmethod
    def get_metabolites(self):
        raise NotImplementedError

    @abstractmethod
    def get_genes(self):
        raise NotImplementedError



class SModel(CBModel,Model):

    def get_reactions(self):
        return self.reactions

    def get_compartmentes(self):
        return self.compartments
    
    def get_metabolites(self):
        return self.metabolites

    def get_genes(self):
        return self.genes    