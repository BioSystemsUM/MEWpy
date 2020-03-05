from enum import Enum


# maps a model to its simulator
# entries take the form:  full_model_class_path -> (simulator_path, simulator_class_name)

    
map_model_simulator ={
    'cobra.core.model.Model'        : ('mewpy.simulation.cobra',    'Simulation'),
    'reframed.core.cbmodel.CBModel' : ('mewpy.simulation.reframed', 'Simulation'),
    'geckopy.gecko.GeckoModel'      : ('mewpy.simulation.cobra',    'GeckoSimulation'),
    'mewpy.model.gecko.GeckoModel'  : ('mewpy.simulation.reframed', 'GeckoSimulation')
}


def get_simulator(model, envcond = None , constraints = None , reference= None):
    """
    Returns a simulator instance for the model
    The simulator instance is model dependent.
    
    This function is invoked by a EA optimization Problem and by evaluation function instances.
    
    arguments:
        model : the model
        envcond (dic) : a dictionary of environmental conditions
        contrainsts (dic) : a dictionary of additional persistent constraints
    """
    name = f"{model.__class__.__module__}.{model.__class__.__name__}"
    try:
        module_name, class_name = map_model_simulator[name]
    except Exception:
        raise ValueError(f"The model [{name}] has no defined simulator.")
    module = __import__(module_name,fromlist=[None])
    class_ = getattr(module, class_name)
    instance = class_(model, envcond = envcond , constraints = constraints, reference= reference)
    return instance



def get_container(model):
    """
    Returns a container for a given model sharing a common interface
    
    :param model 
    """
    try:
        from reframed.core.cbmodel import CBModel
        if isinstance(model, CBModel):
            from mewpy.simulation.reframed import CBModelContainer
            return CBModelContainer(model)
    except:
        pass
    
    try:
        from cobra.core.model import Model
        if isinstance(model, Model):
            from mewpy.simulation.cobra import CobraModelContainer
            return CobraModelContainer(model)
    except:
        pass
    
    raise ValueError("Unrecognized model")


class SimulationMethod(Enum):
    FBA   = 'FBA'
    pFBA  = 'pFBA'
    MOMA  = 'MOMA'
    lMOMA = 'lMOMA'
    ROOM  = 'ROOM'
    NONE  = 'NONE'


class SStatus(Enum):
    """ Enumeration of possible solution status. """
    OPTIMAL = 'Optimal'
    UNKNOWN = 'Unknown'
    SUBOPTIMAL = 'Suboptimal'
    UNBOUNDED = 'Unbounded'
    INFEASIBLE = 'Infeasible'
    INF_OR_UNB = 'Infeasible or Unbounded'
