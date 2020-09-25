from enum import Enum

# Model specific simulators mapping:
# Entries take the form:  full_model_class_path -> (simulator_path, simulator_class_name)
# TODO: use qualified names

map_model_simulator = {
    'geckopy.gecko.GeckoModel': ('mewpy.simulation.cobra', 'GeckoSimulation'),
    'mewpy.model.gecko.GeckoModel': ('mewpy.simulation.reframed', 'GeckoSimulation'),
    'mewpy.model.smoment.SMomentModel': ('mewpy.simulation.reframed', 'GeckoSimulation')
}


def get_simulator(model, envcond=None, constraints=None, reference=None):
    """
    Returns a simulator instance for the model.
    The simulator instance is model dependent.
    Besides able to be used on its own, this function is invoked by EA optimization problems 
    and by evaluation function instances to perform phenotyoe evaluations of candidate solutions.

    
    :param model: the model
    :param dict envcond: A dictionary of environmental conditions.
    :param dict contrainsts: A dictionary of additional persistent constraints.
    :returns: An instance of Simulator
    """
    instance = None
    name = f"{model.__class__.__module__}.{model.__class__.__name__}"
    if name in map_model_simulator:
        module_name, class_name = map_model_simulator[name]
        module = __import__(module_name, fromlist=[None])
        class_ = getattr(module, class_name)
        instance = class_(model, envcond=envcond,
                          constraints=constraints, reference=reference)
    else:
        try:
            from cobra.core.model import Model
            if isinstance(model, Model):
                from .cobra import Simulation
                instance = Simulation(
                    model, envcond=envcond, constraints=constraints, reference=reference)
        except ImportError:
            pass
        if not instance:
            try:
                from reframed.core.cbmodel import CBModel
                if isinstance(model, CBModel):
                    from .reframed import Simulation
                    instance = Simulation(
                        model, envcond=envcond, constraints=constraints, reference=reference)
            except ImportError:
                pass
    if not instance:
        raise ValueError(f"The model <{name}> has no defined simulator.")
    return instance


def get_container(model):
    """
    Returns a container for a given model sharing a common interface.
    A container does not perform any task, it only serves as a basic interface with a phenotype simulator.

    :param model: An instance of a metabolic model.
    :returns: A container.
     
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

    raise ValueError(f"Unrecognized model class: {model.__class__.name}")


class SimulationMethod(Enum):
    
    FBA   = 'FBA'
    pFBA  = 'pFBA'
    MOMA  = 'MOMA'
    lMOMA = 'lMOMA'
    ROOM  = 'ROOM'
    NONE  = 'NONE'
    
    def __eq__(self, other):
        """Overrides equal to enable string name comparison.
        Allows to seamlessly use: 
            SimulationMethod.FBA = SimulationMethod.FBA
            SimulationMethod.FBA = 'FBA'
        without requiring an additional level of comparison (SimulationMethod.FBA.name = 'FBA')
        """
        if isinstance(other,SimulationMethod):
            return super().__eq__(other)
        elif isinstance(other,str):
            return self.name == other
        else:
            return False

    def __hash__(self):
        return hash(self.name)


class SStatus(Enum):
    """ Enumeration of possible solution status. """
    OPTIMAL    = 'Optimal'
    UNKNOWN    = 'Unknown'
    SUBOPTIMAL = 'Suboptimal'
    UNBOUNDED  = 'Unbounded'
    INFEASIBLE = 'Infeasible'
    INF_OR_UNB = 'Infeasible or Unbounded'
