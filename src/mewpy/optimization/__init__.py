from mewpy.utils.constants import EAConstants
from .ea import Solution

engines = dict()

try:
    from .inspyred.ea import EA as InspyredEA
    engines['inspyred'] = InspyredEA
except ImportError:
    pass

try:
    from .jmetal.ea import EA as JMetalEA
    engines['jmetal'] = JMetalEA
except ImportError:
    pass


algorithms ={ 'inspyred' : ['NSGAII'],
              'jmetal': ['NSGAII','SPEA2','NSGAIII']
            }


default_engine = None
preferred_EA = 'NSGAII'

def get_default_engine():

    global default_engine

    if default_engine:
        return default_engine

    engine_order = ['inspyred', 'jmetal']

    for engine in engine_order:
        if engine in list(engines.keys()):
            default_engine = engine
            break

    if not default_engine:
        raise RuntimeError("No EA engine available.")

    return default_engine


def set_default_engine(enginename):
    """ Sets default EA engine.
    
    Parameters:
    
    enginename : (str) engine name (currently available: 'inspyred', 'jmetal')
    """

    global default_engine

    if enginename.lower() in list(engines.keys()):
        default_engine = enginename.lower()
    else:
        raise RuntimeError(f"EA engine {enginename} not available.")


def set_preferred_EA(algorithm):
    "Defines de preferred MOEA."
    global preferred_EA
    global default_engine

    if algorithm in algorithms[get_default_engine()]:
        preferred_EA = algorithm
    else:
        for eng in engines.keys():
            if algorithm in algorithms[eng]:
                preferred_EA = algorithm            
                default_engine = eng
                return
        raise ValueError(f"Algorithm {algorithm} is unavailable.")


def get_preferred_EA():
    global preferred_EA
    return preferred_EA


def get_available_engines():
    return list(engines.keys())




def EA(problem, initial_population=[], max_generations=EAConstants.MAX_GENERATIONS, mp=True, visualizer=False):
    """
    EA running helper

    Parameters:
    
    problem: the optimization problem
    initial_population* (list): the EA initial population
    max_generations (int): the number of iterations of the EA (stopping criteria)
    mp (boolean): if should use multiprocessing
    visualizer (boolean): if the pareto font should be displayed. Requires a graphic environment.
    
    Returns:
    
    An instance of an EA optimizer.
    
    """

    if len(engines) == 0:
        raise RuntimeError('Inspyred or JMetal packages are required.')
    engine = engines[get_default_engine()]
    return engine(problem, initial_population=initial_population, max_generations=max_generations, mp=mp, visualizer=visualizer,algorithm=get_preferred_EA())
