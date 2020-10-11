from mewpy.utils.constants import EAConstants
from .ea import Solution

engines = dict()

try:
    from .inspyred.ea import EA as InspyredEA
    engines['inspyred'] = InspyredEA
except ImportError:
    print("inspyred not available")

try:
    from .jmetal.ea import EA as JMetalEA
    engines['jmetal'] = JMetalEA
except ImportError:
    print("jmetal not available")


algorithms ={ 'inspyred' : ['SA','GA','NSGAII'],
              'jmetal'   : ['SA','GA','NSGAII','SPEA2','NSGAIII']
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
   
    :param str enginename: Optimization engine (currently available: 'inspyred', 'jmetal')
    """

    global default_engine

    if enginename.lower() in list(engines.keys()):
        default_engine = enginename.lower()
    else:
        raise RuntimeError(f"EA engine {enginename} not available.")


def set_preferred_EA(algorithm):
    """Defines de preferred MOEA.
    
    :param str algorithm: The name of the preferred algorithm.
    """
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
    """
    :returns: The name of the preferred MOEA.
    """
    global preferred_EA
    return preferred_EA


def get_available_engines():
    """
    :returns: The list of available engines.
    """
    return list(engines.keys())


def get_available_algorithms():
    """
    :returns: The list of available MOEAs.
    """
    algs = []
    for engine in engines.keys():
        algs.extend(algorithms[engine])
    return list(set(algs))



def EA(problem, initial_population=[], max_generations=EAConstants.MAX_GENERATIONS, mp=True, visualizer=False, algorithm = None ):
    """
    EA running helper. Returns an instance of the EA that reflects the global user configuration settings such as preferred engine and algorithm.

    :param problem: The optimization problem.
    :param list initial_population: The EA initial population.
    :param int max_generations: The number of iterations of the EA (stopping criteria).
    :param bool mp: If multiprocessing should be used. 
    :param bool visualizer: If the pareto font should be displayed. Requires a graphic environment.
    :returns: An instance of an EA optimizer.
    
    """
    if len(engines) == 0:
        raise RuntimeError('Inspyred or JMetal packages are required')

    if algorithm is None or algorithm not in get_available_algorithms():
        algorithm=get_preferred_EA()

    engs = [ k for k,v in algorithms.items() if algorithm in v ]
     
    if get_default_engine() in engs:
        engine = engines[get_default_engine()]
    else:
        engine = engines[engs[0]]
     
    return engine(problem, initial_population=initial_population, max_generations=max_generations, mp=mp, visualizer=visualizer,algorithm=algorithm)
