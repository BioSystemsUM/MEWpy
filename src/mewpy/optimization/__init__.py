from mewpy.utils.constants import EAConstants


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


default_engine = None


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


def get_available_engines():
    return list(engines.keys())


def set_default_engine(enginename):
    """ Sets default EA engine.
    Arguments:
        enginename : (str) engine name (currently available: 'inspyred', 'jmetal')
    """

    global default_engine

    if enginename.lower() in list(engines.keys()):
        default_engine = enginename.lower()
    else:
        raise RuntimeError(f"EA engine {enginename} not available.")


def EA(problem, initial_population=[], max_generations=EAConstants.MAX_GENERATIONS, mp=True, visualizer=False):
    if len(engines) == 0:
        raise RuntimeError('Inspyred or JMetal packages are required.')
    engine = engines[get_default_engine()]
    return engine(problem, initial_population=initial_population, max_generations=max_generations, mp=mp, visualizer=visualizer)
