from ..util.utilities import Singleton


class EAGlobalSettings(Singleton):
    POPULATION_SIZE = 100


def get_default_population_size():
    return EAGlobalSettings.POPULATION_SIZE


def set_default_population_size(size: int):
    EAGlobalSettings.POPULATION_SIZE = size
