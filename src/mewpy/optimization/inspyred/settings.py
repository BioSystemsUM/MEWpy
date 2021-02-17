
from .operators import (uniform_crossover_OU, grow_mutation_OU, shrink_mutation,
                        single_mutation_OU, uniform_crossover_KO, grow_mutation_KO,
                        single_mutation_KO)

from ..settings import get_default_population_size


def get_population_size():
    size = get_default_population_size()
    return size


OU = {
    'variators': [uniform_crossover_OU,
                  grow_mutation_OU,
                  shrink_mutation,
                  single_mutation_OU]
}

KO = {
    'variators': [uniform_crossover_KO,
                  grow_mutation_KO,
                  shrink_mutation,
                  single_mutation_KO]
}

PARAMETERS = {'gs_mutation_rate': 0.1,
              'mutation_rate': 0.1,
              'crossover_rate': 0.9,
              'tournament_size': 7,
              }
