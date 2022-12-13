# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
##############################################################################
Settings for inspyred

Authors: Vitor Pereira
##############################################################################
"""
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
