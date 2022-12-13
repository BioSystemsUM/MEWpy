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
Global settings for Heuristic Optimizations

Author: Vitor Pereira
##############################################################################
"""
from ..util.utilities import Singleton


class EAGlobalSettings(Singleton):
    POPULATION_SIZE = 100


def get_default_population_size():
    return EAGlobalSettings.POPULATION_SIZE


def set_default_population_size(size: int):
    EAGlobalSettings.POPULATION_SIZE = size
