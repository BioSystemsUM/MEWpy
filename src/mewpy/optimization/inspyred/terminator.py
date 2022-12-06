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
Terminators for inspyred

Authors: Vitor Pereira
##############################################################################
"""

def generation_termination(population, num_generations, num_evaluations, args):
    """Return True if the number of generations meets or exceeds a maximum.

    This function compares the number of generations with a specified 
    maximum. It returns True if the maximum is met or exceeded.

    .. Arguments:
       population -- the population of Individuals
       num_generations -- the number of elapsed generations
       num_evaluations -- the number of candidate solution evaluations
       args -- a dictionary of keyword arguments

    Optional keyword arguments in args:

    - *max_generations* -- the maximum generations (default 1) 

    """
    max_gen = args.get('max_generations', 1)
    if isinstance(max_gen, tuple):
        max_generations = max_gen[0]
    else:
        max_generations = max_gen
    return num_generations >= max_generations
