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
Authors: Vitor Pereira
##############################################################################
"""
# Amino acid MW (Da) retrieved from https://modlamp.org/
aa_weights = {'A': 89.093, 'C': 121.158, 'D': 133.103, 'E': 147.129, 'F': 165.189, 'G': 75.067,
              'H': 155.155, 'I': 131.173, 'K': 146.188, 'L': 131.173, 'M': 149.211, 'N': 132.118,
              'P': 115.131, 'Q': 146.145, 'R': 174.20, 'S': 105.093, 'T': 119.119, 'V': 117.146,
              'W': 204.225, 'Y': 181.189}


def calculate_MW(seq, amide=False):
    """Method to calculate the molecular weight [g/mol] of every sequence in the attribute :py:attr:`sequences`.

    :param (str) seq: amino acid sequence 
    :param (boolean) amide: whether the sequences are C-terminally amidated (subtracts 0.95 from the MW).
    :return: array of descriptor values in the attribute :py:attr:`descriptor`
    
    """
    mw = [aa_weights[aa] for aa in seq]
    # sum over AA MW and subtract H20 MW for every
    mw = round(sum(mw) - 18.015 * (len(seq) - 1), 2)
    # if the sequence is amidated, subtract 0.98 from calculated MW (OH - NH2)
    if amide:
        mw -= 0.98
    return mw
