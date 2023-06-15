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
Cofactor Swap optimization Problem

Author: Vitor Pereira 
##############################################################################
"""
from mewpy.problems.problem import AbstractKOProblem
from mewpy.util.constants import COFACTORS

from copy import deepcopy
from typing import Union, TYPE_CHECKING, List

if TYPE_CHECKING:
    from cobra.core import Model
    from reframed.core.cbmodel import CBModel
    from mewpy.optimization.evaluation import EvaluationFunction

SWAPS = [[COFACTORS['NAD'], COFACTORS['NADH']],
         [COFACTORS['NADP'], COFACTORS['NADPH']]]



class CofactorSwapProblem(AbstractKOProblem):

    RX_SUFIX = '_SWAP'

    def __init__(self, model: Union["Model", "CBModel"],
                 fevaluation: List["EvaluationFunction"] = None,
                 compartments:List[str] = ['c'],
                 **kwargs):
        """
        Optimize co-factor swapping
        
        Implements a search for reactions that when swapped improve the given objectives.
        The approach:

        - finds reactions that have all the targeted co-factor pairs e.g. (nad_c -> nadp_c, nadh_c -> nadph_c)

        - adds reactions that have the co-factors swapped and then by a search algorithm switching one off in favor of the
        other

        References
        ----------
        .. [1] King, Zachary A., and Adam M. Feist. "Optimizing Cofactor Specificity of Oxidoreductase Enzymes for the
        Generation of Microbial Production Strains - OptSwap." Industrial Biotechnology 9, no. 4 (August 1,
        2013): 236-46. - doi:10.1089/ind.2013.0005.

        Args:
            model (Union[&quot;Model&quot;, &quot;CBModel&quot;]): _description_
            fevaluation (List[&quot;EvaluationFunction&quot;], optional): _description_. Defaults to None.
        """
        super().__init__(deepcopy(model), fevaluation, **kwargs)
        self.swaps = None
        self.compartments = compartments
        self.rx_swap = dict()
        
        
    def _build_target_list(self):
        # identify the metabolites
        _swaps = []
        for c in self.compartments:
            for [f1,f2] in SWAPS:
                a = self.simulator.metabolite_by_formula(f1,c)
                b = self.simulator.metabolite_by_formula(f2,c)
                if a and b: _swaps.append([a,b])
        swaps = tuple(_swaps)
        # search reactions
        def _search(mets):
            p = all(mets.get(m, False) for m in swaps[0])
            # remove biomasses
            q = all(mets.get(m, False) for m in swaps[1])
            return (p or q) and not (p and q)
        rxns = []
        for rx_id in self.simulator.reactions:
            mets = self.simulator.get_reaction(rx_id).stoichiometry
            if _search(mets):
                rxns.append(rx_id)
                
        # add reactions with swapped cofactors
        dswap = dict(zip(*swaps))
        a,b = tuple(swaps[0])
        
        for rx_id in rxns:
            rx = self.simulator.get_reaction(rx_id)
            st = rx.stoichiometry.copy()
            if a not in st or b not in st:
                continue
            st[dswap[a]]=st.pop(a)
            st[dswap[b]]=st.pop(b)
            
            new_id = rx_id+self.RX_SUFIX
            self.simulator.add_reaction(
                new_id,
                name=rx.name+" SWAP",
                stoichiometry=st,
                lb=0,
                ub=0,
                gpr=rx.gpr,
                annotations=rx.annotations)
            self.rx_swap[rx_id] = new_id
        # define the modification target list
        # the list of reactions with swapped alternatives
        self.simulator.solver = None
        self._trg_list = list(self.rx_swap.keys())
        
    def solution_to_constraints(self, candidate):
        constraints = {}
        for rx in candidate:
            # ko the original reaction
            constraints[rx]=(0,0)
            # define the bound for the cofactor swapped reaction
            lb,ub = self.simulator.get_reaction_bounds(rx)
            constraints[self.rx_swap[rx]]= (lb,ub)
        return constraints