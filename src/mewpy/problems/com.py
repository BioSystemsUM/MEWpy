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
Optimization Problems for Community models

Author: Vitor Pereira
##############################################################################
"""
from copy import deepcopy
from warnings import warn
from collections import OrderedDict
from .problem import AbstractKOProblem
from mewpy.model import CommunityModel
from mewpy.simulation import get_simulator
from mewpy.simulation.simulation import Simulator
from typing import TYPE_CHECKING, List, Union

if TYPE_CHECKING:
    from cobra.core import Model
    from reframed.core.cbmodel import CBModel


class CommunityKOProblem(AbstractKOProblem):
    """
    Community Knockout Optimization Problem.

    :param models: A list of metabolic models.
    :param list fevaluation: A list of callable EvaluationFunctions.

    Optional parameters:

    :param OrderedDict envcond: Environmental conditions.
    :param OrderedDict constraints: Additional constraints to be applied to the model.
    :param int candidate_min_size: The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
    :param int candidate_max_size: The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
    :param list target: List of modification target reactions.
    :param list non_target: List of non target reactions. Not considered if a target list is provided.
    :param float scalefactor: A scaling factor to be used in the LP formulation.
    
    """

    def __init__(self, models: List[Union[Simulator, 'Model', 'CBModel']],
                 fevaluation=[],
                 copy_models: bool = False,
                 **kwargs):

        self.organisms = OrderedDict()
        self.model_ids = list({model.id for model in models})
        self.flavor = kwargs.get('flavor', 'reframed')

        if len(self.model_ids) < len(models):
            warn("Model ids are not unique, repeated models will be discarded.")

        for model in models:
            m = get_simulator(model)
            self.organisms[m.id] = deepcopy(m) if copy_models else m

        self.cmodel = CommunityModel(list(self.organisms.values()), flavor=self.flavor)
        model = self.cmodel.merged_model

        super(CommunityKOProblem, self).__init__(
            model, fevaluation=fevaluation, **kwargs)


    def _build_target_list(self):
        """Target organims, that is, organisms that may be removed from the community.
        """
        target = set(self.model_ids)
        if self.non_target is not None:
            target = target - set(self.non_target)
        self._trg_list = list(target)

    def ext_reactions(self, organism):
        org_mets = set([v for k, v in self.cmodel.metabolite_map.items() if k[0] == organism])
        rxns = [v for k, v in self.cmodel.reaction_map.items() if k[0] == organism]
        ext_mets = set([m for m in self.simulator.metabolites if self.simulator.get_metabolite(m).compartment == 'ext'])
        res = []
        for rxn in rxns:
            st = self.simulator.get_reaction(rxn).stoichiometry
            sub = set([k for k, v in st.items() if v < 0])
            prod = set([k for k, v in st.items() if v > 0])
            if (len(sub-ext_mets) == 0 and len(prod-org_mets) == 0):
                res.append(rxn)
        return res

    def solution_to_constraints(self, candidate):
        ko_organisms = list(candidate.keys())
        constraints = dict()
        for org_id in ko_organisms:
            uptakes = self.ext_reactions(org_id)
            for r_id in uptakes:
                constraints[r_id] = 0
            # ko biomass
            constraints[self.cmodel.organisms_biomass[org_id]] = 0
        return constraints


