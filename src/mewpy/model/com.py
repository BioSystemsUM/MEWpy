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
Compartementalized community Model inspired in REFRAMED.
Can build community models from models loaded from any toolbox
for which there is a Simulator implementation.

Author: Vitor Pereira
##############################################################################
"""
from mewpy.simulation import get_simulator
from mewpy.util.parsing import build_tree, Boolean
from mewpy.util import AttrDict
from copy import deepcopy
from warnings import warn
from numpy import inf


class CommunityModel:

    def __init__(self, models: list, copy_models=False, merge_biomass=False, flavor='reframed'):
        """
        Community Model.

        :param models: A list of metabolic models.
        :param list fevaluation: A list of callable EvaluationFunctions.

        Optional parameters:
        :param bool copy_models: if the models are to be copied, default True.
        :param str flavor: use 'cobrapy' or 'reframed. Default 'reframed'.
        """

        self.organisms = AttrDict()
        self.model_ids = list({model.id for model in models})
        self.flavor = flavor
        self.organisms_biomass = None
        self.organisms_biomass_metabolite = None
        self.biomass = None

        self.reaction_map = None
        self.metabolite_map = None
        self.merge_biomass = merge_biomass
        if merge_biomass:
            self.organisms_coefficient = {}

        if len(self.model_ids) < len(models):
            warn("Model ids are not unique, repeated models will be discarded.")

        for model in models:
            m = get_simulator(model)
            if not m.objective:
                raise ValueError(f"Model {m.id} has no objective")
            self.organisms[m.id] = deepcopy(m) if copy_models else m
            if merge_biomass:
                self.organisms_coefficient[m.id] = 1

        sid = ' '.join(sorted(self.model_ids))
        if flavor == 'reframed':
            from reframed.core.cbmodel import CBModel
            model = CBModel(sid)
        else:
            from cobra.core.model import Model
            model = Model(sid)

        self.comm_model = get_simulator(model)
        self._merge_models()

    def get_community_model(self):
        """Returns a Simulator for the merged model"""
        return self.comm_model

    def size(self):
        return len(self.organisms)

    @property
    def merged_model(self):
        """ Returns a community model (COBRApy or REFRAMED)
        To make compatible with REFRAMED"""
        if self.comm_model is None:
            self._merged_model()
        return self.comm_model.model

    def _merge_models(self):
        """Merges the models. 
        """
        old_ext_comps = []
        ext_mets = []
        self.organisms_biomass = {}
        self.reaction_map = {}
        self.metabolite_map = {}
        if self.merge_biomass:
            self.organisms_biomass_metabolite = {}

        # default IDs
        ext_comp_id = "e"

        comm_growth = "community_growth"

        # create external compartment
        self.comm_model.add_compartment(ext_comp_id, "extracellular environment", external=True)

        # community biomass
        if not self.merge_biomass:
            biomass_id = "community_biomass"
            self.comm_model.add_metabolite(biomass_id, name="Total community biomass", compartment=ext_comp_id)

        # add each organism
        for org_id, model in self.organisms.items():

            def rename(old_id):
                return f"{old_id}_{org_id}"

            def rename_gene(old_id, organism=True):
                if model._g_prefix == self.comm_model._g_prefix:
                    _id = old_id
                else:
                    _id = self.comm_model._g_prefix+old_id[len(model._g_prefix):]
                return rename(_id) if organism else _id

            def rename_met(old_id, organism=True):
                if model._m_prefix == self.comm_model._m_prefix:
                    _id = old_id
                else:
                    _id = self.comm_model._m_prefix+old_id[len(model._m_prefix):]
                return rename(_id) if organism else _id

            def rename_rxn(old_id, organism=True):
                if model._r_prefix == self.comm_model._r_prefix:
                    _id = old_id
                else:
                    _id = self.comm_model._r_prefix+old_id[len(model._r_prefix):]
                return rename(_id) if organism else _id

            # add internal compartments
            for c_id in model.compartments:
                comp = model.get_compartment(c_id)
                if comp.external:
                    old_ext_comps.append(c_id)
                self.comm_model.add_compartment(rename(c_id), name=f"{comp.name} ({org_id})")

            # add metabolites
            for m_id in model.metabolites:
                met = model.get_metabolite(m_id)

                new_id = rename_met(m_id)
                self.comm_model.add_metabolite(new_id,
                                               formula=met.formula,
                                               name=met.name,
                                               compartment=rename(met.compartment)
                                               )
                self.metabolite_map[(org_id, m_id)] = new_id

                if met.compartment in old_ext_comps and m_id not in self.comm_model.metabolites:  # if is external but was not added yet
                    new_mid = rename_met(m_id, False)
                    self.comm_model.add_metabolite(new_mid,
                                                   formula=met.formula,
                                                   name=met.name,
                                                   compartment=ext_comp_id)
                    ext_mets.append(new_mid)

            # add genes
            if self.flavor == 'reframed':
                for g_id in model.genes:
                    gene = model.get_gene(g_id)
                    new_id = rename_gene(g_id)
                    self.comm_model.add_gene(new_id, gene.name)

            # add reactions
            ex_rxns = model.get_exchange_reactions()

            for r_id in model.reactions:
                rxn = model.get_reaction(r_id)

                new_id = rename_rxn(r_id)

                if r_id in ex_rxns:
                    mets = list(rxn.stoichiometry.keys())
                    # this condition should not be necessary...
                    if len(mets) == 1 and rename_met(mets[0], False) in ext_mets:
                        new_stoichiometry = {rename_met(mets[0], False): -1,
                                             rename_met(mets[0]): 1
                                             }
                        self.comm_model.add_reaction(new_id,
                                                     name=rxn.name,
                                                     stoichiometry=new_stoichiometry,
                                                     lb=-inf,
                                                     ub=inf,
                                                     reaction_type='TRP')
                        self.reaction_map[(org_id, r_id)] = new_id
                else:
                    new_stoichiometry = {
                        rename_met(m_id): coeff
                        for m_id, coeff in rxn.stoichiometry.items()
                    }
                    if r_id in [x for x, v in model.objective.items() if v > 0]:
                        if self.merge_biomass:
                            met_id = rename_met(r_id)
                            self.comm_model.add_metabolite(met_id,
                                                           name=f"BIOMASS {org_id}",
                                                           compartment=ext_comp_id)
                            new_stoichiometry[met_id] = 1
                            self.organisms_biomass_metabolite[org_id] = met_id
                        else:
                            new_stoichiometry[biomass_id] = 1

                        self.organisms_biomass[org_id] = new_id

                    if rxn.gpr:
                        t = build_tree(rxn.gpr, Boolean)
                        ren = {x: rename_gene(x) for x in t.get_operands()}
                        new_gpr = t.replace(ren).to_infix()
                    else:
                        new_gpr = rxn.gpr

                    self.comm_model.add_reaction(new_id,
                                                 name=rxn.name,
                                                 stoichiometry=new_stoichiometry,
                                                 lb=rxn.lb,
                                                 ub=rxn.ub,
                                                 gpr=new_gpr,
                                                 annotations=rxn.annotations)

                    self.reaction_map[(org_id, r_id)] = new_id

        # Add exchange reactions
        for m_id in ext_mets:
            m = m_id[len(self.comm_model._m_prefix):] if m_id.startswith(self.comm_model._m_prefix) else m_id
            r_id = f"{self.comm_model._r_prefix}EX_{m}"
            self.comm_model.add_reaction(r_id, name=r_id, stoichiometry={m_id: -1}, lb=-inf, ub=inf, reaction_type="EX")

        if self.merge_biomass:
            biomass_stoichiometry = {met: self.organisms_coefficient[org_id]
                                     for org_id, met in self.organisms_biomass_metabolite.items()
                                     }
        else:
            biomass_stoichiometry = {biomass_id: -1}

        self.comm_model.add_reaction(comm_growth, name="Community growth rate",
                                     stoichiometry=biomass_stoichiometry,
                                     lb=0, ub=inf, reaction_type='SINK')

        self.comm_model.objective = comm_growth
        self.biomass = comm_growth
        return self.comm_model

    def copy(self, copy_models=False, flavor=None):
        models = [m.model for m in self.organisms.values()]
        f = flavor if flavor is not None else self.flavor
        return CommunityModel(models, copy_models=copy_models, flavor=f)
