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
Compartmentalized community model.
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

from typing import Dict, List, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from mewpy.simulation import Simulator
    from cobra.core import Model
    from reframed.core.cbmodel import CBModel


class CommunityModel:
    
    EXT_COMP = "e"
    GROWTH_ID = "community_growth"
    
    def __init__(self, 
                 models:List[Union["Simulator","Model","CBModel"]], 
                 abundances:List[float]= None, 
                 merge_biomasses:bool=False, 
                 copy_models:bool=False,
                 add_compartments=False,
                 flavor:str='reframed'):
        """Community Model.

        :param models: A list of metabolic models.
        :param abundances: A list of relative abundances for each model.
            Default None.
        :param merge_biomasses: If a biomass equation is to be build requiring
            each organism to grow in acordance to a relative abundance.
            Default False.
            If no abundance list is provided all organism will have equal abundance.
        :param add_compartments: If each organism external compartment is to be added
            to the community model. Default False.
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
        self.gene_map = None
        
        self._reverse_map = None
        self._merge_biomasses = True if abundances is not None else merge_biomasses
        self._add_compartments = add_compartments
        
        if len(self.model_ids) < len(models):
            warn("Model ids are not unique, repeated models will be discarded.")

        for model in models:
            m = get_simulator(model)
            if not m.objective:
                raise ValueError(f"Model {m.id} has no objective")
            self.organisms[m.id] = deepcopy(m) if copy_models else m
        
        if self._merge_biomasses:
            if abundances and len(abundances)==len(self.organisms):
                self.organisms_abundance =dict(zip(self.organisms.keys(),abundances))
            else:     
                self.organisms_abundance = {org_id:1 for org_id in self.organisms.keys()}
        
        self._comm_model=None
    
    def init_model(self):
        sid = ' '.join(sorted(self.model_ids))
        if self.flavor == 'reframed':
            from reframed.core.cbmodel import CBModel
            model = CBModel(sid)
        else:
            from cobra.core.model import Model
            model = Model(sid)
        self._comm_model = get_simulator(model)
                                                                   
     
    def clear(self):
        self.organisms_biomass = None
        self.organisms_biomass_metabolite = None
        self.biomass = None
        self.reaction_map = None
        self.metabolite_map = None
        self.gene_map = None
        self._reverse_map = None
        self._comm_model = None
   
    @property
    def add_compartments(self):
        return self._add_compartments
    
    @add_compartments.setter
    def add_compartments(self,value:bool):
        if self._add_compartments == value:
            pass
        else:
            self._add_compartments = value
            self.clear()
            
    @property
    def merge_biomasses(self):
        return self._merge_biomasses
    
    @merge_biomasses.setter
    def merge_biomasses(self,value:bool):
        if self._merge_biomasses == value:
            pass
        else:
            self._merge_biomasses = value
            self.clear()
        
    @property
    def reverse_map(self):
        if self._reverse_map is not None:
            return self._reverse_map
        else:
            self._reverse_map = dict()
            self._reverse_map.update({v:k for k,v in self.reaction_map.items()})
            self._reverse_map.update({v:k for k,v in self.gene_map.items()})
            
    def set_abundance(self,abundances:Dict[str,float],rebuild=False):
        if not self._merge_biomasses:
            raise ValueError("The community model has no merged biomass equation")
        self.organisms_abundance.update(abundances)
        if any([x<0 for x in abundances.values()]):
            raise ValueError("All abundance value need to be non negative.")
        if sum(list(abundances.values()))==0:
            raise ValueError("At leat one organism need to have a positive abundance.")
        # update the biomass equation
        if rebuild:
            self.clear()
            self._merge_models()
        else:
            comm_growth = CommunityModel.GROWTH_ID
            biomass_stoichiometry = {met: -self.organisms_abundance[org_id]
                                     for org_id, met in self.organisms_biomass_metabolite.items()
                                     if self.organisms_abundance[org_id]>0
                                     }
            self._comm_model.add_reaction(comm_growth,
                                         name="Community growth rate",
                                         stoichiometry=biomass_stoichiometry,
                                         lb=0, ub=inf, reaction_type='SINK')
            self._comm_model.objective = comm_growth
            self._comm_model.solver=None

    def get_community_model(self):
        """Returns a Simulator for the merged model"""
        return self.merged_model

    def size(self):
        return len(self.organisms)

    @property
    def merged_model(self):
        """ Returns a community model (COBRApy or REFRAMED)"""
        if self._comm_model is None:
            self._merge_models()
        return self._comm_model

    def _merge_models(self):
        """Merges the models."""
        
        self.init_model()
        
        old_ext_comps = []
        ext_mets = []
        self.organisms_biomass = {}
        self.reaction_map = {}
        self.metabolite_map = {}
        self.gene_map = {}
        self._reverse_map = None
        
        if self._merge_biomasses:
            self.organisms_biomass_metabolite = {}

        # default IDs
        ext_comp_id = CommunityModel.EXT_COMP
        comm_growth = CommunityModel.GROWTH_ID

        # create external compartment
        self._comm_model.add_compartment(ext_comp_id, 
                                        "extracellular environment",
                                        external=True)

        # community biomass
        if not self._merge_biomasses:
            biomass_id = "community_biomass"
            self._comm_model.add_metabolite(biomass_id,
                                           name="Total community biomass",
                                           compartment=ext_comp_id)

        # add each organism
        for org_id, model in self.organisms.items():

            def rename(old_id):
                return f"{old_id}_{org_id}"

            def r_gene(old_id, organism=True):
                if model._g_prefix == self._comm_model._g_prefix:
                    _id = old_id
                else:
                    _id = self._comm_model._g_prefix+old_id[len(model._g_prefix):]
                return rename(_id) if organism else _id

            def r_met(old_id, organism=True):
                if model._m_prefix == self._comm_model._m_prefix:
                    _id = old_id
                else:
                    _id = self._comm_model._m_prefix+old_id[len(model._m_prefix):]
                return rename(_id) if organism else _id

            def r_rxn(old_id, organism=True):
                if model._r_prefix == self._comm_model._r_prefix:
                    _id = old_id
                else:
                    _id = self._comm_model._r_prefix+old_id[len(model._r_prefix):]
                return rename(_id) if organism else _id

            # add internal compartments
            for c_id in model.compartments:
                comp = model.get_compartment(c_id)
                if comp.external:
                    old_ext_comps.append(c_id)
                    if not self._add_compartments:
                        continue    
                self._comm_model.add_compartment(rename(c_id), name=f"{comp.name} ({org_id})")
                
            # add metabolites
            for m_id in model.metabolites:
                met = model.get_metabolite(m_id)
                if met.compartment not in old_ext_comps or self._add_compartments:
                    new_mid = r_met(m_id)
                    self._comm_model.add_metabolite(new_mid,
                                                formula=met.formula,
                                                name=met.name,
                                                compartment=rename(met.compartment)
                                                )
                    self.metabolite_map[(org_id, m_id)] = new_mid
                    
                    
                    

                if met.compartment in old_ext_comps and r_met(m_id, False) not in self._comm_model.metabolites:
                    new_mid = r_met(m_id, False)
                    self._comm_model.add_metabolite(new_mid,
                                                   formula=met.formula,
                                                   name=met.name,
                                                   compartment=ext_comp_id)
                    ext_mets.append(new_mid)
            
            # add genes
            for g_id in model.genes:
                new_id = r_gene(g_id)
                self.gene_map[(org_id,g_id)] = new_id
                if self.flavor == 'reframed':
                    gene = model.get_gene(g_id)
                    self._comm_model.add_gene(new_id, gene.name)
                    
            # add reactions
            ex_rxns = model.get_exchange_reactions()
            
            for r_id in model.reactions:
                rxn = model.get_reaction(r_id)
                new_id = r_rxn(r_id)

                if r_id in ex_rxns:
                    mets = list(rxn.stoichiometry.keys())
                    
                    if self._add_compartments and r_met(mets[0], False) in ext_mets:
                        new_stoichiometry = {r_met(mets[0]): -1,
                                            r_met(mets[0],False): 1
                                            }
                        self._comm_model.add_reaction(new_id,
                                                    name=rxn.name,
                                                    stoichiometry=new_stoichiometry,
                                                    lb=-inf,
                                                    ub=inf,
                                                    reaction_type='TRP')
                        self.reaction_map[(org_id, r_id)] = new_id
                    
                        
                    elif (len(mets) == 1 
                          and r_met(mets[0]) in self._comm_model.metabolites):
                        # some models (e.g. AGORA models) have sink reactions (for biomass) 
                        new_stoichiometry = {r_met(mets[0]): -1}
                        self._comm_model.add_reaction(new_id,
                                                    name=rxn.name,
                                                    stoichiometry=new_stoichiometry,
                                                    lb=0,
                                                    ub=inf,
                                                    reaction_type='SINK')
                        self.reaction_map[(org_id, r_id)] = new_id
                    
                else:
                    if self._add_compartments:
                        new_stoichiometry = {
                            r_met(m_id): coeff
                            for m_id, coeff in rxn.stoichiometry.items()
                            }
                    else:
                        new_stoichiometry = {
                            r_met( m_id, False) if r_met( m_id, False) in ext_mets 
                            else r_met(m_id): coeff
                            for m_id, coeff in rxn.stoichiometry.items()
                            }
                    # assumes that the models' objective is the biomass    
                    if r_id in [x for x, v in model.objective.items() if v > 0]:
                        if self._merge_biomasses:
                            met_id = r_met('Biomass')
                            self._comm_model.add_metabolite(
                                met_id,
                                name=f"Biomass {org_id}",
                                compartment=ext_comp_id)
                            
                            new_stoichiometry[met_id] = 1
                            self.organisms_biomass_metabolite[org_id] = met_id
                            
                            # add biomass sink reaction
                            self._comm_model.add_reaction(
                                r_rxn('Sink_biomass'),
                                name=f"Sink Biomass {org_id}",
                                stoichiometry={met_id:-1},
                                lb=0,
                                ub=inf,
                                reaction_type='SINK')
                            
                        else:
                            new_stoichiometry[biomass_id] = 1

                        self.organisms_biomass[org_id] = new_id

                    if rxn.gpr:
                        t = build_tree(rxn.gpr, Boolean)
                        ren = {x: r_gene(x) for x in t.get_operands()}
                        new_gpr = t.replace(ren).to_infix()
                    else:
                        new_gpr = rxn.gpr

                    self._comm_model.add_reaction(new_id,
                                                 name=rxn.name,
                                                 stoichiometry=new_stoichiometry,
                                                 lb=rxn.lb,
                                                 ub=rxn.ub,
                                                 gpr=new_gpr,
                                                 annotations=rxn.annotations)

                    self.reaction_map[(org_id, r_id)] = new_id

        # Add exchange reactions
        for m_id in ext_mets:
            m = m_id[len(self._comm_model._m_prefix):] if m_id.startswith(self._comm_model._m_prefix) else m_id
            r_id = f"{self._comm_model._r_prefix}EX_{m}"
            self._comm_model.add_reaction(r_id, name=r_id, stoichiometry={m_id: -1}, lb=-inf, ub=inf, reaction_type="EX")

        if self._merge_biomasses:
            # if the biomasses are to be merged add
            # a new product to each organism biomass  
            biomass_stoichiometry = {met: -1*self.organisms_abundance[org_id]
                                     for org_id, met in self.organisms_biomass_metabolite.items()
                                    }
        else:
            biomass_stoichiometry = {biomass_id: -1}

        self._comm_model.add_reaction(comm_growth, name="Community growth rate",
                                     stoichiometry=biomass_stoichiometry,
                                     lb=0, ub=inf, reaction_type='SINK')

        self._comm_model.objective = comm_growth
        self._comm_model.biomass_reaction = comm_growth
        self.biomass = comm_growth
        return self._comm_model

    def copy(self, copy_models=False, flavor=None):
        models = [m.model for m in self.organisms.values()]
        f = flavor if flavor is not None else self.flavor
        return CommunityModel(models, copy_models=copy_models, flavor=f)
