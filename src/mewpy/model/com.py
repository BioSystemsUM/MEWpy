from mewpy.simulation import get_simulator
from mewpy.simulation.simulation import Simulator
from mewpy.util import AttrDict
from copy import deepcopy
from warnings import warn
from numpy import inf


class CommunityModel:
    """
    Community Model.

    :param models: A list of metabolic models.
    :param list fevaluation: A list of callable EvaluationFunctions.

    Optional parameters:
    :param bool copy_models: if the models are to be copied, default True.
    :param str flavor: use 'cobrapy' or 'reframed. Default 'cobrapy'. 

    """

    def __init__(self, models: list, copy_models=True, flavor='cobrapy'):
        self.organisms = AttrDict()
        self.model_ids = list({model.id for model in models})

        if len(self.model_ids) < len(models):
            warn("Model ids are not unique, repeated models will be discarded.")

        for model in models:
            m = model if isinstance(model, Simulator) else get_simulator(model)
            self.organisms[m.id] = deepcopy(m) if copy_models else m

        sid = ' '.join(sorted(self.model_ids))
        if flavor=='reframed':
            from reframed.core.cbmodel import CBModel
            model = CBModel(sid)
            self._met_prefix = 'M_'
            self._rxn_prefix = 'R_'
        else:
            from cobra.core.model import Model
            model = Model(sid)
            self._met_prefix = ''
            self._rxn_prefix = ''

        self.comm_model = get_simulator(model) 
        self._merge_models()        


    def get_community_model(self):
        return self.model


    def _merge_models(self):
        """Merges the models. 
        """
        old_ext_comps = []
        ext_mets = []

        # default IDs
        ext_comp_id = "ext"
        biomass_id = "community_biomass"
        comm_growth = "community_growth"

        # create external compartment
        self.comm_model.add_compartment(ext_comp_id, "extracellular environment", external=True)

        # community biomass
        self.comm_model.add_metabolite(biomass_id, name="Total community biomass",compartment=ext_comp_id)

        self.comm_model.add_reaction(comm_growth, name="Community growth rate",
                                     reversible=False, stoichiometry={biomass_id: -1},
                                     lb=0, ub=inf, objective=1)

        # add each organism

        for org_id, model in self.organisms.items():

            def rename(old_id):
                return f"{old_id}_{org_id}"

            # add internal compartments

            for c_id in model.compartments:
                comp = model.get_compartment(c_id)
                if comp.external:
                    old_ext_comps.append(c_id)
                else:
                    self.comm_model.add_compartment(rename(c_id), name=comp.name)

            # add metabolites

            for m_id in model.metabolites:
                met = model.get_metabolite(m_id)
                if met.compartment not in old_ext_comps:  # if is internal
                    new_id = rename(m_id)
                    self.comm_model.add_metabolite(new_id,
                                              formula=met.formula,
                                              name=met.name,
                                              compartment= rename(met.compartment)
                                             )

                elif m_id not in self.comm_model.metabolites:  # if is external but was not added yet
                    self.comm_model.add_metabolite(m_id,
                                              formula=met.formula,
                                              name=met.name,
                                              compartment=ext_comp_id)
                    ext_mets.append(m_id)

            # add genes

            for g_id in model.genes:
                gene = model.get_gene(g_id)
                new_id = rename(g_id)
                self.comm_model.add_gene(new_id, gene.name)

            # add internal reactions
            ex_rxns = model.get_exchange_reactions()
            for r_id in model.reactions:
                rxn = model.get_reaction(r_id)
                if r_id in ex_rxns:
                    continue

                new_id = rename(r_id)
                new_stoichiometry = {
                    m_id if m_id in ext_mets else rename(m_id): coeff
                    for m_id, coeff in rxn.stoichiometry.items()
                }

                if r_id in [x for x, v in model.objective.items() if v > 0]:
                    new_stoichiometry[biomass_id] = 1

        
                new_gpr = rxn.gpr
                self.comm_model.add_reaction(new_id,
                                             name=rxn.name,
                                             stoichiometry=new_stoichiometry,
                                             lb=rxn.lb,
                                             ub=rxn.ub,
                                             gpr_association=new_gpr)

        # Add exchange reactions

        for m_id in ext_mets:
            r_id = f"{self._rxn_prefix}EX_{m_id[len(self._met_prefix):]}" 
            self.comm_model.add_reaction(r_id, reversible=True, stoichiometry={m_id: -1})
        return self.comm_model


