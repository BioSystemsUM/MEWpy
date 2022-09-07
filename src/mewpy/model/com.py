from mewpy.simulation import get_simulator
from mewpy.simulation.simulation import Simulator
from mewpy.util.parsing import build_tree, Boolean
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

    def __init__(self, models: list, copy_models=False, flavor='cobrapy'):
        self.organisms = AttrDict()
        self.model_ids = list({model.id for model in models})
        self.flavor = flavor
        self.biomasses = None
        self.biomass = None
        self.reaction_map = None
        self.metabolite_map = None

        if len(self.model_ids) < len(models):
            warn("Model ids are not unique, repeated models will be discarded.")

        for model in models:
            m = model if isinstance(model, Simulator) else get_simulator(model)
            self.organisms[m.id] = deepcopy(m) if copy_models else m

        sid = ' '.join(sorted(self.model_ids))
        if flavor == 'reframed':
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
        return self.comm_model

    def size(self):
        return len(self.organisms)

    @property
    def merged_model(self):
        """ To make compatible with REFRAMED"""
        if self.comm_model is None:
            self._merged_model()
        return self.comm_model.model

    def _merge_models(self):
        """Merges the models. 
        """
        old_ext_comps = []
        ext_mets = []
        self.biomasses = {}
        self.reaction_map = {}
        self.metabolite_map = {}

        # default IDs
        ext_comp_id = "ext"
        biomass_id = "community_biomass"
        comm_growth = "community_growth"

        # create external compartment
        self.comm_model.add_compartment(ext_comp_id, "extracellular environment", external=True)

        # community biomass
        self.comm_model.add_metabolite(biomass_id, name="Total community biomass", compartment=ext_comp_id)

        self.comm_model.add_reaction(comm_growth, name="Community growth rate",
                                     stoichiometry={biomass_id: -1},
                                     lb=0, ub=inf)

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
                                                   compartment=rename(met.compartment)
                                                   )
                    self.metabolite_map[(org_id,m_id)] = new_id

                elif m_id not in self.comm_model.metabolites:  # if is external but was not added yet
                    self.comm_model.add_metabolite(m_id,
                                                   formula=met.formula,
                                                   name=met.name,
                                                   compartment=ext_comp_id)
                    ext_mets.append(m_id)

            # add genes
            if self.flavor == 'reframed':
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
                    self.biomasses[org_id]=new_id

                if rxn.gpr:
                    t = build_tree(rxn.gpr, Boolean)
                    ren = {x:rename(x) for x in t.get_operands()}
                    new_gpr = t.replace(ren).to_infix().replace(' | ', ' or ').replace(' & ',' and ')
                else:
                    new_gpr = rxn.gpr
                self.comm_model.add_reaction(new_id,
                                             name=rxn.name,
                                             stoichiometry=new_stoichiometry,
                                             lb=rxn.lb,
                                             ub=rxn.ub,
                                             gpr=new_gpr)
                self.reaction_map[(org_id, r_id)] = new_id
        # Add exchange reactions

        for m_id in ext_mets:
            m = m_id[len(self._met_prefix):] if m_id.startswith(self._met_prefix) else m_id
            r_id = f"{self._rxn_prefix}EX_{m}"
            # self.comm_model.add_reaction(r_id, name=f'{m_id} exchange', stoichiometry={m_id: -1}, lb=-1000, ub=inf)
            self.comm_model.add_reaction(r_id, name=r_id, stoichiometry={m_id: -1}, lb=-inf, ub=inf)

        self.comm_model.objective = comm_growth
        self.biomass = comm_growth
        return self.comm_model
