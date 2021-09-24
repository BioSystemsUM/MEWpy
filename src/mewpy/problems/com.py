import os
from collections import OrderedDict
from time import time
from math import inf
import hashlib
from copy import deepcopy

from reframed.io.sbml import load_cbmodel, parse_gpr_rule
from reframed.core.model import AttrOrderedDict, ReactionType, Compartment, Metabolite
from reframed.core.cbmodel import CBModel, CBReaction, Gene

from .problem import AbstractKOProblem
from mewpy.simulation import get_simulator
from mewpy.simulation.simulation import Simulator


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

    def __init__(self, models: list, fevaluation=[], copy_models=True, **kwargs):
        super(CommunityKOProblem, self).__init__(
            None, fevaluation=fevaluation, **kwargs)
        self.organisms = AttrOrderedDict()
        self.model_ids = list({model.id for model in models})

        if len(self.model_ids) < len(models):
            warn("Model ids are not unique, repeated models will be discarded.")

        for model in models:

            m = model if isinstance(model, Simulator) else get_simulator(model)
            self.organisms[m.id] = deepcopy(m) if copy_models else m

    def merge_models(self, models: list = None):
        """Merges the models. Adapted from REFRAMED 

        :param models: [description], defaults to None.
        :type models: list, optional
        :return: [description]
        :rtype: [type]
        """
        if models is None:
            mlist = self.model_ids
        else:
            mlist = models
        sid = ' '.join(sorted(mlist))
        id = hashlib.md5(sid.encode())
        # TODO: use the id to identify communities already evaluated

        comm_model = CBModel(id)
        old_ext_comps = []
        ext_mets = []

        # default IDs
        ext_comp_id = "ext"
        biomass_id = "community_biomass"
        comm_growth = "community_growth"

        # create external compartment

        comp = Compartment(
            ext_comp_id, "extracellular environment", external=True)
        comm_model.add_compartment(comp)

        # community biomass

        met = Metabolite(biomass_id, "Total community biomass", ext_comp_id)
        comm_model.add_metabolite(met)

        rxn = CBReaction(comm_growth, name="Community growth rate",
                         reversible=False, stoichiometry={biomass_id: -1},
                         lb=0, ub=inf, objective=1)

        comm_model.add_reaction(rxn)

        # add each organism

        for org_id, model in self.organisms.items():

            def rename(old_id):
                return f"{old_id}_{org_id}"

            # add internal compartments

            for c_id in model.compartments:
                comp = model.get_compartment(c_id)
                if comp['external']:
                    old_ext_comps.append(c_id)
                else:
                    new_comp = Compartment(rename(c_id), comp['name'])
                    comm_model.add_compartment(new_comp)

            # add metabolites

            for m_id in model.metabolites:
                met = model.get_metabolite(m_id)
                if met['compartment'] not in old_ext_comps:  # if is internal
                    new_id = rename(m_id)
                    new_met = Metabolite(
                        new_id, met['name'], rename(met['compartment']))
                    new_met.metadata['FORMULA'] = met['formula']
                    comm_model.add_metabolite(new_met)

                elif m_id not in comm_model.metabolites:  # if is external but was not added yet
                    new_met = Metabolite(m_id, met['name'], ext_comp_id)
                    new_met.metadata['FORMULA'] = met['formula']
                    comm_model.add_metabolite(new_met)
                    ext_mets.append(new_met.id)

            # add genes

            for g_id in model.genes:
                gene = model.get_gene(g_id)
                new_id = rename(g_id)
                new_gene = Gene(new_id, gene['name'])
                comm_model.add_gene(new_gene)

            # add internal reactions
            ex_rxns = model.get_exchange_reactions()
            for r_id in model.reactions:
                rxn = model.get_reaction(r_id)
                if r_id in ex_rxns:
                    continue

                new_id = rename(r_id)
                new_stoichiometry = {
                    m_id if m_id in ext_mets else rename(m_id): coeff
                    for m_id, coeff in rxn['stoichiometry'].items()
                }

                if r_id in [x for x, v in model.objective.items() if v > 0]:
                    new_stoichiometry[biomass_id] = 1

                if rxn['gpr'] is None:
                    new_gpr = None
                else:
                    new_gpr = parse_gpr_rule(rxn['gpr'])

                new_rxn = CBReaction(
                    new_id,
                    name=rxn['name'],
                    stoichiometry=new_stoichiometry,
                    lb=rxn['lb'],
                    ub=rxn['ub'],
                    gpr_association=new_gpr
                )

                comm_model.add_reaction(new_rxn)

        # Add exchange reactions

        for m_id in ext_mets:
            r_id = f"R_EX_{m_id[2:]}" if m_id.startswith(
                "M_") else f"R_EX_{m_id}"
            rxn = CBReaction(r_id, reversible=True, stoichiometry={m_id: -1},
                             reaction_type=ReactionType.EXCHANGE)
            comm_model.add_reaction(rxn)
        return comm_model

    def _build_target_list(self):
        """Target organims, i.e., organisms that may be removed from the community.
        """
        print("Building modification target list.")
        target = set(self.model_ids)
        if self.non_target is not None:
            target = target - set(self.non_target)
        self._trg_list = list(target)

    def solution_to_constraints(self, candidate):
        """Returns a model not a list of constraints.

        :param candidate: [description]
        :return: [description]
        """
        ko_organisms = list(candidate.keys())
        models = [x for x in self.model_ids if x not in ko_organisms]
        cmodel = self.merge_models(models)
        return cmodel

    def evaluate_solution(self, solution, decode=True):
        """
        Evaluates a single solution, a community.

        :param solution: The solution to be evaluated (a community model).
        :param decode: If the solution needs to be decoded (convert a list of model ids to a community model).
        :returns: A list of fitness values.
        """
        decoded = {}
        # decoded constraints
        if decode:
            decoded = self.decode(solution)
            cmodel = self.solution_to_constraints(decoded)
        else:
            cmodel = solution
        try:
            p = []
            for method in self.methods:
                sim = get_simulator(cmodel, self.environmental_conditions)
                simulation_result = self.simulator.simulate(method=method, scalefactor=self.scalefactor)
                simulation_results[method] = simulation_result
            # apply the evaluation function(s)
            for f in self.fevaluation:
                v = f(simulation_results, decoded,
                      scalefactor=self.scalefactor)
                p.append(v)
        except Exception as e:
            p = []
            for f in self.fevaluation:
                p.append(f.worst_fitness)
        return p
