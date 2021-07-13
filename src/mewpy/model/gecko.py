import copy
import logging
import os
import re
from collections import OrderedDict
from itertools import chain
from math import isinf

import numpy as np
import pandas as pd
from reframed.core.cbmodel import CBModel, CBReaction
from reframed.core.model import AttrOrderedDict, Metabolite
from reframed.io.sbml import load_cbmodel
from six import iteritems, string_types

from ..util.constants import ModelConstants


class ModelList(object):
    '''Auxilary class to load predifined ecYeast7 models
    '''

    def __init__(self):
        self.DATA_FILES = os.path.join(os.path.dirname(__file__), 'data')
        self.PROTEINS_FILE = os.path.join(self.DATA_FILES, 'proteins.txt')
        self.model_files = dict(
            (re.findall(r'_(.*).xml', f)[0], f) for f in os.listdir(self.DATA_FILES) if f.endswith('.xml'))
        self.models = {}

    def __getitem__(self, item):
        """Get a bundled GECKO model.

        :param item: basestring. Either 'single-pool' for the single-protein pool ecYeastGEM model or
        'multi-pool' for individually modeled protein pools.

        """
        model = None
        try:
            model = load_cbmodel(item)
            file_name = item
        except:
            pass
        if model is None:
            try:
                file_name = self.model_files[item]
                if file_name not in self.models:
                    model = load_cbmodel(os.path.join(os.path.dirname(__file__), 'data/{}'.format(file_name)))
            except KeyError:
                raise KeyError('model name must be one of {}'.format(', '.join(list(self.model_files))))

        self.simplify_model(model)
        self.models[file_name] = model
        return model

    def simplify_model(self, model):
        met_copy = AttrOrderedDict()
        for key, val in model.metabolites.items():
            k = key.replace('__91__', '_').replace('__93__', '')
            val.id = k
            met_copy[k] = val
        model.metabolites = met_copy
        for rxn in model.reactions.keys():
            stoi = model.reactions[rxn].stoichiometry
            nstoi = OrderedDict()
            for key, v in stoi.items():
                k = key.replace('__91__', '_').replace('__93__', '')
                nstoi[k] = v
            model.reactions[rxn].stoichiometry = nstoi
        for rxn in model.reactions.values():
            if isinf(rxn.ub):
                rxn.ub = ModelConstants.REACTION_UPPER_BOUND

    def protein_properties(self):
        return pd.read_csv(self.PROTEINS_FILE, index_col=0)


class GeckoModel(CBModel):
    """Class for representing GECKO models.
    Addapted from the original (https://github.com/SysBioChalmers/GECKO/) to reframed.

    Implement a model class for Genome-scale model to account for Enzyme Constraints, using Kinetics and Omics [1]_.

    :param model: str, A CBModel to apply protein constraints to. Can be 'single-pool' for the bundled ecYeast7 \
        model using only a single pool for all proteins, or 'multi-pool' for the model that has separate pools \
        for all measured proteins.
    :param protein_properties: pd.DataFrame, A data frame that defined molecular weight (g/mol) 'mw', for 'uniprot' \
        proteins and their average 'abundance' in ppm.
    :param sigma: float, The parameter adjusting how much of a protein pool can take part in reactions. Fitted \
        parameter, default is optimized for chemostat experiment in [1]_.
    :param gam: float, The growth associated maintenance cost in mmol / gDW. Default fitted for yeast 8.1.3.
    :param amino_acid_polymerization_cost: float, The cost for turning amino-acids in proteins in mmol / g. \
        Default taken from [2]_.
    :param carbohydrate_polymerization_cost: float, The cost for turning monosaccharides in polysaccharides in \
        mmol / g. Default taken from [2]_.
    :param c_base: float, The carbohydrate content at dilution rate 0.1 / h. Default taken from yeast 8.1.3.
    :param biomass_reaction_id: str, The identifier for the biomass reaction.
    :param protein_reaction_id: str, The identifier for the protein reaction.
    :param carbohydrate_reaction_id: str, The identifier for the carbohydrate reaction.
    :param protein_pool_exchange_id: str, The identifier of the protein pool exchange reaction.
    :param common_protein_pool_id: str, The identifier of the metabolite representing the common protein pool.

    References:

    .. [1] Benjamin J. Sanchez, Cheng Zhang, Avlant Nilsson, Petri-Jaan Lahtvee, Eduard J. Kerkhoven, Jens Nielsen (
       2017). Improving the phenotype predictions of a yeast genome-scale metabolic model by incorporating enzymatic
       constraints. [Molecular Systems Biology, 13(8): 935, http://www.dx.doi.org/10.15252/msb.20167411

       [2] J. Förster, I. Famili, B. Ø. Palsson and J. Nielsen, Genome Res., 2003, 244–253.

    """

    def __init__(self, model, protein_properties=None,
                 sigma=0.46, c_base=0.3855, gam=36.6, amino_acid_polymerization_cost=37.7,
                 carbohydrate_polymerization_cost=12.8, biomass_reaction_id='r_4041',
                 protein_reaction_id='r_4047', carbohydrate_reaction_id='r_4048',
                 protein_pool_exchange_id='prot_pool_exchange', common_protein_pool_id='prot_pool',
                 reaction_prefix=''):

        # load predifined models
        model_list = ModelList()
        if isinstance(model, string_types):
            model = model_list.__getitem__(model)
        elif isinstance(model, CBModel):
            model_list.simplify_model(model)
        else:
            raise ValueError('Model should be a string denomination or a CBModel instance')

        super(GeckoModel, self).__init__(model.id)

        # import CBModel's data
        self.compartments = copy.deepcopy(model.compartments)
        self.metabolites = copy.deepcopy(model.metabolites)
        self.reactions = copy.deepcopy(model.reactions)
        self.genes = copy.deepcopy(model.genes)

        # biomass reaction id (str)
        if biomass_reaction_id not in self.reactions:
            try:
                self.biomass_reaction = model.biomass_reaction
            except:
                self.biomass_reaction = None

        # protein reaction id (CBReaction)
        try:
            self.protein_reaction = self.reactions[protein_reaction_id]
        except KeyError:
            logging.warning(f"Protein reaction {protein_reaction_id} is not in the model")

        # carbohydrate reaction id (CBReaction)
        try:
            self.carbohydrate_reaction = self.reactions[carbohydrate_reaction_id]
        except KeyError:
            logging.warning(f"Carbohydrate reaction {carbohydrate_reaction_id} is not in the model")

        # protein properties DataFrame
        if protein_properties:
            self.protein_properties = protein_properties
        else:
            model_list = ModelList()
            self.protein_properties = model_list.protein_properties()

        # Metabolite of common protein pool
        try:
            self.common_protein_pool = self.metabolites[common_protein_pool_id]
        except KeyError:
            self.common_protein_pool = Metabolite(common_protein_pool_id)
            self.metabolites[common_protein_pool_id] = self.common_protein_pool

        # Reaction identified as protein pool exchange
        if protein_pool_exchange_id in self.reactions.keys():
            self.protein_pool_exchange = self.reactions[protein_pool_exchange_id]
        else:
            self.protein_pool_exchange = CBReaction(protein_pool_exchange_id)
            self.protein_pool_exchange.stoichiometry.update({self.common_protein_pool.id: 1.})
            self.protein_pool_exchange.set_flux_bounds(0, ModelConstants.REACTION_UPPER_BOUND)
            self.add_reaction(self.protein_pool_exchange)

        # multi-pool
        s_protein_exchange = "^" + reaction_prefix + "prot_(.*)_exchange$"
        self.protein_exchange_re = re.compile(s_protein_exchange)
        s_pool_protein_exchange = "^" + reaction_prefix + "draw_prot_(.*)$"
        self.pool_protein_exchange_re = re.compile(s_pool_protein_exchange)
        self.concentrations = pd.Series(np.nan, index=self.proteins)

        # passed parameters
        self.gam = gam
        self.amino_acid_polymerization_cost = amino_acid_polymerization_cost
        self.carbohydrate_polymerization_cost = carbohydrate_polymerization_cost
        self.c_base = c_base
        self.sigma_saturation_factor = sigma

        # initialize remaining variables
        self.measured_ggdw = None
        self.fp_fraction_protein = None
        self.fc_carbohydrate_content = None
        self.p_total = None
        self.c_total = None
        self.p_base = None
        self.fn_mass_fraction_unmeasured_matched = None
        self.fs_matched_adjusted = None
        self.p_measured = None
        self.f_mass_fraction_measured_matched_to_total = None
        self.fm_mass_fraction_matched = None
        self.uniprot = {}
        # protein reverse reactions mapping
        self._protein_rev_reactions = None

    def fraction_to_ggdw(self, fraction):
        """Convert protein measurements in mass fraction of total to g protein / g DW.

        :param fraction: pd.Series, Data of protein measurements which are absolute quantitative fractions \
            of the total amount of these measured proteins. Normalized to sum == 1.
        :returns: pd.Series, g protein / g DW for the measured proteins

        """
        # measurements should be quantitative fractions of the total measured proteins, normalized to unit-length
        fraction = fraction / fraction.sum()
        fraction_measured = self.protein_properties.loc[list(fraction.index), 'abundance'].sum()
        p_measured = self.p_total * fraction_measured
        return fraction.apply(lambda x: x * p_measured)

    def limit_proteins(self, fractions=None, ggdw=None, p_total=0.448, p_base=0.46):
        """Apply proteomics measurements to model.

        Apply measurements in the form of fractions of total of the measured proteins, or directly as g / gDW. Must
        supply exactly one of `fractions` or `ggdw`.

        :param fractions: pd.Series, Protein abundances in fraction of total (normalized to sum to 1).  \
            Ignored if `ggdw` is also supplied.
        :param ggdw: pd.Series, Protein abundances in g / gDW
        :param p_total: float, measured total protein fraction in cell in g protein / g DW. \
            Should be measured for each experiment, the default here is taken from [2]_.
        :param p_base: float, protein content at dilution rate 0.1 / h in g protein / g DW. \
            Default taken from yeast 8.1.3.

        References:

        .. [2] Benjamin J. Sanchez, Cheng Zhang, Avlant Nilsson, Petri-Jaan Lahtvee, Eduard J. Kerkhoven, Jens Nielsen (
           2017). Improving the phenotype predictions of a yeast genome-scale metabolic model by incorporating enzymatic
           constraints. [Molecular Systems Biology, 13(8): 935, http://www.dx.doi.org/10.15252/msb.20167411

        """
        self.p_total = p_total
        self.p_base = p_base
        self.c_total = self.c_base + self.p_base - self.p_total
        self.fp_fraction_protein = self.p_total / self.p_base
        self.fc_carbohydrate_content = self.c_total / self.c_base
        self.measured_ggdw = self.fraction_to_ggdw(fractions) if ggdw is None else ggdw

        # * section 2.5
        # 1. define mmmol_gdw as ub for measured proteins (multi-pool)

        for protein_id, value in iteritems(self.measured_ggdw):
            try:
                mmol_gdw = value / (self.protein_properties.loc[protein_id, 'mw'] / 1000)
                rxn = self.reactions['prot_{}_exchange'.format(protein_id)]
                self.uniprot[rxn.id] = protein_id
            except KeyError:
                pass
            else:
                self.concentrations[protein_id] = value
                rxn.set_flux_bounds(0, mmol_gdw)

        # 2. p_measured is aggregate mass of all matched proteins
        self.p_measured = self.concentrations.sum()
        # 3. fm, mass fraction of measured proteins in the model over total
        self.fm_mass_fraction_matched = self.p_measured / self.p_total
        # 4. mass fraction of unmeasured proteins in the model over all proteins not matched to model
        self.fn_mass_fraction_unmeasured_matched = (
            self.protein_properties.loc[list(self.unmeasured_proteins)].prod(axis=1).sum() /
            self.protein_properties.prod(axis=1).sum()
        )
        self.f_mass_fraction_measured_matched_to_total = (
            self.fn_mass_fraction_unmeasured_matched / (1 - self.fm_mass_fraction_matched))
        # 5. constrain unmeasured proteins by common pool
        self.constrain_pool()
        self.adjust_biomass_composition()

    def constrain_pool(self):
        """Constrain the draw reactions for the unmeasured (common protein pool) proteins.

        Proteins without their own protein pool are collectively constrained by the common protein pool. Remove
        protein pools for all proteins that don't have measurements, along with corresponding draw reactions,
        and add these to the common protein pool and reaction.

        """
        new_reactions = []
        to_remove = []
        # * section 2.5.1
        # 1. and 2. introduce `prot_pool` and exchange reaction done in __init__
        # 3. limiting total usage with the unmeasured amount of protein
        # looks like the matlab code:
        # self.fs_matched_adjusted = ((self.p_total - self.p_measured) / self.p_base *
        #                             self.f_mass_fraction_measured_matched_to_total *
        #                             self.sigma_saturation_factor)
        # but this gives results more like reported:
        self.fs_matched_adjusted = ((self.p_total - self.p_measured) *
                                    self.f_mass_fraction_measured_matched_to_total *
                                    self.sigma_saturation_factor)
        self.protein_pool_exchange.set_flux_bounds(0, self.fs_matched_adjusted)

        # 4. Remove other enzyme usage reactions and replace with pool exchange reactions
        average_mmw = self.protein_properties['mw'].mean() / 1000.
        for protein_id in self.unmeasured_proteins:
            prt = 'prot_{}_exchange'.format(protein_id)
            if prt in self.reactions:
                to_remove.append(prt)
            draw_reaction_id = 'draw_prot_{}'.format(protein_id)
            if draw_reaction_id not in self.reactions.keys():
                draw_rxn = CBReaction(draw_reaction_id)
                # defines bounds
                draw_rxn.set_flux_bounds(0, 1000)
                self.uniprot[draw_rxn.id] = protein_id
                protein_pool = self.metabolites['prot_{}_c'.format(protein_id)]
                try:
                    mmw = self.protein_properties.loc[protein_id, 'mw'] / 1000.
                except KeyError:
                    mmw = average_mmw
                metabolites = {self.common_protein_pool.id: -mmw, protein_pool.id: 1}
                draw_rxn.stoichiometry.update(metabolites)
                new_reactions.append(draw_rxn)
        for rxn in new_reactions:
            self.add_reaction(rxn)
        self.remove_reactions(to_remove)

    def adjust_biomass_composition(self):
        """Adjust the biomass composition.

        After changing the protein and carbohydrate content based on measurements, adjust the corresponding
        coefficients of the biomass reaction.

        """
        for met in self.protein_reaction.stoichiometry:
            is_prot = 'protein' in self.metabolites[met].name
            if not is_prot:
                coefficient = self.fp_fraction_protein * self.protein_reaction.stoichiometry[met]
                self.protein_reaction.stoichiometry[met] = coefficient

        for met in self.carbohydrate_reaction.stoichiometry:
            is_carb = 'carbohydrate' in self.metabolites[met].name
            if not is_carb:
                coefficient = self.fc_carbohydrate_content * self.carbohydrate_reaction.stoichiometry[met]
                self.carbohydrate_reaction.stoichiometry[met] = coefficient

        for met in self.reactions[self.biomass_reaction].stoichiometry:
            sign = -1 if self.reactions[self.biomass_reaction].stoichiometry[met] < 0 else 1
            is_atp = 'ATP' in self.metabolites[met].name
            is_adp = 'ADP' in self.metabolites[met].name
            is_h2o = 'H2O' in self.metabolites[met].name
            is_h = 'H+' in self.metabolites[met].name
            is_p = 'phosphate' in self.metabolites[met].name
            if is_atp or is_adp or is_h2o or is_h or is_p:
                coefficient = sign * (self.gam +
                                      self.amino_acid_polymerization_cost * self.p_total +
                                      self.carbohydrate_polymerization_cost * self.c_total)
                self.reactions[self.biomass_reaction].stoichiometry[met] = coefficient

    def adjust_pool_bounds(self, min_objective=0.05, inplace=False, tolerance=1e-9):
        """Adjust protein pool bounds minimally to make model feasible.

        Bounds from measurements can make the model non-viable or even infeasible. Adjust these minimally by minimizing
        the positive deviation from the measured values.

        :param min_objective: float, The minimum value of for the ojective for calling the model viable.
        :param inplace: bool, Apply the adjustments to the model.
        :param tolerance: float, Minimum non-zero value. Solver specific value.
        :returns: pd.DataFrame, Data frame with the series 'original' bounds and the new 'adjusted' bound, \
            and the optimized 'addition'.

        """
        from reframed.solvers import solver_instance
        solver = solver_instance(self)
        solver.add_constraint(
            'constraint_objective', self.get_objective, sense='>', rhs=min_objective, update=False)
        for pool in self.individual_protein_exchanges:
            solver.add_variable('pool_diff_' + pool.id, lb=0, update=False)
            solver.add_variable('measured_bound_' + pool.id,
                                lb=pool.upper_bound, ub=pool.upper_bound, update=False)
        solver.update()
        solution = solver.solve()
        # with self.model as model:
        #    problem = model.problem
        #    constraint_objective = problem.Constraint(model.objective.expression, name='constraint_objective',
        #                                              lb=min_objective)
        #    to_add = [constraint_objective]
        #    new_objective = S.Zero
        #    for pool in model.individual_protein_exchanges:
        #        ub_diff = problem.Variable('pool_diff_' + pool.id, lb=0, ub=None)
        #        current_ub = problem.Variable('measured_bound_' + pool.id, lb=pool.upper_bound, ub=pool.upper_bound)
        #        constraint = problem.Constraint(pool.forward_variable - current_ub - ub_diff, ub=0,
        #                                        name='pool_ub_' + pool.id)
        #        to_add.extend([ub_diff, current_ub, constraint])
        #        new_objective += ub_diff
        #        pool.bounds = 0, 1000.
        #    model.add_cons_vars(to_add)
        #    model.objective = problem.Objective(new_objective, direction='min')
        #    model.slim_optimize(error_value=None)
        #    primal_values = model.solver.primal_values
        # adjustments = [(pool.id, primal_values['pool_diff_' + pool.id], pool.upper_bound)
        #               for pool in model.individual_protein_exchanges
        #               if primal_values['pool_diff_' + pool.id] > tolerance]
        # result = pd.DataFrame(adjustments, columns=['reaction', 'addition', 'original'])
        # result['adjusted'] = result['addition'] + result['original']
        # if inplace:
        #    for adj in result.itertuples():
        #        model.reactions.get_by_id(adj.reaction).upper_bound = adj.adjusted
        # return result

    @property
    def measured_proteins(self):
        """Get the identifiers of the measured proteins.

        :returns: frozenset, The identifiers for the unmeasured proteins.

        """
        return frozenset(self.concentrations[self.concentrations.notnull()].index)

    @property
    def unmeasured_proteins(self):
        """Get the identifiers of the proteins .

        :returns: frozenset, The protein identifiers for the measured proteins.

        """
        return frozenset(self.concentrations[self.concentrations.isnull()].index)

    @property
    def proteins(self):
        """Get all proteins.

        :returns: frozenset, The set of all proteins identifiers.

        """
        return self.individual_proteins.union(self.pool_proteins)

    @property
    def individual_proteins(self):
        """Get the identifiers for the proteins with their individual abundance pool.

        :returns: frozenset, The set of proteins that have a defined separate pool exchange reaction.

        """
        return frozenset(chain.from_iterable(re.findall(self.protein_exchange_re, rxn)
                                             for rxn in self.protein_exchanges))

    @property
    def pool_proteins(self):
        """Get proteins modeled by common protein pool.

        :returns: frozenset, The set of proteins that have a defined draw reaction.

        """
        return frozenset(chain.from_iterable(re.findall(self.pool_protein_exchange_re, rxn)
                                             for rxn in self.protein_exchanges))

    @property
    def individual_protein_exchanges(self):
        """Individual protein-exchange reactions.

        :returns: frozenset, Set of protein exchange reactions with individual pools

        """
        return (frozenset(rxn for rxn in self.reactions.keys()
                          if re.match(self.protein_exchange_re, rxn)) - {self.protein_pool_exchange.id})

    @property
    def pool_protein_exchanges(self):
        """Protein-exchange reactions by single pool.

        :returns: frozenset, Set of protein exchange reactions for single pool reactions.

        """
        return (frozenset(rxn for rxn in self.reactions
                          if re.match(self.pool_protein_exchange_re, rxn)) -
                {self.protein_pool_exchange})

    @property
    def protein_exchanges(self):
        """Protein-exchange reactions.
            fs_matched_adjusted

        :returns: frozenset, Set of protein exchange reactions (individual and common protein pool reactions).

        """
        return self.individual_protein_exchanges.union(self.pool_protein_exchanges)

    @property
    def protein_rev_reactions(self):
        """
        Pairs of reverse reactions associated with a protein

        :returns: dictionaty, A dictionary which identifies for each protein (key) the list of reversible \
        reactions pairs.

        """
        if not self._protein_rev_reactions:
            proteins = self.proteins
            reactions = self.reactions
            in_sub = {}
            for p in proteins:
                sub = []
                for r_id, rxn in reactions.items():
                    lsub = rxn.get_substrates()
                    for m in lsub:
                        if p in m:
                            sub.append(r_id)
                in_sub[p] = sub
            pairs = {}
            for k, s in in_sub.items():
                revs = [r for r in s if '_REV' in r]
                if len(revs) > 0:
                    for r in revs:
                        lrx = [a for a in s if r.replace('_REV', '') == a]
                        lrx.append(r)
                        if len(lrx) == 2:
                            if k in pairs.keys():
                                pairs[k].append((lrx[0], lrx[1]))
                            else:
                                pairs[k] = [(lrx[0], lrx[1])]
            self._protein_rev_reactions = pairs
        return self._protein_rev_reactions
