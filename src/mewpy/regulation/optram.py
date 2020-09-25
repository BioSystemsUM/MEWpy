from mewpy.problems.problem import AbstractOUProblem
from mewpy.utils.parsing import Boolean, GeneEvaluator, build_tree
from mewpy.utils.constants import EAConstants
from mewpy.regulation.variable import RegulatoryVariable
from collections import OrderedDict
import numpy as np
import pandas as pd
import math


def load_optram(gene_filename, tf_filename, matrix_filename, gene_prefix=''):
    """
    Loads a OptRAM regulatory model from csv files:

        gene_file columns: "Name;id"
        tf_file columns: "Name;Expression"
        matrix_file: genes as rows and tfs as columns
    """
    df_genes = pd.read_csv(gene_filename)
    df_TFs = pd.read_csv(tf_filename)
    mat = pd.read_csv(matrix_filename, header=None)

    genes = OrderedDict()
    for index, row in df_genes.iterrows():
        genes[gene_prefix+row['Name']
              ] = RegGene(gene_prefix+row['Name'], index, row['id'])
    tfs = OrderedDict()
    for index, row in df_TFs.iterrows():
        tf = TF(row['Name'], index, row['Expression'])
        tfs[row['Name']] = tf

    model = OptRAMRegModel(genes, tfs, mat)
    return model


class RegGene(RegulatoryVariable):
    """Genes included in the regulatory model

       args:
        name (str): the gene identifier
        row (int): the associated row in the regulatory model
        id (int): OptRAM ID
        cbm_name (str): the gene corresponding name in the constraint base model (G_XXXXX)

    """

    def __init__(self, name, row, optramid, aliases=None):
        super().__init__(name=name, aliases=aliases)
        self.row = row
        self.optramid = optramid


class TF(RegulatoryVariable):
    """Transcription factor

    args:
        name (str): the TF identifier
        column (int): the associated column in the regulatory model
        expression (str): the TF expression (By default 1)
    """

    def __init__(self, name, column, expression=1):
        super().__init__(name=name)
        self.column = column
        self.expression = expression


class OptRAMRegModel:
    def __init__(self, genes, tfs, regnet):
        """
        OptRAM regulatory network model.

        :param genes: (dic) A dictionary of Gene objects.
        :param tfs: (dic) A dictionary of transcription factors (TF) objects.
        :param regnet: (DataFrame) A panda dataframe containing a matrix of Genes TFs coefficients.
        """
        self.genes = genes
        self.tfs = tfs
        # panda (genes x TFs)
        # TODO: maybe use a list of lists instead of a panda dataframe
        self.regnet = regnet
        self.tf_expression = [0]*len(self.tfs)
        for _, tf in self.tfs.items():
            self.tf_expression[tf.column] = tf.expression


class OptRamProblem(AbstractOUProblem):

    def __init__(self, model, fevaluation, regmodel, **kwargs):
        """
        EA OptRam problem
        """
        super(OptRamProblem, self).__init__(model, fevaluation, **kwargs)
        self.regmodel = regmodel
        # ignore user defined target list
        self._trg_list = None
        # GPR operators
        self._operators = None
        # Reset default OU levels to OptRAM levels if none are provided
        self.levels = kwargs.get('levels', EAConstants.OPTRAM_LEVELS)

    def _build_target_list(self):
        """ The EA target list is the combination [mGene]+[TFs]
        """
        self._trg_list = []
        self._trg_list.extend(list(self.regmodel.genes.keys()))
        self._trg_list.extend(list(self.regmodel.tfs.keys()))

    def target_exp(self, candidate):

        mgenes_p = {}
        # TFs expression vector
        tf_exp_v = self.regmodel.tf_expression.copy()
        # keeps track of the TF whose expression is altered by the candidate
        tf_altered = [False] * len(tf_exp_v)
        for idx, lv_idx in candidate:
            try:
                target = self.target_list[idx]
                lv = self.levels[lv_idx]

                if idx >= len(self.regmodel.genes):
                    # TF expression level encoded in the candidate
                    tf = self.regmodel.tfs[target]
                    c = tf.column
                    tf_exp_v[c] *= lv
                    tf_altered[c] = True
                else:
                    # gene expression level encoded in the candidate
                    #
                    mgenes_p[target] = lv
            except IndexError:
                raise IndexError("Index out of range")

        # update gene expression level with TFs
        # gene expression is p=2^sum(coeff*log2 tfexpr)

        log_tfexpr = [math.log2(x) for x in tf_exp_v]

        for g in self.regmodel.genes:
            g_idx = self.regmodel.genes[g].row
            coeff = abs(self.regmodel.regnet.iloc[g_idx].to_numpy())
            # only genes with altered TFs have their expression level updated
            if np.dot(coeff != 0, tf_altered):
                p = 2 ** np.dot(coeff, log_tfexpr)
                mgenes_p[g] = p
        return mgenes_p

    def decode(self, candidate):
        gr_constraints = OrderedDict()
        mgenes_p = self.target_exp(candidate)
        # Evaluate gpr.
        if not self._operators:
            self._operators = (lambda x, y: min(x, y), lambda x, y: max(x, y))
        evaluator = GeneEvaluator(
            mgenes_p, self._operators[0], self._operators[1])
        for rxn_id in self.simulator.reactions:
            if self.simulator.get_gpr(rxn_id):
                gpr = str(self.simulator.get_gpr(rxn_id))
                tree = build_tree(gpr, Boolean)
                # Apply the operators to obtain a level for the reaction.
                # If a gene in the gpr has no associated level,
                # its factor is 1 (see mewpy.util.parsing.GeneEvaluator)
                lv = tree.evaluate(evaluator.f_operand, evaluator.f_operator)
                # Adds the reaction constraint.
                rev_rxn = self.simul_context.reverse_reaction(rxn_id)
                # Skips if the reverse reaction was already processed.
                if rev_rxn and rev_rxn in gr_constraints.keys():
                    continue
                elif lv < 0:
                    raise ValueError("All UO levels should be positive")
                elif lv == 1:
                    # No constraints are added.
                    continue
                else:
                    gr_constraints.update(
                        self.reaction_constraints(rxn_id, lv))
        return gr_constraints
