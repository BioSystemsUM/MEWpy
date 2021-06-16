import numpy as np
import pandas as pd
from itertools import combinations
from ..simulation import get_simulator
from ..simulation.simulation import Simulator
from ..util.parsing import Boolean, GeneEvaluator, build_tree


class ExpressionSet(object):

    def __init__(self, identifiers: list, conditions: list,
                 expression: np.array, p_values: np.array = None):
        """Expression set. The expression values are a numpy array with shape
        (len(identifiers) x len(conditions)).

        Args:
            identifiers (list): Gene or Proteins identifiers
            conditions (list): Time, experiment,... identifiers.
            expression (np.array): expression values.
            p_values (np.array, optional): p-values. Defaults to None.
        """
        n = len(identifiers)
        m = len(conditions)
        if expression.shape != (n, m):
            raise ValueError(
                f"The shape of the expression {expression.shape} does not "
                f"match the expression and conditions sizes ({n},{m})")

        self._identifiers = identifiers
        self._identifier_index = {iden: idx for idx, iden in
                                  enumerate(identifiers)}
        self._conditions = [str(x) for x in conditions]
        self._condition_index = {cond: idx for idx, cond in
                                 enumerate(self._conditions)}
        self._expression = expression
        self._p_values = p_values

    def shape(self):
        """Returns:
            (tuple): the Expression dataset shape
        """
        return self._expression.shape

    def __getitem__(self, item):
        """
        Index the ExpressionSet.
        """
        return self._expression.__getitem__(item)

    def get_condition(self, condition=None, **kwargs):
        """

        Args:
            condition ([type], optional): [description]. Defaults to None.

        Returns:
            [type]: [description]
        """
        if condition is None:
            values = self[:, :]
        elif isinstance(condition, int):
            values = self[:, condition]
        elif isinstance(condition, str):
            values = self[:, self._condition_index[condition]]

        # format
        format = kwargs.get('format', None)
        if format and condition:
            if format == 'list':
                return values.tolist()
            elif format == 'dict':
                return dict(zip(self._identifiers, values.tolist()))
            else:
                return values
        else:
            return values

    @classmethod
    def from_data_frame(cls, data_frame):
        """Read expression data from a pandas.DataFrame.

        Args:
            data_frame (Dataframe): The expression Dataframe

        Returns:
            ExpressionSet: the expression dataset from the dataframe.
        """
        columns = data_frame.columns.tolist()
        conditions = [c for c in columns if "p-value" not in c]
        p_value_keys = [c for c in columns if "p-value" in c]
        if p_value_keys:
            p_values = data_frame[p_value_keys].values
        else:
            p_values = None

        expression = data_frame[conditions].values
        identifiers = data_frame.index.tolist()
        return ExpressionSet(identifiers, conditions, expression, p_values)

    @classmethod
    def from_csv(cls, file_path: str, **kwargs):
        """Read expression data from a comma separated values (csv) file.

        Args:
            file_path (str): the csv file path.

        Returns:
            ExpressionSet: the expression dataset from the csv file.
        """
        data = pd.read_csv(file_path, **kwargs)
        return cls.from_data_frame(data)

    @property
    def data_frame(self):
        """Build a pandas.DataFrame from the ExpressionProfile.
        Columns headers are conditions and
        line indexes identifiers (genes/proteins)
        """

        if self._p_values is None:
            expression = self._expression
            conditions = self._conditions
        else:
            expression = np.concatenate((self._expression, self.p_values),
                                        axis=1)
            conditions = self._conditions + self.p_value_columns

        return pd.DataFrame(expression,
                            index=self._identifiers,
                            columns=conditions)

    @property
    def p_value_columns(self):
        """ Generate the p-value column names."""
        return [f"{c[0]} {c[1]} p-value"
                for c in combinations(self._conditions, 2)]

    @property
    def p_values(self):
        """Returns the numpy array of p-values.

        Raises:
            ValueError: [description]
        """
        if not self._p_values.all():
            raise ValueError("No p-values defined.")
        else:
            return self._p_values

    @p_values.setter
    def p_values(self, p_values: np.array):
        """Sets p-values

        Args:
            p_values (np.array): [description]

        Raises:
            ValueError: [description]
        """
        if p_values is not None:
            if p_values.shape[1] != len(self.p_value_columns):
                raise ValueError("p-values do not cover all conditions")

        self._p_values = p_values

    @p_values.deleter
    def p_values(self):
        """Delete p_values."""
        self._p_values = None

    def differences(self, p_value=0.005):
        """Calculate the differences based on the MADE method.

        Args:
            p_value (float, optional): [description]. Defaults to 0.005.

        Returns:
            dict: A dictionary of differences
        """

        diff = {}
        for idx, iden in enumerate(self._identifiers):
            diff[iden] = []
            for i in range(1, len(self._conditions)):
                start, end = self._expression[idx, i-1: i+1]
                p_val = self.p_values[idx, i-1]
                if p_val <= p_value:
                    if start < end:
                        diff[iden].append(+1)
                    elif start > end:
                        diff[iden].append(-1)
                    else:
                        diff[iden].append(0)
                else:
                    diff[iden].append(0)
        return diff

    def minmax(self, condition=None):
        """ Return the min and max values for the specified condition.

        Args:
            condition (str): str or int or None, optional (default None)
            The condition to obtain the min and max values for.

        Returns
        -------
        tuple of (min, max)

        """
        values = self.get_condition(condition)
        return np.amin(values), np.amax(values)


def gene_to_reaction_expression(model, gene_exp, and_func=min, or_func=max):
    """Process reaction level from GPRs

    Args:
        model: A model or a MEWpy Simulation
        gene_exp (list): expresion values
        and_func ([type], optional): Function for AND. Defaults to min.
        or_func ([type], optional): Function for OR. Defaults to max.

    Returns:
        dict: Reaction levels
    """
    if isinstance(model, Simulator):
        sim = model
    else:
        sim = get_simulator(model)

    rxn_exp = {}
    evaluator = GeneEvaluator(gene_exp, and_func, or_func)
    for rxn_id in sim.reactions:
        gpr = sim.get_gpr(rxn_id)
        if gpr:
            tree = build_tree(gpr, Boolean)
            lv = tree.evaluate(evaluator.f_operand, evaluator.f_operator)
            rxn_exp[rxn_id] = lv
    return rxn_exp


class Preprocessing(object):
    """Formulation and implementation of preprocessing decisions.
        (A) Types of gene mapping methods
        (B) Types of thresholding approaches (global and local).
        (C) Formulation of combinations of number of states (Global, Local)
        (D) Decisions about the order in which thresholding and gene mapping
        are performed.
        For Order 1, gene expression is converted to reaction activity followed
        by thresholding of reaction activity;
        For Order 2, thresholding ofgene expression is followed by its
        conversion to reaction activity.

        [1]Anne Richelle,Chintan Joshi,Nathan E. Lewis, Assessing key decisions
        for transcriptomic data integration in biochemical networks, PLOS, 2019
        https://doi.org/10.1371/journal.pcbi.1007185
    """

    def __init__(self, model: Simulator, data: ExpressionSet, **kwargs):
        """[summary]

        Args:
            model ([type]): [description]
            data ([type]): [description]
        """
        self.model = model
        self.data = data
        self._conf = kwargs

    def reactions_expression(self, condition, and_func=None, or_func=None):
        exp = self.data.get_condition(condition, format='dict')
        and_func = self._conf.get(
            'and_func', min) if and_func is None else and_func
        or_func = self._conf.get(
            'or_func', max) if or_func is None else or_func
        rxn_exp = gene_to_reaction_expression(
            self.model, exp, and_func, or_func)
        return rxn_exp

    def percentile(self, condition=None, cutoff=0.25):
        """Processes a percentil threshold and returns the respective
        reaction coefficients, ie, a dictionary of reaction:coeff

        Args:
            condition ([type], optional): [description]. Defaults to None.
            cutoff (float, optional): [description]. Defaults to 0.25.

        Returns:
            dict, float: the coefficients and threshold
        """
        rxn_exp = self.reactions_expression(condition)
        threshold = np.percentile(list(rxn_exp.values()), cutoff)
        coeffs = {r_id: threshold-val for r_id,
                  val in rxn_exp.items() if val < threshold}
        return coeffs, threshold
