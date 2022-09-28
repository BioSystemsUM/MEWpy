from typing import Tuple

import numpy as np
import pandas as pd
from itertools import combinations

from pandas._typing import FilePathOrBuffer

from ..simulation import get_simulator
from ..simulation.simulation import Simulator
from ..util.parsing import Boolean, GeneEvaluator, build_tree


class ExpressionSet:

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

        if isinstance(condition, int):
            values = self[:, condition]
        elif isinstance(condition, str):
            values = self[:, self._condition_index[condition]]
        else:
            values = self[:, :]

        # format
        form = kwargs.get('format', 'dict')
        if form and condition is not None:
            if form == 'list':
                return values.tolist()
            elif form == 'dict':
                return dict(zip(self._identifiers, values.tolist()))
            else:
                return values
        else:
            return values

    @classmethod
    def from_dataframe(cls, data_frame):
        """Read expression data from a pandas.DataFrame.

        Args:
            data_frame (Dataframe): The expression Dataframe

        Returns:
            ExpressionSet: the expression dataset from the dataframe.
        """
        # pandas columns might be integers only wich will cause problems in the p-value verification step
        columns = [str(x) for x in data_frame.columns]
        data_frame.columns = columns

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
    def from_csv(cls, file_path: FilePathOrBuffer, **kwargs):
        """Read expression data from a comma separated values (csv) file.

        Args:
            file_path (str): the csv file path.

        Returns:
            ExpressionSet: the expression dataset from the csv file.
        """
        data = pd.read_csv(file_path, **kwargs)
        return cls.from_dataframe(data)

    @property
    def dataframe(self):
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
                start, end = self._expression[idx, i - 1: i + 1]
                p_val = self.p_values[idx, i - 1]
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

    def apply(self, function: None):
        """Apply a function to all expression values.

        :param function: the unary function to be applyied. Default log base 2.
        :type function: callable
        """
        if function is None:
            import math
            def function(x): return math.log(x, 2)
        f = np.vectorize(function)
        self._expression = f(self._expression)

    def quantile_pipeline(self,
                          missing_values: float = None,
                          n_neighbors: int = 5,
                          weights: str = 'uniform',
                          metric: str = 'nan_euclidean',
                          n_quantiles: int = None,
                          q: float = 0.33) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Quantile preprocessing pipeline. It performs the following steps:
            1. KNN imputation of missing values
            2. Quantile transformation
            3. Quantile binarization
        :param missing_values: The placeholder for the missing values.
        All occurrences of missing_values will be imputed.
        :param n_neighbors: Number of neighboring samples to use for imputation.
        :param weights: Weight function used in prediction. Possible values:
            - 'uniform': uniform weights. All points in each neighborhood are weighted equally.
            - 'distance': weight points by the inverse of their distance.
            in this case, closer neighbors of a query point
        :param metric: Metric used to compute the distance between samples. The default metric is nan_euclidean,
        which is the euclidean distance ignoring missing values. Consult sklearn documentation for more information.
        :param n_quantiles: Number of quantiles to be computed.
        It corresponds to the number of landmarks used to discretize
        :param q: Quantile to compute
        :return: Quantile preprocessed expression matrix, quantile expression binarized matrix
        """
        expression = knn_imputation(self._expression,
                                    missing_values=missing_values,
                                    n_neighbors=n_neighbors,
                                    weights=weights,
                                    metric=metric)
        expression = quantile_transformation(expression, n_quantiles=n_quantiles)

        binary_expression = quantile_binarization(expression, q=q)

        expression = pd.DataFrame(expression, self._identifiers, self._conditions)
        binary_expression = pd.DataFrame(binary_expression, self._identifiers, self._conditions)
        return expression, binary_expression


# ----------------------------------------------------------------------------------------------------------------------
# Preprocessing using KNNImputer and Quantile transformation/binarization
# ----------------------------------------------------------------------------------------------------------------------
def knn_imputation(expression: np.ndarray,
                   missing_values: float = None,
                   n_neighbors: int = 5,
                   weights: str = 'uniform',
                   metric: str = 'nan_euclidean') -> np.ndarray:
    """
    KNN imputation of missing values in the expression matrix. It uses the scikit-learn KNNImputer (Consult sklearn
    documentation for more information).
    The default metric is nan_euclidean, which is the euclidean distance ignoring missing values.

    :param expression: Expression matrix
    :param missing_values: The placeholder for the missing values. All occurrences of missing_values will be imputed.
    :param n_neighbors: Number of neighboring samples to use for imputation.
    :param weights: Weight function used in prediction. Possible values:
        - 'uniform': uniform weights. All points in each neighborhood are weighted equally.
        - 'distance': weight points by the inverse of their distance. in this case, closer neighbors of a query point
    :param metric: Metric used to compute the distance between samples. The default metric is nan_euclidean, which is
    the euclidean distance ignoring missing values. Consult sklearn documentation for more information.
    :return: Imputed expression matrix
    """
    try:
        # noinspection PyPackageRequirements
        from sklearn.impute import KNNImputer
    except ImportError:
        raise ImportError('The package scikit-learn is not installed. '
                          'To preprocess gene expression data, please install scikit-learn (pip install scikit-learn).')

    if missing_values is None:
        missing_values = np.nan

    if n_neighbors is None:
        n_neighbors = 5

    if weights is None:
        weights = 'uniform'

    if metric is None:
        metric = 'nan_euclidean'

    imputation = KNNImputer(missing_values=missing_values,
                            n_neighbors=n_neighbors,
                            weights=weights,
                            metric=metric)

    return imputation.fit_transform(expression)


def quantile_transformation(expression: np.ndarray, n_quantiles: int = None) -> np.ndarray:
    """
    Quantile transformation of the expression matrix. It uses the scikit-learn quantile_transform (Consult sklearn
    documentation for more information).
    :param expression: Expression matrix
    :param n_quantiles: Number of quantiles to be computed. It corresponds to the number of landmarks used to discretize
    :return: Quantile transformed expression matrix
    """
    try:
        # noinspection PyPackageRequirements
        from sklearn.preprocessing import quantile_transform
    except ImportError:
        raise ImportError('The package scikit-learn is not installed. '
                          'To preprocess gene expression data, please install scikit-learn (pip install scikit-learn).')

    if n_quantiles is None:
        n_quantiles = expression.shape[1]

    return quantile_transform(expression, n_quantiles=n_quantiles, axis=0, random_state=0)


def quantile_binarization(expression: np.ndarray, q: float = 0.33) -> np.ndarray:
    """
    It computes the q-th quantile of the expression matrix using np.quantile (consult numpy documentation for more
    information). Then, it binarizes the expression matrix using the threshold computed.
    :param expression: Expression matrix
    :param q: Quantile to compute
    :return: Binarized expression matrix
    """
    threshold = np.quantile(expression, q)

    threshold_mask = expression >= threshold
    expression[threshold_mask] = 1
    expression[~threshold_mask] = 0
    return expression


def gene_to_reaction_expression(model, gene_exp, and_func=min, or_func=max):
    """Process reaction level from GPRs

    Args:
        model: A model or a MEWpy Simulation
        gene_exp (dict): gene identifiers and expression values
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
    evaluator = GeneEvaluator(gene_exp, and_func, or_func, unexpressed_value=None)
    for rxn_id in sim.reactions:
        gpr = sim.get_gpr(rxn_id)
        if gpr:
            tree = build_tree(gpr, Boolean)
            op_set = tree.get_operands().intersection(set(gene_exp.keys()))
            if len(op_set) == 0:
                lv = None
            else:
                lv = tree.evaluate(evaluator.f_operand, evaluator.f_operator)
            rxn_exp[rxn_id] = lv
    return rxn_exp


class Preprocessing:
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
            model (Simulator): [description]
            data (ExpressionSet): [description]
            and_func (function): (optional)
            or_func (function): (optional)
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
        # Removes None if maybe is none to evaluate GPRs
        res = {k: v for k, v in rxn_exp.items() if v is not None}
        return res

    def percentile(self, condition=None, cutoff=25):
        """Processes a percentil threshold and returns the respective
        reaction coefficients, ie, a dictionary of reaction:coeff

        Args:
            condition ([type], optional): [description]. Defaults to None.
            cutoff (int, optional): [description]. Defaults to 25.

        Returns:
            dict, float: the coefficients and threshold
        """
        if type(cutoff) is tuple:
            coef = []
            thre = []
            for cut in cutoff:
                rxn_exp = self.reactions_expression(condition)
                threshold = np.percentile(list(rxn_exp.values()), cut)
                coeffs = {r_id: threshold - val for r_id,
                                                    val in rxn_exp.items() if val < threshold}
                coef.append(coeffs)
                thre.append(threshold)
            coeffs = tuple(coef)
            threshold = tuple(thre)
        else:
            rxn_exp = self.reactions_expression(condition)
            threshold = np.percentile(list(rxn_exp.values()), cutoff)
            coeffs = {r_id: threshold - val for r_id,
                                                val in rxn_exp.items() if val < threshold}
        return coeffs, threshold
