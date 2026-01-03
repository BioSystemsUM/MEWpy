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
Expression set for omics data.

Author: Vitor Pereira
Contributors: Paulo Carvalhais
              Fernando Cruz
##############################################################################
"""

from itertools import combinations
from typing import Callable, Optional, Tuple, Union

import numpy as np
import pandas as pd

from mewpy.simulation import Simulator, get_simulator
from mewpy.util.parsing import Boolean, GeneEvaluator, build_tree


class ExpressionSet:

    def __init__(self, identifiers: list, conditions: list, expression: np.array, p_values: np.array = None):
        """Expression set. The expression values are a numpy array with shape
        (len(identifiers) x len(conditions)).

        Args:
            identifiers (list): Gene or Proteins identifiers
            conditions (list): Time, experiment,... identifiers.
            expression (np.array): expression values.
            p_values (np.array, optional): p-values. Defaults to None.

        Raises:
            ValueError: If identifiers or conditions are empty.
            ValueError: If identifiers or conditions contain duplicates.
            ValueError: If expression shape doesn't match identifiers/conditions.
            ValueError: If p_values shape is invalid.
        """
        # Validate non-empty
        if not identifiers:
            raise ValueError("Identifiers cannot be empty")
        if not conditions:
            raise ValueError("Conditions cannot be empty")

        # Check for duplicates
        if len(identifiers) != len(set(identifiers)):
            raise ValueError("Duplicate identifiers found")

        # Convert conditions to strings and check for duplicates
        str_conditions = [str(x) for x in conditions]
        if len(str_conditions) != len(set(str_conditions)):
            raise ValueError("Duplicate conditions found")

        # Validate expression shape
        n = len(identifiers)
        m = len(str_conditions)
        if expression.shape != (n, m):
            raise ValueError(
                f"The shape of the expression {expression.shape} does not "
                f"match the identifiers and conditions sizes ({n},{m})"
            )

        # Validate p_values shape if provided
        if p_values is not None:
            # p_values should have shape (n, C(m, 2)) where C(m, 2) is number of condition pairs
            expected_p_cols = len(list(combinations(str_conditions, 2)))
            if p_values.shape != (n, expected_p_cols):
                raise ValueError(
                    f"p_values shape {p_values.shape} doesn't match expected "
                    f"({n}, {expected_p_cols}) for {m} conditions"
                )

        self._identifiers = identifiers
        self._identifier_index = {iden: idx for idx, iden in enumerate(identifiers)}
        self._conditions = str_conditions
        self._condition_index = {cond: idx for idx, cond in enumerate(self._conditions)}
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

    def get_condition(self, condition: Union[int, str] = None, format: str = None):
        """Retrieves the omics data for a specific condition.

        :param condition: Condition identifier (int index or str name).
                         If None, returns all data.
        :type condition: Union[int, str], optional
        :param format: Output format: "dict", "list", or None for numpy array.
                      Format is only applied for single conditions.
        :type format: str, optional
        :return: Expression values in requested format
        :raises ValueError: If condition identifier is not found
        """
        # Get values based on condition type
        if isinstance(condition, int):
            if condition < 0 or condition >= len(self._conditions):
                raise ValueError(
                    f"Condition index {condition} out of range. " f"Valid range: 0-{len(self._conditions)-1}"
                )
            values = self[:, condition]
        elif isinstance(condition, str):
            if condition not in self._condition_index:
                raise ValueError(f"Unknown condition: '{condition}'. " f"Available conditions: {self._conditions}")
            values = self[:, self._condition_index[condition]]
        else:
            # Return all data
            values = self[:, :]

        # Apply format conversion only for single conditions
        if format and condition is not None:
            if format == "list":
                return values.tolist()
            elif format == "dict":
                return dict(zip(self._identifiers, values.tolist()))

        return values

    @classmethod
    def from_dataframe(cls, data_frame, duplicates='suffix'):
        """Read expression data from a pandas.DataFrame.

        Args:
            data_frame (Dataframe): The expression Dataframe
            duplicates (str): How to handle duplicate identifiers. Options:
                - 'suffix' (default): Keep all rows, rename duplicates with numeric suffixes (_2, _3, etc.)
                - 'error': Raise ValueError if duplicates found
                - 'first': Keep first occurrence of each duplicate
                - 'last': Keep last occurrence of each duplicate
                - 'mean': Average values for duplicate identifiers
                - 'sum': Sum values for duplicate identifiers

        Returns:
            ExpressionSet: the expression dataset from the dataframe.

        Raises:
            ValueError: If duplicates parameter has invalid value
        """
        import warnings

        # Validate duplicates parameter
        valid_strategies = {'error', 'first', 'last', 'mean', 'sum', 'suffix'}
        if duplicates not in valid_strategies:
            raise ValueError(
                f"Invalid duplicates parameter: '{duplicates}'. "
                f"Must be one of: {valid_strategies}"
            )

        # Handle duplicate identifiers based on strategy
        if data_frame.index.duplicated().any():
            if duplicates == 'error':
                # Keep original error behavior
                pass  # Will be caught by ExpressionSet.__init__
            elif duplicates == 'first':
                data_frame = data_frame[~data_frame.index.duplicated(keep='first')]
            elif duplicates == 'last':
                data_frame = data_frame[~data_frame.index.duplicated(keep='last')]
            elif duplicates == 'mean':
                data_frame = data_frame.groupby(data_frame.index).mean()
            elif duplicates == 'sum':
                data_frame = data_frame.groupby(data_frame.index).sum()
            elif duplicates == 'suffix':
                # Count duplicate identifiers and warn user
                duplicate_mask = data_frame.index.duplicated(keep=False)
                n_duplicates = duplicate_mask.sum()
                unique_duplicates = data_frame.index[duplicate_mask].unique()

                warnings.warn(
                    f"Found {n_duplicates} duplicate rows for {len(unique_duplicates)} unique identifiers. "
                    f"Renaming with numeric suffixes (_2, _3, etc.). "
                    f"Duplicate identifiers: {list(unique_duplicates[:10])}"
                    + (f" and {len(unique_duplicates)-10} more..." if len(unique_duplicates) > 10 else ""),
                    UserWarning
                )

                # Create new index with suffixes for duplicates
                new_index = []
                seen = {}
                for idx in data_frame.index:
                    if idx in seen:
                        seen[idx] += 1
                        new_index.append(f"{idx}_{seen[idx]}")
                    else:
                        seen[idx] = 1
                        new_index.append(idx)

                data_frame.index = new_index

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
    def from_csv(cls, file_path, duplicates='suffix', **kwargs):
        """Read expression data from a comma separated values (csv) file.

        Args:
            file_path (str): the csv file path.
            duplicates (str): How to handle duplicate identifiers. Options:
                - 'suffix' (default): Keep all rows, rename duplicates with numeric suffixes (_2, _3, etc.)
                - 'error': Raise ValueError if duplicates found
                - 'first': Keep first occurrence of each duplicate
                - 'last': Keep last occurrence of each duplicate
                - 'mean': Average values for duplicate identifiers
                - 'sum': Sum values for duplicate identifiers
            **kwargs: Additional arguments passed to pandas.read_csv()

        Returns:
            ExpressionSet: the expression dataset from the csv file.

        Raises:
            ValueError: If duplicates parameter has invalid value
        """
        data = pd.read_csv(file_path, **kwargs)
        return cls.from_dataframe(data, duplicates=duplicates)

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
            expression = np.concatenate((self._expression, self.p_values), axis=1)
            conditions = self._conditions + self.p_value_columns

        return pd.DataFrame(expression, index=self._identifiers, columns=conditions)

    @property
    def p_value_columns(self):
        """Generate the p-value column names."""
        return [f"{c[0]} {c[1]} p-value" for c in combinations(self._conditions, 2)]

    @property
    def p_values(self):
        """Returns the numpy array of p-values.

        Raises:
            ValueError: If p-values are not defined.
        """
        if self._p_values is None:
            raise ValueError("No p-values defined.")
        return self._p_values

    @p_values.setter
    def p_values(self, p_values: np.array):
        """Sets p-values array.

        Args:
            p_values (np.array): Numpy array of p-values with shape (n_identifiers, n_condition_pairs).

        Raises:
            ValueError: If p_values shape doesn't match expected dimensions for all condition pairs.
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
            p_value (float, optional): Significance threshold for p-values. Defaults to 0.005.

        Returns:
            dict: A dictionary of differences
        """

        diff = {}
        for idx, iden in enumerate(self._identifiers):
            diff[iden] = []
            for i in range(1, len(self._conditions)):
                start, end = self._expression[idx, i - 1 : i + 1]
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
        """Return the min and max values for the specified condition.

        Args:
            condition (str): str or int or None, optional (default None)
            The condition to obtain the min and max values for.

        Returns
        -------
        tuple of (min, max)

        """
        values = self.get_condition(condition)
        return np.amin(values), np.amax(values)

    def apply(self, function: Optional[Callable[[float], float]] = None):
        """Apply a function to all expression values.

        :param function: Unary function to apply to each element. Defaults to log base 2.
        :type function: Optional[Callable[[float], float]]
        """
        if function is None:
            import math

            def function(x):
                return math.log(x, 2)

        f = np.vectorize(function)
        self._expression = f(self._expression)

    def quantile_pipeline(
        self,
        missing_values: float = None,
        n_neighbors: int = 5,
        weights: str = "uniform",
        metric: str = "nan_euclidean",
        n_quantiles: int = None,
        q: float = 0.33,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
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
        expression = knn_imputation(
            self._expression, missing_values=missing_values, n_neighbors=n_neighbors, weights=weights, metric=metric
        )
        expression = quantile_transformation(expression, n_quantiles=n_quantiles)

        binary_expression = quantile_binarization(expression, q=q)

        expression = pd.DataFrame(expression, self._identifiers, self._conditions)
        binary_expression = pd.DataFrame(binary_expression, self._identifiers, self._conditions)

        return expression, binary_expression


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
    For Order 2, thresholding of gene expression is followed by its
    conversion to reaction activity.

    [1]Anne Richelle,Chintan Joshi,Nathan E. Lewis, Assessing key decisions
    for transcriptomic data integration in biochemical networks, PLOS, 2019
    https://doi.org/10.1371/journal.pcbi.1007185
    """

    def __init__(self, model: Simulator, data: ExpressionSet, **kwargs):
        """Initialize Preprocessing with model and expression data.

        Args:
            model (Simulator): Metabolic model simulator instance.
            data (ExpressionSet): Gene expression data set.
            and_func (function): (optional) Function for AND operation in GPR evaluation.
            or_func (function): (optional) Function for OR operation in GPR evaluation.
        """
        self.model = model
        self.data = data
        self._conf = kwargs

    def reactions_expression(self, condition, and_func=None, or_func=None):
        exp = self.data.get_condition(condition, format="dict")
        and_func = self._conf.get("and_func", min) if and_func is None else and_func
        or_func = self._conf.get("or_func", max) if or_func is None else or_func
        rxn_exp = gene_to_reaction_expression(self.model, exp, and_func, or_func)
        # Removes None if maybe is none to evaluate GPRs
        res = {k: v for k, v in rxn_exp.items() if v is not None}
        return res

    def percentile(self, condition=None, cutoff=25):
        """Processes a percentile threshold and returns the respective
        reaction coefficients, ie, a dictionary of reaction:coeff

        Args:
            condition: The condition identifier. Defaults to None.
            cutoff: Percentile cutoff(s) - int or tuple of ints. Defaults to 25.

        Returns:
            tuple: (coefficients, threshold) where coefficients is dict or tuple of dicts,
                   and threshold is float or tuple of floats
        """
        # Compute reaction expression once (optimization for tuple cutoffs)
        rxn_exp = self.reactions_expression(condition)
        rxn_values = list(rxn_exp.values())

        if isinstance(cutoff, tuple):
            coef = []
            thre = []
            for cut in cutoff:
                threshold = np.percentile(rxn_values, cut)
                coeffs = {r_id: threshold - val for r_id, val in rxn_exp.items() if val < threshold}
                coef.append(coeffs)
                thre.append(threshold)
            coeffs = tuple(coef)
            threshold = tuple(thre)
        else:
            threshold = np.percentile(rxn_values, cutoff)
            coeffs = {r_id: threshold - val for r_id, val in rxn_exp.items() if val < threshold}
        return coeffs, threshold


# ----------------------------------------------------------------------------------------------------------------------
# Preprocessing using KNNImputer and Quantile transformation/binarization
# ----------------------------------------------------------------------------------------------------------------------
def knn_imputation(
    expression: np.ndarray,
    missing_values: float = np.nan,
    n_neighbors: int = 5,
    weights: str = "uniform",
    metric: str = "nan_euclidean",
) -> np.ndarray:
    """
    KNN imputation of missing values in the expression matrix using scikit-learn KNNImputer.

    The default metric is nan_euclidean, which is the euclidean distance ignoring missing values.

    :param expression: Expression matrix
    :param missing_values: Placeholder for missing values to impute. Defaults to np.nan.
    :param n_neighbors: Number of neighboring samples to use for imputation. Defaults to 5.
    :param weights: Weight function used in prediction. Possible values:
        - 'uniform': uniform weights. All points in each neighborhood are weighted equally.
        - 'distance': weight points by the inverse of their distance.
    :param metric: Metric to compute distance between samples. Defaults to 'nan_euclidean'.
    :return: Imputed expression matrix
    """
    try:
        # noinspection PyPackageRequirements
        from sklearn.impute import KNNImputer
    except ImportError:
        raise ImportError(
            "The package scikit-learn is not installed. "
            "To preprocess gene expression data, please install scikit-learn (pip install scikit-learn)."
        )

    imputation = KNNImputer(missing_values=missing_values, n_neighbors=n_neighbors, weights=weights, metric=metric)

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
        raise ImportError(
            "The package scikit-learn is not installed. "
            "To preprocess gene expression data, please install scikit-learn (pip install scikit-learn)."
        )

    if n_quantiles is None:
        n_quantiles = expression.shape[1]

    return quantile_transform(expression, n_quantiles=n_quantiles, axis=0, random_state=0)


def quantile_binarization(expression: np.ndarray, q: float = 0.33) -> np.ndarray:
    """
    Computes the q-th quantile of the expression matrix and binarizes it using the threshold.

    The input array is NOT modified - a new binarized array is returned.

    :param expression: Expression matrix (will not be modified)
    :param q: Quantile to compute (default: 0.33)
    :return: New binarized expression matrix (0s and 1s)
    """
    threshold = np.quantile(expression, q)

    # Create a copy to avoid mutating input
    binary_expression = expression.copy()

    threshold_mask = binary_expression >= threshold
    binary_expression[threshold_mask] = 1
    binary_expression[~threshold_mask] = 0
    return binary_expression
