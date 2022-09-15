import warnings
from typing import Tuple, Dict, Union, TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel

try:

    # noinspection PyPackageRequirements
    from sklearn.impute import KNNImputer
    # noinspection PyPackageRequirements
    from sklearn.preprocessing import quantile_transform

except ImportError as exc:

    warnings.warn('The package scikit-learn is not installed. '
                  'To use PROM analysis, please install scikit-learn (pip install scikit-learn).')

    raise exc

try:

    # noinspection PyPackageRequirements
    from scipy.stats import ks_2samp

except ImportError as exc:

    warnings.warn('The package scipy is not installed. '
                  'To use PROM analysis, please install scipy (pip install scipy).')

    raise exc


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


def quantile_preprocessing_pipeline(expression: pd.DataFrame,
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
    :param expression: Expression matrix
    :param missing_values: The placeholder for the missing values. All occurrences of missing_values will be imputed.
    :param n_neighbors: Number of neighboring samples to use for imputation.
    :param weights: Weight function used in prediction. Possible values:
        - 'uniform': uniform weights. All points in each neighborhood are weighted equally.
        - 'distance': weight points by the inverse of their distance. in this case, closer neighbors of a query point
    :param metric: Metric used to compute the distance between samples. The default metric is nan_euclidean, which is
    the euclidean distance ignoring missing values. Consult sklearn documentation for more information.
    :param n_quantiles: Number of quantiles to be computed. It corresponds to the number of landmarks used to discretize
    :param q: Quantile to compute
    :return: Quantile preprocessed expression matrix, quantile expression binarized matrix
    """
    index = list(expression.index)
    columns = list(expression.columns)
    expression = expression.to_numpy()

    expression = knn_imputation(expression,
                                missing_values=missing_values,
                                n_neighbors=n_neighbors,
                                weights=weights,
                                metric=metric)
    expression = quantile_transformation(expression, n_quantiles=n_quantiles)

    binary_expression = quantile_binarization(expression, q=q)

    return pd.DataFrame(expression, index, columns), pd.DataFrame(binary_expression, index, columns)


def target_regulator_interaction_probability(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                                             expression: pd.DataFrame,
                                             binary_expression: pd.DataFrame) -> Tuple[Dict[str, float],
                                                                                       Dict[str, float]]:
    """
    It computes the conditional probability of a target gene being active when the regulator is inactive.
    It uses the following formula:
        P(target = 1 | regulator = 0) = count(target = 1, regulator = 0) / # samples
    This probability is computed for each combination of target-regulator.
    This method is used in PROM analysis.

    :param model: an integrated Metabolic-Regulatory model aka a MEW model
    :param expression: Quantile preprocessed expression matrix
    :param binary_expression: Quantile preprocessed expression matrix binarized
    :return: Dictionary with the conditional probability of a target gene being active when the regulator is inactive,
    Dictionary with missed interactions
    """
    missed_interactions = {}
    interactions_probabilities = {}

    for interaction in model.yield_interactions():

        target = interaction.target

        if not interaction.regulators or target.id not in expression.index:
            missed_interactions[target.id] = 1
            interactions_probabilities[target.id] = 1
            continue

        target_expression = expression.loc[target.id]
        target_binary = binary_expression.loc[target.id]

        for regulator in interaction.yield_regulators():

            if regulator.id not in expression.index:
                missed_interactions[target.id + regulator.id] = 1
                interactions_probabilities[target.id + regulator.id] = 1
                continue

            regulator_binary = binary_expression.loc[regulator.id]

            target_expression_1_regulator = target_expression[regulator_binary == 1]
            target_expression_0_regulator = target_expression[regulator_binary == 0]

            if len(target_expression_1_regulator) == 0 and len(target_expression_0_regulator) == 0:
                missed_interactions[target.id + regulator.id] = 1
                interactions_probabilities[target.id + regulator.id] = 1
                continue

            _, p_val = ks_2samp(target_expression_1_regulator, target_expression_0_regulator)
            if p_val < 0.05:
                target_binary_0_regulator = target_binary[regulator_binary == 0]

                probability = sum(target_binary_0_regulator) / len(target_binary_0_regulator)

                interactions_probabilities[target.id + regulator.id] = probability
                missed_interactions[target.id + regulator.id] = 0

            else:
                missed_interactions[target.id + regulator.id] = 1
                interactions_probabilities[target.id + regulator.id] = 1

    if sum(missed_interactions.values()) > (0.75 * len(missed_interactions)):
        warnings.warn('Binarization threshold should be changed', Warning, stacklevel=2)

    return interactions_probabilities, missed_interactions
