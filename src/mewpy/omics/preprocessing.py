import warnings
from typing import Tuple, Dict, Union, TYPE_CHECKING, List

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from mewpy.mew.variables import Gene, Target
    from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel

try:

    # noinspection PyPackageRequirements
    from sklearn.impute import KNNImputer
    # noinspection PyPackageRequirements
    from sklearn.preprocessing import quantile_transform
    # noinspection PyPackageRequirements
    from sklearn.linear_model import LinearRegression

except ImportError as exc:

    warnings.warn('The package scikit-learn is not installed. '
                  'To preprocess gene expression data, please install scikit-learn (pip install scikit-learn).')

    raise exc

try:

    # noinspection PyPackageRequirements
    from scipy.stats import ks_2samp

except ImportError as exc:

    warnings.warn('The package scipy is not installed. '
                  'To preprocess gene expression data, please install scipy (pip install scipy).')

    raise exc


# ----------------------------------------------------------------------------------------------------------------------
# Preprocessing using KNNImputer and Quantile transformation/binarization
# Useful for PROM method
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


# ----------------------------------------------------------------------------------------------------------------------
# Probability of Target-Regulator interactions
# Useful for PROM method
# ----------------------------------------------------------------------------------------------------------------------
def target_regulator_interaction_probability(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                                             expression: pd.DataFrame,
                                             binary_expression: pd.DataFrame) -> Tuple[Dict[Tuple[str, str], float],
                                                                                       Dict[Tuple[str, str], float]]:
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
            missed_interactions[(target.id, target.id)] = 1
            interactions_probabilities[(target.id, target.id)] = 1
            continue

        target_expression = expression.loc[target.id]
        target_binary = binary_expression.loc[target.id]

        for regulator in interaction.yield_regulators():

            if regulator.id not in expression.index:
                missed_interactions[(target.id, regulator.id)] = 1
                interactions_probabilities[(target.id, regulator.id)] = 1
                continue

            regulator_binary = binary_expression.loc[regulator.id]

            target_expression_1_regulator = target_expression[regulator_binary == 1]
            target_expression_0_regulator = target_expression[regulator_binary == 0]

            if len(target_expression_1_regulator) == 0 and len(target_expression_0_regulator) == 0:
                missed_interactions[(target.id, regulator.id)] = 1
                interactions_probabilities[(target.id, regulator.id)] = 1
                continue

            _, p_val = ks_2samp(target_expression_1_regulator, target_expression_0_regulator)
            if p_val < 0.05:
                target_binary_0_regulator = target_binary[regulator_binary == 0]

                probability = sum(target_binary_0_regulator) / len(target_binary_0_regulator)

                interactions_probabilities[(target.id, regulator.id)] = probability
                missed_interactions[(target.id, regulator.id)] = 0

            else:
                missed_interactions[(target.id, regulator.id)] = 1
                interactions_probabilities[(target.id, regulator.id)] = 1

    return interactions_probabilities, missed_interactions


# ----------------------------------------------------------------------------------------------------------------------
# Preprocessing using LinearRegression to predict genes expression from regulators co-expression
# Useful for CoRegFlux method
# ----------------------------------------------------------------------------------------------------------------------
def _get_target_regulators(gene: Union['Gene', 'Target'] = None) -> List[str]:
    """
    It returns the list of regulators of a target gene
    :param gene: Target gene
    :return: List of regulators of the target gene
    """
    if gene is None:
        return []

    if gene.is_target():
        return [regulator.id for regulator in gene.yield_regulators()]

    return []


def _filter_influence_and_expression(interactions: Dict[str, List[str]],
                                     influence: pd.DataFrame,
                                     expression: pd.DataFrame,
                                     experiments: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    It filters influence, expression and experiments matrices to keep only the targets and their regulators.
    :param interactions: Dictionary with the interactions between targets and regulators
    :param influence: Influence matrix
    :param expression: Expression matrix
    :param experiments: Experiments matrix
    :return: Filtered influence matrix, filtered expression matrix, filtered experiments matrix
    """
    targets = pd.Index(set(interactions.keys()))
    regulators = pd.Index(set([regulator for regulators in interactions.values() for regulator in regulators]))

    # filter the expression matrix for the target genes only
    expression = expression.loc[expression.index.intersection(targets)].copy()
    influence = influence.loc[influence.index.intersection(regulators)].copy()
    experiments = experiments.loc[experiments.index.intersection(regulators)].copy()
    return influence, expression, experiments


def _predict_experiment(interactions: Dict[str, List[str]],
                        influence: pd.DataFrame,
                        expression: pd.DataFrame,
                        experiment: pd.Series) -> pd.Series:
    predictions = {}

    for target, regulators in interactions.items():

        if not regulators:
            predictions[target] = np.nan
            continue

        if target not in expression.index:
            predictions[target] = np.nan
            continue

        if not set(regulators).issubset(influence.index):
            predictions[target] = np.nan
            continue

        if not set(regulators).issubset(experiment.index):
            predictions[target] = np.nan
            continue

        # a linear regression model is trained for
        # y = expression of the target gene for all samples
        # x1 = influence score of the regulator 1 in the train data set
        # x2 = influence score of the regulator 2 in the train data set
        # x3 ...

        x = influence.loc[regulators].transpose().to_numpy()
        y = expression.loc[target].to_numpy()

        regressor = LinearRegression()
        regressor.fit(x, y)

        # the expression of the target gene is predicted for the experiment
        x_pred = experiment.loc[regulators].to_frame().transpose().to_numpy()
        predictions[target] = regressor.predict(x_pred)[0]

    return pd.Series(predictions)


def predict_gene_expression(model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                            influence: pd.DataFrame,
                            expression: pd.DataFrame,
                            experiments: pd.DataFrame) -> pd.DataFrame:
    """
    It predicts the expression of genes in the experiments set using the co-expression of regulators
    in the expression and influence datasets.
    Adapted from CoRegFlux docs:
        - A GEM model containing GPRs
        - A TRN network containing target genes and the co-activators and co-repressors of each target gene
        - GEM genes and TRN targets must match
        - An influence dataset containing the influence scores of the regulators. Influence score is similar to a
        correlation score between the expression of the regulator and the expression of the target gene.
        Influence dataset format: (rows: regulators, columns: samples). Influence scores are calculated using the
        CoRegNet algorithm in the gene expression dataset.
        - A gene expression dataset containing the expression of the genes. Also called the training dataset.
        The gene expression dataset format: (rows: genes, columns: samples)
        - An experiments dataset containing the influence scores of the regulators in the experiments. These experiments
        are not used for training the linear regression model (the gene state predictor).
        These experiments are condition-specific environmental or genetic conditions that can be used
        to perform phenotype simulations.
        Experiments dataset format: (rows: regulators, columns: experiments/conditions)

    The result is a matrix of predicted gene expression values for each experiment.
    :param model: an integrated Metabolic-Regulatory model aka a MEW model
    :param influence: Influence scores of the regulators in the train data set
    :param expression: Expression of the genes in the train data set
    :param experiments: Influence scores of the regulators in the test data set
    :return: Predicted expression of the genes in the test data set
    """
    # Filtering only the gene expression and influences data of metabolic genes available in the model
    interactions = {target.id: _get_target_regulators(target) for target in model.yield_targets()}
    influence, expression, experiments = _filter_influence_and_expression(interactions=interactions,
                                                                          influence=influence,
                                                                          expression=expression,
                                                                          experiments=experiments)

    predictions = []
    for column in experiments.columns:
        experiment_prediction = _predict_experiment(interactions=interactions,
                                                    influence=influence,
                                                    expression=expression,
                                                    experiment=experiments[column])
        predictions.append(experiment_prediction)

    predictions = pd.concat(predictions, axis=1)
    predictions.columns = experiments.columns
    return predictions.dropna()
