import warnings
from typing import Union, Dict, TYPE_CHECKING

import numpy as np
import pandas as pd

from mewpy.mew.lp import MetabolicLinearizer, Notification
from mewpy.mew.solution import ModelSolution, KOSolution
from mewpy.solvers.solution import Solution
from mewpy.solvers.solver import Solver
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.mew.variables import Regulator
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


class PROM(MetabolicLinearizer):

    def __init__(self,
                 model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                 solver: Union[str, Solver] = None,
                 build: bool = True,
                 attach: bool = True):

        """
        The Probabilistic Regulation of Metabolism (PROM) algorithm predicts the growth phenotype and the flux response
        after transcriptional perturbation, given a metabolic and regulatory network.
        PROM introduces probabilities to represent gene states and gene-transcription factor interactions.

        For more detail consult: https://doi.org/10.1073/pnas.1005139107
        :param model: The metabolic and regulatory model to be simulated.
        :param solver: The solver to be used. If None, a new instance will be created from the default solver.
        :param build: If True, the linear problem will be built upon initialization.
        If False, the linear problem can be built later by calling the build() method.
        :param attach: If True, the linear problem will be attached to the model.
        """
        # prom requires regular fba simulations
        from .fba import FBA
        self._fba = FBA(model=model, solver=None, build=False, attach=False)

        super().__init__(model=model,
                         solver=solver,
                         build=build,
                         attach=attach)

    @property
    def fba(self):
        return self._fba

    def build(self):
        """
        It builds the PROM problem. It also builds a regular FBA problem to be used for the growth prediction.
        The linear problem is then loaded into the solver.
        :return:
        """
        self.fba.build()

        if self._variables or self._constraints:
            self.clean()

    def notification(self, notification: Notification):
        """
        It handles the notifications received from the model and variables.
        The notifications are used to update the linear problem. The linear problem is then loaded into the solver.
        A notification contains a message and a payload. The message carries the content of the changes and the payload
        carries the information about the changes.

        :param notification: A Notification instance.
        :return:
        """
        if notification.content_type == 'reactions' and notification.action == 'add':

            return self.build()

        elif notification.content_type == 'reactions' and notification.action == 'remove':

            return self.build()

        elif notification.content_type == 'metabolites' and notification.action == 'add':

            return self.build()

        elif notification.content_type == 'metabolites' and notification.action == 'remove':

            return self.build()

        elif notification.content_type == 'coefficients' and notification.action == 'set':

            return self.build()

        else:
            return super(PROM, self).notification(notification)

    def fba_optimize(self, constraints=None, objective=None, minimize=None, get_values=True):

        if not constraints:
            constraints = {}

        old_linear_obj = self.fba._linear_objective.copy()
        old_quadratic_obj = self.fba._quadratic_objective.copy()
        old_sense = self.fba.minimize

        tol_bounds = {}
        for rxn in self.model.yield_reactions():

            if rxn.lower_bound == rxn.upper_bound:
                tol_bounds[rxn.id] = (rxn.lower_bound - ModelConstants.TOLERANCE, rxn.upper_bound)

        constraints = {**tol_bounds, **constraints}

        sol = self.fba.optimize(objective=objective,
                                minimize=minimize,
                                constraints=constraints,
                                to_solver=True,
                                get_values=get_values)

        self.fba.set_objective(linear=old_linear_obj, quadratic=old_quadratic_obj, minimize=old_sense)

        if not get_values:

            if sol.fobj is None:
                return 0

            return sol.fobj

        else:

            if not sol.values:
                return {rxn: 0 for rxn in self.model.reactions}, 0, sol.status

            values = {key: val if abs(val) > ModelConstants.TOLERANCE else 0
                      for key, val in sol.values.items()}

            return values, sol.fobj, sol.status

    @staticmethod
    def knn_imputation_transformation(expression,
                                      missing_values=None,
                                      n_neighbors=5,
                                      weights='uniform',
                                      metric='nan_euclidean'):

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

    @staticmethod
    def quantile_transformation(expression,
                                n_quantiles=None):

        if n_quantiles is None:
            n_quantiles = expression.shape[1]

        return quantile_transform(expression, n_quantiles=n_quantiles, axis=0, random_state=0)

    @staticmethod
    def quantile_thresholding(expression, q=0.33):

        threshold = np.quantile(expression, q)

        return threshold

    @staticmethod
    def binarize(expression, threshold):

        threshold_mask = expression >= threshold
        expression[threshold_mask] = 1
        expression[~threshold_mask] = 0

        return expression

    @staticmethod
    def expression_preprocessing(expression,
                                 missing_values=None,
                                 n_neighbors=5,
                                 weights='uniform',
                                 metric='nan_euclidean',
                                 n_quantiles=None,
                                 q=0.33):

        expression_to_process = np.array(expression.values)

        expression_to_process = PROM.knn_imputation_transformation(expression_to_process,
                                                                   missing_values=missing_values,
                                                                   n_neighbors=n_neighbors,
                                                                   weights=weights,
                                                                   metric=metric)
        expression_to_process = PROM.quantile_transformation(expression_to_process, n_quantiles=n_quantiles)

        # binary_expression
        threshold = PROM.quantile_thresholding(expression_to_process, q=q)

        expression = pd.DataFrame(data=expression_to_process,
                                  index=expression.index,
                                  columns=expression.columns)

        binary_expression = expression.copy()

        binary_expression = PROM.binarize(expression=binary_expression, threshold=threshold)

        return expression, binary_expression

    def calculate_probabilities(self, expression, binary_expression):

        missed_interactions = {}
        interactions_probabilities = {}

        # PROM is based on multiple interactions for the same target. In PROM, interactions are considered to be
        # one pair target-regulator. However, a given target can be regulated by several regulators, which means
        # that can exist multiple interactions for the same target. Thus, we have compiled this multiple
        # interactions into the regulatory events.
        # Then, here the rational of PROM is reconstructed again, by creating composed keys: target + regulator
        for interaction in self.model.yield_interactions():

            target = interaction.target
            regulators = interaction.regulators

            if not regulators or target.id not in expression.index:

                missed_interactions[target.id] = 1
                interactions_probabilities[target.id] = 1

            else:

                target_expression = expression.loc[target.id]
                target_binary = binary_expression.loc[target.id]

                for regulator in regulators.values():

                    if regulator.id not in expression.index:
                        missed_interactions[target.id + regulator.id] = 1
                        interactions_probabilities[target.id + regulator.id] = 1

                    else:

                        regulator_binary = binary_expression.loc[regulator.id]

                        target_expression_1_regulator = target_expression[regulator_binary == 1]
                        target_expression_0_regulator = target_expression[regulator_binary == 0]

                        if len(target_expression_1_regulator) > 0 and len(target_expression_0_regulator) > 0:

                            _, p_val = ks_2samp(target_expression_1_regulator, target_expression_0_regulator)

                            if p_val < 0.05:
                                target_binary_0_regulator = target_binary[regulator_binary == 0]

                                probability = sum(target_binary_0_regulator) / len(target_binary_0_regulator)

                                interactions_probabilities[target.id + regulator.id] = probability
                                missed_interactions[target.id + regulator.id] = 0

                            else:
                                missed_interactions[target.id + regulator.id] = 1
                                interactions_probabilities[target.id + regulator.id] = 1

                        else:
                            missed_interactions[target.id + regulator.id] = 1
                            interactions_probabilities[target.id + regulator.id] = 1

        if sum(missed_interactions.values()) > (0.75 * len(missed_interactions)):
            warnings.warn('Binarization threshold should be changed', Warning, stacklevel=2)

        return interactions_probabilities, missed_interactions

    def _max_rates(self):
        # wt-type reference
        reference, _, _ = self.fba_optimize()
        reference_constraints = {key: (reference[key], reference[key])
                                 for key in self._linear_objective}

        rates = {}

        # fva of the reaction at fraction of 1 (for wild-type growth rate)
        for reaction in self.model.reactions:

            min_rxn = self.fba_optimize(constraints=reference_constraints,
                                        objective={reaction: 1},
                                        minimize=True,
                                        get_values=False)

            max_rxn = self.fba_optimize(constraints=reference_constraints,
                                        objective={reaction: 1},
                                        minimize=False,
                                        get_values=False)

            reference_rate = reference[reaction]

            if reference_rate < 0:
                value = min((min_rxn, max_rxn, reference_rate))

            elif reference_rate > 0:
                value = max((min_rxn, max_rxn, reference_rate))

            else:
                value = max((abs(min_rxn), abs(max_rxn), abs(reference_rate)))

            if abs(value) < ModelConstants.TOLERANCE:
                value = 0.0

            rates[reaction] = value

        return rates

    def _optimize_ko(self,
                     regulator: 'Regulator',
                     probabilities: Dict[str, float],
                     reference: Dict[str, float],
                     max_rates: Dict[str, float]):

        prom_constraints = {reaction.id: reaction.bounds for reaction in self.model.yield_reactions()}
        state = {gene: 1 for gene in self.model.genes.keys()}

        # if the regulator to be ko is a metabolic gene, the associated reactions are ko too
        # prom constraints of the associated reactions are set to threshold
        if regulator.is_gene():

            for reaction in regulator.reactions.keys():
                prom_constraints[reaction] = (-ModelConstants.TOLERANCE, ModelConstants.TOLERANCE)

        # finds the target genes of the deleted regulator.
        # finds the reactions associated with these target genes.
        # The reactions' bounds might be changed next, but for now the flag is set to False
        target_reactions = {}
        for target in regulator.yield_targets():

            if target.is_gene():
                # after the regulator ko iteration, this is reset
                state[target.id] = 0

                target_reactions.update({reaction.id: reaction for reaction in target.yield_reactions()})

        # GPR evaluation of each reaction previously found, but using the changed gene_state.
        # If the GPR is evaluated to zero, the reaction bounds will be changed in the future.
        # For that, the reactions dictionary flags must be updated to True.
        inactive_reactions = {}
        for reaction in target_reactions.values():

            if reaction.gpr.is_none:
                continue

            if reaction.gpr.evaluate(values=state):
                continue

            inactive_reactions[reaction.id] = reaction

        # for each target regulated by the regulator
        for target in regulator.yield_targets():

            if not target.is_gene():
                continue

            # composed key for interactions_probabilities
            target_regulator = target.id + regulator.id

            if target_regulator not in probabilities:
                continue

            interaction_probability = probabilities[target_regulator]

            # for each reaction associated with this single target
            for reaction in target.yield_reactions():

                # if the gpr has been evaluated previously to zero,
                # it means that the metabolic genes regulated by this regulator can affect the state of the
                # reaction. Thus, the reaction bounds can be changed using PROM probability.
                # Nevertheless, it is only useful to do that if the probability is inferior to 1, otherwise
                # nothing is changed
                if reaction.id not in inactive_reactions:
                    continue

                if interaction_probability >= 1:
                    continue

                # reaction old bounds
                rxn_lb, rxn_ub = tuple(prom_constraints[reaction.id])

                # probability flux is the upper or lower bound that this reaction can take
                # when the regulator is KO. This is calculated as follows:
                # interaction probability times the reaction maximum limit (determined by fva)
                probability_flux = max_rates[reaction.id] * interaction_probability

                # wild-type flux value for this reaction
                wt_flux = reference[reaction.id]

                # update flux bounds according to probability flux
                if wt_flux < 0:

                    rxn_lb = max((reaction.lower_bound, probability_flux, rxn_lb))
                    rxn_lb = min((rxn_lb, -ModelConstants.TOLERANCE))

                elif wt_flux > 0:

                    rxn_ub = min((reaction.upper_bound, probability_flux, rxn_ub))
                    rxn_ub = max((rxn_ub, ModelConstants.TOLERANCE))

                else:

                    # if it is zero, the reaction is not changed, so that reactions are not activated
                    # by PROM. Only reaction ko is forced by PROM.

                    continue

                prom_constraints[reaction.id] = (rxn_lb, rxn_ub)

        fba_values, fba_objective_value, fba_status = self.fba_optimize(constraints=prom_constraints,
                                                                        minimize=False)

        return ModelSolution(method=self.method,
                             x=fba_values,
                             objective_value=fba_objective_value,
                             status=fba_status,
                             objective_direction='maximize',
                             model=self.model)

    def _optimize(self, probabilities: Dict[str, Union[float, int]]):
        # wild-type reference
        reference, _, _ = self.fba_optimize()

        # max and min fluxes of the reactions
        max_rates = self._max_rates()

        kos = {}

        # knock out of regulators
        for regulator in self.model.yield_regulators():
            ko_solution = self._optimize_ko(regulator=regulator,
                                            probabilities=probabilities,
                                            reference=reference,
                                            max_rates=max_rates)
            kos[regulator.id] = ko_solution
        return kos

    def optimize(self,
                 expression=None,
                 objective=None,
                 minimize=None,
                 constraints=None,
                 to_solver=False,
                 get_values=True,
                 shadow_prices=False,
                 reduced_costs=False,
                 pool_size=0,
                 pool_gap=None) -> Union[KOSolution, Solution]:

        """

        :param expression: Gene expression data frame where rows/index must be set with the gene, target and regulator
        identifiers of the model
        :param objective:
        :param minimize:
        :param constraints:
        :param to_solver:
        :param get_values:
        :param shadow_prices:
        :param reduced_costs:
        :param pool_size:
        :param pool_gap:
        :return:
        """
        expression, binary_expression = self.expression_preprocessing(expression=expression)
        probabilities, _ = self.calculate_probabilities(expression=expression,
                                                        binary_expression=binary_expression)
        kos = self._optimize(probabilities=probabilities)
        return KOSolution(kos)
