from typing import Union, Dict, TYPE_CHECKING, Any, Sequence, Tuple

import pandas as pd

from mewpy.germ.analysis import FBA
from mewpy.germ.solution import ModelSolution, KOSolution
from mewpy.solvers.solution import Solution, Status
from mewpy.solvers.solver import Solver
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.germ.variables import Regulator, Gene, Target
    from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel


def _run_and_decode_solver(lp,
                           additional_constraints: Dict[str, Tuple[float, float]] = None,
                           **kwargs):
    if not additional_constraints:
        additional_constraints = {}

    if not kwargs:
        kwargs = {}

    if 'constraints' in kwargs:
        kwargs['constraints'].update(additional_constraints)

    solution = lp.solver.solve(**kwargs)
    if solution.status == Status.OPTIMAL:
        return solution.fobj
    else:
        return


class PROM(FBA):

    def __init__(self,
                 model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                 solver: Union[str, Solver] = None,
                 build: bool = False,
                 attach: bool = False):
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
        super().__init__(model=model, solver=solver, build=build, attach=attach)

    def _build(self):
        """
        It builds the PROM problem. It also builds a regular FBA problem to be used for the growth prediction.
        The linear problem is then loaded into the solver.
        :return:
        """
        self._build_mass_constraints()
        self._linear_objective = {var.id: value for var, value in self.model.objective.items()}
        self._minimize = False

    def _max_rates(self, solver_kwargs: Dict[str, Any]):
        # wt-type reference
        reference = self.solver.solve(**solver_kwargs)
        if reference.status != Status.OPTIMAL:
            raise RuntimeError('The solver did not find an optimal solution for the wild-type conditions.')
        reference = reference.values.copy()
        reference_constraints = {key: (reference[key] * 0.99, reference[key])
                                 for key in self._linear_objective}

        # fva of the reaction at fraction of 0.99 (for wild-type growth rate)
        rates = {}
        for reaction in self.model.reactions:
            min_rxn = _run_and_decode_solver(self,
                                             additional_constraints=reference_constraints,
                                             **{**solver_kwargs,
                                                'get_values': False,
                                                'linear': {reaction: 1},
                                                'minimize': True})
            max_rxn = _run_and_decode_solver(self,
                                             additional_constraints=reference_constraints,
                                             **{**solver_kwargs,
                                                'get_values': False,
                                                'linear': {reaction: 1},
                                                'minimize': False})

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
                     probabilities: Dict[Tuple[str, str], float],
                     regulator: Union['Gene', 'Regulator'],
                     reference: Dict[str, float],
                     max_rates: Dict[str, float],
                     to_solver: bool = False,
                     solver_kwargs: Dict[str, Any] = None):
        solver_constrains = solver_kwargs.get('constraints', {})

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

            target: Union['Target', 'Gene']

            # composed key for interactions_probabilities
            target_regulator = (target.id, regulator.id)

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

        solution = self.solver.solve(**{**solver_kwargs,
                                        'get_values': True,
                                        'constraints': {**solver_constrains, **prom_constraints}})

        if to_solver:
            return solution

        minimize = solver_kwargs.get('minimize', self._minimize)
        return ModelSolution.from_solver(method=self.method, solution=solution, model=self.model,
                                         minimize=minimize)

    def _optimize(self,
                  initial_state: Dict[Tuple[str, str], float] = None,
                  regulators: Sequence[Union['Gene', 'Regulator']] = None,
                  to_solver: bool = False,
                  solver_kwargs: Dict[str, Any] = None) -> Union[Dict[str, Solution], Dict[str, ModelSolution]]:
        # wild-type reference
        solver_kwargs['get_values'] = True
        reference = self.solver.solve(**solver_kwargs)
        if reference.status != Status.OPTIMAL:
            raise RuntimeError('The solver did not find an optimal solution for the wild-type conditions.')
        reference = reference.values.copy()

        # max and min fluxes of the reactions
        max_rates = self._max_rates(solver_kwargs=solver_kwargs)

        # a single regulator knockout
        if len(regulators) == 1:
            return self._optimize_ko(probabilities=initial_state,
                                     regulator=regulators[0],
                                     reference=reference,
                                     max_rates=max_rates,
                                     to_solver=to_solver,
                                     solver_kwargs=solver_kwargs)

        # multiple regulator knockouts
        kos = {}
        for regulator in regulators:
            ko_solution = self._optimize_ko(probabilities=initial_state,
                                            regulator=regulator,
                                            reference=reference,
                                            max_rates=max_rates,
                                            to_solver=to_solver,
                                            solver_kwargs=solver_kwargs)
            kos[regulator.id] = ko_solution
        return kos

    def optimize(self,
                 initial_state: Dict[Tuple[str, str], float] = None,
                 regulators: Union[str, Sequence['str']] = None,
                 to_solver: bool = False,
                 solver_kwargs: Dict[str, Any] = None) -> Union[KOSolution, Dict[str, Solution]]:
        """
        It solves the PROM linear problem. The linear problem is solved using the solver interface.

        The optimize method allows setting temporary changes to the linear problem. The changes are
        applied to the linear problem reverted to the original state afterward.
        Objective, constraints and solver parameters can be set temporarily.
        :param initial_state: dictionary with the probabilities of
        the interactions between the regulators and the targets.
        :param regulators: list of regulators to be knocked out. If None, all regulators are knocked out.
        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False.
        Otherwise, a ModelSolution is returned.
        :param solver_kwargs: Solver parameters to be set temporarily.
            - linear: A dictionary of linear coefficients to be set temporarily. The keys are the variable names
            and the values are the coefficients. Default: None
            - quadratic: A dictionary of quadratic coefficients to be set temporarily. The keys are tuples of
            variable names and the values are the coefficients. Default: None
            - minimize: Whether to minimize the objective. Default: False
            - constraints: A dictionary with the constraints bounds. The keys are the constraint ids and the values
            are tuples with the lower and upper bounds. Default: None
            - get_values: Whether to retrieve the solution values. Default: True
            - shadow_prices: Whether to retrieve the shadow prices. Default: False
            - reduced_costs: Whether to retrieve the reduced costs. Default: False
            - pool_size: The size of the solution pool. Default: 0
            - pool_gap: The gap between the best solution and the worst solution in the pool. Default: None
        :return: A KOSolution instance or a list of SolverSolution instance if to_solver is True.
        """
        # build solver if out of sync
        if not self.synchronized:
            self.build()

        if not initial_state:
            initial_state = {}

        if not regulators:
            regulators = list(self.model.yield_regulators())
        else:
            if isinstance(regulators, str):
                regulators = [regulators]

            regulators = [self.model.get(regulator) for regulator in regulators]

        if not solver_kwargs:
            solver_kwargs = {}

        # concrete optimize
        solutions = self._optimize(initial_state=initial_state,
                                   regulators=regulators,
                                   to_solver=to_solver,
                                   solver_kwargs=solver_kwargs)

        if to_solver:
            return solutions

        return KOSolution(solutions)


# ----------------------------------------------------------------------------------------------------------------------
# Probability of Target-Regulator interactions
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

    :param model: an integrated Metabolic-Regulatory model aka a GERM model
    :param expression: Quantile preprocessed expression matrix
    :param binary_expression: Quantile preprocessed expression matrix binarized
    :return: Dictionary with the conditional probability of a target gene being active when the regulator is inactive,
    Dictionary with missed interactions
    """
    try:
        # noinspection PyPackageRequirements
        from scipy.stats import ks_2samp
    except ImportError:
        raise ImportError('The package scipy is not installed. '
                          'To compute the probability of target-regulator interactions, please install scipy '
                          '(pip install scipy).')
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
