from typing import Union, Dict, Iterable, Tuple, TYPE_CHECKING, Sequence
from io import TextIOWrapper

from cobra import Model as Cobra_Model
from reframed import CBModel as Reframed_Model
from mewpy.io import Reader, Engines, read_model
from .problem import AbstractKOProblem


if TYPE_CHECKING:
    from mewpy.optimization import EvaluationFunction
    from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel


def load_optorf(regulatory_model: Union[str, TextIOWrapper, Reader],
                metabolic_model: Union[str, TextIOWrapper, Cobra_Model, Reframed_Model, Reader],
                config: dict = None,
                warnings: bool = False):
    """
    The standard method to load an OptORF problem.
    A OptORF problem is a KO strain optimization problem
    that can be created from an integrated Metabolic-Regulatory model.

    The integrated model can be assembled from a metabolic model in the following formats:
        - SBML
        - JSON
        - Cobrapy model
        - Reframed model

    and a regulatory model in the following formats:
        - SBML-qual
        - JSON
        - CSV regulatory network

    :param regulatory_model: the regulatory model or a reader for parsing it
    :param metabolic_model: the metabolic model or a reader for parsing it
    :param config: A configuration dictionary with information for parsing the models
    :param warnings: if True, warnings will be printed
    :return: an integrated model for OptORF
    """
    if not config:
        config = {}

    if not isinstance(regulatory_model, Reader):

        if isinstance(regulatory_model, str):

            file_name = regulatory_model

        elif isinstance(regulatory_model, TextIOWrapper):

            file_name = regulatory_model.name

        else:
            raise ImportError('Invalid file type')

        engine = Engines.BooleanRegulatoryCSV

        if file_name.endswith('.xml') or file_name.endswith('.sbml'):
            engine = Engines.RegulatorySBML

        regulatory_model = Reader(engine=engine,
                                  io=regulatory_model,
                                  **config)

    if not isinstance(metabolic_model, Reader):

        if isinstance(metabolic_model, str):

            engine = Engines.MetabolicSBML

        elif isinstance(metabolic_model, TextIOWrapper):

            engine = Engines.MetabolicSBML

        elif isinstance(metabolic_model, Cobra_Model):

            engine = Engines.CobraModel

        elif isinstance(metabolic_model, Reframed_Model):

            engine = Engines.ReframedModel

        else:
            raise ImportError('Invalid file type')

        metabolic_model = Reader(engine=engine,
                                 io=metabolic_model,
                                 **config)

    return read_model(regulatory_model, metabolic_model, warnings=warnings)


class OptORFProblem(AbstractKOProblem):

    def __init__(self,
                 model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                 fevaluation: Sequence['EvaluationFunction'],
                 initial_state: Dict[str, float] = None,
                 **kwargs):

        """
        OptORF problem using the RFBA implementation analysis.
        The OptORF approach is based on gene and regulator deletion to identify optimization strategies.

        It currently supports KOProblem-based decoding

        For more details consult: https://doi.org/10.1186/1752-0509-4-53
        """
        if isinstance(model, (Cobra_Model, Reframed_Model)):
            raise ValueError(f'OptORF is not available for a model of type {type(model)}.'
                             f'Please use load_optorf() to retrieve an integrated GERM model')

        super(OptORFProblem, self).__init__(model, fevaluation, **kwargs)

        # ignore user defined target list
        self._trg_list = None

        if not initial_state:
            initial_state = {}
        self._initial_state = initial_state.copy()

    def _build_target_list(self):
        """
        Build the target list and initial state for the OptORF problem.
        The target list is the list of regulators in the regulatory layer.
        The initial state is the state of the model before any gene deletion.
        :return:
        """
        # Target list is the combination of genes and regulators available into the mewpy integrated model

        regulators = [regulator.id for regulator in self.model.yield_regulators()
                      if not regulator.is_reaction() and not regulator.is_metabolite()]
        targets = [target.id for target in self.model.yield_targets()
                   if not target.is_reaction() and not target.is_metabolite()]
        genes = [gene.id for gene in self.model.yield_genes()
                 if not gene.is_reaction() and not gene.is_metabolite()]

        self._trg_list = list(set.union(set(regulators), set(targets), set(genes)))

    def _decode_initial_state(self, state: Dict[str, float] = None) -> Dict[str, float]:
        """
        Method responsible for retrieving the initial state of the model.
        The initial state is the state of all regulators found in the Metabolic-Regulatory model.
        :param state: the initial state of the model
        :return: dict of regulatory/metabolic variable keys (regulators) and a value of 0 or 1
        """
        if not state:
            state = {}

        if not self.model.is_regulatory():
            return state

        initial_state = {}
        for regulator in self.model.yield_regulators():
            if regulator.id in state:
                initial_state[regulator.id] = state[regulator.id]

            elif regulator.is_metabolite() and regulator.exchange_reaction:
                if regulator.exchange_reaction.id in state:
                    initial_state[regulator.id] = state[regulator.exchange_reaction.id]

                else:
                    initial_state[regulator.id] = abs(regulator.exchange_reaction.lower_bound)

            else:
                initial_state[regulator.id] = max(regulator.coefficients)

        return initial_state

    def _decode_regulatory_state(self, state: Dict[str, float]) -> Dict[str, float]:
        """
        It solves the boolean regulatory network for a given specific state.
        It also updates all targets having a valid regulatory interaction associated with it for the resulting state

        :param state: dict of regulatory variable keys (regulators) and a value of 0, 1 or float
        (reactions and metabolites predicates)
        :return: dict of target keys and a value of the resulting state
        """
        if not self.model.is_regulatory():
            return {}

        # Solving regulatory model synchronously for the regulators according to the initial state
        # Targets are associated with a single regulatory interaction
        result = {}
        for interaction in self.model.yield_interactions():

            # solving regulators state only
            if not interaction.target.is_regulator():
                continue

            # an interaction can have multiple regulatory events, namely one for 0 and another for 1
            for coefficient, event in interaction.regulatory_events.items():
                if event.is_none:
                    continue

                eval_result = event.evaluate(values=state)
                if eval_result:
                    result[interaction.target.id] = coefficient
                else:
                    result[interaction.target.id] = 0.0
        return result

    def _decode_metabolic_state(self, state: Dict[str, float]) -> Dict[str, float]:
        """
        It solves the boolean regulatory network for a given specific state.
        It also updates all targets having a valid regulatory interaction associated with it for the resulting state

        :param state: dict of regulatory variable keys (regulators) and a value of 0, 1 or float
        (reactions and metabolites predicates)
        :return: dict of target keys and a value of the resulting state
        """
        if not self.model.is_regulatory():
            return {}

        # Solving the whole regulatory model synchronously, as asynchronously would take too much time
        # Targets are associated with a single regulatory interaction
        result = {}
        for interaction in self.model.yield_interactions():

            # an interaction can have multiple regulatory events, namely one for 0 and another for 1
            for coefficient, event in interaction.regulatory_events.items():
                if event.is_none:
                    continue

                eval_result = event.evaluate(values=state)
                if eval_result:
                    result[interaction.target.id] = coefficient
                else:
                    result[interaction.target.id] = 0.0
        return result

    def _decode_constraints(self, state: Dict[str, float]) -> Dict[str, Tuple[float, float]]:
        """
        Method responsible for decoding the RFBA metabolic state, namely the state of all metabolic genes associated
        at least with one reaction in the GPRs rule.

        :param state: dict of regulatory/metabolic variable keys (metabolic target) and a value of 0 or 1
        :return: dict of constraints of the resulting metabolic state
        """
        # This method retrieves the constraints associated with a given metabolic/regulatory state

        if not self.model.is_metabolic():
            return {}

        constraints = {}
        for rxn in self.model.yield_reactions():

            if rxn.gpr.is_none:
                continue

            res = rxn.gpr.evaluate(values=state)

            if not res:
                constraints[rxn.id] = (0.0, 0.0)

        return constraints

    def decode(self, candidate: Iterable[int]) -> Dict[str, Tuple[float, float]]:
        """
        OptORF candidate decode.
        From a candidate aka a set of regulators to be knocked out, it runs a synchronous boolean simulation of the
        regulatory layer. The regulatory state is used to infer the metabolic state, namely metabolic reactions affected
        by the KOs. A second synchronous boolean simulation of the GPRs is performed with the regulatory state.
        Finally, the metabolic state aka constraints are decoded from the inactive reactions.

        :param candidate: a candidate (a set of regulators to be knocked out)
        :return: a metabolic state (a set of constraints)
        """
        # the perturbation state to be used to decode the initial state for simulation
        candidate_state = {self.target_list[idx]: 0 for idx in candidate}

        # the initial state of the model
        initial_state = {**self._initial_state, **candidate_state}
        initial_state = self._decode_initial_state(state=initial_state)

        # Regulatory state from a synchronous boolean simulation
        # noinspection PyTypeChecker
        regulatory_state = self._decode_regulatory_state(state=initial_state)

        # After a simulation of the regulators outputs, the state of the targets are retrieved now
        metabolic_state = self._decode_metabolic_state(state={**initial_state, **regulatory_state})

        regulatory_constraints = self._decode_constraints(metabolic_state)
        return regulatory_constraints
