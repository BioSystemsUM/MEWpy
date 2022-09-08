from typing import Union, Dict, Iterable, Tuple
from io import TextIOWrapper

from cobra import Model as Cobra_Model
from reframed import CBModel as Reframed_Model
from ..io import Reader, Engines, read_model
from .problem import AbstractKOProblem


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

        engine = Engines.RegulatoryCSV

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

            engine = Engines.CobrapyModel

        elif isinstance(metabolic_model, Reframed_Model):

            engine = Engines.ReframedModel

        else:
            raise ImportError('Invalid file type')

        metabolic_model = Reader(engine=engine,
                                 io=metabolic_model,
                                 **config)

    return read_model(regulatory_model, metabolic_model, warnings=warnings)


class OptORFProblem(AbstractKOProblem):

    def __init__(self, model, fevaluation, **kwargs):

        """
        OptORF problem using the RFBA implementation analysis.
        The OptORF approach is based on gene and regulator deletion to identify optimization strategies.

        It currently supports KOProblem-based decoding

        For more details consult: https://doi.org/10.1186/1752-0509-4-53

        """
        if isinstance(model, (Cobra_Model, Reframed_Model)):

            raise ValueError(f'OptORF is not available for a model of type {type(model)}.'
                             f'Please use load_optorf() to retrieve an integrated mewpy model')

        super(OptORFProblem, self).__init__(model, fevaluation, **kwargs)

        # ignore user defined target list
        self._trg_list = None

        self._initial_state = {}

    def _build_target_list(self):
        """
        Build the target list and initial state for the OptORF problem.
        The target list is the list of regulators in the regulatory layer.
        The initial state is the state of the model before any gene deletion.
        :return:
        """
        # Target list is the combination of genes and regulators available into the mewpy integrated model

        regulators = [regulator.id for regulator in self.model.yield_regulators()
                      if not regulator.is_reaction() and not regulator.is_metabolite() and not regulator.is_gene()]

        self._trg_list = [gene for gene in self.model.genes] + regulators

        self._initial_state = {trg: 1 for trg in self._trg_list}

    def _sim_bool(self, state: Dict[str, int]) -> Dict[str, int]:
        """
        Simulate the model with a boolean state.
        :param state: a boolean state
        :return: the simulated state
        """
        # see RFBA implementation for details

        result = {}

        for interaction in self.model.yield_interactions():

            is_empty = True
            active_coefficient = None

            for coefficient, regulatory_event in interaction.regulatory_events.items():

                if not regulatory_event.is_none:

                    is_empty = False

                    regulatory_event_eval = regulatory_event.evaluate(values=state)

                    if regulatory_event_eval:
                        active_coefficient = coefficient

                        break

            if is_empty:
                target_val = interaction.target.coefficient.default_coefficient

            elif active_coefficient is not None:
                target_val = active_coefficient

            else:
                target_val = 0.0

            result[interaction.target.id] = target_val

        return result

    def _decode_metabolic_state(self, state: Dict[str, int]) -> Dict[str, Tuple[float, float]]:
        """
        Decode a metabolic state from the regulatory state.
        :param state: a regulatory state
        :return: a metabolic state (reactions constraints)
        """
        # see RFBA implementation for details

        constraints = {}

        # gpr over model.reactions
        for rxn in self.model.yield_reactions():

            if rxn.gpr.is_none:
                continue

            else:

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
        initial_state = self._initial_state.copy()

        for t_idx in candidate:

            initial_state[self.target_list[t_idx]] = 0

        # simulation of the regulatory network to determine the affected metabolic state
        state = self._sim_bool(initial_state)

        # find the affected reactions
        constraints = self._decode_metabolic_state(state)

        return constraints
