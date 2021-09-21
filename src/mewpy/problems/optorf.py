from typing import Union
from io import TextIOWrapper

try:
    # noinspection PyPackageRequirements
    from cobra import Model as Cobra_Model

except ImportError:
    Cobra_Model = str

try:
    # noinspection PyPackageRequirements
    from reframed import CBModel as Reframed_Model

except ImportError:

    Reframed_Model = str

from ..io import Reader, Engines, read_model
from .problem import AbstractKOProblem


# TODO: should it be in io?
def load_optorf(regulatory_model: Union[str, TextIOWrapper, Reader],
                metabolic_model: Union[str, TextIOWrapper, Cobra_Model, Reframed_Model, Reader],
                config: dict = None,
                warnings: bool = False):
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

        # TODO: this is the main disadvantage of the new model integration.
        #  Integrated models (metabolic-regulatory) are only available with conversion to mewpy models.
        #  However, one can easily converted a reframed or cobra model to mewpy with the IO read that I implemented.
        if isinstance(model, (Cobra_Model, Reframed_Model)):

            raise ValueError(f'OptORF is not available for a model of type {type(model)}.'
                             f'Please use load_optorf() to retrieve an integrated mewpy model')

        super(OptORFProblem, self).__init__(model, fevaluation, **kwargs)

        # ignore user defined target list
        self._trg_list = None

        self._initial_state = {}

    def _build_target_list(self):

        # Target list is the combination of genes and regulators available into the mewpy integrated model

        regulators = [regulator.id for regulator in self.model.yield_regulators()
                      if not regulator.is_reaction() and not regulator.is_metabolite() and not regulator.is_gene()]

        self._trg_list = [gene for gene in self.model.genes] + regulators

        self._initial_state = {trg: 1 for trg in self._trg_list}

    def _sim_bool(self, state):

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
                target_val = interaction.target.coefficient.active_coefficient

            elif active_coefficient is not None:
                target_val = active_coefficient

            else:
                target_val = 0.0

            result[interaction.target.id] = target_val

        return result

    def _decode_metabolic_state(self, state):

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

    def decode(self, candidate):

        """

        OptORF candidate decode

        :param candidate:
        :return:
        """

        for t_idx in candidate:

            self._initial_state[self.target_list[t_idx]] = 0

        # simulation of the regulatory network to determine the affected metabolic state
        state = self._sim_bool(self._initial_state)

        # find the affected reactions
        constraints = self._decode_metabolic_state(state)

        return constraints
