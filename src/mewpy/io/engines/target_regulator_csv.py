import os
from functools import partial
from typing import Union, TYPE_CHECKING

import numpy as np
import pandas as pd

from mewpy.io.dto import VariableRecord, DataTransferObject, FunctionTerm
from mewpy.germ.algebra import Expression, Symbol
from mewpy.germ.models import RegulatoryModel
from .engine import Engine
from .engines_utils import csv_warning

if TYPE_CHECKING:
    from mewpy.germ.models import RegulatoryModel, Model, MetabolicModel


class TargetRegulatorRegulatoryCSV(Engine):

    def __init__(self, io, config, model=None):
        """
        Engine for CSV and TXT Target-Regulator interactions based csv
        """
        super().__init__(io, config, model)

    @property
    def model_type(self):
        return 'regulatory'

    @property
    def model(self):

        if self._model is None:
            identifier = self.get_identifier()

            return RegulatoryModel(identifier=identifier)

        return self._model

    def build_data_frame(self):

        sep = self.config.get('sep', ',')
        target_col = self.config.get('target_col', 0)
        regulator_col = self.config.get('regulator_col', 1)
        header = self.config.get('header', None)
        filter_nan = self.config.get('filter_nan', False)

        names = {target_col: 'targets', regulator_col: 'regulator'}

        try:
            df = pd.read_csv(self.io, sep=sep, header=header)

        except BaseException as exc:
            self.clean()
            self.close()
            raise exc

        cols = []

        for j, col in enumerate(df.columns):
            if j in names:
                cols.append(names[j])
            else:
                del df[col]

        df.columns = cols
        df.index = df.loc[:, 'targets']

        if filter_nan:

            df = df.dropna(subset=['regulator'])

        else:

            df = df.replace(np.nan, '', regex=True)

        self.dto.data_frame = df

    def get_identifier(self):

        if os.path.exists(self.io):
            _, identifier = os.path.split(self.io)
            return os.path.splitext(identifier)[0]

    def open(self, mode='r'):

        self._dto = DataTransferObject()

        if not os.path.exists(self.io):
            raise OSError(f'{self.io} is not a valid input. Provide the path or file handler')

        self.dto.id = self.get_identifier()

    def parse(self):

        if self.dto is None:
            raise OSError('File is not open')

        if self.dto.id is None:
            raise OSError('File is not open')

        # -----------------------------------------------------------------------------
        # CSV/TXT to pandas dataframe
        # -----------------------------------------------------------------------------

        self.build_data_frame()

        targets = self.dto.data_frame.index.unique()

        for target in targets:
            # -----------------------------------------------------------------------------
            # Target
            # -----------------------------------------------------------------------------
            target_id = target.replace(' ', '')

            target_aliases = self.dto.data_frame.loc[target, self.dto.aliases_columns]

            target_record = VariableRecord(id=target_id,
                                           name=target_id,
                                           aliases=set(target_aliases))

            self.variables[target_id].add('target')

            self.dto.targets[target_id] = target_record

            # -----------------------------------------------------------------------------
            # Regulators and Function terms
            # -----------------------------------------------------------------------------
            regulators_mask = self.dto.data_frame.index == target
            regulators = self.dto.data_frame.loc[regulators_mask, 'regulator'].unique()

            regulator_records = {}
            function_terms = {}
            for i, regulator in enumerate(regulators):
                self.variables[regulator].add('regulator')

                regulator_record = VariableRecord(id=regulator,
                                                  name=regulator,
                                                  aliases={regulator})

                regulator_records[regulator] = regulator_record

                self.dto.regulators[regulator] = regulator_record

                function_terms[str(i)] = FunctionTerm(id=str(i), symbolic=Symbol(value=regulator), coefficient=i)

            # -----------------------------------------------------------------------------
            # Interaction
            # -----------------------------------------------------------------------------

            interaction_id = f'{target_id}_interaction'

            interaction_record = VariableRecord(id=interaction_id,
                                                name=interaction_id,
                                                aliases={target_id},
                                                target=target_record,
                                                function_terms=function_terms,
                                                regulators=regulator_records)

            self.dto.interactions[interaction_id] = interaction_record

            self.variables[interaction_id].add('interaction')

    def read(self,
             model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
             variables=None):

        if not model:
            model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = self.model

        if not variables:
            variables = self.variables

        if self.dto.id:
            model._id = self.dto.id

        processed_regulators = set()

        for interaction_id, interaction_record in self.dto.interactions.items():

            target_record = interaction_record.target

            target, warning = target_record.to_variable(model=model,
                                                        types=variables.get(target_record.id, {'target'}),
                                                        name=target_record.name,
                                                        aliases=target_record.aliases)

            if warning:
                self.warnings.append(partial(csv_warning, warning))

            regulators_records = interaction_record.regulators

            regulators = {}

            for regulator_id, regulator_record in regulators_records.items():

                regulator, warning = regulator_record.to_variable(model=model,
                                                                  types=variables.get(regulator_id, {'regulator'}),
                                                                  name=regulator_record.name,
                                                                  aliases=regulator_record.aliases)

                if warning:
                    self.warnings.append(partial(csv_warning, warning))

                regulators[regulator_id] = regulator

                processed_regulators.add(regulator_id)

            regulatory_events = {}

            for func_term in interaction_record.function_terms.values():
                expression_regulators = {symbol.name: regulators[symbol.name]
                                         for symbol in func_term.symbolic.atoms(symbols_only=True)}

                regulatory_events[func_term.coefficient] = Expression(symbolic=func_term.symbolic,
                                                                      variables=expression_regulators)

            interaction, warning = interaction_record.to_variable(model=model,
                                                                  types=variables.get(interaction_id, {'interaction'}),
                                                                  name=interaction_record.name,
                                                                  aliases=interaction_record.aliases,
                                                                  target=target,
                                                                  regulatory_events=regulatory_events)

            if warning:
                self.warnings.append(partial(csv_warning, warning))

            model.add(interaction, comprehensive=True)

        if len(processed_regulators) != len(self.dto.regulators):

            for regulator_id, regulator_record in self.dto.regulators.items():

                if regulator_id not in processed_regulators:

                    regulator, warning = regulator_record.to_variable(model=model,
                                                                      types=variables.get(regulator_id, {'regulator'}),
                                                                      name=regulator_record.name,
                                                                      aliases=regulator_record.aliases)

                    if warning:
                        self.warnings.append(partial(csv_warning, warning))

                    model.add(regulator)

        return model

    def write(self):

        pass

    def close(self):

        if hasattr(self.io, 'close'):
            self.io.close()

    def clean(self):
        self._dto = None
