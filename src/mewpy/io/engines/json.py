import json
import os
from typing import Union, TYPE_CHECKING

from mewpy.germ.models import Model
from mewpy.io.dto import DataTransferObject

from .engine import Engine

if TYPE_CHECKING:
    from mewpy.germ.models import RegulatoryModel, Model, MetabolicModel


class JSON(Engine):
    def __init__(self, io, config, model=None):

        """
        Engine for JSON files
        """
        super().__init__(io, config, model)

    @property
    def model_type(self):
        return 'metabolic'

    @property
    def model(self):

        if self._model is None:
            return Model(identifier='model')

        return self._model

    def open(self,  mode='r'):

        self._dto = DataTransferObject()

        if mode == 'r':
            if isinstance(self.io, str):

                if not os.path.exists(self.io):
                    raise OSError(f'{self.io} is not a valid input. Provide the path or file handler')

                self._io = open(self.io)

        elif mode == 'w':
            pass

        else:
            raise ValueError(f'{mode} mode is not recognized. Try one of the following: r, w')

    def parse(self):

        if self.dto is None:
            raise OSError('Model is not open')

        try:

            self.dto.model = json.load(self.io)

        except BaseException as exc:

            self.close()
            self.clean()

            raise exc

        if 'types' not in self.dto.model:
            self.close()
            self.clean()

            raise OSError(f'{self.io} is not a valid json GERM model')

    def read(self,
             model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
             variables=None):

        return Model.from_dict(self.dto.model, variables=True)

    def write(self):

        dict_model = self.model.to_dict(variables=True)

        if hasattr(self.io, 'close'):

            try:
                json.dump(dict_model, self.io)

            except BaseException as exc:

                self.close()
                self.clean()

                raise exc

        else:
            with open(self.io, 'w') as file_path:
                json.dump(dict_model, file_path)

            self.clean()

    def close(self):

        if hasattr(self.io, 'close'):
            self.io.close()

    def clean(self):
        self._dto = None
