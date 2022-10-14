from abc import ABCMeta, abstractmethod
from collections import defaultdict

from typing import Union, TYPE_CHECKING, Dict, Optional

if TYPE_CHECKING:
    from mewpy.germ.models import RegulatoryModel, MetabolicModel, Model
    from mewpy.io.dto import DataTransferObject
    from io import TextIOWrapper
    from cobra import Model as CobraModel
    from reframed import CBModel as ReframedModel


class Engine(metaclass=ABCMeta):

    def __init__(self,
                 io: Union[str, 'TextIOWrapper', 'CobraModel', 'ReframedModel'],
                 config: dict,
                 model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None):
        """
        The engine interface for reading/writing models. The abstract properties and methods should be fully
        implemented in the concrete engines.

        The concrete engines will open, parse and read specific types to a GERM model,
        such as:
            - regulatory sbml
            - metabolic sbml
            - regulatory csv
            - cbm models objects from external packages
            - json files

        The concrete engines will also open and write specific types from an adequate GERM model,
        such as:
            - regulatory sbml
            - metabolic sbml
            - regulatory csv
            - cbm models objects from external packages
            - regulatory sbml
            -jsons

        The only difference between the builders and the engines is that engines read the model in five steps:
            - open
            - parse
            - read
            - close
            - clean

        Builders, on the other hand, just read and thus performing all operations at once.
        This is important for reading multiple files into a single model.

        For reading files in stages, multiple readers must be created and the director must be used to merge all the
        readers into a single model

        The only difference between the builders and the engines is that engines write the model in four steps:
                    - open
                    - write
                    - close
                    - clean

        Builders, on the other hand, just write and thus performing all operations at once.
        This is important for writing multiple files out of a single model.

        For writing files in stages, multiple writers must be created and the director must be used to merge all the
        writers into a single model

        :param io: file path, file handler or cbm model object from cobrapy or reframed
        :param config: dictionary of multiple configurations to be used when reading
        :param model: in case of writing, a GERM model must be provided
        """
        self._io = io
        self._config = config
        self._model = model
        self._dto = None
        self._variables = defaultdict(set)
        self._warnings = []

    @property
    def io(self) -> Union[str, 'TextIOWrapper', 'CobraModel', 'ReframedModel']:
        return self._io

    @property
    def config(self) -> dict:
        return self._config

    @property
    def model(self) -> Union['Model', 'MetabolicModel', 'RegulatoryModel']:
        return self._model

    @property
    @abstractmethod
    def model_type(self) -> Optional[str]:
        return

    @property
    def dto(self) -> 'DataTransferObject':
        return self._dto

    @property
    def variables(self) -> Dict[str, set]:
        return self._variables

    @property
    def warnings(self) -> list:
        return self._warnings

    @abstractmethod
    def open(self, mode: str = 'r'):
        """
        Open makes the preparation for reading or writing

        :param mode: One of the following: r - read; w - write;
        :return:
        """
        pass

    @abstractmethod
    def parse(self):
        """
        Parses a given file to a DataTransferObject

        :return:
        """
        pass

    @abstractmethod
    def read(self,
             model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
             variables: dict = None) -> Union['Model', 'MetabolicModel', 'RegulatoryModel']:
        """
        Reads a model into a GERM model. If a model is provided, the read method will increment further variables
        If a variables dictionary is provided, multi-type variables can be built together with the ones available in
        the model

        Reading is performed from the middle DataTransferObject

        :param model: A valid GERM model
        :param variables: A dictionary of variables already built to be updated during reading
        :return: GERM model
        """
        pass

    @abstractmethod
    def write(self):
        """
        Writes a GERM model.
        :return:
        """

        pass

    @abstractmethod
    def close(self):
        """
        Closes the engine.
        :return:
        """
        pass

    @abstractmethod
    def clean(self):
        """
        Cleans the engine.
        :return:
        """
        pass
