from pathlib import Path
from typing import Type, Union, TYPE_CHECKING

from .builder import Builder
from .engines import Engines

if TYPE_CHECKING:
    from mewpy.io.engines.engine import Engine
    from io import TextIOWrapper
    from cobra import Model as Cobra_Model
    from reframed import CBModel as Reframed_Model


class Reader(Builder):
    """
    The Reader is just a wrapper for an engine. It just wraps configurations and provides the simple read API
    Yet, this object is the user interface to read files, file handlers, cbm models from cobrapy and reframed

    Wrapping is accomplished using the composition pattern

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
    """
    def __init__(self,
                 engine: Union[Type['Engine'], Engines],
                 io: Union[str, Path, 'TextIOWrapper', 'Cobra_Model', 'Reframed_Model'],
                 sep: str = ',',
                 id_col: int = 0,
                 target_col: int = 0,
                 regulator_col: int = 1,
                 rule_col: int = 1,
                 co_activating_col: int = 1,
                 co_repressing_col: int = 2,
                 aliases_cols: Union[int, tuple, list] = None,
                 header: Union[None, int] = None,
                 filter_nan: bool = False,
                 config: dict = None):

        """
        Read a given file type (or model type) into a GERM model (metabolic, regulatory, regulatory-metabolic).

        It currently supports xml, sbml, json, csv and txt. It also supports parsing of cobrapy and reframed models.

        An engine listed in the Engines enumerator of the mewpy.io.engines module must be provided. This engine
        determines the model type.

        :param engine: A valid Engine listed in the Engines enumerator of the mewpy.io.engines.
        It will handle file parsing and model building
        :param io: A valid string path or IO is acceptable.
        Alternatively, it can also be provided a cobrapy and reframed model
        :param sep: column separator
        :param id_col: identifier column
        :param target_col: target column
        :param regulator_col: regulator column
        :param rule_col: regulatory rule column
        :param co_activating_col: activating regulators column. Regulators must be separated by spaces
        :param co_repressing_col: repressing regulators' column. Regulators must be separated by spaces
        :param aliases_cols: aliases multiple columns
        :param header: If there is a header, please provide the integer for the header row
        :param filter_nan: filter regulatory interactions without explicit regulatory rule
        :param config: dictionary with additional configurations
        """
        if not engine:
            raise ValueError('Nothing to read. Please provide an engine')

        if not io:
            raise ValueError('Nothing to read. Please provide a path, file handler or model')

        if isinstance(io, Path):
            io = str(io)

        if not isinstance(engine, Engines):

            old_engine = engine

            engine = Engines.get(old_engine)

            if engine is None:
                raise ValueError(f'{old_engine} is not supported. See available engines at {Engines}')

        engine = engine.value

        if not aliases_cols:
            aliases_cols = []

        if not config:
            config = {}

        params_config = dict(sep=sep,
                             id_col=id_col,
                             target_col=target_col,
                             rule_col=rule_col,
                             co_activating_col=co_activating_col,
                             co_repressing_col=co_repressing_col,
                             aliases_cols=aliases_cols,
                             header=header,
                             filter_nan=filter_nan)

        config.update(params_config)

        super(Reader, self).__init__()

        self._engine = engine(io=io, config=config)

    # reading with contexts

    def __enter__(self):

        self._context = True

        self.engine.open()
        self.engine.parse()

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._context = False

        self.engine.close()
        self.engine.clean()

    # to understand this read method, please consult the director, builder and reader init stubs
    def read(self, warnings: bool = True):
        """
        Reading a GERM model, namely metabolic, regulatory or both encoded into one file type
        (e.g. sbml, csv, cobrapy, reframed, json, etc).
        The file is closed upon reading or failure

        :param warnings: Whether to launch warnings found during reading
        """

        if self._context:

            model = self.engine.read()

        else:
            self.engine.open()
            self.engine.parse()

            model = self.engine.read()

            self.engine.close()
            self.engine.clean()

        if warnings:

            for warning in self.warnings:
                warning()

        model.clean_history()

        return model
