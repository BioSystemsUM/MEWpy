from pathlib import Path
from typing import Type, Union, TYPE_CHECKING

from .builder import Builder
from .engines import Engines

if TYPE_CHECKING:
    from .engines.engine import Engine
    from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel
    from cobra import Model as CobraModel
    from reframed import CBModel as ReframedModel


class Writer(Builder):
    """
    The Writer is just a wrapper for an engine. It just wraps configurations and provides the simple API write
    Yet, this object is the user interface to write files, file handlers, cbm models from cobrapy and reframed

    Wrapping is accomplished using the composition pattern

    The only difference between the builders and the engines is that engines write the model in four steps:
                - open
                - write
                - close
                - clean

    Builders, on the other hand, just write and thus performing all operations at once.
    This is important for writing multiple files out of a single model.

    For writing files in stages, multiple writers must be created and the director must be used to merge all the
    writers into a single model
    """
    def __init__(self,
                 engine: Union[Type['Engine'], Engines],
                 io: Union[str, Path, 'CobraModel', 'ReframedModel'],
                 model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                 config: dict = None):

        """
        Write a given model type into a file of any type.

        It currently supports xml, sbml and json. It also supports writing of cobrapy and reframed models.

        An engine listed in the Engines enumerator of the mewpy.io.engines module must be provided. This engine
        determines the file type.

        :param engine: A valid Engine listed in the Engines enumerator of the mewpy.io.engines.
        It will handle file writing
        :param io: A valid string path or IO is acceptable.
        Alternatively, it can also be provided a cobrapy and reframed model that will be filled
        :param model: A valid metabolic, regulatory or both GERM model
        :param config: Dictionary with additional configurations

        """
        if not engine:
            raise ValueError('Nothing to write. Please provide an engine')

        if not io:
            raise ValueError('Nothing to write. Please provide a path, file handler or model')

        if not model:
            raise ValueError('Nothing to write. Please provide a GERM model')

        if isinstance(io, Path):
            io = str(io)

        if not isinstance(engine, Engines):

            old_engine = engine

            engine = Engines.get(old_engine)

            if engine is None:
                raise ValueError(f'{old_engine} is not supported. See available engines at {Engines}')

        engine = engine.value

        if not config:
            config = {}

        super(Writer, self).__init__()

        self._engine = engine(io=io, config=config, model=model)

    def __enter__(self):

        self._context = True

        self.engine.open()

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):

        self._context = False

        self.engine.close()
        self.engine.clean()

    # to understand this read method consult the director, builder and reader init stubs
    def write(self, warnings: bool = True):
        """
        Writing a metabolic, regulatory or GERM model to a file type
        (e.g. sbml, json, cobrapy, reframed, etc).
        The file is closed upon writing or failure

        :param warnings: Whether to launch warnings found during reading
        """

        if self._context:

            self.engine.write()

        else:
            self.engine.open()

            self.engine.write()

            self.engine.close()
            self.engine.clean()

        if warnings:

            for warning in self.warnings:
                warning()
