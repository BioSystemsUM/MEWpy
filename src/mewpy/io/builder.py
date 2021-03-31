from typing import Union, TYPE_CHECKING

if TYPE_CHECKING:

    from io import TextIOWrapper

    from mewpy.io.engines.engine import Engine

    try:
        # noinspection PyPackageRequirements
        from cobra import Model as Cobra_Model

    except ImportError:
        Cobra_Model = str

    try:
        # noinspection PyPackageRequirements
        from reframed import ReframedModel as Reframed_Model

    except ImportError:
        Reframed_Model = str


class Builder:

    def __init__(self):
        """
        Builder base for an engine wrapper of any kind.
        That is, the builder does not provide anything by itself.
        Instead, this class is to be used by the reader and writer, which are just wrappers of engines.

        # Wrapping is accomplished using the composition pattern

        Engines, on the other hand, do the hard work of opening, parsing and reading/writing a specific file or model.
        Basically, engines are, actually, the builders of models or files

        """
        self._engine: 'Engine' = None

        # useful for the context manager and the read or write methods
        self._context = False

    @property
    def engine(self) -> 'Engine':
        return self._engine

    @property
    def io(self) -> Union[str, 'TextIOWrapper', 'Cobra_Model', 'Reframed_Model']:
        return self.engine.io

    @property
    def config(self) -> dict:
        return self.engine.config

    @property
    def warnings(self) -> tuple:
        return tuple(self.engine.warnings)
