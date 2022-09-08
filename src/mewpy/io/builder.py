from typing import Union, TYPE_CHECKING

if TYPE_CHECKING:

    from io import TextIOWrapper

    from mewpy.io.engines.engine import Engine
    from cobra import Model as Cobra_Model
    from reframed import CBModel as Reframed_Model


class Builder:

    def __init__(self):
        """
        Builder base for an engine wrapper of any kind.
        That is, the builder does not provide anything by itself.
        Instead, this class is to be used by the reader and writer, which are just wrappers of engines.

        Wrapping is accomplished using the composition pattern

        Engines, on the other hand, do the hard work of opening, parsing and reading/writing a specific file or model.
        Basically, engines are, actually, the builders of models or files

        """
        self._engine = None

        # useful for the context manager and the read or write methods
        self._context = False

    @property
    def engine(self) -> 'Engine':
        """
        Returns the engine associated with this builder
        :return: engine
        """
        return self._engine

    @property
    def io(self) -> Union[str, 'TextIOWrapper', 'Cobra_Model', 'Reframed_Model']:
        """
        Returns the IO associated with this builder
        :return: IO
        """
        return self.engine.io

    @property
    def config(self) -> dict:
        """
        Returns the configuration associated with this builder
        :return: configuration
        """
        return self.engine.config

    @property
    def warnings(self) -> tuple:
        """
        Returns the warnings associated with this builder
        :return: warnings
        """
        return tuple(self.engine.warnings)
