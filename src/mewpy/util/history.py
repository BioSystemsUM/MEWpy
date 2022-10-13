from collections.abc import Callable
from typing import TYPE_CHECKING, Union
from functools import partial, wraps

import pandas as pd

if TYPE_CHECKING:
    from mewpy.germ.models import Model
    from mewpy.germ.variables import Variable


class HistoryManager:

    def __init__(self):

        self._history = []
        self._undo_able_commands = []
        self._temp_stack = []
        self._redo_able_commands = []

    def __str__(self):
        return f'History: {len(self._undo_able_commands)} undos and {len(self._redo_able_commands)} redos'

    @property
    def history(self):

        return pd.DataFrame(data=self._history, columns=['method', 'args', 'kwargs', 'object'])

    @property
    def undo_able_commands(self):
        return self._undo_able_commands

    @property
    def redo_able_commands(self):
        return self._redo_able_commands

    def _do(self, undo=True):

        if undo:
            method = self.undo_able_commands.pop()
            redo_command = self._temp_stack.pop()
            self.redo_able_commands.append(redo_command)

        else:
            method = self.redo_able_commands.pop()

        method()

    def undo(self) -> None:
        self._do(undo=True)

    def redo(self) -> None:
        self._do(undo=False)

    def reset(self) -> None:

        while len(self.undo_able_commands) > 0:
            self.undo()

    def restore(self) -> None:

        while len(self.redo_able_commands) > 0:
            self.redo()

    def __call__(self, *args, **kwargs) -> None:

        return self.queue_command(*args, **kwargs)

    def queue_command(self,
                      undo_func: Callable,
                      func: Callable,
                      undo_args: tuple = None,
                      undo_kwargs: dict = None,
                      args: tuple = None,
                      kwargs: dict = None,
                      obj: 'Model' = None) -> None:

        if not undo_args:
            undo_args = ()

        if not undo_kwargs:
            undo_kwargs = {}

        self.undo_able_commands.append(partial(undo_func, *undo_args, **undo_kwargs))

        if not args:
            args = ()

        if not kwargs:
            kwargs = {}

        self._temp_stack.append(partial(func, *args, **kwargs))

        self._history.append((func.__name__, str(args), str(kwargs), str(obj)))


def recorder(func: Callable):

    @wraps(func)
    def wrapper(self: Union['Model', 'Variable'], value):

        history = self.history

        old_value = getattr(self, func.__name__)

        if old_value != value:

            history.queue_command(undo_func=func,
                                  undo_args=(self, old_value),
                                  func=func,
                                  args=(self, value))

        func(self, value)

    return wrapper
