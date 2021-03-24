from abc import ABCMeta, abstractmethod


class ModelSolutionInterface(metaclass=ABCMeta):

    """
    Solution interface
    """

    @property
    @abstractmethod
    def objective_direction(self):
        return

    @property
    @abstractmethod
    def method(self):
        return

    @property
    @abstractmethod
    def model(self):
        return

    @property
    @abstractmethod
    def objective(self):
        return

    @property
    @abstractmethod
    def objective_value(self):
        return

    @property
    @abstractmethod
    def status(self):
        return

    @property
    @abstractmethod
    def x(self):
        return

    @property
    @abstractmethod
    def simulator(self):
        return

    @abstractmethod
    def to_summary(self):
        return

    @abstractmethod
    def to_series(self):
        return

    @abstractmethod
    def to_frame(self):
        return

    @abstractmethod
    def __repr__(self):
        return

    @abstractmethod
    def __str__(self):
        return
