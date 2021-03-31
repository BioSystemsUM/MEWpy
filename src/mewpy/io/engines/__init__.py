from enum import Enum

from .csv import RegulatoryCSV, CoExpressionCSV, TargetRegulatorCSV
from .cobra import CobrapyModel, ReframedModel
from .json import JSON
from .sbml import RegulatorySBML, MetabolicSBML


class Engines(Enum):
    """
    List of all engines to read/write files/models

    """

    RegulatoryCSV = RegulatoryCSV
    CoExpressionCSV = CoExpressionCSV
    TargetRegulatorCSV = TargetRegulatorCSV
    RegulatorySBML = RegulatorySBML
    MetabolicSBML = MetabolicSBML
    CobrapyModel = CobrapyModel
    ReframedModel = ReframedModel
    JSON = JSON

    @classmethod
    def has_engine(cls, engine):
        try:
            return cls[engine.__name__]

        except KeyError:

            return False

    @classmethod
    def get(cls, engine, default=None):
        try:
            return cls[engine.__name__]

        except KeyError:

            return default
