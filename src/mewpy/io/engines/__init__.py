from enum import Enum

from .boolean_csv import BooleanRegulatoryCSV
from . co_expression_csv import CoExpressionRegulatoryCSV
from .target_regulator_csv import TargetRegulatorRegulatoryCSV
from .cobra_model import CobraModel
from .reframed_model import ReframedModel
from .json import JSON
from .metabolic_sbml import MetabolicSBML
from .regulatory_sbml import RegulatorySBML


class Engines(Enum):
    """
    List of all engines to read/write files/models

    """

    BooleanRegulatoryCSV = BooleanRegulatoryCSV
    CoExpressionRegulatoryCSV = CoExpressionRegulatoryCSV
    TargetRegulatorRegulatoryCSV = TargetRegulatorRegulatoryCSV
    RegulatorySBML = RegulatorySBML
    MetabolicSBML = MetabolicSBML
    CobraModel = CobraModel
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
