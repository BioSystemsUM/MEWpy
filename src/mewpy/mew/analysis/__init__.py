from enum import Enum

from .fba import FBA
from .pfba import pFBA
from .rfba import RFBA
from .srfba import SRFBA
from .prom import PROM, target_regulator_interaction_probability
from .coregflux import CoRegFlux, predict_gene_expression
from .metabolic_analysis import slim_fba, slim_pfba, fva, single_gene_deletion, single_reaction_deletion
from .regulatory_analysis import regulatory_truth_table
from .integrated_analysis import (slim_rfba, slim_srfba, slim_prom, slim_coregflux,
                                  ifva, isingle_regulator_deletion, isingle_reaction_deletion, isingle_gene_deletion,
                                  find_conflicts)


class Analysis(Enum):
    """
    Enumeration of the available analysis methods.
    """
    FBA = FBA
    pFBA = pFBA
    RFBA = RFBA
    SRFBA = SRFBA
    PROM = PROM
    CoRegFlux = CoRegFlux

    @classmethod
    def has_analysis(cls, analysis: str) -> bool:
        """
        Check if the analysis is available.

        :param analysis: The name of the analysis
        :return: True if the analysis is available, False otherwise
        """
        return analysis in cls.__members__

    @classmethod
    def get(cls, analysis: str, default=FBA) -> 'Analysis':
        """
        Get the analysis class.

        :param analysis: the analysis name
        :param default: the default analysis to return if the analysis is not available
        :return: the analysis class
        """
        try:
            return cls[analysis]

        except KeyError:

            return default
