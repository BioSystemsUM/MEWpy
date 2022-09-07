from enum import Enum

from .fba import FBA, milpFBA, pFBA
from .integrated_analysis import (slim_rfba, slim_srfba, ifva, isingle_regulator_deletion, isingle_reaction_deletion,
                                  isingle_gene_deletion)
from .metabolic_analysis import slim_fba, slim_pfba, slim_milp_fba, fva, single_gene_deletion, single_reaction_deletion
from .milp_bool import milpBool
from .regulatory_analysis import slim_sim_bool, slim_milp_bool, single_regulator_deletion, regulatory_events
from .rfba import RFBA
from .sim_bool import SimBool
from .srfba import SRFBA


class Analysis(Enum):
    """
    Enumeration of the available analysis methods.
    """
    FBA = FBA
    milpFBA = milpFBA
    pFBA = pFBA
    RFBA = RFBA
    SRFBA = SRFBA
    SynchronousBoolean = SimBool
    BooleanMILP = milpBool

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
