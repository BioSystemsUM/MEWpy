from enum import Enum

from .fba import FBA, milpFBA, pFBA
from .rfba import RFBA
from .srfba import SRFBA
from .sim_bool import SimBool
from .milp_bool import milpBool

from .metabolic_analysis import slim_fba, slim_pfba, slim_milp_fba, fva, single_gene_deletion, single_reaction_deletion
from .regulatory_analysis import slim_sim_bool, slim_milp_bool, single_regulator_deletion, regulatory_events
from .integrated_analysis import (slim_rfba, slim_srfba, ifva, isingle_regulator_deletion, isingle_reaction_deletion,
                                  isingle_gene_deletion)


class Analysis(Enum):
    FBA = FBA
    milpFBA = milpFBA
    pFBA = pFBA
    RFBA = RFBA
    SRFBA = SRFBA
    SynchronousBoolean = SimBool
    BooleanMILP = milpBool

    @classmethod
    def has_analysis(cls, analysis):
        try:
            return cls[analysis]

        except KeyError:

            return False

    @classmethod
    def get(cls, analysis, default=FBA):
        try:
            return cls[analysis]

        except KeyError:

            return default
