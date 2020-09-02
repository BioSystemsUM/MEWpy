"""
Implements a crossmodel simulation of solutions
"""
from mewpy.simulation import get_simulator
from mewpy.utils.utilities import Parser
from mewpy.optimization import Solution
import pandas as pd

map_problem_notation={"GeckoRKOProblem":"Protein" ,
                      "GeckoROUProblem":"Protein",
                      "GKOProblem"     :"Systematic",
                      "GOUProblem"     :"Systematic",
                      "OptRamProblem"  :"Systematic"
                     }

class NotationTranslator:
    def __init__(self, database, from_notation, to_notation,sep=';'):
        if isinstance(database,pd.DataFrame):
            self.db = database
        elif isinstance(database,str):
            self.db = pd.read_csv(database,sep=sep)
        self.from_notation = from_notation
        self.to_notation = to_notation

    def columns(self):
        return self.db.columns.tolist()

    def translate(self, value):
        """
        Translate a single value
        """
        if self.from_notation == self.to_notation:
            return value 

        l = self.db.loc[self.db[self.from_notation] == value][self.to_notation].tolist()
        if l and len(l)>0:
            return l[0]
        else:
            raise ValueError(f'Value {value} not found.')


    def translate_representation(self,representation,source_prefix="",destination_prefix=""):
        """
        Translates constraints.
        constraints are defined as a dictionary (OU) or a list (KO)
        """
        
        p = len(source_prefix)
        if isinstance(representation,list):
            return [destination_prefix+self.translate(value[p:]) for value in representation]
        else:
            return {destination_prefix+self.translate(value[p:]):level for value,level in representation.items()}

    def get_list(self,name):
        return self.db[name].tolist()



