"""
Implements a crossmodel simulation of solutions
"""
import pandas as pd

map_problem_notation = {"GeckoRKOProblem": "Protein",
                        "GeckoROUProblem": "Protein",
                        "GKOProblem": "Systematic",
                        "GOUProblem": "Systematic",
                        "OptRamProblem": "Systematic"
                        }


class NotationTranslator:
    def __init__(self, database, from_notation, to_notation, sep=';'):
        if isinstance(database, pd.DataFrame):
            self.db = database
        elif isinstance(database, str):
            self.db = pd.read_csv(database, sep=sep)
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

        la = self.db.loc[self.db[self.from_notation] == value][self.to_notation].tolist()
        if la and len(la) > 0:
            return la[0]
        else:
            raise ValueError(f'Value {value} not found.')

    def translate_representation(self, representation, source_prefix="", destination_prefix=""):
        """
        Translates constraints.
        constraints are defined as a dictionary (OU) or a list (KO)
        """

        p = len(source_prefix)
        if isinstance(representation, list):
            return [destination_prefix+self.translate(value[p:]) for value in representation]
        else:
            return {destination_prefix+self.translate(value[p:]): level for value, level in representation.items()}

    def get_list(self, name):
        return self.db[name].tolist()
