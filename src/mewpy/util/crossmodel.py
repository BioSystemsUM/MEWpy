"""
Implements a crossmodel simulation of solutions
"""
import pandas as pd


class NotationTranslator:
    def __init__(self, database, from_notation, to_notation, admissible=None, sep=';'):
        if isinstance(database, pd.DataFrame):
            self.db = database
        elif isinstance(database, str):
            self.db = pd.read_csv(database, sep=sep)
        self.from_notation = from_notation
        self.to_notation = to_notation
        self.admissible = admissible

    def columns(self):
        return self.db.columns.tolist()

    def translate(self, value):
        """
        Translate a single value
        """

        la = self.db.loc[self.db[self.from_notation].str.contains(value, na=False)][self.from_notation].tolist()
        lb = None
        if la and len(la) > 0:
            if len(la) == 1:
                lb = self.db.loc[self.db[self.from_notation] == la[0]][self.to_notation].tolist()
            else:
                for x in la:
                    tokens = x.split(' ')
                    if value in tokens:
                        lb = self.db.loc[self.db[self.from_notation] == x][self.to_notation].tolist()
                        break
            if not lb:
                raise ValueError(f'Value {value} not found.')
            if len(lb) > 1:
                raise ValueError(f'More than a value {value} found. {lb}')
            s = lb[0]
            tokens = s.split(' ')
            if tokens and len(tokens) == 1:
                return tokens[0]
            else:
                idx = 0
                while idx < len(tokens):
                    res = tokens[idx]
                    if res in self.admissible:
                        return res
                    idx += 1
                raise ValueError(f'Value {value} correspondences not in admissible.')
        else:
            raise ValueError(f'Value {value} not found.')

    def translate_representation(self, representation, source_prefix="", destination_prefix=""):
        """
        Translates constraints.
        constraints are defined as a dictionary (OU) or a list (KO)
        """

        p = len(source_prefix)
        if isinstance(representation, list):
            return [destination_prefix + self.translate(value[p:]) for value in representation]
        else:
            return {destination_prefix + self.translate(value[p:]): level for value, level in representation.items()}

    def get_list(self, name):
        return self.db[name].tolist()
