import os
import pandas as pd
import numpy as np
from mewpy.utils.parsing import BOOLEAN_SPECIAL_CHARS


def read_tabular_regulatory_model(file_or_filepath, sep=None, id_col=None, rule_col=None, aliases_cols=None, header=None,
                                  special_chars=None, compartments_ids=None):

    if not sep:
        sep = ','

    if not id_col:
        id_col = 0

    if not rule_col:
        rule_col = 1

    if not isinstance(file_or_filepath, str):
        raise TypeError("File or file path must be a string")

    if not isinstance(sep, str):
        raise TypeError("Separator must be a string")

    if not isinstance(id_col, int):
        raise TypeError(
            "The column of genes identifiers must be its index in the file")

    if not isinstance(rule_col, int):
        raise TypeError(
            "The column of regulatory rules must be its index in the file")

    if aliases_cols:

        if not isinstance(aliases_cols, list):
            raise TypeError(
                "The aliases columns must be a list of all its indexes in the file")

        if not isinstance(aliases_cols[0], int):
            raise TypeError(
                "The aliases columns must be a list of all its indexes in the file")

    else:
        aliases_cols = []

    if not isinstance(header, int) and header is not None:
        raise TypeError("The header must be a int value or None")

    if special_chars:

        if not isinstance(special_chars, list):
            raise TypeError(
                "The special chars must be a list of all chars to be parsed out from the file")

        if not isinstance(special_chars[0], str):
            raise TypeError(
                "The special chars must be a list of all chars to be parsed out from the file")

        for char in special_chars:

            if char not in BOOLEAN_SPECIAL_CHARS:

                new_char = char.replace(')','_').replace('(','_').replace(' ','_')
                BOOLEAN_SPECIAL_CHARS[char] = '_' + new_char + '_'

    if compartments_ids:

        if not isinstance(compartments_ids, list):
            raise TypeError("The compartments ids must be a list of all compartments identifiers to be parsed out from "
                            "the file")

        if not isinstance(compartments_ids[0], str):
            raise TypeError("The compartments ids must be a list of all compartments identifiers to be parsed out from "
                            "the file")

        for char in special_chars:

            if char not in BOOLEAN_SPECIAL_CHARS:

                new_char = char.replace(')', '_').replace('(', '_').replace(' ', '_')
                BOOLEAN_SPECIAL_CHARS[char] = '_' + new_char + '_'

    names = {id_col: 'ids', rule_col: 'rules'}
    names.update({j: 'aliases_' + str(i + 1)
                  for i, j in enumerate(aliases_cols)})

    csv = pd.read_csv(file_or_filepath, sep, header=header)

    csv = csv.replace(np.nan, '', regex=True)

    cols = []

    for j, col in enumerate(csv.columns):
        if j in names:
            cols.append(names[j])
        else:
            del csv[col]

    csv.columns = cols
    csv.index = csv.loc[:, 'ids']

    return csv

def read_tabular_aliases(file_or_filepath, sep=None, id_col=None, aliases_cols=None, header=None):

    if not sep:
        sep = ','

    if not id_col:
        id_col = 0

    if not aliases_cols:
        aliases_cols = [1]

    if not isinstance(file_or_filepath, str):
        raise TypeError("File or file path must be a string")

    if not isinstance(sep, str):
        raise TypeError("Separator must be a string")

    if not isinstance(id_col, int):
        raise TypeError("The column of genes identifiers must be its index in the file")

    if not isinstance(aliases_cols, list):
        raise TypeError("The aliases columns must be a list of all its indexes in the file")

    if not isinstance(aliases_cols[0], int):
        raise TypeError("The aliases columns must be a list of all its indexes in the file")

    if not isinstance(header, int) and header is not None:
        raise TypeError("The header must be a int value or None")

    names = {id_col: 'ids'}
    names.update({j: 'aliases_' + str(i + 1) for i, j in enumerate(aliases_cols)})

    csv = pd.read_csv(file_or_filepath, sep, header=header)

    csv = csv.replace(np.nan, '', regex=True)

    cols = []

    for j, col in enumerate(csv.columns):
        if j in names:
            cols.append(names[j])
        else:
            del csv[col]

    csv.columns = cols
    csv.index = csv.loc[:, 'ids']

    return csv

if __name__ == '__main__':

    file = "iMC1010v4.csv"
    os.chdir('../../../examples/models/')

    df = read_tabular_regulatory_model(
        file, sep=';', id_col=2, rule_col=3, header=None, special_chars=['X'])
    print(df.columns)
    print(df.loc['XylE', 'rules'])
    print(type(df.loc['XylE', 'rules']))
