import pandas as pd


def read_gene_expression_dataset(path,
                                 sep=',',
                                 gene_col=0,
                                 header=None,
                                 filter_nan=False):
    if not path:
        raise ValueError(f'{path} is not a valid path')

    if not sep:
        sep = ','

    if gene_col is None:
        gene_col = 0

    df = pd.read_csv(path, sep=sep, header=header)

    df.index = df.iloc[:, gene_col]

    df = df.drop(df.columns[gene_col], axis=1)

    if filter_nan:
        return df.dropna()

    return df


def read_coregflux_influence_matrix(path,
                                    sep=',',
                                    gene_col=0,
                                    header=None,
                                    filter_nan=False):
    if not path:
        raise ValueError(f'{path} is not a valid path')

    if not sep:
        sep = ','

    if gene_col is None:
        gene_col = 0

    df = pd.read_csv(path, sep=sep, header=header)

    df.index = df.iloc[:, gene_col]

    df = df.drop(df.columns[gene_col], axis=1)

    if filter_nan:
        return df.dropna()

    return df
