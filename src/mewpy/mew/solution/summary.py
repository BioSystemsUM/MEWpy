import pandas as pd


class Summary:

    def __init__(self,
                 inputs: pd.DataFrame = None,
                 outputs: pd.DataFrame = None,
                 objective: pd.DataFrame = None,
                 df: pd.DataFrame = None,
                 metabolic: pd.DataFrame = None,
                 regulatory: pd.DataFrame = None):
        """
        A summary of a ModelSolution

        :param inputs: the inputs of the model
        :param outputs: the outputs of the model
        :param objective: the objective of the model
        :param df: the data frame of the summary
        :param metabolic: the metabolic summary
        :param regulatory: the regulatory summary
        :return:
        """
        if inputs is None:
            inputs = pd.DataFrame()

        if outputs is None:
            outputs = pd.DataFrame()

        if objective is None:
            objective = pd.DataFrame()

        if df is None:
            df = pd.DataFrame()

        if metabolic is None:
            metabolic = pd.DataFrame()

        if regulatory is None:
            regulatory = pd.DataFrame()

        self.inputs = inputs
        self.outputs = outputs
        self.objective = objective
        self.df = df
        self.metabolic = metabolic
        self.regulatory = regulatory

    def _repr_html_(self):
        """
        It returns a html representation of the linear problem
        :return:
        """
        return self.df.to_html()
