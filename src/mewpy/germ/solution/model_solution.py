from typing import Union, TYPE_CHECKING, Dict, Optional, Tuple, Any

import pandas as pd

from .summary import Summary
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.solvers import Solution
    from mewpy.germ.variables import Variable
    from mewpy.germ.lp import LinearProblem
    from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel


class ModelSolution:
    """
    ModelSolution can be used to retrieve the results of a simulation using metabolic, regulatory
    or integrated analysis.
    It contains the main information about the solution:
        - method
        - variables values (x)
        - objective value
        - status
        - objective direction
        - reduced costs
        - shadow prices

    It also contains the model and the simulator used to obtain the solution.

    All solution attributes can be retrieved directly from a ModelSolution object.
    Alternatively, three export methods are available:
        - to_series: returns a pandas Series with the fluxes of the reactions and regulatory variables
        - to_frame: returns a pandas DataFrame with the fluxes of the reactions, regulatory variables
        and environmental conditions
        - to_summary: returns a pandas DataFrame with the input and output fluxes of the reactions, regulatory variables
         and the environmental conditions
    """
    def __init__(self,
                 method: str,
                 x: Dict[str, float],
                 objective_value: float,
                 status: str,
                 objective_direction: str = 'maximize',
                 reduced_costs: Dict[str, float] = None,
                 shadow_prices: Dict[str, float] = None,
                 model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                 simulator: 'LinearProblem' = None,
                 tol: float = ModelConstants.TOLERANCE):
        """
        ModelSolution can be used to retrieve the results of a simulation using metabolic, regulatory
        or integrated analysis.
        It contains the main information about the solution:
            - method
            - variables values (x)
            - objective value
            - status
            - objective direction
            - reduced costs
            - shadow prices

        It also contains the model and the simulator used to obtain the solution.

        All solution attributes can be retrieved directly from a ModelSolution object.
        Alternatively, three export methods are available:
            - to_series: returns a pandas Series with the fluxes of the reactions and regulatory variables
            - to_frame: returns a pandas DataFrame with the fluxes of the reactions, regulatory variables
            and environmental conditions
            - to_summary: returns a pandas DataFrame with the input and output fluxes of the reactions,
            regulatory variables and the environmental conditions

        :param method: analysis method used to obtain the solution
        :param x: dictionary with the variables values
        :param objective_value: objective value obtained in the simulation
        :param status: the status of the solution obtained from the solver
        :param objective_direction: the direction of the objective function. Default is 'maximize'
        :param reduced_costs: The reduced costs of the variables. Default is None
        :param shadow_prices: The shadow prices of the constraints. Default is None
        :param model: The model used to obtain the solution. Default is None
        :param simulator: The simulator used to obtain the solution. Default is None
        :param tol: The tolerance used to round the solution. Default is 1e-6
        """
        if not method:
            method = 'fba'

        if not x:
            x = {}

        if not objective_value:
            objective_value = 0.0

        if not status:
            status = 'feasible'

        if not objective_direction:
            objective_direction = 'maximize'

        if not reduced_costs:
            reduced_costs = {}

        if not shadow_prices:
            shadow_prices = {}

        self._method = method
        self._x = x
        self._objective_value = objective_value
        self._status = status
        self._objective_direction = objective_direction
        self._reduced_costs = reduced_costs
        self._shadow_prices = shadow_prices
        self._model = model
        self._simulator = simulator
        self.tol = tol

    # ---------------------------------
    # Buil-in
    # ---------------------------------
    def __str__(self):
        return f'{self.method} Solution\n Objective value: {self.objective_value}\n Status: {self.status}'

    def __repr__(self):
        return self.__str__()

    def _repr_html_(self):
        """
        It returns a html representation of the linear problem
        :return:
        """

        return f"""
        <table>
            <tr>
                <td>Method</td>
                <td>{self.method}</td>
            </tr>
            <tr>
                <td>Model</td>
                <td>{self.model}</td>
            </tr>
            <tr>
                <th>Objective</th>
                <td>{self.objective}</td>
            </tr>
            <tr>
                <th>Objective value</th>
                <td>{self.objective_value}</td>
            </tr>
            <tr>
                <th>Status</th>
                <td>{self.status}</td>
            </tr>
        </table>
        """

    @staticmethod
    def _filter_mid_term_variables(x: Dict[str, float],
                                   model: Union['Model', 'MetabolicModel', 'RegulatoryModel']) -> Dict[str, float]:
        """
        It filters out midterm variables from the solution.
        Internal use only.
        :param x: The solution dictionary
        :param model: The model used to obtain the solution
        :return: A dictionary with the solution without mid-term variables
        """
        if model is None:
            return x

        new_x = {}

        for variable, value in x.items():

            model_variable = model.get(variable, None)

            if model_variable is not None:
                new_x[variable] = value

        return new_x

    @property
    def objective_direction(self) -> str:
        """
        The direction of the objective function.
        :return: The direction of the objective function
        """
        return self._objective_direction

    @property
    def method(self) -> str:
        """
        The method used to obtain the solution.
        :return: The method used to obtain the solution
        """
        return self._method

    @property
    def model(self) -> Union['Model', 'MetabolicModel', 'RegulatoryModel']:
        """
        The model used to obtain the solution.
        :return: The model used to obtain the solution
        """
        return self._model

    @property
    def objective(self) -> Optional[str]:
        """
        A string representation of the objective function of the model.
        :return: A string representation of the objective function of the model
        """
        if self.model:

            if hasattr(self.model, 'objective'):
                return ' + '.join([obj.id for obj in self.model.objective])

        return

    @property
    def objective_value(self) -> float:
        """
        The objective value obtained in the simulation.
        :return: The objective value obtained in the simulation
        """
        return self._objective_value

    @property
    def simulator(self) -> 'LinearProblem':
        """
        The simulator used to obtain the solution.
        :return: A LinearProblem-like object defining the simulator used to obtain the solution
        """
        return self._simulator

    @property
    def status(self) -> str:
        """
        The status of the solution obtained from the solver.
        :return: The status of the solution obtained from the solver
        """
        return self._status

    @property
    def x(self) -> Dict[str, float]:
        """
        The variables' values obtained in the solution of the linear problem.
        :return: A dictionary with the variables values
        """
        return self._filter_mid_term_variables(x=self._x, model=self.model)

    @property
    def shadow_prices(self):
        """
        The shadow prices of the variables. How much of each variable would be needed to
        increase the objective value.
        :return: A dictionary with the shadow prices of the constraints
        """
        return self._filter_mid_term_variables(x=self._shadow_prices, model=self.model)

    @property
    def reduced_costs(self):
        """
        The reduced costs of the variables. The objective value to increase
        to assume a positive value in the optimal solution.
        :return: A dictionary with the reduced costs of the variables
        """
        return self._filter_mid_term_variables(x=self._reduced_costs, model=self.model)

    def _get_variable_info(self, variable: 'Variable') -> Tuple[Any, str, Optional[float]]:
        """
        It returns the information about a variable in the solution.
        Internal use only.
        :param variable:
        :return:
        """
        identifier = variable.id

        x = self._x.get(variable.id, None)

        v_type = ', '.join(variable.types)

        return identifier, v_type, x

    def _metabolic_frame(self) -> pd.DataFrame:
        """
        It returns a pandas DataFrame with the fluxes of the reactions and environmental conditions.
        Internal use only.
        :return: A pandas DataFrame with the fluxes of the reactions and environmental conditions
        """

        results = {}

        for variable in self.model.yield_reactions():
            _id, v_type, x = self._get_variable_info(variable)

            results[_id] = (variable.id, v_type, x)

        return pd.DataFrame.from_dict(results,
                                      orient='index',
                                      columns=['reaction', 'variable type', 'flux'])

    def _regulatory_frame(self) -> pd.DataFrame:
        """
        It returns a pandas DataFrame with the regulatory variables.
        Internal use only.
        :return: A pandas DataFrame with the regulatory variables
        """
        results = {}

        for variable in self.model.yield_regulators():
            _id, v_type, x = self._get_variable_info(variable)

            results[_id] = (variable.id, v_type, x)

        for variable in self.model.yield_targets():
            _id, v_type, x = self._get_variable_info(variable)

            results[_id] = (variable.id, v_type, x)

        return pd.DataFrame.from_dict(results,
                                      orient='index',
                                      columns=['regulatory variable',
                                               'variable type',
                                               'expression coefficient'])

    def _metabolic_environmental_conditions_frame(self) -> pd.DataFrame:
        """
        It returns a pandas DataFrame with the environmental conditions.
        Internal use only.
        :return: A pandas DataFrame with the environmental conditions
        """
        results = {}

        for variable in self.model.yield_exchanges():

            if variable.is_reaction():
                _id, v_type, x = self._get_variable_info(variable)

                metabolite = next(iter(variable.metabolites.keys()))

                lb, ub = variable.bounds

                results[_id] = (_id, v_type, metabolite, lb, ub, x)

        return pd.DataFrame.from_dict(results,
                                      orient='index',
                                      columns=['exchange', 'variable type', 'metabolite',
                                               'lower bound', 'upper bound',
                                               'flux'])

    def _regulatory_environmental_conditions_frame(self) -> pd.DataFrame:
        """
        It returns a pandas DataFrame with the environmental conditions.
        Internal use only.
        :return: A pandas DataFrame with the environmental conditions
        """
        results = {}

        for variable in self.model.yield_environmental_stimuli():

            if variable.is_regulator():
                _id, v_type, x = self._get_variable_info(variable)

                lb, ub = variable.coefficients

                results[_id] = (_id, v_type, lb, ub, x)

        return pd.DataFrame.from_dict(results,
                                      orient='index',
                                      columns=['regulatory variable', 'variable type',
                                               'minimum coefficient', 'maximum coefficient',
                                               'expression coefficient'])

    def _regulatory_inputs_frame(self) -> pd.DataFrame:
        """
        It returns a pandas DataFrame with the inputs of the regulatory model.
        Internal use only.
        :return: A pandas DataFrame with the inputs of the regulatory model
        """
        results = {}

        for variable in self.model.yield_environmental_stimuli():

            _id, v_type, x = self._get_variable_info(variable)

            if variable.is_regulator():
                results[_id] = (_id, v_type, x)

        return pd.DataFrame.from_dict(results,
                                      orient='index',
                                      columns=['regulator', 'variable type', 'expression coefficient'])

    def _regulatory_outputs_frame(self) -> pd.DataFrame:
        """
        It returns a pandas DataFrame with the outputs of the regulatory model.
        Internal use only.
        :return: A pandas DataFrame with the outputs of the regulatory model
        """
        results = {}

        for variable in self.model.yield_targets():
            _id, v_type, x = self._get_variable_info(variable)

            results[_id] = (_id, v_type, x)

        return pd.DataFrame.from_dict(results,
                                      orient='index',
                                      columns=['target', 'variable type', 'expression coefficient'])

    def _metabolic_inputs_frame(self) -> pd.DataFrame:
        """
        It returns a pandas DataFrame with the inputs of the metabolic model.
        Internal use only.
        :return: A pandas DataFrame with the inputs of the metabolic model
        """
        results = {}

        for variable in self.model.yield_exchanges():

            _id, v_type, x = self._get_variable_info(variable)

            if variable.is_reaction() and x < -self.tol:
                metabolite = next(iter(variable.metabolites.keys()))

                results[_id] = (_id, v_type, metabolite, x)

        return pd.DataFrame.from_dict(results,
                                      orient='index',
                                      columns=['reaction', 'variable type', 'metabolite', 'flux'])

    def _metabolic_outputs_frame(self) -> pd.DataFrame:
        """
        It returns a pandas DataFrame with the outputs of the metabolic model.
        Internal use only.
        :return: A pandas DataFrame with the outputs of the metabolic model
        """
        results = {}

        for variable in self.model.yield_exchanges():

            _id, v_type, x = self._get_variable_info(variable)

            if variable.is_reaction() and x > self.tol:
                metabolite = next(iter(variable.metabolites.keys()))

                results[_id] = (_id, v_type, metabolite, x)

        return pd.DataFrame.from_dict(results,
                                      orient='index',
                                      columns=['reaction', 'variable type', 'metabolite', 'flux'])

    def _metabolic_summary_frame(self) -> pd.DataFrame:
        """
        It returns a pandas DataFrame with the summary of the metabolic model.
        Internal use only.
        :return: A pandas DataFrame with the summary of the metabolic model
        """
        results = {}

        for variable in self.model.yield_exchanges():

            _id, v_type, x = self._get_variable_info(variable)

            if variable.is_reaction() and x < -self.tol:
                metabolite = next(iter(variable.metabolites.keys()))

                results[_id] = (_id, v_type, metabolite, 'input', x)

            if variable.is_reaction() and x > self.tol:
                metabolite = next(iter(variable.metabolites.keys()))

                results[_id] = (_id, v_type, metabolite, 'output', x)

        return pd.DataFrame.from_dict(results,
                                      orient='index',
                                      columns=['reaction', 'variable type', 'metabolite', 'role', 'flux'])

    def _regulatory_summary_frame(self) -> pd.DataFrame:
        """
        It returns a pandas DataFrame with the summary of the regulatory model.
        Internal use only.
        :return: A pandas DataFrame with the summary of the regulatory model
        """
        results = {}

        for variable in self.model.yield_environmental_stimuli():

            _id, v_type, x = self._get_variable_info(variable)

            if variable.is_regulator():
                results[_id] = (_id, v_type, 'input', x)

        for variable in self.model.yield_targets():
            _id, v_type, x = self._get_variable_info(variable)

            results[_id] = (_id, v_type, 'output', x)

        return pd.DataFrame.from_dict(results,
                                      orient='index',
                                      columns=['regulatory variable', 'variable type', 'role',
                                               'expression coefficient'])

    def _objective_frame(self) -> pd.DataFrame:
        """
        It returns a pandas DataFrame with the objective value.
        Internal use only.
        :return: A pandas DataFrame with the objective value
        """
        return pd.DataFrame([[self.objective_value, self.objective_direction]],
                            index=[self.objective],
                            columns=['value', 'direction'])

    def _get_frame(self, to: str, dimension: str) -> pd.DataFrame:
        """
        It returns a pandas DataFrame with the results of the model.
        Internal use only.
        :param to: The type of the results to be returned
        :param dimension: The dimension of the results to be returned
        :return: A pandas DataFrame with the results of the model
        """
        dim_to = f'{dimension}_{to}'

        if dim_to == 'regulatory_environmental_conditions':

            return self._regulatory_environmental_conditions_frame()

        elif dim_to == 'metabolic_environmental_conditions':

            return self._metabolic_environmental_conditions_frame()

        elif dim_to == 'regulatory_frame':

            return self._regulatory_frame()

        elif dim_to == 'metabolic_frame':

            return self._metabolic_frame()

        elif dim_to == 'metabolic_inputs':

            return self._metabolic_inputs_frame()

        elif dim_to == 'metabolic_outputs':

            return self._metabolic_outputs_frame()

        elif dim_to == 'regulatory_inputs':

            return self._regulatory_inputs_frame()

        elif dim_to == 'regulatory_outputs':

            return self._regulatory_outputs_frame()

        elif dim_to == 'objective_objective':

            return self._objective_frame()

        elif dim_to == 'regulatory_summary':

            return self._regulatory_summary_frame()

        elif dim_to == 'metabolic_summary':

            return self._metabolic_summary_frame()

        else:
            return pd.DataFrame()

    def to_frame(self, dimensions: Tuple[str, ...] = None) -> pd.DataFrame:
        """
        It returns a  pandas DataFrame having the summary results of the simulation.

        Example:
        >>> from mewpy.germ.analysis import FBA
        >>> from mewpy.io import read_sbml
        >>> model = read_sbml('e_coli_core.xml')
        >>> fba = FBA(model)
        >>> solution = fba.optimize()
        >>> solution.to_frame()

        :param dimensions: The dimensions of the results to be returned. If None, all dimensions are returned.
        possible values are: 'regulatory', 'metabolic', 'objective'
        :return: A pandas DataFrame having the summary results of the model
        """
        if not dimensions:
            dimensions = self.model.types

        frames = [self._get_frame(to='environmental_conditions', dimension=dimension)
                  for dimension in dimensions]

        if frames:
            return pd.concat(frames,
                             axis=1,
                             join='outer',
                             keys=dimensions)

        return pd.DataFrame()

    def to_series(self) -> pd.Series:
        """
        It returns a pandas Series with the values of the linear problem variables.
        Namely, the X variables can include reaction fluxes, metabolite concentrations, gene variables,
        and regulatory variables.
        :return: A pandas Series with the values of the linear problem variables
        """
        return pd.Series(self.x)

    def _summary_builder(self, dimensions: Tuple[str, ...]) -> Summary:
        """
        It returns a namedtuple with pandas DataFrame having all results of the model.
        Internal use only.
        :param dimensions: The dimensions of the results to be returned
        :return: A Summary with pandas DataFrame having all results of the model
        """
        frames = {}

        inputs_frames = [self._get_frame(to='inputs', dimension=dimension)
                         for dimension in dimensions]
        if inputs_frames:
            inputs_frame = pd.concat(inputs_frames,
                                     axis=1,
                                     join='outer',
                                     keys=dimensions)
        else:
            inputs_frame = pd.DataFrame()
        frames['inputs'] = inputs_frame

        outputs_frames = [self._get_frame(to='outputs', dimension=dimension)
                          for dimension in dimensions]
        if outputs_frames:
            outputs_frame = pd.concat(outputs_frames,
                                      axis=1,
                                      join='outer',
                                      keys=dimensions)
        else:
            outputs_frame = pd.DataFrame()
        frames['outputs'] = outputs_frame

        objective_frame = self._get_frame(to='objective', dimension='objective')
        frames['objective'] = objective_frame

        summary_frames = [self._get_frame(to='summary', dimension=dimension)
                          for dimension in dimensions]
        if summary_frames:
            frame = pd.concat(summary_frames,
                              axis=1,
                              join='outer',
                              keys=dimensions)

        else:
            frame = pd.DataFrame()
        frames['df'] = frame

        for i, dimension in enumerate(dimensions):
            frames[dimension] = summary_frames[i]

        return Summary(**frames)

    def to_summary(self, dimensions: Tuple[str, ...] = None) -> Summary:
        """
        It returns a summary with pandas DataFrame of the simulation.

        According to the dimensions, the Summary will have the following attributes:
            - inputs: A pandas DataFrame with the inputs of the metabolic and regulatory model
            - outputs: A pandas DataFrame with the outputs of the metabolic and regulatory model
            - objective: A pandas DataFrame with the objective value
            - metabolic: A pandas DataFrame with the summary of the metabolic model
            - regulatory: A pandas DataFrame with the summary of the regulatory model
            - df: A pandas DataFrame with the summary of the metabolic and regulatory models

        Example:
        >>> from mewpy.germ.analysis import FBA
        >>> from mewpy.io import read_sbml
        >>> model = read_sbml('e_coli_core.xml')
        >>> fba = FBA(model)
        >>> solution = fba.optimize()
        >>> solution.to_summary()

        :param dimensions: The dimensions of the results to be returned. If None, all dimensions are returned.
        possible values are: 'regulatory', 'metabolic', 'objective'
        :return: A namedtuple with pandas DataFrame having all results of the model
        """
        if not dimensions:
            dimensions = self.model.types

        return self._summary_builder(dimensions=dimensions)

    @classmethod
    def from_solver(cls,
                    method: str,
                    solution: 'Solution',
                    **kwargs) -> 'ModelSolution':
        """
        It returns a solution object from a solver solution.
        :param method: The method used to solve the problem
        :param solution: The solution object returned by the solver
        :param kwargs: Additional arguments
        :return: A new ModelSolution object
        """
        minimize = kwargs.pop('minimize', False)
        if minimize:
            objective_direction = 'minimize'
        else:
            objective_direction = 'maximize'

        return cls(method=method,
                   x=solution.values,
                   objective_value=solution.fobj,
                   objective_direction=objective_direction,
                   status=solution.status.value.lower(),
                   reduced_costs=solution.reduced_costs,
                   shadow_prices=solution.shadow_prices,
                   **kwargs)
