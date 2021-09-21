from collections import namedtuple
from typing import Type, Union, TYPE_CHECKING

# TODO: this module largely depends on pandas dataframes. Should it be set as package requirement?
# noinspection PyPackageRequirements
from pandas import Series, DataFrame, concat

from mewpy.util.constants import ModelConstants
from .solution import ModelSolutionInterface

if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


# TODO: methods stubs and type hinting
class PolymorphicSolution:

    @classmethod
    def from_cobamp(cls: Type['ModelSolution'],
                    method,
                    solution,
                    **kwargs):
        return cls(method=method,
                   x=solution.var_values(),
                   objective_value=solution.objective_value(),
                   status=solution.status(),
                   **kwargs)

    @classmethod
    def from_solver(cls: Type['ModelSolution'],
                    method,
                    solution,
                    **kwargs):
        return cls(method=method,
                   x=solution.values,
                   objective_value=solution.fobj,
                   status=solution.status.value.lower(),
                   reduced_costs=solution.reduced_costs,
                   shadow_prices=solution.shadow_prices,
                   **kwargs)


# TODO: methods stubs and type hinting
class ModelSolution(ModelSolutionInterface, PolymorphicSolution):

    def __init__(self,
                 method,
                 x,
                 objective_value,
                 status,
                 objective_direction='maximize',
                 reduced_costs=None,
                 shadow_prices=None,
                 model=None,
                 simulator=None,
                 tol=ModelConstants.TOLERANCE):

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

    @staticmethod
    def _filter_mid_term_variables(x, model):

        if model is None:
            return x

        new_x = {}

        for variable, value in x.items():

            model_variable = model.get(variable, None)

            if model_variable is not None:
                new_x[variable] = value

        return new_x

    @property
    def objective_direction(self):
        return self._objective_direction

    @property
    def method(self):
        return self._method

    @property
    def model(self) -> Union['Model', 'MetabolicModel', 'RegulatoryModel']:
        return self._model

    @property
    def objective(self):

        if self.model:

            if hasattr(self.model, 'objective'):
                return ' + '.join([obj.id for obj in self.model.objective])

        return None

    @property
    def objective_value(self):
        return self._objective_value

    @property
    def simulator(self):
        return self._simulator

    @property
    def status(self):
        return self._status

    @property
    def x(self):
        return self._filter_mid_term_variables(x=self._x, model=self.model)

    @property
    def shadow_prices(self):
        return self._filter_mid_term_variables(x=self._shadow_prices, model=self.model)

    @property
    def reduced_costs(self):
        return self._filter_mid_term_variables(x=self._reduced_costs, model=self.model)

    def _get_variable_info(self, variable):

        identifier = variable.id

        x = self._x.get(variable.id, None)

        v_type = ', '.join(variable.types)

        return identifier, v_type, x

    def _metabolic_frame(self):

        results = {}

        for variable in self.model.yield_reactions():
            _id, v_type, x = self._get_variable_info(variable)

            results[_id] = (variable.id, v_type, x)

        return DataFrame.from_dict(results,
                                   orient='index',
                                   columns=['reaction', 'variable type', 'flux'])

    def _regulatory_frame(self):

        results = {}

        for variable in self.model.yield_regulators():
            _id, v_type, x = self._get_variable_info(variable)

            results[_id] = (variable.id, v_type, x)

        for variable in self.model.yield_targets():
            _id, v_type, x = self._get_variable_info(variable)

            results[_id] = (variable.id, v_type, x)

        return DataFrame.from_dict(results,
                                   orient='index',
                                   columns=['regulatory variable',
                                            'variable type',
                                            'expression coefficient'])

    def _metabolic_environmental_conditions_frame(self):

        results = {}

        for variable in self.model.yield_exchanges():

            if variable.is_reaction():
                _id, v_type, x = self._get_variable_info(variable)

                metabolite = next(iter(variable.metabolites.keys()))

                lb, ub = variable.bounds

                results[_id] = (_id, v_type, metabolite, lb, ub, x)

        return DataFrame.from_dict(results,
                                   orient='index',
                                   columns=['exchange', 'variable type', 'metabolite',
                                            'lower bound', 'upper bound',
                                            'flux'])

    def _regulatory_environmental_conditions_frame(self):

        results = {}

        for variable in self.model.yield_environmental_stimuli():

            if variable.is_regulator():
                _id, v_type, x = self._get_variable_info(variable)

                lb, ub = variable.coefficient.bounds

                results[_id] = (_id, v_type, lb, ub, x)

        return DataFrame.from_dict(results,
                                   orient='index',
                                   columns=['regulatory variable', 'variable type',
                                            'minimum coefficient', 'maximum coefficient',
                                            'expression coefficient'])

    def _regulatory_inputs_frame(self):

        results = {}

        for variable in self.model.yield_environmental_stimuli():

            _id, v_type, x = self._get_variable_info(variable)

            if variable.is_regulator():
                results[_id] = (_id, v_type, x)

        return DataFrame.from_dict(results,
                                   orient='index',
                                   columns=['regulator', 'variable type', 'expression coefficient'])

    def _regulatory_outputs_frame(self):

        results = {}

        for variable in self.model.yield_targets():
            _id, v_type, x = self._get_variable_info(variable)

            results[_id] = (_id, v_type, x)

        return DataFrame.from_dict(results,
                                   orient='index',
                                   columns=['target', 'variable type', 'expression coefficient'])

    def _metabolic_inputs_frame(self):

        results = {}

        for variable in self.model.yield_exchanges():

            _id, v_type, x = self._get_variable_info(variable)

            if variable.is_reaction() and x < -self.tol:
                metabolite = next(iter(variable.metabolites.keys()))

                results[_id] = (_id, v_type, metabolite, x)

        return DataFrame.from_dict(results,
                                   orient='index',
                                   columns=['reaction', 'variable type', 'metabolite', 'flux'])

    def _metabolic_outputs_frame(self):

        results = {}

        for variable in self.model.yield_exchanges():

            _id, v_type, x = self._get_variable_info(variable)

            if variable.is_reaction() and x > self.tol:
                metabolite = next(iter(variable.metabolites.keys()))

                results[_id] = (_id, v_type, metabolite, x)

        return DataFrame.from_dict(results,
                                   orient='index',
                                   columns=['reaction', 'variable type', 'metabolite', 'flux'])

    def _metabolic_summary_frame(self):

        results = {}

        for variable in self.model.yield_exchanges():

            _id, v_type, x = self._get_variable_info(variable)

            if variable.is_reaction() and x < -self.tol:
                metabolite = next(iter(variable.metabolites.keys()))

                results[_id] = (_id, v_type, metabolite, 'input', x)

            if variable.is_reaction() and x > self.tol:
                metabolite = next(iter(variable.metabolites.keys()))

                results[_id] = (_id, v_type, metabolite, 'output', x)

        return DataFrame.from_dict(results,
                                   orient='index',
                                   columns=['reaction', 'variable type', 'metabolite', 'role', 'flux'])

    def _regulatory_summary_frame(self):

        results = {}

        for variable in self.model.yield_environmental_stimuli():

            _id, v_type, x = self._get_variable_info(variable)

            if variable.is_regulator():
                results[_id] = (_id, v_type, 'input', x)

        for variable in self.model.yield_targets():
            _id, v_type, x = self._get_variable_info(variable)

            results[_id] = (_id, v_type, 'output', x)

        return DataFrame.from_dict(results,
                                   orient='index',
                                   columns=['regulatory variable', 'variable type', 'role', 'expression coefficient'])

    def _objective_frame(self):

        return DataFrame([[self.objective_value, self.objective_direction]],
                         index=[self.objective],
                         columns=['value', 'direction'])

    def _get_frame(self, to, dimension):

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
            return DataFrame()

    def _frame_builder(self, to, dimensions):

        fields = tuple(dimensions) + ('frame',)

        Frame = namedtuple('Frame', fields)

        frames = [self._get_frame(to=to, dimension=dimension)
                  for dimension in dimensions]

        if frames:
            frame = concat(frames,
                           axis=1,
                           join='outer',
                           keys=dimensions)

        else:
            frame = DataFrame()

        frames.append(frame)

        # noinspection PyArgumentList
        return Frame(*frames)

    def to_frame(self, dimensions=None, environmental_conditions=False):

        if not dimensions:
            dimensions = self.model.types

        if environmental_conditions:
            return self._frame_builder(to='environmental_conditions', dimensions=dimensions)

        return self._frame_builder(to='frame', dimensions=dimensions)

    def to_series(self):

        return Series(self.x)

    def _summary_builder(self, dimensions):

        fields = ('inputs', 'outputs', 'objective', 'frame') + tuple(dimensions)

        Summary = namedtuple('Summary', fields)

        frames = []

        inputs_frames = [self._get_frame(to='inputs', dimension=dimension)
                         for dimension in dimensions]

        if inputs_frames:
            inputs_frame = concat(inputs_frames,
                                  axis=1,
                                  join='outer',
                                  keys=dimensions)

        else:
            inputs_frame = DataFrame()

        frames.append(inputs_frame)

        outputs_frames = [self._get_frame(to='outputs', dimension=dimension)
                          for dimension in dimensions]

        if outputs_frames:
            outputs_frame = concat(outputs_frames,
                                   axis=1,
                                   join='outer',
                                   keys=dimensions)

        else:
            outputs_frame = DataFrame()

        frames.append(outputs_frame)

        objective_frame = self._get_frame(to='objective', dimension='objective')

        frames.append(objective_frame)

        summary_frames = [self._get_frame(to='summary', dimension=dimension)
                          for dimension in dimensions]

        if summary_frames:
            frame = concat(summary_frames,
                           axis=1,
                           join='outer',
                           keys=dimensions)

        else:
            frame = DataFrame()

        frames.append(frame)
        frames.extend(summary_frames)

        # noinspection PyArgumentList
        return Summary(*frames)

    def to_summary(self, dimensions=None):

        if not dimensions:
            dimensions = self.model.types

        return self._summary_builder(dimensions=dimensions)

    def __repr__(self):
        return f'{self.method} Solution'

    def __str__(self):
        return f'{self.method} {self.status} solution: {self.objective_value}'
