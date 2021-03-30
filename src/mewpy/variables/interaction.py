from typing import Any, TYPE_CHECKING, Dict, Union, Generator, Tuple

# TODO: this module depends on pandas dataframes. Should it be set as package requirement?
try:
    # noinspection PyPackageRequirements
    from pandas import concat, DataFrame

except ImportError:

    concat = False
    DataFrame = False

from mewpy.algebra import Expression, parse_expression
from mewpy.lp import Notification
from mewpy.util.utilities import generator
from mewpy.util.serialization import serialize
from mewpy.util.history import recorder
from mewpy.io.engines.engines_utils import expression_warning
from .variable import Variable, variables_from_symbolic


# Preventing circular dependencies that only happen due to type checking
if TYPE_CHECKING:
    from .regulator import Regulator
    from .target import Target
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


# TODO: methods stubs
class Interaction(Variable, variable_type='interaction', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier: Any,
                 target: 'Target' = None,
                 regulatory_events: Dict[Union[float, int], Expression] = None,
                 **kwargs):

        """
        A regulatory interaction is regularly associated with a target and the regulatory
        events modelling the coefficients of this target variable.

        The regulatory events determine the multiple coefficients that the target variable can take
        depending on the logic and conditions encoded into an expression object.

        The expression must hold the multiple variables/regulators that participate in the regulatory event/expression.

        The set of regulators is inferred from the regulatory events.
        Also, the regulatory events are used to generate a regulatory table,
        namely the possible coefficients of the target for the active (or not) coefficients of the regulators.

        :param identifier: the interaction identifier
        :param target: the target variable that is regulated by the regulators/expressions logic
        :param regulatory_events: A dictionary comprehending coefficient-expression pairs.
        The expression when evaluated to true results into the corresponding coefficient
        that must be associated with the target coefficients

        """

        if not target:
            target = None

        if not regulatory_events:
            regulatory_events = {0.0: Expression()}

        self._regulatory_events = {}
        self._target = None

        super().__init__(identifier,
                         **kwargs)

        # it handles addition of a target. It fulfills the correct containers/attributes
        self.add_target(target, history=False)

        # it handles setting a regulatory event. It fulfills the correct containers/attributes for the regulators
        # associated with the expression og a given event
        for coefficient, expression in regulatory_events.items():
            self.add_regulatory_event(coefficient=coefficient, expression=expression, history=False)

        # There is an alternative constructor for building an interaction with a stringify-like expression rule.
        # Parsing takes the hard job to infer the correct relation between regulators

    # -----------------------------------------------------------------------------
    # Variable type manager
    # -----------------------------------------------------------------------------

    @property
    def types(self):

        # noinspection PyUnresolvedReferences
        _types = {Interaction.variable_type}

        _types.update(super(Interaction, self).types)

        return _types

    # -----------------------------------------------------------------------------
    # Built-in
    # -----------------------------------------------------------------------------

    def __str__(self):

        if self.target:
            string_repr = f'{self.target.id}: '

        else:
            string_repr = ''

        for expression in self.yield_expressions():
            string_repr += expression.to_string()

        return string_repr

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------

    @serialize('regulatory_events', 'regulatory_events', '_regulatory_events')
    @property
    def regulatory_events(self) -> Dict[Union[float, int], Expression]:
        return self._regulatory_events.copy()

    @serialize('target', 'target', '_target')
    @property
    def target(self) -> 'Target':
        return self._target

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------

    @regulatory_events.setter
    @recorder
    def regulatory_events(self, value: Dict[Union[float, int], Expression]):

        for coefficient, expression in self.regulatory_events.items():
            self.remove_regulatory_event(coefficient=coefficient, expression=expression,
                                         remove_orphans=True, history=False)

        for coefficient, expression in value.items():
            self.add_regulatory_event(coefficient=coefficient, expression=expression, history=False)

    @target.setter
    @recorder
    def target(self, value: 'Target'):

        self.remove_target(remove_from_model=True, history=False)
        self.add_target(value, history=False)

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------

    @property
    def regulators(self) -> Dict[str, 'Regulator']:

        return {reg_id: regulator
                for expression in self.yield_expressions()
                for reg_id, regulator in expression.variables.items()}

    def _compute_regulatory_table(self, active_states=True):

        if concat is False:
            raise RuntimeError('pandas must be installed to compute regulatory tables')

        truth_tables = []

        for coefficient, expression in self._regulatory_events.items():
            df = expression.truth_table(active_states=active_states, coefficient=coefficient)

            df.index = [self.target.id if self.target else 'result'] * df.shape[0]

            truth_tables.append(df)

        # noinspection PyCallingNonCallable
        regulatory_events = concat(truth_tables)

        regulatory_events = regulatory_events[['result'] + [col for col in regulatory_events.columns
                                                            if col != 'result']]

        return regulatory_events

    @property
    def regulatory_table(self) -> DataFrame:

        return self._compute_regulatory_table(active_states=True)

    def full_regulatory_table(self) -> DataFrame:

        return self._compute_regulatory_table(active_states=False)

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------

    def yield_regulatory_events(self) -> Generator[Tuple[Union[float, int], Expression], None, None]:

        return ((coefficient, expression) for coefficient, expression in self._regulatory_events.items())

    def yield_expressions(self) -> Generator[Expression, None, None]:

        return generator(self._regulatory_events)

    def yield_regulators(self) -> Generator['Regulator', None, None]:

        return generator(self.regulators)

    # -----------------------------------------------------------------------------
    # Polymorphic constructors
    # -----------------------------------------------------------------------------

    @classmethod
    def from_string(cls,
                    identifier: Any,
                    stringify_rule: str,
                    target: 'Target',
                    coefficient: Union[float, int] = 1.0,
                    model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                    **kwargs) -> 'Interaction':

        try:

            symbolic = parse_expression(stringify_rule)

            regulators = variables_from_symbolic(symbolic=symbolic, types=('regulator',), model=model)

            expression = Expression(symbolic=symbolic, variables=regulators)

        except SyntaxError as exc:

            expression_warning(f'{stringify_rule} cannot be parsed')

            raise exc

        instance = cls(identifier=identifier,
                       target=target,
                       regulatory_events={coefficient: expression},
                       model=model,
                       **kwargs)

        return instance

    @classmethod
    def from_expression(cls,
                        identifier: Any,
                        expression: Expression,
                        target: 'Target',
                        coefficient: Union[float, int] = 1.0,
                        model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                        **kwargs) -> 'Interaction':

        if not isinstance(expression, Expression):
            raise TypeError(f'expression must be an {Expression} object')

        instance = cls(identifier=identifier,
                       target=target,
                       regulatory_events={coefficient: expression},
                       model=model,
                       **kwargs)

        return instance

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------

    def update(self,
               regulatory_events: Dict[Union[float, int], Expression] = None,
               target: 'Target' = None,
               **kwargs):

        super(Interaction, self).update(**kwargs)

        if target is not None:
            self.target: 'Target' = target

        if regulatory_events is not None:
            self.regulatory_events: Dict[Union[float, int], Expression] = regulatory_events

    def add_target(self, target: 'Target', history=True):

        if history:
            self.history.queue_command(undo_func=self.remove_target,
                                       undo_kwargs={'remove_from_model': True,
                                                    'history': False},
                                       func=self.add_target,
                                       kwargs={'target': target,
                                               'history': history})

        if self.target and not target:
            return self.remove_target(remove_from_model=True, history=False)

        elif not self.target and not target:
            return

        else:
            self._target = target
            self.target._interaction = self

            if self.model:
                self.model.add((self.target,), 'target', comprehensive=False, history=False)

                notification = Notification(content=(self,),
                                            content_type='interactions',
                                            action='add')

                self.model.notify(notification)

    def remove_target(self, remove_from_model=True, history=True):

        if history:
            self.history.queue_command(undo_func=self.add_target,
                                       undo_kwargs={'target': self.target,
                                                    'history': False},
                                       func=self.remove_target,
                                       kwargs={'remove_from_model': remove_from_model,
                                               'history': history})

        if self.target:
            target = self.target

            self.target._interaction = None
            self._target = None

            if self.model:

                if remove_from_model:
                    self.model.remove((target,), 'target', remove_orphans=False, history=False)

                notification = Notification(content=(self,),
                                            content_type='interactions',
                                            action='add')

                self.model.notify(notification)

    def add_regulatory_event(self,
                             coefficient: Union[float, int],
                             expression: Expression,
                             history=True):

        if not isinstance(expression, Expression):
            raise TypeError(f'expression must be an {Expression} object. '
                            f'To set None, provide an empty Expression()')

        if history:
            self.history.queue_command(undo_func=self.remove_regulatory_event,
                                       undo_kwargs={'coefficient': coefficient,
                                                    'expression': expression,
                                                    'remove_orphans': True,
                                                    'history': False},
                                       func=self.add_regulatory_event,
                                       kwargs={'coefficient': coefficient,
                                               'expression': expression,
                                               'history': history})

        to_add = []

        for regulator in expression.variables.values():
            regulator.update(interactions={self.id: self},
                             model=self.model)

            if self.model:
                if regulator.id not in self.model.regulators:
                    to_add.append(regulator)

        self._regulatory_events[coefficient] = expression

        if self.model:
            self.model.add(to_add, 'regulator', comprehensive=False, history=False)

            notification = Notification(content=(self,),
                                        content_type='interactions',
                                        action='add')

            self.model.notify(notification)

    def remove_regulatory_event(self,
                                coefficient: Union[float, int],
                                expression: Expression,
                                remove_orphans=True,
                                history=True):

        if not isinstance(expression, Expression):
            raise TypeError(f'expression must be an {Expression} object. '
                            f'To set None, provide an empty Expression()')

        if history:
            self.history.queue_command(undo_func=self.add_regulatory_event,
                                       undo_kwargs={'coefficient': coefficient,
                                                    'expression': expression,
                                                    'history': False},
                                       func=self.remove_regulatory_event,
                                       kwargs={'coefficient': coefficient,
                                               'expression': expression,
                                               'remove_orphans': remove_orphans,
                                               'history': history})

        to_remove = []

        for regulator in expression.variables.values():

            # noinspection PyProtectedMember
            del regulator._interactions[self.id]

            # noinspection PyProtectedMember
            if not regulator._interactions:
                to_remove.append(regulator)

        self._regulatory_events.pop(coefficient)

        if self.model:

            if remove_orphans:
                self.model.remove(to_remove, 'regulator', remove_orphans=False, history=False)

            notification = Notification(content=(self,),
                                        content_type='interactions',
                                        action='add')

            self.model.notify(notification)
