from typing import Any, TYPE_CHECKING, Dict, Union, Generator, Tuple

import pandas as pd

from mewpy.germ.algebra import Expression, parse_expression
from mewpy.util.utilities import generator
from mewpy.germ.models.serialization import serialize
from mewpy.util.history import recorder
from mewpy.io.engines.engines_utils import expression_warning
from .variable import Variable, variables_from_symbolic


# Preventing circular dependencies that only happen due to type checking
if TYPE_CHECKING:
    from .regulator import Regulator
    from .target import Target
    from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel


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
        namely the possible coefficients of the target for the default (or not) coefficients of the regulators.

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
            target_str = f'{self.target.id} || '

        else:
            target_str = ''

        expression_str = ' + '.join([f'{coefficient} = {expression.to_string()}'
                                     for coefficient, expression in self.regulatory_events.items()
                                     if not expression.is_none])

        return target_str + expression_str

    def _interaction_to_html(self):
        """
        It returns a html dict representation.
        """
        html_dict = {'Target': self.target,
                     'Regulators': ', '.join(self.regulators),
                     'Regulatory events': ', '.join([f'{coefficient} = {expression.to_string()}'
                                                     for coefficient, expression in self.regulatory_events.items()
                                                     if not expression.is_none])}
        return html_dict

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------
    @serialize('regulatory_events', 'regulatory_events', '_regulatory_events')
    @property
    def regulatory_events(self) -> Dict[Union[float, int], Expression]:
        """
        A dictionary comprehending coefficient-expression key-value pairs.
        The expression when evaluated to true results into the corresponding coefficient that must be associated
        with the target coefficients
        """
        return self._regulatory_events.copy()

    @serialize('target', 'target', '_target')
    @property
    def target(self) -> 'Target':
        """
        The target variable that is regulated by the regulators/expressions logic
        """
        return self._target

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------
    @regulatory_events.setter
    @recorder
    def regulatory_events(self, value: Dict[Union[float, int], Expression]):
        """
        Setting a dictionary comprehending coefficient-expression key-value pairs
        to replace the current regulatory events.
        This setter adds and removes regulatory events.
        It also handles the linear problems associated to this interaction
        """
        if not value:
            value = {0.0: Expression()}

        for coefficient, expression in self.regulatory_events.items():
            self.remove_regulatory_event(coefficient=coefficient,
                                         remove_orphans=True, history=False)

        for coefficient, expression in value.items():
            self.add_regulatory_event(coefficient=coefficient, expression=expression, history=False)

    @target.setter
    @recorder
    def target(self, value: 'Target'):
        """
        Setting a new Target variable.
        This setter adds and removes target variable.
        It also handles the linear problems associated to this interaction
        """
        self.remove_target(remove_from_model=True, history=False)
        self.add_target(value, history=False)

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------
    @property
    def regulators(self) -> Dict[str, 'Regulator']:
        """
        It returns a dictionary with the regulators associated to this interaction
        """
        return {reg_id: regulator
                for expression in self.yield_expressions()
                for reg_id, regulator in expression.variables.items()}

    @property
    def regulatory_truth_table(self) -> pd.DataFrame:
        """
        It calculates the truth table for this interaction based on the regulatory events.
        The truth table is composed by the combination of values taken by empty, numeric and symbolic variables
        available in the regulatory events of this interaction.
        That is, the regulatory table is a pandas DataFrame with the possible coefficients
        of the target variable for the default (or not) coefficients of the regulators.

        :return: It returns a pandas DataFrame with the truth table
        """
        truth_tables = []

        for coefficient, expression in self._regulatory_events.items():
            df = expression.truth_table(strategy='max', coefficient=coefficient)

            df.index = [self.target.id if self.target else 'result'] * df.shape[0]

            truth_tables.append(df)

        regulatory_events = pd.concat(truth_tables)

        return regulatory_events

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------
    def yield_regulatory_events(self) -> Generator[Tuple[Union[float, int], Expression], None, None]:
        """
        It yields all regulatory events
        """
        return ((coefficient, expression) for coefficient, expression in self._regulatory_events.items())

    def yield_expressions(self) -> Generator[Expression, None, None]:
        """
        It yields all expressions
        """
        return generator(self._regulatory_events)

    def yield_regulators(self) -> Generator['Regulator', None, None]:
        """
        It yields all regulators
        """
        return generator(self.regulators)

    # -----------------------------------------------------------------------------
    # Polymorphic constructors
    # -----------------------------------------------------------------------------
    @classmethod
    def from_string(cls,
                    identifier: Any,
                    rule: str,
                    target: 'Target',
                    coefficient: Union[float, int] = 1.0,
                    model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
                    **kwargs) -> 'Interaction':
        """
        A regulatory interaction is regularly associated with a target and the regulatory
        events modelling the coefficients of this target variable.

        A regulatory interaction can be assembled from a string object encoding an algebra expression

        :param identifier: The interaction identifier
        :param rule: A string object encoding an algebra expression
        :param target: The target variable that is regulated by the regulators/expressions logic
        :param coefficient: the coefficient that should be taken by the target
        whether this expression is evaluated positively
        :param model: If a model is provided, this interaction will be associated to the model
        :param kwargs: Additional arguments
        :return: A new interaction object
        """
        try:

            symbolic = parse_expression(rule)

            regulators = variables_from_symbolic(symbolic=symbolic, types=('regulator',), model=model)

            expression = Expression(symbolic=symbolic, variables=regulators)

        except SyntaxError as exc:

            expression_warning(f'{rule} cannot be parsed')

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

        """
        A regulatory interaction is regularly associated with a target and the regulatory
        events modelling the coefficients of this target variable.

        A regulatory interaction can be assembled from an expression object encoding an algebra expression

        :param identifier: The interaction identifier
        :param expression: An expression object encoding an algebra expression
        :param target: The target variable that is regulated by the regulators/expressions logic
        :param coefficient: the coefficient that should be taken by the target
        whether this expression is evaluated positively
        :param model: If a model is provided, this interaction will be associated to the model
        :param kwargs: Additional arguments
        :return: A new interaction object
        """

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
        """
        It performs an update operation to this interaction.
        The update operation is similar to a dictionary update.

        It also handles the linear problems associated to this interaction

        Note that, some update operations are not registered in history.
        It is strongly advisable to use update outside history context manager

        :param target: the target variable that is regulated by the regulators/expressions logic
        :param regulatory_events: A dictionary comprehending coefficient-expression pairs.
        The expression when evaluated to true results into the corresponding coefficient
        that must be associated with the target coefficients
        :param kwargs: Other arguments for the base variable, such as identifier, name, etc
        :return:
        """
        super(Interaction, self).update(**kwargs)

        if target is not None:
            self.target: 'Target' = target

        if regulatory_events is not None:
            self.regulatory_events: Dict[Union[float, int], Expression] = regulatory_events

    def add_target(self, target: 'Target', history=True):
        """
        It adds a new target to this regulatory interaction.
        If a target is already associated with this interaction, the current target will be removed and replaced by the
        new target.

        It also handles the linear problems associated to this interaction

        This operation can be reverted using the model history
        :param target: the target variable that is regulated by the regulators/expressions logic
        :param history: Whether to register this operation in the model history
        """
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
                self.model.add(self.target, comprehensive=False, history=False)

                self.model.notify()

    def remove_target(self, remove_from_model=True, history=True):
        """
        It removes the current target of this regulatory interaction.

        It also handles the linear problems associated to this interaction

        This operation can be reverted using the model history
        :param remove_from_model: Whether this operation should be reflected in the model too.
        It is usually advisable to remove a target from the interaction and model at the same time
        :param history: Whether to register this operation in the model history
        """
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
                    self.model.remove(target, remove_orphans=False, history=False)

                self.model.notify()

    def add_regulatory_event(self,
                             coefficient: Union[float, int],
                             expression: Expression,
                             history=True):
        """
        It adds a new regulatory event to this interaction.
        If a regulatory event with the same coefficient is already associated with this interaction, the current
        regulatory event will be removed and replaced by the new regulatory event.

        It also handles the linear problems associated to this interaction

        This operation can be reverted using the model history
        :param coefficient: the coefficient that should be taken by the target for this expression
        :param expression: An expression object encoding an algebra expression
        :param history: Whether to register this operation in the model history
        :return: None
        """

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
            self.model.add(*to_add, comprehensive=False, history=False)

            self.model.notify()

    def remove_regulatory_event(self,
                                coefficient: Union[float, int],
                                remove_orphans=True,
                                history=True):
        """
        It removes a regulatory event from this interaction.
        If the expression is not associated with the coefficient, nothing will happen.

        It also handles the linear problems associated to this interaction

        This operation can be reverted using the model history
        :param coefficient: the coefficient that should be taken by the target for this expression
        :param remove_orphans: Whether to remove orphaned regulators from the model
        :param history: Whether to register this operation in the model history
        :return:
        """
        expression = self.regulatory_events.get(coefficient)
        if expression is None:
            raise ValueError(f'No regulatory event associated with coefficient {coefficient}')

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
                self.model.remove(*to_remove, remove_orphans=False, history=False)

            self.model.notify()
