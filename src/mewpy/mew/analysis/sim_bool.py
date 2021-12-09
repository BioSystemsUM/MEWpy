from typing import Union

from mewpy.solvers.solver import Solver
from mewpy.mew.solution import ModelSolution
from mewpy.model import Model, MetabolicModel, RegulatoryModel


# TODO: type hinting and documentation
class SimBool:

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = True,
                 attach: bool = True):
        """
        Simulation of a Boolean network (SimBool) for a regulatory model.

        :param model: a mewpy Model, MetabolicModel, RegulatoryModel or all. The model is used to retrieve
        variables and constraints to the linear problem

        :param solver: A Solver, CplexSolver, GurobiSolver or OptLangSolver instance.
        Alternatively, the name of the solver is also accepted.
        The solver interface will be used to load and solve a linear problem in a given solver.
        If none, a new solver is instantiated. An instantiated solver may be used but it will be overwritten
        if build is true.

        :param build: Whether to build the linear problem upon instantiation. Default: False
        :param attach: Whether to attach the linear problem to the model upon instantiation. Default: False
        """

        self.model = model

    def build(self):

        pass

    def optimize(self, initial_state=None, to_dict=False):

        res = {}

        # solving regulatory rule by regulatory rule (synchronously though, as asynchronously would take too much time)
        # Each target has one and only one regulatory interaction
        # The regulatory model object only supports this strategy, though this can easily be altered. If so,
        # this solve method must be changed too

        if not initial_state:
            initial_state = {}

        else:
            initial_state = initial_state.copy()

        if self.model.is_regulatory():

            reg_state = {reg.id: reg.coefficient.active_coefficient
                         for reg in self.model.yield_regulators()}

            # update with the most recent state
            reg_state.update(initial_state)

            # interaction over model.interactions.

            for interaction in self.model.yield_interactions():

                # regulatory event over interaction.regulatory events.
                # an interaction can have multiple regulatory events.
                # Regularly, there is one event for each possible coefficient of the target variable.
                # If the regulatory event is evaluated to True or 1, this is the outcome coefficient of the target
                # By default the outcome of the target is a zero coefficient

                target_val = 0

                for coefficient, regulatory_event in interaction.regulatory_events.items():

                    reg_event_val = regulatory_event.evaluate(values=reg_state)

                    if reg_event_val:
                        target_val = coefficient

                res[interaction.target.id] = target_val

        if to_dict:
            return res

        return ModelSolution(method='SynchronousBoolean', x=res, objective_value=0.0, status='Feasible')
