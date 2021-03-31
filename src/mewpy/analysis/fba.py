from typing import Union

from mewpy.model import Model, MetabolicModel, RegulatoryModel
from mewpy.solvers.solution import Solution, Status
from mewpy.solvers.solver import VarType, Solver
from mewpy.lp import Notification, ConstraintContainer, VariableContainer, MetabolicLinearizer, GPRLinearizer
from mewpy.solution import ModelSolution


# TODO: missing documentation and typing
class FBA(MetabolicLinearizer):

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = True,
                 attach: bool = False):

        """
        Flux Balance Analysis (FBA) of a metabolic model. Regular implementation of a FBA for a metabolic model.

        For more details consult: https://dx.doi.org/10.1038%2Fnbt.1614

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

        super().__init__(model=model,
                         solver=solver,
                         build=build,
                         attach=attach)

    def build(self):

        if self._variables or self._constraints:
            self.clean()

        if self.model.is_metabolic():
            self.metabolite_reaction_lookup(reactions=self.model.yield_reactions(),
                                            metabolites=self.model.yield_metabolites())

            self.constraints_from_gprs(inplace=True)

            self.update()

            self.set_objective(linear={var.id: value for var, value in self.model.objective.items()},
                               minimize=False)

    def notification(self, notification: Notification):

        if notification.content_type == 'coefficients' and notification.action == 'set':

            # if the notification is a coefficient object, the bounds must be processed accordingly
            # by each method implementation. In the case of FBA, the reactions bounds are inferred
            # from the coefficients and set regularly.
            # Nevertheless, if the coefficients of a gene are being set, the gprs associated with this gene
            # must be evaluated and then the resulting constraints must be set to the reactions

            variable = notification.content.id
            lb, ub = notification.content.bounds
            active_coefficient = notification.content.active_coefficient

            if variable in self._cols:

                # the variable is a reaction, and bounds must be set normally

                return super(FBA, self).set_bounds(variable=variable, lb=lb, ub=ub)

            else:

                # the variable is a gene, so new constraints must be inferred
                # from the gene coefficients and gprs evaluation

                return self.constraints_from_gene(variable=variable, lb=active_coefficient, ub=ub, inplace=True)

        else:
            return super(FBA, self).notification(notification)

    def optimize(self,
                 objective=None,
                 minimize=False,
                 constraints=None,
                 to_solver=False,
                 get_values=True,
                 shadow_prices=False,
                 reduced_costs=False,
                 pool_size=0,
                 pool_gap=None) -> Union[ModelSolution, Solution]:

        self._set_objective(linear=objective, minimize=minimize)

        if constraints:
            constraints = constraints.copy()

        solver_solution = self._solver.solve(linear=objective,
                                             quadratic=None,
                                             minimize=minimize,
                                             model=None,
                                             constraints=constraints,
                                             get_values=get_values,
                                             shadow_prices=shadow_prices,
                                             reduced_costs=reduced_costs,
                                             pool_size=pool_size,
                                             pool_gap=pool_gap)

        if to_solver:
            return solver_solution

        if minimize:
            sense = 'minimize'

        else:
            sense = 'maximize'

        return ModelSolution.from_solver(method=self.method, solution=solver_solution, model=self.model,
                                         objective_direction=sense)

    def constraints_from_gprs(self, inplace=False):

        rxns, lbs, ubs = [], [], []

        gene_state = {gene.id: gene.coefficient.active_coefficient for gene in self.model.yield_genes()}

        if self.model.is_metabolic():

            # gprs over reactions
            for rxn in self.model.yield_reactions():

                if rxn.gpr.is_none:
                    continue

                else:

                    res = rxn.gpr(values=gene_state)

                    if not res:

                        variable = rxn.id
                        lb = 0.0
                        ub = 0.0

                        if inplace:
                            super(FBA, self).set_bounds(variable=variable, lb=lb, ub=ub)

                        rxns.append(rxn.id)
                        lbs.append(0.0)
                        ubs.append(0.0)

        return rxns, lbs, ubs

    def constraints_from_gene(self, variable, lb, ub, inplace=False):

        rxns, rxn_lbs, rxn_ubs = [], [], []

        gene_state = {gene.id: gene.coefficient.active_coefficient for gene in self.model.yield_genes()}
        gene_state[variable] = lb

        # Either the gene is not available in the model or it is not a gene. The later can occur if a multi type model
        # has an fba simulator attached but a the bounds of a regulator are changed
        model_var = self.model.get(variable, None)

        if model_var is None:
            return

        if not model_var.is_gene():
            return

        for rxn in model_var.yield_reactions():

            if rxn.gpr.is_none:

                continue

            else:

                res = rxn.gpr(values=gene_state)

                if not res:
                    rxn_lb = 0.0
                    rxn_ub = 0.0

                else:
                    rxn_lb = rxn.lower_bound
                    rxn_ub = rxn.upper_bound

                rxn_id = rxn.id
                rxns.append(rxn_id)
                rxn_lbs.append(rxn_lb)
                rxn_ubs.append(rxn_ub)

                if inplace:
                    super(FBA, self).set_bounds(variable=rxn_id, lb=rxn_lb, ub=rxn_ub)

        return rxns, rxn_lbs, rxn_ubs


class milpFBA(MetabolicLinearizer, GPRLinearizer):

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = True,
                 attach: bool = False):

        """
        Mixed-Integer Flux Balance Analysis (FBA) of a metabolic model.
        Regular implementation of a FBA for a metabolic model.
        Additionally, GPRs are linearized into mixed-integer constraints, so that genes are variables of the linear
        problem to be solved.

        Check the mewpy.lp.linearizers module for more detail regarding GPRs linearization

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

        super().__init__(model=model,
                         solver=solver,
                         build=build,
                         attach=attach)

    def build(self):

        if self._variables or self._constraints:
            self.clean()

        if self.model.is_metabolic():
            self.metabolite_reaction_lookup(reactions=self.model.yield_reactions(),
                                            metabolites=self.model.yield_metabolites())

            self.add_gprs(self.model.yield_reactions())

            self.update()

            self.set_objective(linear={var.id: value for var, value in self.model.objective.items()},
                               minimize=False)

    def notification(self, notification: Notification):

        if notification.content_type == 'coefficients' and notification.action == 'set':

            # if the notification is a coefficient object, the whole gpr must be build again,
            # as lower and upper bounds of variables are also encoded into the constraints

            self.update_milp_constraints(notification.content.id)

            lb, ub = notification.content.bounds

            return super(milpFBA, self).set_bounds(notification.content.id, lb=lb, ub=ub)

        else:
            return super(milpFBA, self).notification(notification)

    def update_milp_constraints(self, reaction):

        variable = self.model.get(reaction)

        if variable.is_reaction():
            self.add_gprs([variable])

    def _optimize_with_constraints(self,
                                   objective=None,
                                   minimize=False,
                                   constraints=None,
                                   get_values=True,
                                   shadow_prices=False,
                                   reduced_costs=False,
                                   pool_size=0,
                                   pool_gap=None):

        constraints = constraints.copy()

        old_simulators = self.model.simulators.copy()

        self.model._simulators = [self]

        with self.model:

            for variable, bounds in constraints.items():
                variable_obj = self.model.get(variable)

                # sending notification. Notification processing will update variables and constraints accordingly
                if isinstance(bounds, (tuple, list, set)):
                    variable_obj.coefficient.coefficients = bounds

                else:
                    variable_obj.coefficient.coefficients = (bounds,)

            # updating the last modifications
            self.update()

            # solving for the last modifications
            solver_solution = self._solver.solve(linear=objective,
                                                 quadratic=None,
                                                 minimize=minimize,
                                                 model=None,
                                                 constraints=None,
                                                 get_values=get_values,
                                                 shadow_prices=shadow_prices,
                                                 reduced_costs=reduced_costs,
                                                 pool_size=pool_size,
                                                 pool_gap=pool_gap)

        # context exit. restore. update
        self.update()

        self.model._simulators = old_simulators

        return solver_solution

    def optimize(self,
                 objective=None,
                 minimize=False,
                 constraints=None,
                 to_solver=False,
                 get_values=True,
                 shadow_prices=False,
                 reduced_costs=False,
                 pool_size=0,
                 pool_gap=None):

        self._set_objective(linear=objective, minimize=minimize)

        if constraints:
            solver_solution = self._optimize_with_constraints(objective=objective,
                                                              minimize=minimize,
                                                              constraints=constraints,
                                                              get_values=get_values,
                                                              shadow_prices=shadow_prices,
                                                              reduced_costs=reduced_costs,
                                                              pool_size=pool_size,
                                                              pool_gap=pool_gap)

        else:
            solver_solution = self._solver.solve(linear=objective,
                                                 quadratic=None,
                                                 minimize=minimize,
                                                 model=None,
                                                 constraints=None,
                                                 get_values=get_values,
                                                 shadow_prices=shadow_prices,
                                                 reduced_costs=reduced_costs,
                                                 pool_size=pool_size,
                                                 pool_gap=pool_gap)

        if to_solver:
            return solver_solution

        if minimize:
            sense = 'minimize'

        else:
            sense = 'maximize'

        return ModelSolution.from_solver(method=self.method, solution=solver_solution, model=self.model,
                                         objective_direction=sense)


class pFBA(FBA):

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = True,
                 attach: bool = False):

        """
        Parsimonious Flux Balance Analysis (FBA) of a metabolic model.
        Regular implementation of a pFBA for a metabolic model.

        This pFBA implementation was heavily inspired by pFBA implementation of reframed python package. Take a look at
        the source: https://github.com/cdanielmachado/reframed and https://reframed.readthedocs.io/en/latest/

        For more details consult: https://doi.org/10.1038/msb.2010.47

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

        super().__init__(model=model,
                         solver=solver,
                         build=build,
                         attach=attach)

    def build(self):

        super(pFBA, self).build()

    def optimize(self,
                 objective=None,
                 minimize=False,
                 fraction=None,
                 constraints=None,
                 reactions=None,
                 to_solver=False,
                 get_values=True,
                 shadow_prices=False,
                 reduced_costs=False,
                 pool_size=0,
                 pool_gap=None) -> Union[ModelSolution, Solution]:

        self._set_objective(linear=objective, minimize=minimize)

        if constraints:
            constraints = constraints.copy()

        sol = self._solver.solve(linear=objective,
                                 quadratic=None,
                                 minimize=minimize,
                                 model=None,
                                 constraints=constraints,
                                 get_values=False,
                                 shadow_prices=False,
                                 reduced_costs=False,
                                 pool_size=0,
                                 pool_gap=None)

        if sol.status != Status.OPTIMAL:
            return sol

        constraint = ConstraintContainer(name='pfba_constraints', coefs=[], lbs=[], ubs=[])
        variable = VariableContainer(name='pfba_variables', sub_variables=[], lbs=[], ubs=[], variables_type=[])

        constraint.coefs.append(self._linear_objective)
        if fraction is None:

            constraint.lbs.append(sol.fobj)
            constraint.ubs.append(sol.fobj)

        else:

            constraint.lbs.append(sol.fobj * fraction)
            constraint.ubs.append(sol.fobj)

        if not reactions:
            reactions = self.model.reactions.keys()

        pfba_objective = {}

        for rxn in reactions:

            rxn_obj = self.model.get(rxn)

            if rxn_obj.reversibility:

                rxn_forward = f'{rxn}_forward'
                rxn_reverse = f'{rxn}_reverse'

                variable.sub_variables.extend([rxn_forward, rxn_reverse])
                variable.lbs.extend([0, 0])
                variable.ubs.extend([rxn_obj.upper_bound, rxn_obj.upper_bound])
                variable.variables_type.extend([VarType.CONTINUOUS, VarType.CONTINUOUS])

                constraint.coefs.extend([{rxn: -1, rxn_forward: 1},
                                         {rxn: 1, rxn_reverse: 1}])
                constraint.lbs.extend([0, 0])
                constraint.ubs.extend([rxn_obj.upper_bound, rxn_obj.upper_bound])

                pfba_objective[rxn_forward] = 1
                pfba_objective[rxn_reverse] = 1

            else:

                pfba_objective[rxn_obj.id] = 1

        self.add((variable, constraint))
        self.update()

        old_linear_objective = self._linear_objective.copy()
        old_quadratic_objective = self._quadratic_objective.copy()
        old_sense = self._minimize

        self._set_objective(linear=pfba_objective, minimize=True)

        solver_solution = self._solver.solve(linear=pfba_objective,
                                             quadratic=None,
                                             minimize=True,
                                             model=None,
                                             constraints=constraints,
                                             get_values=get_values,
                                             shadow_prices=shadow_prices,
                                             reduced_costs=reduced_costs,
                                             pool_size=pool_size,
                                             pool_gap=pool_gap)

        self.remove((variable, constraint))
        self.update()
        self.set_objective(linear=old_linear_objective, quadratic=old_quadratic_objective, minimize=old_sense)

        if to_solver:
            return solver_solution

        if minimize:
            sense = 'minimize'

        else:
            sense = 'maximize'

        return ModelSolution.from_solver(method=self.method, solution=solver_solution, model=self.model,
                                         objective_direction=sense)
