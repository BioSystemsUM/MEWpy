from abc import abstractmethod
from typing import Union, TYPE_CHECKING, Tuple, Dict, Any

from numpy import zeros

from mewpy.mew.solution import ModelSolution
from mewpy.solvers.solution import Solution
from mewpy.solvers.solver import Solver
from .linear_containers import ConstraintContainer, VariableContainer
from .linear_utils import LinkedList, Node, get_solver_instance

if TYPE_CHECKING:
    from mewpy.mew.models import Model, MetabolicModel, RegulatoryModel


class LinearProblem:

    def __init__(self,
                 model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                 solver: Union[str, Solver] = None,
                 build: bool = False,
                 attach: bool = False):
        """
        Linear programing base implementation. A mewpy model is converted into a linear problem using reframed/mewpy
        solver interface. Both CPLEX and Gurobi solvers are currently supported. Other solvers may also be supported
        using an additional OptLang solver interface. However, CPLEX and Gurobi are recommended for certain problems.

        A linear problem is linked with a given model via asynchronous updates.
        That is, alterations to the model are sent to all attached simulators via notification objects.
        Notifications are processed accordingly by all linear problems attached to the model.
        Each implementation of a linear problem (e.g. FBA, RFBA, SRFBA, etc) is responsible
        for processing the notifications in the correct way.

        A linear problem has one and only one solver object.
        Alterations to a linear problem are promptly forced in the solver by building a new solver instance.
        Alternatively, one can impose temporary constraints during problem optimization
        (see the method for further details)

        Notes for developers:
        A linear problem object is an observer (observer pattern) of a mew model.
        A notification with model changes is sent to all observers (linear problems).
        The linear problem implementation specific for each method processes the notification accordingly.
        Finally, when the linear problem is updated,
        all variables and constraints added to the linear problem are implemented and kept in sync with the solver
        This can avoid consecutive building of the solver, namely a lazy loading

        :param model: a mewpy Model, MetabolicModel, RegulatoryModel or all. The model is used to retrieve
        variables and constraints to the linear problem
        :param solver: A Solver, CplexSolver, GurobiSolver or OptLangSolver instance.
        Alternatively, the name of the solver is also accepted.
        The solver interface will be used to load and solve a linear problem in a given solver.
        If none, a new solver is instantiated. An instantiated solver may be used,
        but it will be overwritten if build is true.
        :param build: Whether to build the linear problem upon instantiation. Default: False
        :param attach: Whether to attach the linear problem to the model upon instantiation. Default: False
        """
        if not model:
            raise ValueError('A valid model must be provided')

        self._model = model

        # this is useful to restore the linear problem to the point of init
        self._initial_solver = solver
        self._solver = None
        self._synchronized = False

        # Simulator index engine uses a linked list with a built-in dictionary. This allows fast access to the index
        # of a given variable or constraint.
        # Note that, some variables or constraints can comprise multiple rows or columns,
        # so that the job of keeping in track of all indexes of a given variable/constraint is actually way
        # harder than it seems for simple linear problems (e.g. fba)
        self._cols = LinkedList()
        self._rows = LinkedList()

        # one to one indexing of all variables
        self._sub_cols = LinkedList()

        # Holding constraints and variables objects
        self._constraints = {}
        self._variables = {}

        # the objective is a dict variable_id: coefficient
        self._linear_objective = {}
        self._quadratic_objective = {}
        self._minimize = True

        if build:
            self.build()

        if attach:
            self.model.attach(self)

    # -----------------------------------------------------------------------------
    # Built-in
    # -----------------------------------------------------------------------------
    def __str__(self):
        return f"{self.method} for {self.model.id}"

    def __repr__(self):
        return self.__str__()

    def _repr_html_(self):
        """
        It returns a html representation of the linear problem
        :return:
        """
        if self.solver:
            solver = self.solver.__class__.__name__
        else:
            solver = 'None'

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
                <th>Variables</th>
                <td>{len(self.variables)}</td>
            </tr>
            <tr>
                <th>Constraints</th>
                <td>{len(self.constraints)}</td>
            </tr>
            <tr>
                <th>Objective</th>
                <td>{self.objective}</td>
            </tr>
            <tr>
                <th>Solver</th>
                <td>{solver}</td>
            </tr>
            <tr>
                <th>Synchronized</th>
                <td>{self.synchronized}</td>
            </tr>
        </table>
        """

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------
    @property
    def method(self) -> str:
        """
        Name of the method implementation to build and solve the linear problem
        :return: the name of the class
        """
        return self.__class__.__name__

    @property
    def model(self) -> Union['Model', 'MetabolicModel', 'RegulatoryModel']:
        """
        Mew model of this simulator
        :return: a mewpy Model, MetabolicModel, RegulatoryModel or all
        """
        return self._model

    @property
    def solver(self) -> Solver:
        """
        mewpy solver instance for this linear problem. It contains an interface for the concrete solver
        :return: A Solver, CplexSolver, GurobiSolver or OptLangSolver instance
        """
        return self._solver

    @property
    def synchronized(self) -> bool:
        """
        Whether the linear problem is synchronized with the model
        :return:
        """
        return self._synchronized

    @property
    def constraints(self) -> Dict[str, ConstraintContainer]:
        """
        A copy of the constraints' container.
        This container holds all ConstraintContainer objects for this linear problem.
        Note that, a constraint container can hold several constraints/rows
        :return: copy of the constraints dictionary
        """
        return self._constraints.copy()

    @property
    def variables(self) -> Dict[str, VariableContainer]:
        """
        A copy of the variables' container.
        This container holds all VariableContainer objects for this linear problem.
        Note that, a variable container can hold several variables/columns
        :return: copy of the variables dictionary
        """
        return self._variables.copy()

    @property
    def objective(self) -> Dict[Union[str, Tuple[str, str]], Union[float, int]]:
        """
        A copy of the objective dictionary. Keys are either variable identifiers or tuple of variable identifiers.
        Values are the corresponding coefficients
        Note that, linear and quadratic objectives can be encoded in the objective dictionary.
        See the set_objective method for further detail
        :return: copy of the objective dictionary
        """
        return {**self._linear_objective, **self._quadratic_objective}

    @property
    def minimize(self) -> bool:
        """
        The linear problem objective sense/direction
        :return: a boolean whether the linear problem objective sense/direction is minimization
        """
        return bool(self._minimize)

    # -----------------------------------------------------------------------------
    # Dynamic attributes
    # -----------------------------------------------------------------------------
    @property
    def matrix(self):
        """
        The linear problem matrix
        :return: a matrix as numpy array
        """
        return self._get_matrix()

    @property
    def bounds(self):
        """
        The linear problem bounds
        :return: bounds as list of tuples
        """
        return self.get_bounds(as_list=True)

    @property
    def b_bounds(self):
        """
        The linear problem b bounds (constraints bounds)
        :return: b bounds as list of tuples
        """
        return self.get_bounds(b_bounds=True, as_list=True)

    @property
    def shape(self):
        """
        The linear problem shape
        :return: a tuple with the number of rows and columns
        """
        return int(len(self._rows)), int(len(self._cols))

    # -----------------------------------------------------------------------------
    # MEWpy solver
    # -----------------------------------------------------------------------------
    def build_solver(self, variables: bool = True, constraints: bool = True, objective: bool = True):
        """
        It creates a new solver instance and adds the current state (variables, constraints) of the linear problem
        to the solver.
        :param variables: Whether to add variables to the solver. Default: True
        :param constraints: Whether to add constraints to the solver. Default: True
        :param objective: Whether to add the objective to the solver. Default: True
        :return:
        """
        if variables or constraints:
            self._solver = get_solver_instance(self._initial_solver)

        if variables:
            for variable in self._variables.values():

                # Using mewpy/reframed solver interface ...
                for name, (lb, ub, var_type) in variable.items():
                    self.solver.add_variable(var_id=name, lb=lb, ub=ub, vartype=var_type, update=False)

            self.solver.update()

        if constraints:
            for i, constraint in enumerate(self._constraints.values()):

                # Using mewpy/reframed solver interface ...
                for j, (coef, lb, ub) in constraint.items():

                    cnt_id = str(i + j)

                    if lb == ub:
                        rhs = lb
                        self.solver.add_constraint(constr_id=cnt_id, lhs=coef, sense='=', rhs=rhs, update=False)

                    else:
                        cnt_id_f = f'{cnt_id}_forward'
                        rhs = lb
                        self.solver.add_constraint(constr_id=cnt_id_f, lhs=coef, sense='>', rhs=rhs, update=False)

                        cnt_id_r = f'{cnt_id}_reverse'
                        rhs = ub
                        self.solver.add_constraint(constr_id=cnt_id_r, lhs=coef, sense='<', rhs=rhs, update=False)

            self.solver.update()

        if objective:
            linear_objective = {}

            if self._linear_objective:

                for k, v in self._linear_objective.items():

                    if k not in self._cols and k not in self._sub_cols:
                        raise ValueError(f'{k} is not a variable of this linear problem')

                    linear_objective[k] = v

            quadratic_objective = {}

            if self._quadratic_objective:

                for (k1, k2), v in self._quadratic_objective.items():

                    if k1 not in self._cols and k1 not in self._sub_cols:
                        raise ValueError(f'{k1} is not a variable of this linear problem')

                    if k2 not in self._cols and k2 not in self._sub_cols:
                        raise ValueError(f'{k2} is not a variable of this linear problem')

                    quadratic_objective[(k1, k2)] = v

            self.solver.set_objective(linear_objective, quadratic_objective, self._minimize)
            self.solver.update()

    # -----------------------------------------------------------------------------
    # Clean
    # -----------------------------------------------------------------------------
    def clean(self):
        """
        It cleans the linear problem object by removing all variables and constraints
        :return:
        """
        self._synchronized = False
        self._solver = get_solver_instance(self._initial_solver)
        self._cols = LinkedList()
        self._rows = LinkedList()
        self._sub_cols = LinkedList()
        self._constraints = {}
        self._variables = {}
        self._linear_objective = {}
        self._quadratic_objective = {}
        self._minimize = True
        return

    # -----------------------------------------------------------------------------
    # Build
    # -----------------------------------------------------------------------------
    @abstractmethod
    def _build(self):
        """
        Abstract method for the concrete build method
        :return:
        """
        pass

    def build(self) -> 'LinearProblem':
        """
        Abstract implementation
        :return:
        """
        # clean first
        self.clean()

        # concrete build
        self._build()

        # build solver
        self.build_solver(variables=True, constraints=True, objective=True)

        # update status
        self._synchronized = True
        return self

    # -----------------------------------------------------------------------------
    # Optimization
    # -----------------------------------------------------------------------------
    @abstractmethod
    def _optimize(self, solver_kwargs: Dict[str, Any] = None, **kwargs) -> Solution:
        """
        Abstract method for the concrete optimization method
        :param solver_kwargs: solver specific keyword arguments
        :param kwargs: keyword arguments
        :return:
        """
        pass

    def optimize(self,
                 to_solver: bool = False,
                 solver_kwargs: Dict[str, Any] = None,
                 **kwargs) -> Union[ModelSolution, Solution]:
        """
        It solves the linear problem. The linear problem is solved using the solver interface.

        The optimize method allows setting temporary changes to the linear problem. The changes are
        applied to the linear problem reverted to the original state afterward.
        Objective, constraints and solver parameters can be set temporarily.

        The solution is returned as a ModelSolution instance, unless to_solver is True. In this case,
        the solution is returned as a SolverSolution instance.

        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False.
        Otherwise, a ModelSolution is returned.
        :param solver_kwargs: Solver parameters to be set temporarily.
            - linear: A dictionary of linear coefficients to be set temporarily. The keys are the variable names
            and the values are the coefficients. Default: None
            - quadratic: A dictionary of quadratic coefficients to be set temporarily. The keys are tuples of
            variable names and the values are the coefficients. Default: None
            - minimize: Whether to minimize the objective. Default: False
            - constraints: A dictionary with the constraints bounds. The keys are the constraint ids and the values
            are tuples with the lower and upper bounds. Default: None
            - get_values: Whether to retrieve the solution values. Default: True
            - shadow_prices: Whether to retrieve the shadow prices. Default: False
            - reduced_costs: Whether to retrieve the reduced costs. Default: False
            - pool_size: The size of the solution pool. Default: 0
            - pool_gap: The gap between the best solution and the worst solution in the pool. Default: None
        :return: A ModelSolution instance or a SolverSolution instance if to_solver is True.
        """
        # build solver if out of sync
        if not self.synchronized:
            self.build()

        if not solver_kwargs:
            solver_kwargs = {}

        # concrete optimize
        solution = self._optimize(solver_kwargs=solver_kwargs, **kwargs)

        if to_solver:
            return solution

        minimize = solver_kwargs.get('minimize', self._minimize)
        return ModelSolution.from_solver(method=self.method, solution=solution, model=self.model,
                                         minimize=minimize)

    # -----------------------------------------------------------------------------
    # Update - Observer interface
    # -----------------------------------------------------------------------------
    def update(self):
        """
        It updates the linear problem object by adding/removing variables and constraints
        Note that linear problems are not updated after each addition/removal of a variable or constraint to the model.
        This is done to avoid unnecessary updates of the solver. Instead, the update is done when the method
        `build` is called. If required, this method is called by the simulation methods (e.g. fba, pfba, etc) before the
        optimization process in the `optimize` method.
        :return:
        """
        self._synchronized = False

    # -----------------------------------------------------------------------------
    # Objective
    # -----------------------------------------------------------------------------
    def set_objective(self,
                      linear: Union[str, Dict[str, Union[float, int]]] = None,
                      quadratic: Dict[Tuple[str, str], Union[float, int]] = None,
                      minimize: bool = True):
        """
        A dictionary of the objective for the linear problem.
        Keys must be variables of the linear problem,
        whereas values must be the corresponding coefficients as int or float.
        It can be changed during optimization

        :param linear: a dictionary of linear coefficients or variable identifier (that is set with a coefficient of 1)
        :param quadratic: a dictionary of quadratic coefficients.
        Note that keys must be a tuple of reaction pairs to be summed up to a quadratic objective function
        :param minimize: whether to solve a minimization problem. This parameter is True by default
        :return:
        """
        if linear is None:
            linear = {}

        if quadratic is None:
            quadratic = {}

        if isinstance(linear, str):
            linear = {linear: 1}

        if not isinstance(linear, dict):
            raise TypeError(f'linear objective must be a dictionary, not {type(linear)}')

        if not isinstance(quadratic, dict):
            raise TypeError(f'quadratic objective must be a dictionary, not {type(quadratic)}')

        if not isinstance(minimize, bool):
            raise TypeError(f'minimize must be a boolean, not {type(minimize)}')

        self._linear_objective = linear
        self._quadratic_objective = quadratic
        self._minimize = minimize
        self.build_solver(objective=True)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations - add/remove variables and constraints
    # -----------------------------------------------------------------------------
    def add_constraints(self, *constraints: ConstraintContainer):
        for constraint in constraints:
            if constraint.name in self._rows:
                # The constraint is replaced, as the linear problem behaves like a set
                # This also mimics the solver interface behavior

                old_constraint = self._constraints[constraint.name]

                self.remove_constraints(old_constraint)

            node = constraint.to_node()

            self._rows.add(node)

            self._update_constraint_coefs(constraint)

            self._constraints[constraint.name] = constraint

    def remove_constraints(self, *constraints: ConstraintContainer):
        for constraint in constraints:
            if constraint.name in self._rows:
                self._rows.pop(constraint.name)
                self._constraints.pop(constraint.name)

    def add_variables(self, *variables: VariableContainer):
        for variable in variables:
            if variable.name in self._cols:
                # The variable is replaced, as the linear problem behaves like a set
                # This also mimics the solver interface behavior

                old_variable = self._variables[variable.name]

                self.remove_variables(old_variable)

            node = variable.to_node()

            self._cols.add(node)

            self._variables[variable.name] = variable

            for sub_variable in variable.keys():
                sub_node = Node(value=sub_variable, length=1)

                self._sub_cols.add(sub_node)

    def remove_variables(self, *variables: VariableContainer):
        for variable in variables:
            if variable.name in self._cols:

                self._cols.pop(variable.name)
                self._variables.pop(variable.name)

                for sub_variable in variable.keys():
                    self._sub_cols.pop(sub_variable)

    def _update_constraint_coefs(self, constraint: ConstraintContainer):

        # some constraint coefficients might have keys that refer to the variable name and not sub-variable name.
        # Since only sub-variable names are added to the solvers,
        # these keys must be updated to the last sub-variable name. Note that, the last sub-variable name is regularly
        # the one that matters, as the initial sub-variables regularly decide the outcome of the last one.

        new_coefs = []

        for coefficient in constraint:

            new_coef = {}

            for var, coef in coefficient.items():

                if var in self._sub_cols:

                    sub_var = var

                else:

                    variable = self._variables[var]

                    sub_var = variable.sub_variables[-1]

                new_coef[sub_var] = coef

            new_coefs.append(new_coef)

        constraint.coefs = new_coefs

    # -----------------------------------------------------------------------------
    # Getters
    # -----------------------------------------------------------------------------
    def index(self, variable=None, constraint=None, as_list=False, as_int=False, default=None):
        """
        It returns the index of a variable or constraint
        :param variable: a variable container
        :param constraint: a constraint container
        :param as_list: a boolean indicating whether the index should be returned as a list
        :param as_int: a boolean indicating whether the index should be returned as an integer
        :param default: a default value to be returned if the variable or constraint is not found
        :return: the index of the variable or constraint
        """
        if variable is None and constraint is None:
            raise ValueError('Please provide a variable or constraint')

        if constraint is not None:

            slc = self._rows.get(constraint)

        else:

            slc = self._cols.get(variable, self._sub_cols.get(variable))

            if slc is None:
                return default

        if as_list:
            return [i for i in range(slc.start, slc.stop)]

        elif as_int:
            return slc.stop - 1

        else:

            return slc

    def _get_matrix(self):

        matrix = zeros(self.shape)

        n = 0
        for cnt in self._constraints.values():

            for coef in cnt.coefs:

                for var, value in coef.items():

                    m = self.index(variable=var, as_int=True)

                    if m is None:
                        continue

                    matrix[n, m] = value

                n += 1

        return matrix

    def _get_b_bounds(self, as_list=False, as_tuples=False):

        if as_list:

            b_bounds = []

            for cnt in self._constraints.values():
                bds = list(zip(cnt.lbs, cnt.ubs))

                b_bounds.extend(bds)

            return b_bounds

        elif as_tuples:

            lbs = []
            ubs = []

            for cnt in self._constraints.values():
                lbs.extend(cnt.lbs)
                ubs.extend(cnt.ubs)

            return tuple(lbs), tuple(ubs)

        else:
            return {key: (cnt.lbs, cnt.ubs) for key, cnt in self._constraints.items()}

    def _get_bounds(self, as_list=False, as_tuples=False):

        if as_list:

            bounds = []

            for var in self._variables.values():
                bds = list(zip(var.lbs, var.ubs))

                bounds.extend(bds)

            return bounds

        elif as_tuples:

            lbs = []
            ubs = []

            for var in self._variables.values():
                lbs.extend(var.lbs)
                ubs.extend(var.ubs)

            return tuple(lbs), tuple(ubs)

        else:
            return {key: (var.lbs, var.ubs) for key, var in self._variables.items()}

    def _get_variable_bounds(self, variable):

        variable = self._variables.get(variable)

        if variable is None:
            lb, ub = self._get_bounds(as_list=True)[self._sub_cols.get(variable)]

            return [lb], [ub]

        return variable.lbs, variable.ubs

    def _get_constraint_bounds(self, constraint):

        constraint = self._constraints.get(constraint)

        return constraint.lbs, constraint.ubs

    def get_bounds(self,
                   variable=None,
                   constraint=None,
                   b_bounds=False,
                   as_list=False,
                   as_tuples=False):
        """
        It returns the bounds of a variable or constraint
        :param variable: a variable container
        :param constraint: a constraint container
        :param b_bounds: a boolean indicating whether the bounds of the constraints should be returned
        :param as_list: a boolean indicating whether the bounds should be returned as a list
        :param as_tuples: a boolean indicating whether the bounds should be returned as a tuple
        :return: the bounds of the variable or constraint
        """
        if variable is not None:

            return self._get_variable_bounds(variable=variable)

        elif constraint is not None:

            return self._get_constraint_bounds(constraint)

        elif b_bounds:

            return self._get_b_bounds(as_list=as_list, as_tuples=as_tuples)

        else:
            return self._get_bounds(as_list=as_list, as_tuples=as_tuples)
