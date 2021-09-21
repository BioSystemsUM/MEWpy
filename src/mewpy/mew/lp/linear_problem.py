from abc import ABCMeta, abstractmethod
from typing import Union, TYPE_CHECKING, Tuple, Dict, Set, List

from numpy import zeros

from mewpy.solvers import get_default_solver, solvers
from mewpy.solvers.solution import Solution
from mewpy.solvers.solver import Solver, VarType
from mewpy.mew.solution import ModelSolution

from .linear_containers import ConstraintContainer, VariableContainer
from .notification import Notification
from .linear_utils import LinkedList, Node, integer_coefficients

try:
    # noinspection PyPackageRequirements
    from cobamp.core.linear_systems import GenericLinearSystem
    # noinspection PyPackageRequirements
    from cobamp.core.optimization import LinearSystemOptimizer

except ImportError:

    GenericLinearSystem = None
    LinearSystemOptimizer = None

if TYPE_CHECKING:
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


class LinearProblemInterface(metaclass=ABCMeta):
    """

    Interface for a Linear Problem. The following attributes and methods must be implemented to set up alternative
    problems

    """

    @property
    @abstractmethod
    def model(self):
        return

    @property
    @abstractmethod
    def constraints(self):
        return

    @property
    @abstractmethod
    def variables(self):
        return

    @property
    @abstractmethod
    def objective(self):
        return

    @property
    @abstractmethod
    def minimize(self):
        return

    @property
    @abstractmethod
    def solver(self):
        return

    @abstractmethod
    def notification(self, notification: Notification):
        pass

    @abstractmethod
    def optimize(self):
        pass

    @abstractmethod
    def build(self):
        pass

    @abstractmethod
    def clean(self):
        pass

    @abstractmethod
    def add(self, items):
        pass

    @abstractmethod
    def remove(self, items):
        pass

    @abstractmethod
    def update(self):
        pass


# TODO: missing documentation and typing
class LinearProblem(LinearProblemInterface):

    def __init__(self,
                 model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                 solver: Union[str, Solver] = None,
                 build: bool = True,
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

        A linear problem has one and only one solver object that connect with each other via synchronous updates.
        That is, alterations to a linear problem are directly reflected in the solver.
        Alternatively, one can impose temporary constraints during problem optimization
        (see the method for further details)

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
        if not model:
            raise ValueError('A valid model must be provided')

        self._model = model

        # this is usefull to restore the linear problem to the point of init
        self._initial_solver = solver

        self._solver = self.mewpy_solver(solver)

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

        # Observer pattern: a linear problem object is an observer of a mewpy model.
        # When the model is changed, a notification with the change is sent to all observers (linear problems).
        # The linear problem implementation specific for each method processes the notification accordingly.
        # Finally, when the linear problem is updated,
        # all variables and constraints added to the linear problem are implemented and kept in sync with the solver
        # This can avoid consecutive building of the solver, namely a lazy loading

        # a dictionary may be too much to store pending changes, but we only want to hold the latest alterations.
        # Thus, old pending constraints or variables are replaced by the recent ones.
        # This allows the last optimization to be a direct reflection of the model latest state
        # Note that, adding and removing items from a set can be an overhead
        self._variables_queue: Dict[str, VariableContainer] = {}
        self._constraints_queue: Dict[str, ConstraintContainer] = {}

        if build:
            self.build()

        if attach:
            self.model.attach(self)

    # -----------------------------------------------------------------------------
    # Correct solver helper
    # -----------------------------------------------------------------------------
    @staticmethod
    def mewpy_solver(solver):

        if solver is None:
            solver_name = get_default_solver()

            SolverType = solvers[solver_name]

            solver = SolverType()

        elif isinstance(solver, str):

            SolverType = solvers.get(solver, None)

            if SolverType is None:
                raise ValueError(f'{solver} is not listed as valid solver. Check the valid solvers: {solvers}')

            solver = SolverType()

        elif isinstance(solver, Solver):

            pass

        else:
            raise ValueError(f'Invalid solver {solver}. Check the valid solvers: {solvers}')

        return solver

    # -----------------------------------------------------------------------------
    # Built-in
    # -----------------------------------------------------------------------------

    def __str__(self):
        return f"{self.method} for {self.model.id}"

    def __repr__(self):
        return f"{self.method}: {self.model.id}"

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

        mewpy model of this simulator

        :return: a mewpy Model, MetabolicModel, RegulatoryModel or all
        """

        return self._model

    @property
    def solver(self) -> Union[Solver]:
        """
        mewpy solver instance for this linear problem. It contains an interface for the concrete solver

        :return: A Solver, CplexSolver, GurobiSolver or OptLangSolver instance
        """

        return self._solver

    @property
    def constraints(self) -> Dict[str, ConstraintContainer]:
        """
        A copy of the constraints container.
        This container holds all ConstraintContainer objects for this linear problem.
        Note that, a constraint container can hold several constraints/rows

        :return: copy of the constraints dictionary
        """

        return self._constraints.copy()

    @property
    def variables(self) -> Dict[str, VariableContainer]:
        """
        A copy of the variables container.
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
        return self._get_matrix()

    @property
    def bounds(self):
        return self.get_bounds(as_list=True)

    @property
    def b_bounds(self):
        return self.get_bounds(b_bounds=True, as_list=True)

    @property
    def shape(self):
        return int(len(self._rows)), int(len(self._cols))

    # -----------------------------------------------------------------------------
    # Observer pattern - Model notification system
    # -----------------------------------------------------------------------------

    def notification(self, notification: Notification):

        if notification.content_type == 'coefficients' and notification.action == 'set':

            variable = notification.content.id
            lb, ub = notification.content.bounds

            return self.set_bounds(variable=variable, lb=lb, ub=ub)

        elif notification.content_type == 'objectives' and notification.action == 'set':

            return self.set_objective(**notification.content)

        return

    # -----------------------------------------------------------------------------
    # Objective
    # -----------------------------------------------------------------------------

    def _set_objective(self,
                       linear: Union[str, Dict[str, Union[float, int]]] = None,
                       quadratic: Dict[Tuple[str, str], Union[float, int]] = None,
                       minimize: bool = True):
        """

        INTERNAL USE ONLY. See set_objective method

        :param linear:
        :param quadratic:
        :param minimize:
        :return:
        """

        # This helper method verify inputs and sets the objective in the linear problem object

        linear_objective = {}

        if linear:

            if isinstance(linear, str):

                linear = {linear: 1.0}

            else:

                for k, _ in linear.items():

                    if k not in self._cols and k not in self._sub_cols:
                        raise ValueError(f'{k} is not a variable of this linear problem')

            linear_objective = linear

        if linear_objective:
            self._linear_objective = linear_objective

        quadratic_objective = {}

        if quadratic:

            for (k1, k2), _ in quadratic.items():

                if k1 not in self._cols and k1 not in self._sub_cols:
                    raise ValueError(f'{k1} is not a variable of this linear problem')

                if k2 not in self._cols and k2 not in self._sub_cols:
                    raise ValueError(f'{k2} is not a variable of this linear problem')

            quadratic_objective = quadratic

        if quadratic_objective:
            self._quadratic_objective = quadratic_objective

        if minimize is not None and minimize != self._minimize:
            self._minimize = minimize

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

        # This helper method verify inputs and sets the objective in the linear problem object
        self._set_objective(linear=linear, quadratic=quadratic, minimize=minimize)

        # If the inputs are correct, the objective is set in the solver
        self._solver.set_objective(linear=linear, quadratic=quadratic, minimize=minimize)

    def set_bounds(self, variable: str, lb: Union[float, int], ub: Union[float, int]):
        """
        Set the bounds of a given linear variable available in the linear problem.
        Note that, only the bounds of atom variables can be set

        :param variable: the variable identifier
        :param lb: lower bound value
        :param ub: upper bound value

        :return:
        """

        # Either the variable is not available in the model or it is not a variable of this lp.
        # The later can occur if a multi type model has for instance an fba simulator attached
        # but a the bounds of a regulator are changed
        if variable not in self._cols and variable not in self._sub_cols:
            return

        old_lb, old_ub = self._get_variable_bounds(variable=variable)

        if len(old_lb) > 1:
            raise ValueError(f'{variable} bounds cannot be set')

        else:
            old_lb = old_lb[0]
            old_ub = old_ub[0]

        lb = lb if lb is not None else old_lb
        ub = ub if ub is not None else old_ub

        if (lb, ub) in integer_coefficients:
            var_type = VarType.INTEGER

        else:
            var_type = VarType.CONTINUOUS

        var = self._variables[variable]
        var.lbs = [lb]
        var.ubs = [ub]
        var.variables_type = [var_type]

        self.add([var])

    # -----------------------------------------------------------------------------
    # Optimization
    # -----------------------------------------------------------------------------

    def optimize(self, *args, **kwargs) -> Union[ModelSolution, Solution]:
        """

        Abstract implementation

        :return: ModelSolution or Solution (from solver) objects
        """

        # The concrete implementation is defined by each simulation method, e.g. fba, pfba, etc

    # -----------------------------------------------------------------------------
    # Operations/Manipulations - build and clean
    # -----------------------------------------------------------------------------

    @abstractmethod
    def build(self):
        """

        Abstract implementation

        :return:
        """

        # The concrete implementation is defined by each simulation method, e.g. fba, pfba, etc

    def build_solver(self):
        """

        It creates an new solver instance and adds the current state (variables, constraints) of the linear problem

        :return:
        """

        self._solver = self.mewpy_solver(self._initial_solver)

        for variable in self._variables.values():

            # Using mewpy/reframed solver interface ...
            for name, (lb, ub, var_type) in variable.items():
                self._solver.add_variable(var_id=name, lb=lb, ub=ub, vartype=var_type, update=False)

        self.solver.update()

        for i, constraint in enumerate(self._constraints.values()):

            # Using mewpy/reframed solver interface ...
            for j, (coef, lb, ub) in constraint.items():

                cnt_id = str(i + j)

                if lb == ub:
                    rhs = lb
                    self._solver.add_constraint(constr_id=cnt_id, lhs=coef, sense='=', rhs=rhs, update=False)

                else:
                    cnt_id_f = f'{cnt_id}_forward'
                    rhs = lb
                    self._solver.add_constraint(constr_id=cnt_id_f, lhs=coef, sense='>', rhs=rhs, update=False)

                    cnt_id_r = f'{cnt_id}_reverse'
                    rhs = ub
                    self._solver.add_constraint(constr_id=cnt_id_r, lhs=coef, sense='<', rhs=rhs, update=False)

        self.solver.update()

        self.set_objective(linear=self._linear_objective,
                           quadratic=self._quadratic_objective,
                           minimize=self._minimize)

    def clean(self):

        self._solver = self.mewpy_solver(self._initial_solver)
        self._cols = LinkedList()
        self._rows = LinkedList()
        self._sub_cols = LinkedList()
        self._constraints = {}
        self._variables = {}
        self._linear_objective = {}
        self._quadratic_objective = {}
        self._minimize = True
        self._variables_queue: Dict[str, VariableContainer] = {}
        self._constraints_queue: Dict[str, ConstraintContainer] = {}

    # -----------------------------------------------------------------------------
    # Operations/Manipulations - add, remove
    # -----------------------------------------------------------------------------

    def stack_container(self, container, addition):

        if isinstance(container, ConstraintContainer):
            self._constraints_queue[container.name] = (container, addition)

        if isinstance(container, VariableContainer):
            self._variables_queue[container.name] = (container, addition)

    def add(self, containers: Union[List[Union[VariableContainer, ConstraintContainer]],
                                    Set[Union[VariableContainer, ConstraintContainer]],
                                    Tuple[Union[VariableContainer, ConstraintContainer]]]):

        for container in containers:
            self.stack_container(container, addition=True)

    def remove(self, containers: Union[List[Union[VariableContainer, ConstraintContainer]],
                                       Set[Union[VariableContainer, ConstraintContainer]],
                                       Tuple[Union[VariableContainer, ConstraintContainer]]]):

        for container in containers:
            self.stack_container(container, addition=False)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations - Update
    # -----------------------------------------------------------------------------

    def update(self):

        build_problem = False

        if self._constraints_queue or self._variables_queue:
            build_problem = True

        build_new_objective = False

        for container, addition in self._variables_queue.values():

            if addition:
                self._add_variable(container)

            else:

                build_new_objective = True

                self._remove_variable(container)

        for container, addition in self._constraints_queue.values():

            if addition:
                self._add_constraint(container)

            else:
                self._remove_constraint(container)

        self._variables_queue: Dict[str, VariableContainer] = {}
        self._constraints_queue: Dict[str, ConstraintContainer] = {}

        if build_new_objective:

            new_linear = {}
            new_quadratic = {}

            for k, val in self._linear_objective.items():

                if k in self._cols or k in self._sub_cols:
                    new_linear[k] = val

            for (k1, k2), val in self._quadratic_objective.items():

                if (k1 in self._cols or k1 in self._sub_cols) and (k2 in self._cols or k2 in self._sub_cols):
                    new_quadratic[(k1, k2)] = val

            self.set_objective(linear=new_linear, quadratic=new_quadratic, minimize=self.minimize)

        # updating the solver with new constraints is horrible!!! It sometimes yields infeasible solutions
        if build_problem:
            self.build_solver()

    # -----------------------------------------------------------------------------
    # Operations/Manipulations - LinearProblem getters
    # -----------------------------------------------------------------------------
    def get_linear_system(self) -> Tuple[GenericLinearSystem, LinearSystemOptimizer]:

        if GenericLinearSystem is None or LinearSystemOptimizer is None:
            raise RuntimeError('CoBAMP is not installed or was not found')

        lb, ub = self.get_bounds(as_tuples=True)
        b_lb, b_ub = self.get_bounds(b_bounds=True, as_tuples=True)

        var_names = list(self._sub_cols.keys())
        var_types = [types.value.lower() for var in self._variables.values() for types in var.variables_type]

        linear_system = GenericLinearSystem(S=self.matrix,
                                            lb=lb,
                                            ub=ub,
                                            b_lb=b_lb,
                                            b_ub=b_ub,
                                            var_names=var_names,
                                            var_types=var_types)

        linear_system.build_problem()
        optimizer = LinearSystemOptimizer(linear_system, build=False)

        if hasattr(self.model, 'objective'):

            objective = zeros(self.shape[1], )

            for k, v in self.model.objective.items():
                objective[self.index(variable=k)] = v

            linear_system.set_objective(objective, False)

        return linear_system, optimizer

    def index(self, variable=None, constraint=None, as_list=False, as_int=False, default=None):

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
                   as_tuples=False
                   ):

        if variable is not None:

            return self._get_variable_bounds(variable=variable)

        elif constraint is not None:

            return self._get_constraint_bounds(constraint)

        elif b_bounds:

            return self._get_b_bounds(as_list=as_list, as_tuples=as_tuples)

        else:
            return self._get_bounds(as_list=as_list, as_tuples=as_tuples)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations - LinearProblem add/remove variables and constraints
    # -----------------------------------------------------------------------------

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

    def _add_constraint(self, constraint: ConstraintContainer):

        if constraint.name in self._rows:
            # The constraint is replaced, as the linear problem behaves like a set
            # This also mimics the solver interface behavior

            old_constraint = self._constraints[constraint.name]

            self._remove_constraint(old_constraint)

        node = constraint.to_node()

        self._rows.add(node)

        self._update_constraint_coefs(constraint)

        self._constraints[constraint.name] = constraint

    def _remove_constraint(self, constraint: ConstraintContainer):

        if constraint.name in self._rows:
            self._rows.pop(constraint.name)
            self._constraints.pop(constraint.name)

    def _add_variable(self, variable: VariableContainer):

        if variable.name in self._cols:
            # The variable is replaced, as the linear problem behaves like a set
            # This also mimics the solver interface behavior

            old_variable = self._variables[variable.name]

            self._remove_variable(old_variable)

        node = variable.to_node()

        self._cols.add(node)

        self._variables[variable.name] = variable

        for sub_variable in variable.keys():
            sub_node = Node(value=sub_variable, length=1)

            self._sub_cols.add(sub_node)

    def _remove_variable(self, variable: VariableContainer):

        if variable.name in self._cols:

            self._cols.pop(variable.name)
            self._variables.pop(variable.name)

            for sub_variable in variable.keys():
                self._sub_cols.pop(sub_variable)
