from abc import ABC, abstractmethod
from mewpy.utils.constants import EAConstants
from collections import OrderedDict


class SolutionInterface(ABC):

    @abstractmethod
    def get_fitness(self):
        """
        returns a list of fitness values
        """
        raise NotImplementedError

    def get_representation(self):
        """
        returns a set representation of the solution
        """
        raise NotImplementedError


class Solution(SolutionInterface):

    def __init__(self, values, fitness, constraints=None):
        """
        EA Solution

        args:
            candidate: a set() representation of the solution 
            fitness:  a list of fitness values
        """
        self.values = values
        self.fitness = fitness
        self.constraints = OrderedDict() if constraints is None else constraints
        self._is_maximize = True

    def get_fitness(self):
        return self.fitness

    def get_representation(self):
        return self.values

    def get_constraints(self):
        return self.constraints

    def __gt__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize=self._is_maximize) == 1
        return False

    def __lt__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize=self._is_maximize) == -1
        return False

    def __ge__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize=self._is_maximize) != -1
        return False

    def __le__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize=self._is_maximize) != 1
        return False

    def __copy__(self):
        import copy
        values = copy.copy(self.values)
        fitness = self.fitness.copy()
        new_solution = Solution(values, fitness)
        return new_solution


class AbstractEA():

    def __init__(self, problem, initial_population=[], max_generations=EAConstants.MAX_GENERATIONS, mp=True, visualizer=False):

        self.problem = problem
        self.initial_population = initial_population
        self.max_generations = max_generations
        self.visualizer = visualizer
        self.mp = mp
        self.final_population = None

    def run(self):
        """ Runs the optimization for the defined problem.
            The number of objectives is defined to be the number of evaluation
            functions in fevalution. If there are more than one objective, 
            NSGAII is used as optimization engine. 
        """

        if self.problem.fevaluation is None or len(self.problem.fevaluation) == 0:
            raise ValueError("At leat one objective should be provided.")

        # builds the target list
        self.problem.pre_process()

        if self.problem.number_of_objectives == 1:
            final_pop = self._run_so()
        else:
            final_pop = self._run_mo()

        self.final_population = self._convertPopulation(final_pop)
        return self.final_population

    @abstractmethod
    def _convertPopulation(self, population):
        raise NotImplementedError

    @abstractmethod
    def _run_so(self):
        raise NotImplementedError

    @abstractmethod
    def _run_mo(self):
        raise NotImplementedError


def dominance_test(solution1, solution2, maximize=True):
    """
    Testes Pareto dominance
    args
        solution1 : The first solution 
        solution2 : The second solution
        maximize (bool): maximization (True) or minimization (False)

    returns 
         1 : if the first solution dominates the second 
        -1 : if the second solution dominates the first
         0 : if non of the solutions dominates the other
    """
    best_is_one = 0
    best_is_two = 0

    values1 = solution1.get_fitness()
    values2 = solution2.get_fitness()

    for i in range(len(values1)):
        value1 = values1[i]
        value2 = values2[i]
        if value1 != value2:
            if value1 < value2:
                best_is_two = 1
            if value1 > value2:
                best_is_one = 1

    if best_is_one > best_is_two:
        if maximize:
            result = 1
        else:
            result = -1
    elif best_is_two > best_is_one:
        if maximize:
            result = -1
        else:
            result = 1
    else:
        result = 0

    return result


def non_dominated_population(population, maximize=True, filter_duplicate=False):
    """
    returns the non dominated solutions from the population.
    """
    #population.sort(reverse = True)
    non_dominated = []
    for i in range(len(population)-1):
        individual = population[i]
        j = 0
        dominates = True
        while j < len(population) and dominates:
            if dominance_test(individual, population[j], maximize=maximize) == -1:
                dominates = False
            else:
                j += 1
        if dominates:
            non_dominated.append(individual)

    if filter_duplicate:
        result = filter_duplicates(non_dominated)
    else:
        result = non_dominated
    return result


def filter_duplicates(population):
    """ Filters equal solutions from a population
    """
    def remove_equal(individual, population):
        to_remove = []
        for i in range(len(population)):
            other = population[i]
            if individual.fitness == other.fitness and set(individual.get_representation()) == set(other.get_representation()):
                to_remove.append(i)
        for i in to_remove:
            del population[i]
        return population

    fitered_list = []
    l = population
    while len(l) > 1:
        individual = l[0]
        fitered_list.append(individual)
        remove_equal(individual, l[1:])
    return fitered_list
