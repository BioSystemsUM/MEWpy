"""JMetalpy operators
"""

import copy
import random
from typing import List

from jmetal.core.operator import Mutation, Crossover
from jmetal.core.solution import Solution

from .problem import KOSolution, OUSolution
from ...util.constants import EAConstants


class ShrinkMutation(Mutation[Solution]):
    """ Shrink mutation. A gene is removed from the solution.

    :param probability: (float), The mutation probability.
    :param min_size: (int) the solution minimum size.

    """

    def __init__(self, probability: float = 0.1, min_size: int = EAConstants.MIN_SOLUTION_SIZE):
        super(ShrinkMutation, self).__init__(probability=probability)
        self.min_size = min_size

    def execute(self, solution: Solution) -> Solution:
        """
        Apply the mutation.

        :param solution: The candidate solution to be mutated.
        :returns: A mutated solution.

        """
        if random.random() <= self.probability and solution.number_of_variables > self.min_size:
            var = copy.copy(solution.variables)
            index = random.randint(0, len(var) - 1)
            del var[index]
            solution.variables = var
            solution.number_of_variables = len(var)
        return solution

    def get_name(self):
        return 'Shrink Mutation'


class GrowMutationKO(Mutation[KOSolution]):
    """ Grow mutation. A gene is added to the solution.

    :param probability: (float), The mutation probability.
    :param min_size: (int) the solution minimum size.

    """

    def __init__(self, probability: float = 0.1, max_size: int = EAConstants.MAX_SOLUTION_SIZE):
        super(GrowMutationKO, self).__init__(probability=probability)
        self.max_size = max_size

    def execute(self, solution: Solution) -> Solution:
        """
        Apply the mutation.

        :param solution: The candidate solution to be mutated.
        :returns: A mutated solution.

        """
        if random.random() <= self.probability and solution.number_of_variables < self.max_size:
            mutant = copy.copy(solution.variables)
            newElem = random.randint(solution.lower_bound, solution.upper_bound)
            while newElem in mutant:
                newElem = random.randint(solution.lower_bound, solution.upper_bound)
            mutant.append(newElem)
            solution.variables = mutant
            solution.number_of_variables = len(mutant)
        return solution

    def get_name(self):
        return 'Grow Mutation KO'


class GrowMutationOU(Mutation[OUSolution]):
    """ Grow mutation. A gene is added to the solution.

    :param probability: (float), The mutation probability.
    :param min_size: (int) the solution minimum size.

    """

    def __init__(self, probability: float = 0.1, max_size: int = EAConstants.MAX_SOLUTION_SIZE):
        super(GrowMutationOU, self).__init__(probability=probability)
        self.max_size = max_size

    def execute(self, solution: Solution) -> Solution:
        """
        Apply the mutation.

        :param solution: The candidate solution to be mutated.
        :returns: A mutated solution.

        """
        if random.random() <= self.probability and solution.number_of_variables < self.max_size:
            mutant = copy.copy(solution.variables)
            idx = random.randint(solution.lower_bound[0], solution.upper_bound[0])
            idxs = [a for (a, b) in mutant]
            while idx in idxs:
                idx = random.randint(solution.lower_bound[0], solution.upper_bound[0])
            lv = random.randint(solution.lower_bound[1], solution.upper_bound[1])
            mutant.append((idx, lv))
            solution.variables = mutant
            solution.number_of_variables = len(mutant)
        return solution

    def get_name(self):
        return 'Grow Mutation OU'


class UniformCrossoverKO(Crossover[KOSolution, KOSolution]):
    """Uniform Crossover for KO solutions

    :param probability: (float) The probability of crossover.
    :param max_size: (int) The solution maximum size.

    """

    def __init__(self, probability: float = 0.1, max_size: int = EAConstants.MAX_SOLUTION_SIZE):
        super(UniformCrossoverKO, self).__init__(probability=probability)
        self.max_size = max_size

    def execute(self, parents: List[KOSolution]) -> List[KOSolution]:
        if len(parents) != 2:
            raise Exception('The number of parents is not two: {}'.format(len(parents)))

        offspring = [copy.deepcopy(parents[0]), copy.deepcopy(parents[1])]

        if random.random() <= self.probability and (
                offspring[0].number_of_variables > 1 or offspring[1].number_of_variables > 1):
            mom = set(copy.copy(offspring[0].variables))
            dad = set(copy.copy(offspring[1].variables))
            intersection = mom & dad
            otherElems = list((mom | dad).difference(intersection))
            child1 = copy.copy(intersection)
            child2 = copy.copy(intersection)

            while len(otherElems) > 0:
                elemPosition = random.randint(0, len(otherElems) - 1) if len(otherElems) > 1 else 0
                if len(child1) == self.max_size or len(child2) == 0:
                    child2.add(otherElems[elemPosition])
                elif len(child2) == self.max_size or len(child1) == 0:
                    child1.add(otherElems[elemPosition])
                else:
                    r = random.random()
                    if r <= 0.5:
                        child1.add(otherElems[elemPosition])
                    else:
                        child2.add(otherElems[elemPosition])
                otherElems.pop(elemPosition)

            offspring[0].variables = list(child1)
            offspring[0].number_of_variables = len(child1)
            offspring[1].variables = list(child2)
            offspring[1].number_of_variables = len(child2)
        return offspring

    def get_number_of_parents(self) -> int:
        return 2

    def get_number_of_children(self) -> int:
        return 2

    def get_name(self):
        return 'Uniform Crossover KO'


class MutationContainer(Mutation[Solution]):
    """A container for the mutation operators.

    :param probability: (float) The probability of applying a mutation.
    :param mutators: (list) The list of mutators.

    """

    def __init__(self, probability: float = 0.5, mutators=[]):
        super(MutationContainer, self).__init__(probability=probability)
        self.mutators = mutators

    def execute(self, solution: Solution) -> Solution:
        # randomly select a mutator and apply it
        if random.random() <= self.probability:
            idx = random.randint(0, len(self.mutators) - 1)
            mutator = self.mutators[idx]
            return mutator.execute(solution)
        else:
            return solution

    def get_name(self):
        return 'Mutation container'


class UniformCrossoverOU(Crossover[OUSolution, OUSolution]):
    """
        Uniform Crossover for OU solutions
    """

    def __init__(self, probability: float = 0.1, max_size: int = EAConstants.MAX_SOLUTION_SIZE):
        super(UniformCrossoverOU, self).__init__(probability=probability)
        self.max_size = max_size

    def execute(self, parents: List[KOSolution]) -> List[KOSolution]:
        if len(parents) != 2:
            raise Exception('The number of parents is not two: {}'.format(len(parents)))

        offspring = [copy.deepcopy(parents[0]), copy.deepcopy(parents[1])]

        if random.random() <= self.probability and (
                offspring[0].number_of_variables > 1 or offspring[1].number_of_variables > 1):
            mom = set(copy.copy(offspring[0].variables))
            dad = set(copy.copy(offspring[1].variables))
            intersection = mom & dad
            otherElems = list((mom | dad).difference(intersection))
            child1 = copy.copy(intersection)
            child2 = copy.copy(intersection)

            while len(otherElems) > 0:
                elemPosition = random.randint(0, len(otherElems) - 1) if len(otherElems) > 1 else 0
                if len(child1) == self.max_size or len(child2) == 0:
                    child2.add(otherElems[elemPosition])
                elif len(child2) == self.max_size or len(child1) == 0:
                    child1.add(otherElems[elemPosition])
                else:
                    r = random.random()
                    if r <= 0.5:
                        child1.add(otherElems[elemPosition])
                    else:
                        child2.add(otherElems[elemPosition])
                otherElems.pop(elemPosition)

            offspring[0].variables = list(child1)
            offspring[0].number_of_variables = len(child1)
            offspring[1].variables = list(child2)
            offspring[1].number_of_variables = len(child2)
        return offspring

    def get_number_of_parents(self) -> int:
        return 2

    def get_number_of_children(self) -> int:
        return 2

    def get_name(self):
        return 'Uniform Crossover OU'


class SingleMutationKO(Mutation[KOSolution]):
    """
    Mutates a single element
    """

    def __init__(self, probability: float = 0.1):
        super(SingleMutationKO, self).__init__(probability=probability)

    def execute(self, solution: Solution) -> Solution:
        if random.random() <= self.probability:
            mutant = copy.copy(solution.variables)
            index = random.randint(0, len(mutant) - 1)
            newElem = random.randint(solution.lower_bound, solution.upper_bound)
            while newElem in mutant:
                newElem = random.randint(solution.lower_bound, solution.upper_bound)
            mutant[index] = newElem
            solution.variables = mutant
        return solution

    def get_name(self):
        return 'Single Mutation KO'


class SingleMutationOU(Mutation[OUSolution]):
    """
    Mutates a single element
    """

    def __init__(self, probability: float = 0.1):
        super(SingleMutationOU, self).__init__(probability=probability)

    def execute(self, solution: Solution) -> Solution:
        if random.random() <= self.probability:
            mutant = copy.copy(solution.variables)
            lix = [i for (i, j) in mutant]
            index = random.randint(0, len(mutant) - 1)
            idx, idy = mutant[index]
            is_mutate_idx = False
            if random.random() > 0.5:
                idx = random.randint(solution.lower_bound[0], solution.upper_bound[0])
                while idx in lix:
                    idx = random.randint(solution.lower_bound[0], solution.upper_bound[0])
                is_mutate_idx = True
            lv = random.randint(solution.lower_bound[1], solution.upper_bound[1])
            if not is_mutate_idx:
                while lv == idy:
                    lv = random.randint(solution.lower_bound[1], solution.upper_bound[1])
            mutant[index] = (idx, lv)
            solution.variables = mutant
        return solution

    def get_name(self):
        return 'Single Mutation KO'


class SingleMutationOULevel(Mutation[OUSolution]):
    """
    Mutates the expression level of a single element
    """

    def __init__(self, probability: float = 0.1):
        super(SingleMutationOULevel, self).__init__(probability=probability)

    def execute(self, solution: Solution) -> Solution:
        if random.random() <= self.probability:
            mutant = copy.copy(solution.variables)
            index = random.randint(0, len(mutant) - 1)
            idx, idy = mutant[index]
            lv = random.randint(solution.lower_bound[1], solution.upper_bound[1])
            while lv == idy:
                lv = random.randint(solution.lower_bound[1], solution.upper_bound[1])
            mutant[index] = (idx, lv)
            solution.variables = mutant
        return solution

    def get_name(self):
        return 'Single Mutation KO'


class GaussianMutation(Mutation[Solution]):
    """
     A Gaussian mutator
    """

    def __init__(self, probability: float = 0.1,
                 gaussian_mutation_rate: float = 0.1,
                 gaussian_mean: float = 0.0,
                 gaussian_std: float = 1.0):
        super(GaussianMutation, self).__init__(probability=probability)
        self.gaussian_mutation_rate = gaussian_mutation_rate
        self.gaussian_mean = gaussian_mean
        self.gaussian_std = gaussian_std

    def execute(self, solution: Solution) -> Solution:
        if random.random() <= self.probability:
            mutant = copy.copy(solution.variables)
            for i, m in enumerate(mutant):
                if random.random() < self.gaussian_mutation_rate:
                    v = m + random.gauss(self.gaussian_mean, self.gaussian_std)
                    while v < solution.lower_bound[i] or v > solution.upper_bound[i]:
                        v = m + random.gauss(self.gaussian_mean, self.gaussian_std)
                    solution.variables[i] = v
        return solution

    def get_name(self):
        return 'Gaussian Mutator'

class UniformCrossover(Crossover[Solution, Solution]):
    """
        Uniform Crossover
    """

    def __init__(self, probability: float = 0.1):
        super(UniformCrossover, self).__init__(probability=probability)

    def execute(self, parents: List[Solution]) -> List[Solution]:
        if len(parents) != 2:
            raise Exception('The number of parents is not two: {}'.format(len(parents)))

        offspring = [copy.deepcopy(parents[0]), copy.deepcopy(parents[1])]

        if random.random() <= self.probability and (
                offspring[0].number_of_variables > 1 or offspring[1].number_of_variables > 1):
            mom = copy.copy(offspring[0].variables)
            dad = copy.copy(offspring[1].variables)
            child1 = []
            child2 = []
            for p in range(len(mom)):
                if random.random() <= 0.5:
                    child1.append(mom[p])
                    child2.append(dad[p])
                else:
                    child1.append(dad[p])
                    child2.append(mom[p])
                
            offspring[0].variables = list(child1)
            offspring[0].number_of_variables = len(child1)
            offspring[1].variables = list(child2)
            offspring[1].number_of_variables = len(child2)
        return offspring

    def get_number_of_parents(self) -> int:
        return 2

    def get_number_of_children(self) -> int:
        return 2

    def get_name(self):
        return 'Uniform Crossover'



class SingleRealMutation(Mutation[Solution]):
    """
    Mutates a single element
    """

    def __init__(self, probability: float = 0.1):
        super(SingleRealMutation, self).__init__(probability=probability)

    def execute(self, solution: Solution) -> Solution:
        if random.random() <= self.probability:
            index = random.randint(0, solution.number_of_variables - 1)
            solution.variables[index] = solution.lower_bound[index] + \
                (solution.upper_bound[index] - solution.lower_bound[index]) * random.random()
        return solution

    def get_name(self):
        return 'Single Real Mutation'


REP_INT = {
    "SHRINK": ShrinkMutation,
    "GROWKO": GrowMutationKO,
    "GROWOU": GrowMutationOU,
    "UCROSSKO": UniformCrossoverKO,
    "UCROSSOU": UniformCrossoverOU,
    "SMUTKO": SingleMutationKO,
    "SMUTOU": SingleMutationOU,
    "SMLEVEL": SingleMutationOULevel
}


def build_ko_operators(problem):
    crossover = UniformCrossoverKO(0.8, problem.candidate_max_size)
    mutators = []

    # add shrink and growth mutation only if max size != min size
    if problem.candidate_max_size != problem.candidate_min_size:
        mutators.append(GrowMutationKO(
            1.0, max_size=problem.candidate_max_size))
        mutators.append(ShrinkMutation(
            1.0, min_size=problem.candidate_min_size))

    mutators.append(SingleMutationKO(1.0))
    mutations = MutationContainer(0.3, mutators=mutators)
    return crossover, mutations


def build_ou_operators(problem):
    crossover = UniformCrossoverOU(0.5, problem.candidate_max_size)
    mutators = []
    _max = problem.candidate_max_size
    _min = problem.candidate_min_size
    _t_size = problem.bounder.upper_bound[0]+1

    # add shrink and growth mutation only if max size != min size
    # and do not add if single ou if  max size == min size == targets size
    if _max != _min:
        mutators.append(GrowMutationOU(
            1.0, max_size=problem.candidate_max_size))
        mutators.append(ShrinkMutation(
            1.0, min_size=problem.candidate_min_size))
        mutators.append(SingleMutationOU(1.0))
    elif _min != _t_size:
        mutators.append(SingleMutationOU(1.0))
    else:
        pass

    mutators.append(SingleMutationOULevel(1.0))
    mutations = MutationContainer(0.3, mutators=mutators)
    return crossover, mutations
