import ast
import os
from collections import OrderedDict
import numpy as np
from .constants import ModelConstants
from ..optimization.ea import Solution, non_dominated_population
from ..optimization.ea import filter_duplicates
from ..simulation import get_simulator, SimulationMethod

# xml escapes
escapes = {'&': '&#38', '<': '&lt', '>': '&gt', "'": '&#39', '"': '&#34'}


def str_constraints(constraints, separator=";"):
    """ Converts a dictionary o constraints into a string representation.
        Constraints may of the form  'key:value' or 'key:(value,value)'
    """
    tokens = []
    for k, v in constraints.items():
        tokens.append(k)
        if isinstance(v, (float, int)):
            tokens.extend([str(v), ""])
        elif isinstance(v, tuple):
            tokens.extend(map(str, list(v)))
        else:
            raise ValueError("Unrecognized value type")
    out_str = separator.join(x for x in tokens)
    return out_str


def xml_constraints(constraints, separator=""):
    """ Converts a dictionary o constraints into a xml string representation.
    Constraints may of the form  'key:value' or 'key:(value,value)'


    :param constrainst: dic, a dictionary of constraints
    :param separator: str, optional
    :return: a xml string representation
    """
    tokens = []
    tokens.append("<constraints>")
    for k, v in constraints.items():
        tokens.append("<constraint><reaction>{}</reaction>".format(k))
        if isinstance(v, (float, int)):
            lb = v
            ub = v
        elif isinstance(v, tuple):
            lb = v[0]
            ub = v[1]
        else:
            raise ValueError("Unrecognized value type")
        tokens.append("<lb>{}</lb>".format(lb))
        tokens.append("<ub>{}</ub>".format(ub))
        tokens.append("</constraint>")
    tokens.append("</constraints>")
    out_str = separator.join(x for x in tokens)
    return out_str


def population_to_csv(problem, candidates, filename, simplify=False, non_dominated=False, sep=';'):
    """ Saves a population of solution candidates to csv
    """
    if non_dominated:
        population = non_dominated_population(
            candidates, problem.is_maximization)
    else:
        population = candidates

    # should probably simplify recursevly until no changes
    if simplify:
        pop = []
        for s in population:
            pop.extend(problem.simplify(s))
        population = pop

    population.sort(reverse=True)
    with open(filename, 'w') as f:
        title = [str(x).replace(sep, " ") for x in problem.fevaluation]
        title.extend([sep, "solution", "\n"])
        f.write(sep.join(x for x in title))
        for individual in population:
            line = []
            if problem.number_of_objectives == 1:
                line.append(str(individual.fitness))
            else:
                line.extend(map(str, individual.fitness))
            line.append("")
            line.append(str(individual.values))
            line.append("")
            line.append(str(individual.constraints))
            line.append("\n")
            f.write(sep.join(x for x in line))


def population_to_xml(problem, candidates, filename, simplify=False, non_dominated=True):
    """ Saves a population of solution candidates to xml
    """

    if non_dominated:
        population = non_dominated_population(
            candidates, problem.is_maximization)
    else:
        population = candidates

    # should recursevly simplify until the initial
    # and simplified populations are the same
    if simplify:
        pop = []
        for s in population:
            pop.extend(problem.simplify(s))
        population = pop

    population.sort(reverse=True)
    with open(filename, 'w') as f:
        f.write("<population>")
        f.write("<ea><objectives>")
        for obj in problem.fevaluation:
            f.write("<objective>{}</objective>".format(str(obj)))
        f.write("</objectives>")
        f.write("</ea><candidates>")
        for individual in population:
            line = []
            line.append("<candidate>")
            if problem.number_of_objectives == 1:
                line.append("<fitness>{}</fitness>".format(individual.fitness))
            else:
                for value in individual.fitness:
                    line.append("<fitness>{}</fitness>".format(value))
            constraints = individual.get_constraints()
            line.append(xml_constraints(constraints))
            line.append("</candidate>")
            f.write("".join(x for x in line))
        f.write("</candidates></population>")


def is_number(s):
    """Easy way to verify if a string is a number
    """
    return s.replace('.', '').replace('-', '').replace('e', '').isdigit()


class Parser:
    def __init__(self, obj_labels):
        """
        Parser to retreive solutions saved as csv.

        :param obj_labels: (list) optimization objective labels.
        """
        self.obj_labels = obj_labels
        self.n_obj = len(obj_labels)
        self.obj_labels.append('Size')
        self.population = []
        self.targets_stats = None

    def __parse_line(self, line, separator, parse_constraints=True):
        line = line.replace(';;', ';')
        tokens = line.split(separator)
        tokens.remove('\n')
        # fitness
        fitness = [abs(float(x)) for x in tokens[:self.n_obj]]
        # encoding
        values = ast.literal_eval(tokens[self.n_obj])

        if parse_constraints:
            # temporary compatibility with previous version
            if tokens[self.n_obj + 1].startswith('{'):
                constrainst = ast.literal_eval(tokens[self.n_obj + 1])
            elif tokens[self.n_obj + 1].startswith('OrderedDict'):
                constrainst = eval(tokens[self.n_obj + 1])
            else:
                constrainst = {}
                i = self.n_obj + 1
                while i < len(tokens) - 2:
                    if is_number(tokens[i + 2]):
                        constrainst[tokens[i]] = (float(tokens[i + 1]), float(tokens[i + 2]))
                        i = i + 3
                    else:
                        constrainst[tokens[i]] = (float(tokens[i + 1]), float(tokens[i + 1]))
                        i = i + 2

        else:
            constrainst = None
        if np.prod(fitness) != 0:
            solution = Solution(values, fitness, constrainst)
            self.population.append(solution)

    def parse_results(self, path, separator=';', non_dominated=True, maximize=True, filter_duplicate=True,
                      parse_constraints=True):
        """
        Parse all csv file in the path.

        Parameters

        :param path: (str) the path to the directory containing the csv files;
        :param separator: (str) the separator used in the csv files;
        :param non_dominated: (boolean) filter dominated solutions;
        :param maximize: (boolean) if the solutions are from a maximization problem;
        :param filter_duplicate: (boolean) remove duplicates solutions.

        """
        for r, _, fl in os.walk(path):
            for file in fl:
                if '.csv' in file:
                    # print('Parsing :',file)
                    with open(os.path.join(r, file), 'r') as f:
                        # skip the first line
                        f.readline()
                        line = f.readline()
                        while line:
                            self.__parse_line(line, separator, parse_constraints=parse_constraints)
                            line = f.readline()
        if non_dominated:
            self.population = non_dominated_population(
                self.population, maximize=maximize, filter_duplicate=filter_duplicate)
        elif filter_duplicate:
            self.population = filter_duplicates(self.population)

    def compute_fluxes(self, model, biomass, product, carbon_source, envcond=None,
                       population=None, minimal_growth=0.0, minimal_product=0.0):
        """
        Computes flux analises for each solution contained in the population.

        :param model: The model.
        :param biomass: (str) Biomass reaction id.
        :param product: (str) Product reaction id.
        :param carbon_source: (str) Carbon source reaction id.
        :param envcond: (dict) Dictionary of environmental conditions. Default None.
        :param population: (list) List of solutions. If none is given the parsed population is considered.
        :param minimal_growth: (float) Solution with growth below "minimal_growth" are discarded. Default 0.0.
        :param minimal_product: (float) Solution with product yield below "minimal_product" are discarded. Default 0.0.

        """
        ModelConstants.RESET_SOLVER = True
        simul = get_simulator(model, envcond=envcond)
        simul.set_objective(biomass)
        reactions = [product, carbon_source]
        if population is None:
            population = self.population
        if not population or len(population) == 0:
            raise RuntimeError("Population is empty")
        res = simul.simulate(
            objective={biomass: 1}, method=SimulationMethod.pFBA)
        self.wt_biomass = res.fluxes[biomass]
        self.wt_products = {}
        # wt_ref = res.fluxes

        for rx in reactions:
            self.wt_products[rx] = res.fluxes[rx]

        filtered_population = []
        data_fitness = []
        # data_constraints = []
        data_fba = []
        data_lmoma = []
        data_fva = []

        n = len(population[0].fitness)
        if n == self.obj_labels:
            labels = self.obj_labels
        else:
            labels = [f"Obj_{i}" for i in range(1, n + 1)]

        for solution in population:
            r = []
            u = []
            m = []

            constraints = {}
            constraints.update(solution.constraints)

            try:
                res = simul.simulate(
                    objective={biomass: 1}, method=SimulationMethod.pFBA, constraints=constraints)
                if res.fluxes:
                    # pFBA
                    biomass_value = res.fluxes[biomass]
                    if biomass_value < minimal_growth or res.fluxes[product] < minimal_product:
                        continue
                    r.append(biomass_value)
                    for rx in reactions:
                        r.append(res.fluxes[rx])
                    # lMOMA
                    res = simul.simulate(objective={
                        biomass: 1}, method=SimulationMethod.lMOMA, constraints=constraints)
                    u.append(res.fluxes[biomass])
                    for rx in reactions:
                        u.append(res.fluxes[rx])
                    # FVA
                    constraints[biomass] = (biomass_value * 0.99, 100000.0)
                    res = simul.simulate(objective={product: 1}, constraints=constraints, maximize=False)
                    if res.fluxes:
                        m.append(res.fluxes[product])
                    else:
                        m.append(0)
                    res = simul.simulate(objective={product: 1}, constraints=constraints, maximize=True)
                    if res.fluxes:
                        m.append(res.fluxes[product])
                    else:
                        m.append(0)

                filtered_population.append(solution)
                data_fitness.append(solution.get_fitness())
                # data_constraints.append(str(solution.constraints))
                data_fba.append(r)
                data_lmoma.append(u)
                data_fva.append(m)
            except Exception:
                continue

        data_values = [str(individual.values) for individual in filtered_population]
        import pandas as pd
        df_values = pd.DataFrame(data_values, columns=["Solution"])
        df_fitness = pd.DataFrame(data_fitness, columns=labels)
        df_fba = pd.DataFrame(data_fba, columns=['Biomass_pFBA', product + '_pFBA', carbon_source + '_pFBA'])
        df_lmoma = pd.DataFrame(data_lmoma, columns=['Biomass_lMOMA', product + '_lMOMA', carbon_source + '_lMOMA'])
        df_fva = pd.DataFrame(data_fva, columns=[product + '_FVAmin_99', product + '_FVAmax_99'])
        # df_constraints = pd.DataFrame(data_constraints, columns=['Constraints'])

        df_all = df_values.join(df_fitness).join(df_fba).join(df_lmoma).join(df_fva)  # .join(df_constraints)

        return df_all, filtered_population

    def reaction_stats(self):
        """
        Reactions statistics
        """
        stats = OrderedDict()
        for solution in self.population:
            rxns = solution.constraints.keys()
            for rx in rxns:
                if rx in stats.keys():
                    stats[rx] += 1
                else:
                    stats[rx] = 1
        self.stats = stats

    def target_stats(self, population=None):
        """
        Target statistics
        """
        stats = OrderedDict()
        if not population:
            population = self.population
        for solution in population:
            targets = list(solution.values.keys())
            for t in targets:
                if t in stats.keys():
                    stats[t] += 1
                else:
                    stats[t] = 1
        import operator
        sorted_d = OrderedDict(sorted(stats.items(),
                                      key=operator.itemgetter(1), reverse=True))

        self.targets_stats = sorted_d
        return sorted_d

    def find_all_solution(self, target, population=None):
        """
        Finds all solutions containing the list of targets

        :param targets: (list) List of target ids.
        :param population: (list) List of solution where to perform the search. If none is provided, the search is \
            performed on the entire population.
        :returns: A list of solutions.
        """
        result = []
        if not isinstance(target, list):
            target = [target]
        if not population:
            population = self.population
        for solution in population:
            if all([v in solution.values.keys() for v in target]):
                result.append(solution)
        return result
