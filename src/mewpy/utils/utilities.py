from mewpy.optimization.ea import non_dominated_population
import copy

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

        parameters:
            constrainst: dic
                A dictionary of constraints
            separator: str, optional
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


def population_to_csv(problem, candidates, filename, simplify=False, non_dominated=True):
    """ Saves a population of solution candidates to csv
    """
    if non_dominated:
        population = non_dominated_population(
            candidates, problem.is_maximization)
    else:
        population = candidates

    population.sort(reverse=True)
    with open(filename, 'w') as f:
        title = [str(x).replace(";", " ") for x in problem.fevaluation]
        title.extend([";", "solution", "\n"])
        f.write(";".join(x for x in title))
        for individual in population:
            line = []
            if problem.number_of_objectives == 1:
                line.append(str(individual.fitness))
            else:
                line.extend(map(str, individual.fitness))
            line.append("")
            if simplify:
                constraints = problem.simplify(individual)
            else:
                constraints = individual.get_constraints()
            line.append(str_constraints(constraints))
            line.append("\n")
            f.write(";".join(x for x in line))


def population_to_xml(problem, candidates, filename, simplify=False, non_dominated=True):
    """ Saves a population of solution candidates to xml
    """

    if non_dominated:
        population = non_dominated_population(
            candidates, problem.is_maximization)
    else:
        population = candidates

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

            if simplify:
                constraints = problem.simplify(individual)
            else:
                constraints = individual.get_constraints()
            line.append(xml_constraints(constraints))
            line.append("</candidate>")
            f.write("".join(x for x in line))
        f.write("</candidates></population>")
