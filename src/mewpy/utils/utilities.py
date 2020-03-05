from inspyred.ec.emo import Pareto
import copy

#xml escapes
escapes = {'&' : '&#38', '<' :'&lt', '>' :'&gt', "'" :'&#39', '"' : '&#34'} 

def str_constraints(constraints, separator = ";"):
    """ Converts a dictionary o constraints into a string representation.
        Constraints may of the form  'key:value' or 'key:(value,value)'
    """
    tokens = []
    for k,v in constraints.items():
        tokens.append(k)
        if isinstance(v,(float,int)):
            tokens.extend([str(v),""])
        elif isinstance(v,tuple):
            tokens.extend(map(str,list(v)))
        else:
            raise ValueError("Unrecognized value type")
    out_str = separator.join(x for x in tokens)
    return out_str



def xml_constraints(constraints, separator = ""):
    """ Converts a dictionary o constraints into a xml string representation.
        Constraints may of the form  'key:value' or 'key:(value,value)'

        parameters:
            constrainst: dic
                A dictionary of constraints
            separator: str, optional
    """
    tokens =[]
    tokens.append("<constraints>")
    for k,v in constraints.items():
        tokens.append("<constraint><reaction>{}</reaction>".format(k))
        if isinstance(v,(float,int)):
            lb = v
            ub = v
        elif isinstance(v,tuple):
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




def population_to_csv(problem, candidates,filename, simplify = False, non_dominated = True):
    """ Saves a population of solution candidates to csv
    """
    if non_dominated:
        population = non_dominated_population(candidates, problem.is_maximization)
    else:
        population = candidates
    
    population.sort(reverse = True)
    with open(filename,'w') as f:
        title = [str(x).replace(";"," ") for x in problem.fevaluation]
        title.extend([";","solution","\n"])
        f.write(";".join(x for x in title))
        for individual in population:
            line = []
            if problem.number_of_objectives == 1:
                line.append(str(individual.fitness))
            else:
                line.extend(map(str,individual.fitness))
            line.append("")    
            if simplify:
                constraints = problem.simplify(individual)
            else:     
                constraints = problem.get_constraints(individual)
            line.append(str_constraints(constraints))
            line.append("\n")
            f.write(";".join(x for x in line))




def population_to_xml(problem,candidates,filename, simplify = False , non_dominated = True):
    """ Saves a population of solution candidates to xml
    """

    if non_dominated:
        population = non_dominated_population(candidates, problem.is_maximization)
    else:
        population = candidates

    population.sort(reverse = True)
    with open(filename,'w') as f:
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
                constraints = problem.get_constraints(individual)
            line.append(xml_constraints(constraints))
            line.append("</candidate>")
            f.write("".join(x for x in line))
        f.write("</candidates></population>")





def dominance_test(solution1, solution2, maximize = True):
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

    if isinstance(solution1.fitness,Pareto):
        values1 = solution1.fitness.values 
        values2 = solution2.fitness.values
    else :
        values1 = [solution1.fitness] 
        values2 = [solution2.fitness]

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




def non_dominated_population(population, maximize = True, filter_duplicate = False ):
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
            if dominance_test(individual,population[j],maximize= maximize) == -1:
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
            other =  population[i]
            if individual.fitness == other.fitness and set(individual.candidate)==set(other.candidate):
                to_remove.append(i)
        for i in to_remove:
            del population[i]
        return population    

    fitered_list = []
    l = population
    while len(l)>1:
        individual = l[0]
        fitered_list.append(individual)
        remove_equal(individual,l[1:])
    return fitered_list    
        


