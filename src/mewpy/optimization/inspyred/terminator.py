
def generation_termination(population, num_generations, num_evaluations, args):
    """Return True if the number of generations meets or exceeds a maximum.

    This function compares the number of generations with a specified 
    maximum. It returns True if the maximum is met or exceeded.

    .. Arguments:
       population -- the population of Individuals
       num_generations -- the number of elapsed generations
       num_evaluations -- the number of candidate solution evaluations
       args -- a dictionary of keyword arguments

    Optional keyword arguments in args:

    - *max_generations* -- the maximum generations (default 1) 

    """
    max_gen = args.get('max_generations', 1)
    if isinstance(max_gen, tuple):
        max_generations = max_gen[0]
    else:
        max_generations = max_gen
    return num_generations >= max_generations
