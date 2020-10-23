

class ModelConstants:
    # Default reactions upper bound (used instead of Inf)
    REACTION_UPPER_BOUND = 10000
    # reset solver
    RESET_SOLVER = False


class EAConstants:
    # Default fitness value for invalid solution on minimization problems
    MAX_THRESHOLD = 100000000.0
    # Default fitness value for invalid solution on maximization problems
    MIN_THRESHOLD = 0.0
    # Default solution minimum size (minimum number of reactions/proteins/genes to KO/OU)
    MIN_SOLUTION_SIZE = 1
    # Default solution maximum size (maximum number of reactions/proteins/genes to KO/OU)
    MAX_SOLUTION_SIZE = 10
    # Default OU levels (Note: 0 levels are interpreted as KO)
    LEVELS = [0, 1/32, 1/16, 1/8, 1/4, 1/2, 2, 4, 8, 16, 32]
    # Default OU levels for OptRAM(Note: instead of 0 use 0.001 as a logaritmic function is used to process expression values)
    OPTRAM_LEVELS = [0.001, 1/32, 1/16, 1/8, 1/4, 1/2, 2, 4, 8, 16, 32]
    # Number of cpus for multiprocessor evaluation of candidates.
    # If NUM_CPUS is non positive, half of the available cpus are used
    NUM_CPUS = -1
    # Maximum number of generations (used as stopping criteria for the EA)
    MAX_GENERATIONS = 100
    # Use probabilistic modofication target list when applicable
    PROB_TARGET = True
    # DEBUG
    DEBUG = False
