class ModelConstants:
    # Default reactions upper bound (used instead of Inf)
    REACTION_UPPER_BOUND = 10000
    # Default reactions lower bound (used instead of -Inf)
    REACTION_LOWER_BOUND = -10000
    # Default tolerance for bounds and coefficients
    TOLERANCE = 1E-10
    # reset solver
    RESET_SOLVER = False
    # Multiprocessing engine. If ray is installed, it is used by default.
    MP_EVALUATOR = 'ray'
    # Solver timeout
    SOLVER_TIMEOUT = 3600


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
    LEVELS = [0, 1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 2, 4, 8, 16, 32]
    # Default OU levels for OptRAM(Note: instead of 0 use 0.001 as a logaritmic function is used
    # to process expression values)
    OPTRAM_LEVELS = [0.001, 1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 2, 4, 8, 16, 32]
    # Number of cpus for multiprocessor evaluation of candidates.
    # If NUM_CPUS is non positive, half of the available cpus are used
    NUM_CPUS = -1
    # Maximum number of generations (used as stopping criteria for the EA)
    MAX_GENERATIONS = 100
    # Use probabilistic modofication target list when applicable
    PROB_TARGET = True
    # DEBUG
    DEBUG = False
    # Dump on interrupt
    KILL_DUMP = False


# calculated from https://github.com/HegemanLab/atomicWeightsDecimal (refer to this repository for credits)
atomic_weights = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
                  'O': 15.9994, 'F': 18.9984032, 'Ne': 2.01797, 'Na': 22.98977, 'Mg': 24.305, 'Al': 26.981538,
                  'Si': 28.0855, 'P': 30.973761, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
                  'Sc': 44.95591, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938049, 'Fe': 55.845,
                  'Co': 58.9332, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64, 'As': 74.9216,
                  'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585, 'Zr': 91.224,
                  'Nb': 92.90638, 'Mo': 95.94, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42, 'Ag': 107.8682,
                  'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.6, 'I': 126.90447, 'Xe': 131.293,
                  'Cs': 132.90545, 'Ba': 137.327, 'La': 138.9055, 'Ce': 140.116, 'Pr': 140.90765, 'Nd': 144.24,
                  'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92534, 'Dy': 162.5, 'Ho': 164.93032,
                  'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04, 'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479,
                  'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.078, 'Au': 196.96655,
                  'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.98038, 'Th': 232.0381, 'Pa': 231.03588,
                  'U': 238.02891, 'Tc': 98.0, 'Pm': 145.0, 'Po': 209.0, 'At': 210.0, 'Rn': 222.0, 'Fr': 223.0,
                  'Ra': 226.0, 'Ac': 227.0, 'Np': 237.0, 'Pu': 244.0, 'Am': 243.0, 'Cm': 247.0, 'Bk': 247.0,
                  'Cf': 251.0, 'Es': 252.0, 'Fm': 257.0, 'Md': 258.0, 'No': 259.0, 'Lr': 262.0, 'Rf': 261.0,
                  'Db': 262.0, 'Sg': 266.0, 'Bh': 264.0, 'Hs': 277.0, 'Mt': 268.0, 'Ds': 281.0, 'Rg': 272.0,
                  'Cn': 285.0, 'Uuq': 289.0, 'Uuh': 292.0}
