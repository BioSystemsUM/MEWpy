__author__ = 'BiSBII CEB University of Minho'
__email__ = 'vpereira@ceb.uminho.pt'
__version__ = '0.1.21'


def info():
    
    print('MEWpy version:',__version__)
    print('Author:',__author__)
    print('Contact:',__email__,'\n')

    from .simulation import __MEWPY_sim_solvers__, get_default_solver
    print('Available LP solvers:',' '.join(__MEWPY_sim_solvers__))
    print('Default LP solver:',get_default_solver(),'\n')
    
    from .solvers import __MEWPY_ode_solvers__, get_default_ode_solver
    print('Available ODE solvers:',' '.join(list(__MEWPY_ode_solvers__.keys())))
    print('Default ODE solver:', get_default_ode_solver(),'\n')
    
    from mewpy.util.utilities import get_all_subclasses
    from mewpy.problems.problem import AbstractProblem
    c = get_all_subclasses(AbstractProblem)
    print('Optimization Problems:',' '.join(sorted([x.__name__ for x in c])))
    print()
    from mewpy.optimization import engines, get_default_engine, get_available_algorithms
    print('Available EA engines:',' '.join(list(engines.keys())))
    print('Default EA engine:', get_default_engine())
    print('Available EAs:', ' '.join(sorted(get_available_algorithms())),'\n')
    