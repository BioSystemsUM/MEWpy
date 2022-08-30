__author__ = 'BiSBII CEB University of Minho'
__email__ = 'vpereira@ceb.uminho.pt'
__version__ = '0.1.18'


def info():
    
    print('MEWpy version:',__version__)
    print('Author:',__author__)
    print('Contact:',__email__,'\n')

    from .simulation import solvers, get_default_solver
    print('Available LP solvers:',' '.join(solvers))
    print('Default LP solver:',get_default_solver(),'\n')
    
    from .solvers import ode_solvers, get_default_ode_solver
    print('Available ODE solvers:',' '.join(list(ode_solvers.keys())))
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
    