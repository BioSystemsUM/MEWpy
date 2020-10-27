## Options



MEWpy makes available a large set of options, some being globally defined in *mewpy.util.constants*.



**Number of processors for parallel solutions evaluation.**

By default, MEWpy uses half of the available treads to run parallel evaluations.  However, a user may define the number of parallel threads by altering the *NUM_CPUS* constant in mewpy.util.constants:

```python
from mewpy.utils.constants import EAConstants
# uses 32 parallel threads
EAConstants.NUM_CPUS = 32
```



**Over and under expression folds.**

Over- and under-expression optimization problems have a set of possible folds globally defined. It is possible, however, to define such folds directly in a problem definition, for example, in a genes' under expression problem:

```python
# the model has already been loaded

# Define the regulation folds for under-expression or deletion.
# 0 for deletion.
levels = [1/8,1/4,1/2,0]

from mewpy.problems import GOUProblem
problem = GOUProblem(model,levels=levels)
```



**Number of modifications.**

The minimum and the maximum number of modifications may be defined directly in the problem definition. For example to allow a maximum of  gene deletions:

```python
from mewpy.problems import GKOProblem
problem = GKOProblem(model,candidate_max_size=6)
```

Likewise the minimum number of modifications may be explecitly defined:

```python
from mewpy.problems import GKOProblem
problem = GKOProblem(model,candidate_min_size=4,candidate_max_size=6)
```

The default minimum number of modifications is 1, while the maximum is 30. When both the minimum and the maximum number of modifications are equal, all solutions will have the same common number of modifications.

**Optimization algorithm.**

MEWpy resorts to the Inspyred and JMetalPy packages to evolve modification solutions. If both packages are installed, MEWpy uses Inpyred by default, running the Non-dominated Sorting Genetic Algorithm (NSGA-II) for multi-objective optimizations and a Genetic Algorithm (GA) for single objective problems. To alter the engine preference to JMetalPy the following must be added to your script:

```python
from mewpy.optimization import set_default_engine
set_default_engine('jmetal')
```

Also, MEWpy allows to define the optimization algorithm when configuring the EA for the addressed optimization problem:

```python
ea = EA(problem, max_generations=ITERATIONS, algorithm='NSGAIII')
```

Note that when using a single objective, only Simulated Annealing (SA) and Genetic Algorithm (GA) are allowed. Any other configuration will be ignored. The same is true when choosing single objective algorithms to solve multi objective problems. To run multi objective problems using a single objective algorithms please refer to the AggregatedSum evaluation function.