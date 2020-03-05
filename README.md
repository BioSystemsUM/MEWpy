MEWpy
============

  

MEWpy is an integrated Metabolic Engineering Workbench, under development, for strain design optimization. Developed in Python, MEWpy offers methods to explore different classes of constraint-based models (CBM). 

* Simulation: allows to simulate different metabolic models considering different phenotype.
  
* Optimization: performs Evolutionary Computation based strain design optimization by knocking out or over/under expressing reactions, genes or enzymes.


*  MEWPy currently supports [REFRAMED](<https://github.com/cdanielmachado/reframed>) and [COBRApy](<https://opencobra.github.io/cobrapy/>) based models, notably:

	* CBM;
	* [GECKO](https://doi.org/10.15252/msb.20167411);
	* [OptRAM](https://doi.org/10.1371/journal.pcbi.1006835).

* Next releases will enable simulation and optimization strategies for [Expression and Thermodynamics Flux](https://doi.org/10.1371/journal.pcbi.1006835) (ETFL) and [Metabolism and Macromolecular Expression](https://doi.org/10.1371/journal.pcbi.1006302) (ME) models, as well as Boolean Regulatory Networks. 

  
  
* The optimization engine relies on either [inspyred](<https://github.com/aarongarrett/inspyred>) or [jMetalPy](<https://github.com/jMetal/jMetalPy>) packages, which are used for creating biologically-inspired computational intelligence algorithms in Python.

  
  
  
  

### Usage Examples

  
  
  

### Instalation

  

1. clone the repository

2. run ``python setup.py install``

  
  
  

MEWPy requires a compatible solver for linear programming problems, with installed Python dependencies installed, from the following list:

  

-  [CPLEX](<https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/>) (preferred)
-  [GUROBI](<http://www.gurobi.com>)
-  [GLPK](<https://www.gnu.org/software/glpk/>)

  
  

### Credits and License

Developed at:

* Centre of Biological Engineering, University of Minho (2019-2020)

  
Released under an Apache License.



