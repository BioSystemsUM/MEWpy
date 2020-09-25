[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI version](https://badge.fury.io/py/mewpy.svg)](https://badge.fury.io/py/mewpy)
[![Documentation Status](https://readthedocs.org/projects/mewpy/badge/?version=latest)](https://mewpy.readthedocs.io/en/latest/?badge=latest)


MEWpy
============

  

MEWpy is an integrated Metabolic Engineering Workbench, under development, for strain design optimization. Developed in Python, MEWpy offers methods to explore different classes of constraint-based models (CBM). 

* Simulation: allows to simulate different metabolic models considering different phenotype.
  
* Optimization: performs Evolutionary Computation based strain design optimization by knocking out (KO) or over/under expressing (OU) reactions, genes or enzymes.

  

  **Metabolic Constraints**

  | Method       | Strategy |
  | ------------ | -------- |
  | Reactions    | KO / OU  |
  | Genes (GPRs) | KO / OU  |

  **Enzymatic Constraints**

  | Method                                                       | Strategy |
  | ------------------------------------------------------------ | -------- |
  | [GECKO](https://doi.org/10.15252/msb.20167411)               | KO / OU  |
  | [sMOMENT](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3329-9) | KO / OU  |

  **Regulatory Constraints**

  | Method                                                 | Strategy |
  | ------------------------------------------------------ | -------- |
  | [OptRAM](https://doi.org/10.1371/journal.pcbi.1006835) | KO / OU  |
  | [OptORF](https://doi.org/10.1186/1752-0509-4-53)       | KO       |



*  MEWPy currently supports [REFRAMED](<https://github.com/cdanielmachado/reframed>) and [COBRApy](<https://opencobra.github.io/cobrapy/>) simulation environments.
*  The optimization engine relies on either [inspyred](<https://github.com/aarongarrett/inspyred>) or [jMetalPy](<https://github.com/jMetal/jMetalPy>) packages, which are used for creating biologically-inspired computational intelligence algorithms in Python.


* Next releases will enable simulation and optimization strategies for [Expression and Thermodynamics Flux](https://doi.org/10.1371/journal.pcbi.1006835) (ETFL) and [Metabolism and Macromolecular Expression](https://doi.org/10.1371/journal.pcbi.1006302) (ME) models. 

  
  
  
  

### Usage Examples

Examples are provided as [jupiter notebooks](examples) and as [python scripts](src/mewpy/unittests).


### Instalation

Installing from Pypi package repository:
  
  ``pip install mewpy``

Installing from github:

1. clone the repository 
  
  ``git clone https://github.com/BioSystemsUM/mewpy.git``

2. run ``python setup.py install``

  

  

MEWPy requires a compatible solver for linear programming problems, with installed Python dependencies installed, from the following list:


-  [CPLEX](<https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/>) (preferred)
-  [GUROBI](<http://www.gurobi.com>)
-  [GLPK](<https://www.gnu.org/software/glpk/>)

  


### Credits and License

Developed at Centre of Biological Engineering, University of Minho (2019-2020)

This project has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement number 814408.  

Released under an Apache License.



