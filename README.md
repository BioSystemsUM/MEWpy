[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI version](https://badge.fury.io/py/mewpy.svg)](https://badge.fury.io/py/mewpy)
[![Documentation Status](https://readthedocs.org/projects/mewpy/badge/?version=latest)](https://mewpy.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/BioSystemsUM/mewpy.svg?branch=master)](https://travis-ci.org/BioSystemsUM/mewpy)

# MEWpy

MEWpy is an integrated Metabolic Engineering Workbench for strain design optimization. It offers methods to explore different classes of constraint-based models (CBM) for:

- Simulation: allows to simulate  steady-state metabolic models, considering different formulations (e.g., GECKO, ETFL) and kinetic models;
- Omics data integration (eFlux, GIMME, iMAT);
- Optimization: performs Evolutionary Computation based strain design optimization by knocking out (KO) or over/under expressing (OU) reactions, genes or enzymes.

MEWPy currently supports [REFRAMED](https://github.com/cdanielmachado/reframed) and [COBRApy](https://opencobra.github.io/cobrapy/) simulation environments. The optimization engine relies on either [inspyred](https://github.com/aarongarrett/inspyred) or [jMetalPy](https://github.com/jMetal/jMetalPy) packages.

### Examples

Examples are provided as [jupiter notebooks](examples) and as [python scripts](examples).

### Documentation

The package documentation is available at [mewpy.readthedocs.io](https://mewpy.readthedocs.io).

### Instalation

Installing from Pypi package repository:

`pip install mewpy`

Installing from github:

1. clone the repository

`git clone https://github.com/BioSystemsUM/mewpy.git -b master`

2. run `python setup.py install`

MEWPy requires a compatible linear programming solver, with installed Python dependencies, from the following list:

- [CPLEX](https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) 
- [GUROBI](http://www.gurobi.com)
- [GLPK](https://www.gnu.org/software/glpk/)

### Cite

Vítor Pereira, Fernando Cruz, Miguel Rocha, MEWpy: a computational strain optimization workbench in Python, Bioinformatics, 2021;, btab013, [https://doi.org/10.1093/bioinformatics/btab013](https://doi.org/10.1093/bioinformatics/btab013)

### Credits and License

Developed at Centre of Biological Engineering, University of Minho (2019- )

This project has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement number 814408.

Released under an Apache License.
