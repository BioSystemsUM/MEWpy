
MEWpy examples
============

This folder contains examples on how MEWpy can be used to perform phenotype simulations and Computational Strain Optimization.

## Jupyter Notebooks:

- The [simulation](01-simulation.ipynb) and [optimization](02-optimization.ipynb) notebooks, illustrate how MEWpy may be used to perform basic analysis tasks, phenotype simulations and optimizations using GSMMs.
- The [kinetic](03-kinetic.ipynb) illustrates how MEWpy may be use to load and work with kinetic models.
- The [ROUProblem](04-ROUproblem.ipynb) example illustrates how MEWpy may be used to identify reaction bounds modifications that favor the production of the aromatic amino acid (AAA) L-tyrosine in yeast.
- The [GOUProblem](05-GOUproblem.ipynb) example aims to increase the production of the same AAA in E. coli by modifying genes expression.
- The [GeckoKOProblem](06-GeckoKOProblem.ipynb), also targets the same goal in yeast by exploiting a GECKO model and by deleting enzymatic
- The [GERM_Models](GERM_Models.ipynb) example illustrates how to work with GERM models in MEWpy.
- The [GERM_Models_analysis](GERM_Models_analysis.ipynb) example illustrates how to work with GERM models analysis in MEWpy.
proteins.

## Python scripts

- [load.py](scripts/): Examples on how to load models and get a model specific phenotype simulator.
- [geneopt.py](scripts/): Examples on how to run CSO thar modify genes expression.
- [gpr_eval.py](scripts/): Examples on how to evaluate GPRs.
- [geckosimulation.py](scripts/): Examples on how to run phenotype simulations on GECKO models.
- [cobra_geckoopt.py](scripts/): Examples on how to run CSOs using the original GECKO model.
- [geckoopt.py](scripts/): Examples on how to run CSOs using the original the REFRAME based GECKO model.
- [eacomparison.py](scripts/): Compares the performance of distincts MOEAs on solving a GECKO OU problem.
- [smoment.py](scripts/): Examples on how to run CSOs using sMOMENT and AUTOPACMEN GECKO like models for E.coli.
- [kinetic.py](scripts/): Example on how to use MEWpy for strain design using a kinetic model.
- [optorf.py](scripts/): OptORF examples.
- [optram.py](scripts/): OptRAM example.
- [germ_models_analysis.py](scripts/): Examples on how to perform GERM Model analysis in MEWpy.

More information can be found in the MEWpy documentation https://mewpy.readthedocs.io.