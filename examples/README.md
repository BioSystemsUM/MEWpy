
MEWpy examples
============


This folder contains examples on how MEWpy can be used to perform phenotype simulations and Computational Strain Optimization. 

## Jupyter Notebooks:

- The firs example, [simulation](simulation.ipynb), illustrates how MEWpy may be used to perform basic analysis tasks  and phenotype simulations with GSMMs.

- The second example, [ROUProblem](ROUproblem.ipynb), illustrates how MEWpy may be used to identify reaction bounds modifications that favor the production of the aromatic amino acid (AAA) L-tyrosine in yeast.

- The third example, [GOUProblem](GOUproblem.ipynb), aims to increase the production of the same AAA in E. coli by modifying genes expression.

- The last example, [GeckoKOProblem](GeckoKOProblem.ipynb), also targets the same goal in yeast by exploiting a GECKO model and by deleting enzymatic proteins.

## Python scripts

- [load.py](scripts/): Examples on how to load models and get a model specific phenotype simulator.
- [geneopt.py](scripts/): Examples on how to run CSO thar modify genes expression.
- [gpr_eval.py](scripts/): Examples om how to evaluate GPRs.
- [geckosimulation.py](scripts/): Examples on how to run phenotype simulations on GECKO models.
- [cobra_geckoopt.py](scripts/): Examples on how to run CSOs using the original GECKO model.
- [geckoopt.py](scripts/): Examples on how to run CSOs using the original the REFRAME based GECKO model.
- [eacomparison.py](scripts/): Compares the performance of distincts MOEAs on solving a GECKO OU problem.
- [smoment.py](scripts/): Examples on how to run CSOs using sMOMENT and AUTOPACMEN GECKO like models for E.coli.
- [optram.py](scripts/): OptRAM example.
- [rfba.py](scripts/): rFBA example.
- [srfba.py](scripts/): srFBA example.
- [optorf.py](scripts/): OptORF examples.


More information can be found in the MEWpy documentation https://mewpy.readthedocs.io.