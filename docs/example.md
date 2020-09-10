# Examples



## Metabolic Constraints



### Simulation



MEWpy enables to use REFRAMED and COBRApy phenotype simulators seamlessly, depending on the loaded model.

Next, an example is presented where FBA and pFBA simulations are run over a model loaded with REFRAMED or COBRApy. 



Loading and FBA simulation

```python
# using REFRAMED
from reframed.io.sbml import load_cbmodel
model = load_cbmodel('iJO1366SL.xml', flavor='cobra')

# or using COBRApy
from cobra.io import read_sbml_model
model = read_sbml_model('iJO1366SL.xml')

# build a phenotype simulator
from mewpy.simulation import get_simulator, SimulationMethod   
simul = get_simulator(model)
simul.summary()

# FBA 
simul.simulate()

# pFBA
simul.simulate(method = SimulationMethod.pFBA)
```



Simulations may include additional metabolic constraints on reaction fluxes, for example by limiting oxygen consumption.



```python
# define reaction fluxes constraints
constraints = {'R_EX_o2_LPAREN_e_RPAREN_': (-5.0, 100000.0)}

# FBA 
simul.simulate(constraints = constraints)

```



### Optimization

 The definition of an optimization problem encompasses specifying the optimization objectives (BPCY, WYIELD, etc ...), environmental conditions and EA configurations. 



```python  
# load the model
from reframed.io.sbml import load_cbmodel
model = load_cbmodel('iJO1366SL.xml', flavor='cobra')

# Define the target
PRODUCT_ID = 'R_EX_tyr_DASH_L_LPAREN_e_RPAREN_'
BIOMASS_ID = 'R_Ec_biomass_iJO1366_core_53p95M'

# environmental conditions
envcond = {'R_EX_o2_LPAREN_e_RPAREN_'  : (-9.66, 100000.0),
           'R_EX_glc_LPAREN_e_RPAREN_' : (-12.5,100000.0)}
}

# Optimization objectives
from mewpy.optimization.evaluation import  BPCY, WYIELD
evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method=SimulationMethod.lMOMA)
evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)

# build a new problem instance
from mewpy.problems import ROUProblem
problem = ROUProblem(model, fevaluation=[
                         evaluator_1, evaluator_2], envcond=envcond)

# run the optimization
from mewpy.optimization import EA
ea = EA(problem, max_generations= 100, visualizer=True)
final_pop = ea.run()

```



Optimizations over genes' expression modifications are run by setting and running the intended problem. Gene deletion optimization problems are defined a GKOProblem while gene over- or under expression optimization using the GOUProblem class.



```python  

# build a new problem instance
from mewpy.problems import GOUProblem
problem = GOUProblem(model, fevaluation=[
                         evaluator_1, evaluator_2], envcond=envcond)

# run the optimization
from mewpy.optimization import EA
ea = EA(problem, max_generations= 100, visualizer=True)
final_pop = ea.run()

```







## Enzymatic Constraints 



MEWpy enables strain optimization using Genome-scale models enhanced with enzymatic (kcat) parameters and enzyme mass constraints:

* MEWpy supports [GECKO](https://doi.org/10.15252/msb.20167411) models, from the original COBRApy based implementation, but also implemented over REFRAMED package.

* MEWpy also supports sMOMENT and GECKO like models obtained from [AutoPACMEN](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3329-9) . 





### Simulation

Loading a GECKO model and performing a FBA simulation

```python
# load the GECKO model
from mewpy.model.gecko import GeckoModel
model = GeckoModel('single-pool')

# or using the original geckopy model based on COBRApy
from from geckopy import GeckoModel
model = GeckoModel('single-pool')

# build a phenotype simulator
from mewpy.simulation import get_simulator, SimulationMethod   
simul = get_simulator(model)
simul.summary()

# FBA 
simul.simulate()

# pFBA
simul.simulate(method = SimulationMethod.pFBA)

```



### Optimization

The optimization API is common to the one defined for the optimization of metabolic constraints (Reactions and Genes),  . MEWpy automatically  selects the phenotype simulator for the loaded model.

```python  
# load the model
from mewpy.model.gecko import GeckoModel
model = GeckoModel('single-pool')

# Define the target
PRODUCT_ID = 'r_1913'
BIOMASS_ID = 'r_2111'

# environmental conditions
envcond = {'r_1714_REV' : (-12.5,100000.0)}


# Optimization objectives
from mewpy.optimization.evaluation import  BPCY, WYIELD
evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method=SimulationMethod.lMOMA)
evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)

# build a new problem instance for enzymatic OU
from mewpy.problems import GeckoROUProblem
problem = GeckoROUProblem(model, fevaluation=[
                         evaluator_1, evaluator_2], envcond=envcond)

# run the optimization
from mewpy.optimization import EA
ea = EA(problem, max_generations= 100, visualizer=True)
final_pop = ea.run()

```







## Regulatory Constraints



MEWpy implements computational strain design optimization with regulatory constraints. Presently two methods are available, [OptRAM](https://doi.org/10.1371/journal.pcbi.1006835) and OptORF.



### OptRAM Example



```python
from mewpy.regulation.optram import OptRamProblem, load_optram
  
# regulatory matrix Genes x TFs   
matrix_file = 'regnet.csv'
# csv file mapping genes names entries in the regulatory matrix 
gene_file = 'mgene.csv'
# csv with TFs expression 
tf_file ='TFnames.csv'

BIOMASS_ID = 'r_2111'
PRODUCT_ID = 'r_1913' #TYR
GLC = 'r_1714'

# build the regulatory network
# add the prefix 'G_' to genes. Only for REFRAMED models
regnet = load_optram(gene_file, tf_file, matrix_file, gene_prefix='')

# load the model
from cobra.io import read_sbml_model
model = read_sbml_model('yeast_7.6-optram.xml')


# define the optimization objectives
from mewpy.simulation import SimulationMethod
from mewpy.optimization.evaluation import BPCY, WYIELD

evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method=SimulationMethod.lMOMA)
evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID, parsimonious=True)

# environmental conditions
envcond = {GLC:(-12.5,10000)}

# instantiate the problem
problem = OptRamProblem(model, [evaluator_1, evaluator_2],
                            regnet, candidate_min_size=1, candidate_max_size=6, envcond = envcond)


from mewpy.optimization import EA
ea = EA(problem, max_generations=100, mp=True)
final_pop = ea.run()
```
