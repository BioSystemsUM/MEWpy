# Phenotype Simulation



## Loading metabolic models

Models can be loaded using REFRAMED or COBRApy . 

```python
# using REFRAMED
from reframed.io.sbml import load_cbmodel
model = load_cbmodel('iJO1366.xml', flavor='cobra')

# using COBRApy
from cobra.io import read_sbml_model
model = read_sbml_model('iJO1366.xml')
```

A simulator object provides a common interface to realize the main phenotype analysis tasks. The *get_simulator* function returns a simulator, a wrapper,  for the provided model. The simulator interface remains the same regardless of how the model was loaded, using REFRAMED or COBRApy. This simplify the use of both environments and ease the management of future changes and deprecation on their APIs.

```python
# build a phenotype simulator
from mewpy.simulation import get_simulator
simul = get_simulator(model)
```

The simulator offers a wide API, and enable to perform basic tasks, such as, list  metabolites, reactions, genes and compartments:

```python
simul.metabolites
```

```
['10fthf_c', '12dgr120_c', '12dgr140_c', '12dgr141_c', '12dgr160_c', '12dgr161_c', '12dgr180_c', '12dgr181_c', '12ppd__R_c', '12ppd__S_c', '13dpg_c', '14dhncoa_c', '14glucan_c', '15dap_c', '1ddecg3p_c', '1hdec9eg3p_c', '1hdecg3p_c', '1odec11eg3p_c', 'progly_p', 
...
'psclys_p', 'pser__L_p', 'ptrc_p', 'pydam_p', 'pydx_p', 'pydxn_p', 'pyr_p', 'quin_p', 'r5p_p', 'rfamp_p', 'rib__D_p', 'rmn_p', 'sbt__D_p', 'sel_p', 'ser__D_p', 'ser__L_p', 'skm_p', 'slnt_p', 'so2_p']
```

```python
simul.reactions
```

```
['EX_cm_e', 'EX_cmp_e', 'EX_co2_e', 'EX_cobalt2_e', 'DM_4crsol_c', 'DM_5drib_c', 'DM_aacald_c', 'DM_amob_c', 'DM_mththf_c', 'EX_colipa_e', 'DM_oxam_c', 'EX_glc__D_e', 'EX_glcn_e', 'BIOMASS_Ec_iJO1366_WT_53p95M', 'EX_glcr_e', 'EX_colipap_e', 'EX_glcur_e', 
...
'UPPRT', 'RNDR1', 'URACPAH', 'URAt2pp_copy1', 'RNDR1b', 'URAt2pp_copy2', 'URAtex', 'RNDR2', 'URDGLYCD', 'UREAtex', 'RNDR2b', 'UREAtpp', 'RNDR3', 'RNDR3b', 'RNDR4', 'RNDR4b', 'RNTR1c2', 'RNTR2c2', 'RNTR3c2']
```

```
simul.genes
```

```
['b1377', 'b0241', 'b0929', 'b2215', 'b0653', 'b0655', 'b0118', 'b1276', 'b4032', 'b3359', 'b3528', 'b4033', 'b4035', 'b4034', 'b4036', 'b4213', 'b4123', 'b2835', 'b4138', 'b4077', 'b3503', 'b1064', 'b1747', 'b1539', 'b1748', 'b1090', 'b3475', 'b2836', 'b2563', 'b3553', 
...
'b2235', 'b2676', 'b2675', 'b4238', 'b4237', 'b3386', 'b4301', 'b2914', 'b0638', 'b1679', 'b1682', 'b1681', 'b1680', 'b1683', 'b1684', 'b0222', 'b1745']
```

```python
simul.compartments
```

```
{'c': 'cytosol', 'e': 'extracellular space', 'p': 'periplasm'}
```

A simulator may also be loaded considering environmental conditions, that will be considered during phenotype simulations. In the next example, glucose consumption is limited to 10 mmol/gDW/h while oxygen is set to unlimited.

```python
# environmental conditions
envcond = {'EX_glc__D_e': (-10.0, 100000.0),
           'EX_o2_e':(-1000,1000)}

simul = get_simulator(model,envcond=envcond)
```

All phenotype simulations will consider the imposed environmental conditions, and as such they only need to be set once. Also, these conditions do not persistently alter the model, which can be reused with a different simulator instance.   



## Phenotype simulation

Phenotype simulations are also run using the simulator instance using  the `simulate` method. 

```python
# FBA 
result = simul.simulate()
# or 
result = simul.simulate(method='FBA')
```

```
objective: 0.9823718127269817
Status: OPTIMAL
```

Flux Balance Analysis (FBA) can be run without identifying any method, or  by passing the 'FBA' as method parameter. Other phenotype simulation methods may also be run using one of the identifiers:

- Flux Balance Analysis: `method = 'FBA'`
- Parsimonious FBA:`method = 'pFBA'`
- Minimization of Metabolic Adjustment:`method = 'MOMA'`
- Linear MOMA: `method = 'lMOMA'` 
- Regulatory on/off minimization of metabolic flux: `method = 'ROOM'`

```python
# pFBA
result = simul.simulate(method = 'pFBA')
```

```
objective: 699.0222751839458
Status: OPTIMAL
```



## Reaction fluxes

The phenotype simulation result object, besides the objective value and solver status, also include reaction fluxes in the form of a dictionary:

```python
result.fluxes 
```

```
OrderedDict([('EX_cm_e', -0.0), ('EX_cmp_e', 0.0), ('EX_co2_e', 19.675222635663296), ('EX_cobalt2_e', -2.455929531817453e-05), ('DM_4crsol_c', 0.0002190689142381168), ('DM_5drib_c', 0.00022103365786357076), ('DM_aacald_c', -0.0), ('DM_amob_c', 1.9647436254539624e-06), ('DM_mththf_c', 0.0004401025721020457), ('EX_colipa_e', -0.0), ('DM_oxam_c', -0.0), ('EX_glc__D_e', -10.0), ('EX_glcn_e', 0.0), 
...
('RNDR2b', 0.0), ('UREAtpp', 0.0), ('RNDR3', 0.0), ('RNDR3b', 0.0), ('RNDR4', 0.0), ('RNDR4b', 0.0), ('RNTR1c2', 0.02570474085181419), ('RNTR2c2', 0.02654073926444485), ('RNTR3c2', 0.02654073926444485)])
```

It is also to possible to retrieve reaction fluxes in the form of a data frame:

```python
result.data_frame
```

```
       Reaction ID       Flux
0          EX_cm_e  -0.000000
1         EX_cmp_e   0.000000
2         EX_co2_e  19.675223
3     EX_cobalt2_e  -0.000025
4      DM_4crsol_c   0.000219
...            ...        ...
2578         RNDR4   0.000000
2579        RNDR4b   0.000000
2580       RNTR1c2   0.025705
2581       RNTR2c2   0.026541
2582       RNTR3c2   0.026541

[2583 rows x 2 columns]
```

Individual reaction flux values can be obtained from the dictionary representation. For example, the *Prephenate dehydratase* reaction flux can be obtained from the previous pFBA simulation using the reaction identifier:

```python
result.fluxes['PPNDH']
0.18199911388486417
```



## Retrieving and setting the model objective

The simulation objective, when running FBA or pFBA phenotype simulations, is, by default, the model objective which can be seen using the simulator.

```python
simul.objective
```

```python
{'BIOMASS_Ec_iJO1366_core_53p95M': 1.0}
```

The simulator may also be used to change the model objective, for example, to optimize the ATP maintenance requirement (ATPM) :

```python
simul.objective = 'ATPM'
# or
simul.objective = {'ATPM':1}
```

The last enables to define objectives as linear expressions where the dictionary values are the expression coefficients.



## Adding additional constraints to phenotype simulations

Simulations may include additional metabolic constraints on reaction fluxes. From the previous pFBA simulation one can observe that the organism does not produce L-tyrosine:

```python
result.fluxes['EX_tyr__L_e']
0.0
```

Additional constraints may be added to the model so that the organism start to produce this aromatic amino acid. We may alter, for example, the *3-dehydroquinate dehydratase* reaction bounds, starting by verifying its initial bounds:

```python
# initial bounds
simul.get_reaction_bounds('DHQTi')
(0.0, 1000.0)
```

```python
# additional constraint
constraints ={'DHQTi', (2.995204503155621, 10000)}
# run a pFBA simulation accounting with the new constraint
result = simul.simulate(method='pFBA',constraints=constraints)

result.fluxes['EX_tyr__L_e']
2.8033869589601976
```

We also need to verify that the organism continues to grow:

```python
res.fluxes['BIOMASS_Ec_iJO1366_core_53p95M']
0.5033009222721073
```

It is also possible to plot the production envelope:

```python
from mewpy.visualization.envelope import plot_flux_envelope
plot_flux_envelope(simul,'BIOMASS_Ec_iJO1366_core_53p95M','EX_tyr__L_e',constraints = constraints)
```

![envelope](C:\Users\vmsap\Python\mewpy\docs\envelope.png)

The `simulate` method includes additional parameters, such as the optimization direction. For a full description please refer to the module documentation. 



## Flux Variability Analysis

The simulator interface also allows to perform Flux Variability Analysis (FVA) for L-tyrosine:

```python
simul.FVA(reactions=['EX_tyr__L_e'])
{'EX_tyr__L_e': [0.0, 0.5789655774986189]}
```

By default, MEWpy sets the model objective fraction to 90%, however this fraction may be altered. For example, one might want to consider a fraction of 10% from optimal growth:

```python
simul.FVA(reactions=['EX_tyr__L_e'],obj_frac=0.1)
{'EX_tyr__L_e': [0.0, 5.152853412084507]}
```

The FVA simulations are run considering the defined environmental conditions. Additional constraints may be added, or changed, such as increasing glucose availability. 

```python
simul.FVA(reactions=['EX_tyr__L_e'],obj_frac=0.1, constraints={'EX_glc__D_e': (-20.0, 100000.0)})
{'EX_tyr__L_e': [0.0, 10.366329645681093]}
```

 COBRApy users may have noticed that this same task would have required many additional coding lines if using the COBRApy API directly.

## Genes and reactions essentially

Gene and reaction essentially tests identify , respectively, the list of genes and reactions whose deletion would prevent the organism to grow. 

```python
simul.essential_reactions
['EX_cobalt2_e', 'DM_4crsol_c', 'DM_5drib_c', 'DM_amob_c', 'DM_mththf_c', 'BIOMASS_Ec_iJO1366_core_53p95M', 'EX_k_e', 'EX_cu2_e', 'EX_meoh_e', 'EX_mg2_e', 'EX_mn2_e', 'EX_mobd_e', 'EX_nh4_e', 'EX_ca2_e', 'EX_ni2_e', 'EX_cl_e', 'EX_pi_e', 'EX_zn2_e', 'EX_so4_e', '5DOAN', 'A5PISO', '3OAR140', 'ACCOAC', '3OAS140', 'ACGK', 'ACGS', 'ACHBS', 'ACLS', 'ACODA', 'ACONTa', 'ACONTb', 'ACOTA', 'AGPAT160', 'AGPAT161', 'AGPR', 'ADCL', 'AHCYSNS', 'ADCS', 'AICART', 'AIRC2', 'AIRC3', 'ALAALAr', 'ALAR', 'ADSK', 'ADSL1r', 'ADSL2r', 'ADSS', 'AMAOTr', 'AMPMS2', 'ASP1DC', 'ASPCT', 'ASPK', 'ANPRT', 'ANS', 'AOXSr2', 'ASPTA', 'CAt6pp', 'ATPPRT', 'BMOCOS', 'BMOGDS1', 'BMOGDS2', 'BPNT', 'CDPMEK', 'BTS5', 'CA2tex', 'CHORM', 'CHORS', 'CHRPL', 'CLt3_2pp', 'CLtex', 'COBALT2tex', 'COBALT2tpp', 'CPMPS', 'CYSS', 'CYSTL', 'APRAUR', 'CS', 'CTPS2', 'CU2tex', 'CU2tpp', 'DAPDC', 'DAPE', 'DDPA', 'DASYN160', 'DASYN161', 'DB4PS', 'DHAD1', 'DHAD2', 'DBTS', 'DHDPRy', 'DHDPS', 'DHFR', 'DHFS', 'DNMPPA', 'DNTPPA', 'DHNPA2r', 'DHORTS', 'DPCOAK', 'DPR', 'DHPPDA2', 'ARGSL', 'DHPS2', 'DHPTDCs2', 'ARGSS', 'DHQS', 'DHQTi', 'DTMPK', 'DXPRIi', 'DXPS', 'DMATT', 'E4PD', 'ASAD', 'EGMEACPR', 'EPMEACPR', 'FCLT', 'FMNAT', 'G1PACT', 'G1SAT', 'G3PD2', 'G5SADs', 'GK1', 'GLNS', 'GLUPRT', 'GLUR', 'GLUTRR', 'GLUTRS', 'GCALDD', 'GF6PTA', 'HBZOPT', 'HCO3E', 'GMPS2', 'HISTD', 'HISTP', 'HMBS', 'GRTT', 'HPPK2', 'HSDy', 'GTPCI', 'GTPCII2', 'HSK', 'HSST', 'HSTPT', 'IPMD', 'IPPMIa', 'IPPMIb', 'IPPS', 'K2L4Aabcpp', 'K2L4Aabctex', 'KARA1', 'KARA2', 'KDOCT2', 'KDOPP', 'KDOPS', 'Ktex', 'ICDHyr', 'ICYSDS', 'IG3PS', 'IGPDH', 'IGPS', 'LPADSS', 'ILETA', 'IMPC', 'LEUTAi', 'MALCOAMT', 'MEOHtex', 'MEOHtrpp', 'MEPCT', 'METAT', 'METS', 'MG2tex', 'MCOATA', 'MCTP1App', 'MECDPDH5', 'MECDPS', 'MOADSUx', 'MOAT', 'MOAT2', 'MOBDabcpp', 'MOBDtex', 'MOCOS', 'MOHMT', 'MPTAT', 'MPTG', 'MNtex', 'MPTS', 'MPTSS', 'NH4tex', 'NH4tpp', 'MTHFR2', 'MTHTHFSs', 'NI2tex', 'NADK', 'NNATr', 'NNDPR', 'NADS1', 'NDPK2', 'NDPK4', 'OPMEACPS', 'ORPT', 'P5CR', 'OCBT', 'PANTS', 'OCTDPS', 'OGMEACPD', 'OGMEACPR', 'OGMEACPS', 'OHPBAT', 'OMCDC', 'PAPPT3', 'OMPDC', 'OPHBDC', 'PDX5PS', 'OPMEACPD', 'OPMEACPR', 'PE160abcpp', 'PE161abcpp', 'PGAMT', 'PERD', 'PHETA1', 'PMDPHT', 'PMEACPE', 'PMPK', 'PNTK', 'PItex', 'PPBNGS', 'PPCDC', 'PPNCL2', 'PPND', 'PPNDH', 'PRAGSr', 'PRAIS', 'PRAIi', 'PRAMPC', 'PRASCSi', 'PRATPP', 'PRFGS', 'PRMICI', 'PSCVT', 'PSD160', 'PSD161', 'SADT2', 'SDPDS', 'SDPTA', 'PSSA160', 'SERAT', 'PSSA161', 'PTPATi', 'SHCHD2', 'SHCHF', 'SHK3Dr', 'SHKK', 'SHSL1', 'SO4tex', 'QULNS', 'THRD_L', 'THRS', 'THZPSN3', 'RBFK', 'RBFSa', 'SULR', 'TMDS', 'RBFSb', 'TMPK', 'TMPPP', 'TDSK', 'THDPS', 'TYRL', 'TYRTA', 'U23GAAT', 'UAAGDS', 'UAGAAT', 'USHD', 'UAGCVT', 'UAGDP', 'UAGPT3', 'UAMAGS', 'UAMAS', 'UAPGR', 'UDCPDP', 'UDCPDPS', 'RHCCE', 'Zn2tex', 'UGMDDS', 'UHGADA', 'UMPK', 'UPP3MT', 'UPP3S', 'UPPDC1']
```

```python
simul.essential_genes
['b3359', 'b2019', 'b1096', 'b1812', 'b0180', 'b3360', 'b1093', 'b2323', 'b0827', 'b3857', 'b4214', 'b0775', 'b0004', 'b0159', 'b3196', 'b1094', 'b3040', 'b2750', 'b1131', 'b1208', 'b4177', 's0001', 'b2599', 'b2600', 'b2329', 'b4039', 'b3018', 'b3958', 'b4006', 'b0522', 'b0523', 'b0185', 'b2316', 'b3256', 'b3255', 'b0914', 'b0783', 'b0781', 'b3959', 'b2818', 'b0720', 'b2780', 'b1091', 'b3957', 'b3201', 'b4261', 'b4262', 'b3200', 'b3199', 'b0774', 'b3994', 'b1263', 'b1264', 'b0776', 'b0103', 'b3774', 'b2838', 'b0414', 'b3809', 'b0175', 'b3041', 'b0778', 'b1098', 'b3648', 'b3960', 'b3771', 'b3172', 'b0173', 'b0420', 'b0031', 'b2478', 'b2315', 'b1288', 'b3433', 'b3058', 'b1062', 'b3177', 'b3389', 'b0131', 'b1693', 'b0421', 'b0029', 'b4245', 'b2574', 'b1415', 'b0928', 'b3729', 'b2763', 'b2764', 'b0475', 'b0182', 'b2312', 'b3967', 'b1210', 'b2400', 'b0777', 'b1092', 'b2507', 'b0025', 'b2515', 'b2746', 'b2747', 'b2942', 'b2153', 'b1277', 'b4040', 'b0784', 'b3633', 'b3730', 'b2020', 'b2022', 'b0134', 'b0009', 'b1069', 'b0154', 'b0785', 'b0826', 'b3805', 'b3941', 'b2103', 'b0142', 'b3608', 'b0003', 'b4013', 'b2021', 'b2615', 'b2530', 'b3807', 'b1740', 'b1136', 'b2025', 'b2023', 'b1262', 'b3770', 'b0639', 'b0109', 'b0073', 'b0072', 'b0071', 'b0074', 'b0918', 'b3198', 'b1215', 'b2751', 'b2752', 'b2472', 'b3607', 'b3187', 'b0907', 'b1281', 'b3843', 'b2311', 'b3642', 'b3368', 'b3939', 'b0386', 'b0133', 'b0087', 'b2762', 'b0052', 'b2564', 'b2320', 'b0915', 'b3176', 'b0166', 'b3992', 'b3990', 'b0423', 'b4407', 'b2827', 'b0417', 'b3993', 'b3412', 'b3974', 'b1260', 'b1261', 'b3991', 'b0369', 'b0179', 'b3639', 'b0085', 'b0181', 'b3189', 'b0090', 'b0088', 'b0091', 'b3972', 'b3850', 'b0174', 'b4005', 'b2499', 'b2026', 'b2476', 'b2557', 'b2024', 'b0086', 'b0096', 'b0908', 'b4160', 'b2585', 'b3804', 'b3997', 'b3634', 'b0524', 'b0750', 'b1662', 'b0415', 'b2687']
```



## Simulation time

The MEWpy framework is a tool for Computational Strain Optimization that resorts to Evolutionary Computation to find metabolic modifications that favor defined optimization goals. As Evolutionary Computation algorithms are very time consuming, the time required by each phenotype simulation greatly reflects in the overall computation time of an optimization. REFRAMED and COBRApy have different architecture and consequently phenotype simulations have different time requirements. Experiments may be run to compare the running time of each of the environments. In the next example with compare both on a pFBA simulation, starting by loading a same Genome-scale Metabolic Model (GSMM) and get a simulator for each of the environments.  

```python
from mewpy.simulation import get_simulator

from cobra.io import read_sbml_model
model_cobra = read_sbml_model('iJO1366.xml')
sim_cobra = get_simulator(model_cobra)
```

```python
from reframed.io.sbml import load_cbmodel
model_reframed = load_cbmodel('iJO1366.xml')
sim_reframed = get_simulator(model_reframed)
```

We set as objective the optimization of biomass production `BIOMASS_Ec_iJO1366_core_53p95M`:

```python
sim_cobra.objective = 'BIOMASS_Ec_iJO1366_core_53p95M'
sim_reframed.objective = 'R_BIOMASS_Ec_iJO1366_core_53p95M'
```

It is important to note that REFRAMED uses prefixes on identifiers: 

- "'R\_"  for reactions;
- "G\_" for genes;
- "M\_" for metabolites.

MEWpy makes available a panoply of utilities such as context timer that may be use to asses the time required by phenotype simulations.  

```
from mewpy.utils.utilities import Timer
```

COBRApy pFBA simulation time:

```python
with Timer():
	sim_cobra.simulate(method='pFBA')
	
objective: 699.0222751839475
Status: OPTIMAL
Elapsed time: 2.885654 seconds
```

REFRAMED pFBA simulation time:

```python
with Timer():
	sim_reframed.simulate(method='pFBA')
	
objective: 698.2600975241428
Status: OPTIMAL
Elapsed time: 0.273235 seconds
```

Those results are obtained with *out of the box* configurations, and although the solver, in this particular case CPLEX, may be fine tuned for each of the environments, REFRAMED continues to run faster mainly due to its simplified architecture, especially when constraints are added. While COBRApy offers a broader API, REFRAMED is faster, which is an important feature in combinatorial Computational Strain Optimization. This partially justifies why MEWpy supports both. Still, a user may choose the preferred phenotype simulation environment. 

