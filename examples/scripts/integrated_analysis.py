import os
from pathlib import Path

from mewpy.mew.analysis import (FBA, pFBA, fva, single_reaction_deletion, single_gene_deletion, SRFBA, RFBA, ifva,
                                isingle_reaction_deletion, isingle_gene_deletion, isingle_regulator_deletion)
from mewpy.io import read_model, Engines, Reader
from mewpy.optimization import EA
from mewpy.optimization.evaluation import BPCY, WYIELD
from mewpy.problems import OptORFProblem
from mewpy.simulation import get_simulator
from mewpy.util.io import population_to_csv


def read():
    """
    Reads the model from the given file.
    :return: the model
    """
    # current directory
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    reg_path = path.joinpath('models', 'regulation')

    # E. coli core constraint-based model directory
    cbm_model_f = str(reg_path.joinpath('e_coli_core.xml'))

    # E. coli core Transcriptional Regulatory Network directory
    reg_model_f = str(reg_path.joinpath('e_coli_core_trn.csv'))

    # reader for the metabolic model
    metabolic_reader = Reader(Engines.MetabolicSBML, cbm_model_f)

    # reader for the regulatory model
    regulatory_reader = Reader(Engines.RegulatoryCSV,
                               reg_model_f,
                               sep=',',
                               id_col=1,
                               rule_col=2,
                               aliases_cols=[0],
                               header=0)

    # reading the integrated metabolic-regulatory model
    model = read_model(metabolic_reader, regulatory_reader)
    return model


def integrated_analysis():
    """
    Performs an integrated analysis of the model.
    :return:
    """
    model = read()

    # Biomass reaction identifier. The model objective function is set to be the biomass reaction, as regular practice.
    _BIOMASS_ID = 'Biomass_Ecoli_core'
    model.objective = {_BIOMASS_ID: 1}

    # composition of integrated metabolic-regulatory model
    print(f'Model types: {model.types}')
    print(f'Model current simulators: {model.simulators}')

    print(f'Model interactions: {len(model.interactions)}')
    print(f'Model targets: {len(model.targets)}')
    print(f'Model regulators: {len(model.regulators)}')
    print(f'Model environmental stimuli: {len(model.environmental_stimuli)}')

    print(f'Model objective: {model.objective}')
    print(f'Model reactions: {len(model.reactions)}')
    print(f'Model metabolites: {len(model.metabolites)}')
    print(f'Model genes: {len(model.genes)}')
    print(f'Model sinks: {len(model.sinks)}')
    print(f'Model demands: {len(model.demands)}')
    print(f'Model exchanges: {len(model.exchanges)}')

    print(f'Model compartments: {model.compartments}')
    print(f'Model external compartment: {model.external_compartment}')

    # glucose-exchange reaction identifier. Glucose is the main carbon source for E. coli.
    # Thus, the glucose exchange reaction bounds are set to -10 and 100000.0
    _GLC = 'EX_glc__D_e'
    model.get(_GLC).bounds = (-10.0, 100000.0)

    # -----------------------------------------
    # Flux Analysis using MEWpy simulator
    # -----------------------------------------

    # The MEWpy simulator can be easily created using get_simulator function
    simulator = get_simulator(model)

    # retrieving essential reactions or genes
    essential_reactions = simulator.essential_reactions()
    essential_genes = simulator.essential_genes()

    # FBA (default method of the simulator)
    sol = simulator.simulate()

    # pFBA (default method of the simulator)
    sol = simulator.simulate(method='pFBA')

    # -----------------------------------------
    # Additional Flux Analysis
    # -----------------------------------------

    # FBA
    simulator = FBA(model)
    sol = simulator.optimize()

    # pFBA
    simulator = pFBA(model)
    sol = simulator.optimize()

    # simulating reaction deletions using pFBA
    reactions_deletion = single_reaction_deletion(model=model,
                                                  method='pfba',
                                                  reactions=list(model.reactions.keys())[0:10])

    # simulating genes deletions using FBA
    genes_deletion = single_gene_deletion(model=model,
                                          method='fba',
                                          genes=list(model.genes.keys())[0:10])

    # simulating FVA
    fva_sol = fva(model=model, method='pfba', fraction=0.9, reactions=list(model.reactions.keys())[0:10])

    # -----------------------------------------
    # End of the Additional Flux Analysis
    # -----------------------------------------

    # simulating Steady-State Regulatory FBA
    simulator = SRFBA(model)
    sol = simulator.optimize()

    # simulating Regulatory FBA
    simulator = RFBA(model)
    # simulating Steady State Regulatory FBA
    sol = simulator.optimize()
    # simulating Dynamic Regulatory FBA
    sol = simulator.optimize(dynamic=True)

    # -----------------------------------------
    # Additional Integrated Flux Analysis
    # -----------------------------------------

    # Integrated reaction deletion using SRFBA
    sol = isingle_reaction_deletion(model, method='srfba', reactions=list(model.reactions.keys())[0:10])

    # Integrated gene deletion using SRFBA
    sol = isingle_gene_deletion(model, method='srfba', genes=list(model.genes.keys())[0:10])

    # Integrated regulator deletion using SRFBA
    sol = isingle_regulator_deletion(model, method='srfba', regulators=list(model.regulators.keys())[0:10])

    # Integrated FVA using SRFBA
    sol = ifva(model, method='srfba', reactions=list(model.reactions.keys())[0:10])

    # -----------------------------------------
    # End of the Additional Integrated Flux Analysis
    # -----------------------------------------


def optorf_optimization():
    """
    OptORF optimization
    :return:
    """
    model = read()

    # Biomass reaction identifier. The model objective function is set to be the biomass reaction, as regular practice.
    _BIOMASS_ID = 'Biomass_Ecoli_core'
    model.objective = {_BIOMASS_ID: 1}

    # glucose-exchange reaction identifier. Glucose is the main carbon source for E. coli.
    # Thus, the glucose exchange reaction bounds are set to -10 and 100000.0
    _GLC = 'EX_glc__D_e'
    model.get(_GLC).bounds = (-10.0, 100000.0)

    # setting each metabolic gene state/coefficient to 1 (active state)
    for gene in model.yield_genes():
        gene.coefficient.coefficients = (1,)

    # Succinate is a compound of interest that can be synthesized by E. coli.
    # In this case, it will be the target product for the optimization
    _SUC = 'EX_succ_e'

    # The evaluators Biomass-Product Coupled Yield and Weighted Yield will be used as evaluation metrics
    # of each set of candidates
    evaluator_1 = BPCY(_BIOMASS_ID, _SUC)
    evaluator_2 = WYIELD(_BIOMASS_ID, _SUC)

    # A OptORFProblem problem is created.
    # The OptORF approach is based on gene and regulator deletion to identify optimization strategies.
    # This problem uses a RFBA-like flux analysis to simulate a model constraints derived
    # from the gene and regulator deletions
    problem = OptORFProblem(model, [evaluator_1, evaluator_2], candidate_max_size=6)

    # setting up the Evolutionary Algorithm runner with multiprocessing and maximum of 10 generations
    ea = EA(problem, max_generations=10, mp=True)
    final_pop = ea.run()

    # saving the results
    filename = f'OPTORF_{_SUC}_KO.csv'
    population_to_csv(problem, final_pop, filename, simplify=False)


if __name__ == '__main__':
    integrated_analysis()
    optorf_optimization()
