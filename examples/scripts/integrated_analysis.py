import os
from pathlib import Path

import numpy as np

from mewpy.io import read_model, Engines, Reader, read_gene_expression_dataset
from mewpy.mew.analysis import (FBA, pFBA, fva, single_reaction_deletion, single_gene_deletion, SRFBA, RFBA, ifva,
                                isingle_reaction_deletion, isingle_gene_deletion, isingle_regulator_deletion, PROM)
from mewpy.mew.analysis.coregflux import CoRegFlux
from mewpy.omics.preprocessing import quantile_preprocessing_pipeline, target_regulator_interaction_probability
from mewpy.simulation import get_simulator


def read_ecoli_core():
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
    regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                               reg_model_f,
                               sep=',',
                               id_col=1,
                               rule_col=2,
                               aliases_cols=[0],
                               header=0)

    # reading the integrated metabolic-regulatory model
    model = read_model(metabolic_reader, regulatory_reader)
    return model


def ecoli_core_integrated_analysis():
    """
    Performs an integrated analysis of the E. coli core integrated model.
    :return:
    """
    model = read_ecoli_core()

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

    # The MEWpy simulator can be easily created using get_simulator function
    simulator = get_simulator(model)

    # retrieving essential reactions or genes
    essential_reactions = simulator.essential_reactions()
    essential_genes = simulator.essential_genes()

    # FBA (default method of the simulator)
    sol = simulator.simulate()

    # pFBA (default method of the simulator)
    sol = simulator.simulate(method='pFBA')

    # FBA
    simulator = FBA(model)
    sol = simulator.optimize()

    # pFBA
    simulator = pFBA(model)
    sol = simulator.optimize()

    # simulating reaction deletions using pFBA
    reactions_deletion = single_reaction_deletion(model=model,
                                                  reactions=list(model.reactions.keys())[0:10])

    # simulating genes deletions using FBA
    genes_deletion = single_gene_deletion(model=model,
                                          genes=list(model.genes.keys())[0:10])

    # simulating FVA
    fva_sol = fva(model=model, fraction=0.9, reactions=list(model.reactions.keys())[0:10])

    # simulating Steady-State Regulatory FBA
    simulator = SRFBA(model)
    sol = simulator.optimize()

    # simulating Regulatory FBA
    simulator = RFBA(model)
    # simulating Steady State Regulatory FBA
    sol = simulator.optimize()
    # simulating Dynamic Regulatory FBA
    sol = simulator.optimize(dynamic=True)

    # Integrated reaction deletion using SRFBA
    sol = isingle_reaction_deletion(model, method='srfba', reactions=list(model.reactions.keys())[0:10])

    # Integrated gene deletion using SRFBA
    sol = isingle_gene_deletion(model, method='srfba', genes=list(model.genes.keys())[0:10])

    # Integrated regulator deletion using SRFBA
    sol = isingle_regulator_deletion(model, method='srfba', regulators=list(model.regulators.keys())[0:10])

    # Integrated FVA using SRFBA
    sol = ifva(model, method='srfba', reactions=list(model.reactions.keys())[0:10])


def iMM904_integrated_analysis():
    """
    Performs an integrated analysis of the iMM904 integrated model.
    :return:
    """
    # current directory
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    reg_path = path.joinpath('models', 'regulation')

    # iMM904 constraint-based model directory
    cbm_model_f = str(reg_path.joinpath('iMM904.xml'))

    # iMM904 Transcriptional Regulatory Network directory
    reg_model_f = str(reg_path.joinpath('iMM904_trn.csv'))

    # reader for the metabolic model
    metabolic_reader = Reader(Engines.MetabolicSBML, cbm_model_f)

    # reader for the regulatory model
    regulatory_reader = Reader(Engines.CoExpressionRegulatoryCSV,
                               reg_model_f,
                               sep=',',
                               target_col=2,
                               co_activating_col=3,
                               co_repressing_col=4,
                               header=0)

    # reading the integrated metabolic-regulatory model
    model = read_model(metabolic_reader, regulatory_reader)

    # Biomass reaction identifier. The model objective function is set to be the biomass reaction, as regular practice.
    _BIOMASS_ID = 'BIOMASS_SC5_notrace'
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

    # make gene expression predictions using CoRegFlux linear regression method
    # It is commented because it takes around 1 minute to run.
    # Results are stored in the file 'iMM904_gene_expression_prediction.csv'

    # expression_file = reg_path.joinpath('iMM904_gene_expression.csv')
    # influence_file = reg_path.joinpath('iMM904_influence.csv')
    # experiments_file = reg_path.joinpath('iMM904_experiments.csv')
    # expression = read_gene_expression_dataset(expression_file,
    #                                           sep=';',
    #                                           gene_col=0,
    #                                           header=0)
    #
    # influence = read_coregflux_influence_matrix(influence_file,
    #                                             sep=';',
    #                                             gene_col=0,
    #                                             header=0)
    #
    # experiments = read_coregflux_influence_matrix(experiments_file,
    #                                               sep=';',
    #                                               gene_col=0,
    #                                               header=0)
    #
    # gene_expression_prediction = predict_gene_expression(model=model, influence=influence, expression=expression,
    #                                                      experiments=experiments)

    gene_expression_prediction_f = reg_path.joinpath('iMM904_gene_expression_prediction.csv')
    gene_expression_prediction = read_gene_expression_dataset(gene_expression_prediction_f,
                                                              sep=',',
                                                              gene_col=0,
                                                              header=0)

    co_reg_flux = CoRegFlux(model).build()

    initial_state = gene_expression_prediction.iloc[:, 0].to_dict()
    metabolites = {'glc__D_e': 16.6, 'etoh_e': 0}
    growth_rate = 0.45
    time_steps = np.arange(1, 20, 1)
    sol = co_reg_flux.optimize(initial_state=initial_state,
                               metabolites=metabolites,
                               growth_rate=growth_rate,
                               dynamic=True,
                               time_steps=time_steps)

def iNJ661_integrated_analysis():
    """
    Performs an integrated analysis of the iNJ661 integrated model.
    :return:
    """
    # current directory
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    reg_path = path.joinpath('models', 'regulation')

    # iJN661 constraint-based model directory
    cbm_model_f = str(reg_path.joinpath('iNJ661.xml'))

    # iJN661 Transcriptional Regulatory Network directory
    reg_model_f = str(reg_path.joinpath('iNJ661_trn.csv'))

    # reader for the metabolic model
    metabolic_reader = Reader(Engines.MetabolicSBML, cbm_model_f)

    # reader for the regulatory model
    regulatory_reader = Reader(Engines.TargetRegulatorRegulatoryCSV,
                               io=reg_model_f,
                               sep=';',
                               target_col=0,
                               regulator_col=1,
                               header=None)

    # reading the integrated metabolic-regulatory model
    model = read_model(metabolic_reader, regulatory_reader)

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

    # computing PROM target-regulator interaction probabilities using quantile preprocessing pipeline
    expression_file = reg_path.joinpath('iNJ661_gene_expression.csv')
    expression = read_gene_expression_dataset(expression_file, sep=';', gene_col=0, header=None)
    quantile_expression, binary_expression = quantile_preprocessing_pipeline(expression)
    initial_state, _ = target_regulator_interaction_probability(model,
                                                                expression=quantile_expression,
                                                                binary_expression=binary_expression)

    prom = PROM(model).build()
    sol = prom.optimize(initial_state=initial_state)


if __name__ == '__main__':
    ecoli_core_integrated_analysis()
    iMM904_integrated_analysis()
    iNJ661_integrated_analysis()

