import os
from pathlib import Path

from mewpy.omics import ExpressionSet

from mewpy.io import read_model, Engines, Reader
from mewpy.germ import *
from mewpy.simulation import get_simulator


def read_ecoli_core():
    """
    Reads the model from the given file.
    :return: the model
    """
    # current directory
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    reg_path = path.joinpath('models', 'germ')

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
                               id_col=0,
                               rule_col=2,
                               aliases_cols=[1],
                               header=0)

    # reading the integrated metabolic-regulatory model
    model = read_model(metabolic_reader, regulatory_reader)
    return model


# noinspection DuplicatedCode
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
    simulator.essential_reactions()
    simulator.essential_genes()

    # FBA (default method of the simulator)
    simulator.simulate()

    # pFBA (default method of the simulator)
    simulator.simulate(method='pFBA')

    # FBA
    simulator = FBA(model)
    simulator.optimize()

    # pFBA
    simulator = pFBA(model)
    simulator.optimize()

    # simulating reaction deletions using pFBA
    single_reaction_deletion(model=model,
                             reactions=list(model.reactions.keys())[0:10])

    # simulating genes deletions using FBA
    single_gene_deletion(model=model,
                         genes=list(model.genes.keys())[0:10])

    # simulating FVA
    fva(model=model, fraction=0.9, reactions=list(model.reactions.keys())[0:10])

    # simulating Steady-State Regulatory FBA
    simulator = SRFBA(model)
    simulator.optimize()

    # simulating Regulatory FBA
    simulator = RFBA(model)
    # simulating Steady State Regulatory FBA
    simulator.optimize()
    # simulating Dynamic Regulatory FBA
    simulator.optimize(dynamic=True)

    # Integrated reaction deletion using SRFBA
    isingle_reaction_deletion(model, method='srfba', reactions=list(model.reactions.keys())[0:10])

    # Integrated gene deletion using SRFBA
    isingle_gene_deletion(model, method='srfba', genes=list(model.genes.keys())[0:10])

    # Integrated regulator deletion using SRFBA
    isingle_regulator_deletion(model, method='srfba', regulators=list(model.regulators.keys())[0:10])

    # Integrated FVA using SRFBA
    ifva(model, method='srfba', reactions=list(model.reactions.keys())[0:10])


def read_imc1010():
    """
    Reads the model from the given file.
    :return: the model
    """
    # current directory
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    reg_path = path.joinpath('models', 'germ')

    # E. coli core constraint-based model directory
    cbm_model_f = str(reg_path.joinpath('iJR904.xml'))

    # E. coli core Transcriptional Regulatory Network directory
    reg_model_f = str(reg_path.joinpath('iMC1010.csv'))

    # reader for the metabolic model
    metabolic_reader = Reader(Engines.MetabolicSBML, cbm_model_f)

    # reader for the regulatory model
    regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                               reg_model_f,
                               sep=',',
                               id_col=0,
                               rule_col=4,
                               aliases_cols=[1, 2, 3],
                               header=0)

    # reading the integrated metabolic-regulatory model
    model = read_model(metabolic_reader, regulatory_reader)
    return model


# noinspection DuplicatedCode
def iMC1010_integrated_analysis():
    """
    Performs an integrated analysis of the E. coli iMC1010 integrated model.
    :return:
    """
    model = read_imc1010()

    # Biomass reaction identifier. The model objective function is set to be the biomass reaction, as regular practice.
    _BIOMASS_ID = 'BiomassEcoli'
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
    _GLC = 'EX_glc_DASH_D_e'
    model.get(_GLC).bounds = (-10.0, 100000.0)

    # The MEWpy simulator can be easily created using get_simulator function
    simulator = get_simulator(model)

    # retrieving essential reactions or genes
    simulator.essential_reactions()
    simulator.essential_genes()

    # FBA (default method of the simulator)
    simulator.simulate()

    # pFBA (default method of the simulator)
    simulator.simulate(method='pFBA')

    # FBA
    simulator = FBA(model)
    sol = simulator.optimize()
    print(sol.objective_value)

    # pFBA
    simulator = pFBA(model)
    sol = simulator.optimize()
    print(sol.objective_value)

    # simulating reaction deletions using pFBA
    single_reaction_deletion(model=model,
                             reactions=list(model.reactions.keys())[0:10])

    # simulating genes deletions using FBA
    single_gene_deletion(model=model,
                         genes=list(model.genes.keys())[0:10])

    # simulating FVA
    fva(model=model, fraction=0.9, reactions=list(model.reactions.keys())[0:10])

    # simulating Steady-State Regulatory FBA
    simulator = SRFBA(model)
    sol = simulator.optimize()
    print(sol.objective_value)

    # simulating Regulatory FBA
    simulator = RFBA(model)
    # simulating Steady State Regulatory FBA
    sol = simulator.optimize()
    print(sol.objective_value)
    # simulating Dynamic Regulatory FBA
    simulator.optimize(dynamic=True)

    # Integrated reaction deletion using SRFBA
    isingle_reaction_deletion(model, method='srfba', reactions=list(model.reactions.keys())[0:10])

    # Integrated gene deletion using SRFBA
    isingle_gene_deletion(model, method='srfba', genes=list(model.genes.keys())[0:10])

    # Integrated regulator deletion using SRFBA
    isingle_regulator_deletion(model, method='srfba', regulators=list(model.regulators.keys())[0:10])

    # Integrated FVA using SRFBA
    ifva(model, method='srfba', reactions=list(model.reactions.keys())[0:10])


# noinspection DuplicatedCode
def iNJ661_integrated_analysis():
    """
    Performs an integrated analysis of the iNJ661 integrated model.
    :return:
    """
    # current directory
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    reg_path = path.joinpath('models', 'germ')

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
    expression = ExpressionSet.from_csv(file_path=expression_file, sep=';', index_col=0, header=None)
    quantile_expression, binary_expression = expression.quantile_pipeline()
    initial_state, _ = target_regulator_interaction_probability(model,
                                                                expression=quantile_expression,
                                                                binary_expression=binary_expression)

    prom_ = PROM(model).build()
    prom_.optimize(initial_state=initial_state)

    slim_prom(model, initial_state=initial_state, regulator='Rv0001')


# noinspection DuplicatedCode
def iMM904_integrated_analysis():
    """
    Performs an integrated analysis of the iMM904 integrated model.
    :return:
    """
    # current directory
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    reg_path = path.joinpath('models', 'germ')

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

    # expression = ExpressionSet.from_csv(reg_path.joinpath('iMM904_gene_expression.csv'),
    #                                     sep=';', index_col=0, header=0).dataframe
    # influence = ExpressionSet.from_csv(reg_path.joinpath('iMM904_influence.csv'),
    #                                    sep=';', index_col=0, header=0).dataframe
    # experiments = ExpressionSet.from_csv(reg_path.joinpath('iMM904_experiments.csv'),
    #                                      sep=';', index_col=0, header=0).dataframe
    #
    # gene_expression_prediction = predict_gene_expression(model=model, influence=influence, expression=expression,
    #                                                      experiments=experiments)

    gene_expression_prediction_f = reg_path.joinpath('iMM904_gene_expression_prediction.csv')
    gene_expression_prediction = ExpressionSet.from_csv(gene_expression_prediction_f,
                                                        sep=',',
                                                        index_col=0,
                                                        header=0).dataframe

    initial_state = list(gene_expression_prediction.to_dict().values())
    co_reg_flux = CoRegFlux(model).build()
    co_reg_flux.optimize(initial_state=initial_state[0])

    slim_coregflux(model, initial_state=initial_state[0])

    co_reg_flux = CoRegFlux(model).build()
    metabolites = {'glc__D_e': 16.6, 'etoh_e': 0}
    growth_rate = 0.45
    # time steps in the dataset
    time_steps = list(range(1, 14))
    co_reg_flux.optimize(initial_state=initial_state,
                         metabolites=metabolites,
                         growth_rate=growth_rate,
                         time_steps=time_steps)


if __name__ == '__main__':
    ecoli_core_integrated_analysis()
    iMC1010_integrated_analysis()
    iNJ661_integrated_analysis()
    iMM904_integrated_analysis()
