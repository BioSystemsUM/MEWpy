from pathlib import Path

from mewpy.io import read_model, Engines, Reader
from mewpy.optimization import EA
from mewpy.optimization.evaluation import BPCY, WYIELD
from mewpy.problems import OptORFProblem
from mewpy.util.io import population_to_csv


def optorf_ec():
    """
    The OptORF approach was used in the E. coli core model to identify metabolic engineering strategies
    for overproduction of succinate.
    :return:
    """
    path = Path(__file__).parent.parent.joinpath('models', 'germ')
    metabolic_reader = Reader(Engines.MetabolicSBML, path.joinpath('e_coli_core.xml'))
    regulatory_reader = Reader(Engines.BooleanRegulatoryCSV, path.joinpath('e_coli_core_trn.csv'),
                               sep=',', id_col=0, rule_col=2, aliases_cols=[1], header=0)
    model = read_model(metabolic_reader, regulatory_reader)

    BIOMASS_ID = 'Biomass_Ecoli_core'
    GLC = 'EX_glc__D_e'
    PRODUCT_ID = 'EX_succ_e'

    # objective
    model.objective = {BIOMASS_ID: 1}
    model.get(GLC).bounds = (-10.0, 100000.0)
    model.get(BIOMASS_ID).lower_bound = 0.1

    evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID)
    evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)
    problem = OptORFProblem(model, [evaluator_1, evaluator_2], candidate_max_size=10)

    ea = EA(problem, max_generations=10, mp=True)
    final_pop = ea.run()

    filename = "OPTORF_{}_KO_{}.csv".format(PRODUCT_ID, "ec")
    population_to_csv(problem, final_pop, filename, simplify=False)


def optorf_imc():
    """
    The OptORF approach was used in the iMC1010 model to identify metabolic engineering strategies
    for overproduction of succinate.
    :return:
    """
    path = Path(__file__).parent.parent.joinpath('models', 'germ')
    metabolic_reader = Reader(Engines.MetabolicSBML, path.joinpath('iJR904.xml'))
    regulatory_reader = Reader(Engines.BooleanRegulatoryCSV, path.joinpath('iMC1010.csv'),
                               sep=',', id_col=0, rule_col=4, aliases_cols=[1, 2, 3], header=0)
    model = read_model(metabolic_reader, regulatory_reader)

    BIOMASS_ID = 'BiomassEcoli'
    GLC = 'EX_glc_DASH_D_e'
    PRODUCT_ID = 'EX_succ_e'

    initial_state = {
        'Stringent': 0.0,
        'high-NAD': 0.0,
        'AGDC': 0.0,
    }

    model.objective = {BIOMASS_ID: 1}
    model.get(GLC).bounds = (-18.5, 0.0)
    model.get(BIOMASS_ID).lower_bound = 0.1

    evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID)
    evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)
    problem = OptORFProblem(model, [evaluator_1, evaluator_2], initial_state=initial_state, candidate_max_size=6)

    ea = EA(problem, max_generations=10, mp=True)
    final_pop = ea.run()

    filename = "OPTORF_{}_KO_{}.csv".format(PRODUCT_ID, "iJR904_srfba")
    population_to_csv(problem, final_pop, filename, simplify=False)


if __name__ == '__main__':
    optorf_ec()
    optorf_imc()
