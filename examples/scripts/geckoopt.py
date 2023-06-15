# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
Author: Vitor Pereira

GECKO E. coli model optimization
"""

from collections import OrderedDict
from time import time

from mewpy.model.gecko import GeckoModel
from mewpy.optimization import EA
from mewpy.optimization.evaluation import BPCY, WYIELD, TargetFlux, ModificationType
from mewpy.problems.gecko import GeckoKOProblem, GeckoOUProblem

ITERATIONS = 10


def ec_gecko_ko(compound, display=False, filename=None):
    """ GECKO enzyme deletion example.
    It runs a multi objective optimization for the increased production of a certain compound on E. Coli.
    The GECKO model is the yeast model companion from the GECKO paper "Improving the phenotype predictions
    of a yeast genome‐scale metabolic model by incorporating enzymatic
    constraints" https://doi.org/10.15252/msb.20167411.
    Runs over the MEWpy implementation.

    :param compound: A target reaction identifier.
    :param display: Prints the best solution.
    :param filename: If given, saves the results as csv to filename.

    """

    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(dir_path, '../../../examples/models/gecko')
    DATA_FILE = os.path.join(PATH, 'eciML1515_batch.xml')
    from reframed.io.sbml import load_cbmodel
    m = load_cbmodel(DATA_FILE)

    model = GeckoModel(m, biomass_reaction_id='R_BIOMASS_Ec_iML1515_core_75p37M',
                       protein_pool_exchange_id='R_prot_pool_exchange', reaction_prefix='R_')
    model.set_objective({'R_BIOMASS_Ec_iML1515_core_75p37M': 1.0})

    envcond = OrderedDict()

    # the evaluation (objective) functions

    evaluator_1 = BPCY("R_BIOMASS_Ec_iML1515_core_75p37M", compound, method='lMOMA')
    evaluator_2 = WYIELD("R_BIOMASS_Ec_iML1515_core_75p37M", compound)
    # The optimization problem
    problem = GeckoKOProblem(model,
                             fevaluation=[evaluator_1, evaluator_2],
                             envcond=envcond,
                             prot_prefix='R_draw_prot_',
                             candidate_max_size=6)

    # A new instance of the EA optimizer
    ea = EA(problem, max_generations=ITERATIONS)
    # runs the optimization
    ea.run()
    # optimization results
    if display:
        print(ea.dataframe())
    # save final population to file
    if filename:
        print("Saving solutions to file")
        ea.dataframe().to_csv(filename)


def ec_gecko_ou(compound, display=False, filename=None):
    """ GECKO enzyme over/under expression example.
    It runs a multi objective optimization for the increased production of a certain compound on E. Coli.
    The GECKO model is the yeast model companion from the GECKO paper "Improving the phenotype predictions
    of a yeast genome‐scale metabolic model by incorporating enzymatic
    constraints" https://doi.org/10.15252/msb.20167411.
    Runs over the MEWpy implementation.

    :param compound: A target reaction identifier.
    :param display: Prints the best solution.
    :param filename: If given, saves the results as csv to filename.

    """
    
    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    PATH = os.path.join(dir_path, '../../examples/models/ec/')
    DATA_FILE = os.path.join(PATH, 'eciML1515_batch.xml')
    print(DATA_FILE)
    from reframed.io.sbml import load_cbmodel
    m = load_cbmodel(DATA_FILE)
    model = GeckoModel(m, biomass_reaction_id='R_BIOMASS_Ec_iML1515_core_75p37M',
                       protein_pool_exchange_id='R_prot_pool_exchange', reaction_prefix='R_')
    model.set_objective({'R_BIOMASS_Ec_iML1515_core_75p37M': 1.0})

    # change protein pool bound (suggested by Leslie)
    model.reactions['R_prot_pool_exchange'].ub = 0.26
    # define environmental consitions
    envcond = OrderedDict()

    # the evaluation (objective) functions
    evaluator_1 = BPCY("R_BIOMASS_Ec_iML1515_core_75p37M", compound, method='lMOMA')
    # FVA MAX is strangely very high... changing the default alpha (0.3) to compensate..
    evaluator_2 = WYIELD("R_BIOMASS_Ec_iML1515_core_75p37M", compound, alpha=0.01)

    evaluator_3 = TargetFlux(compound)

    # The optimization problem
    problem = GeckoOUProblem(model,
                             fevaluation=[evaluator_1, evaluator_2, evaluator_3],
                             envcond=envcond,
                             prot_prefix='R_draw_prot_',
                             candidate_max_size=30)

    # A new instance of the EA optimizer
    ea = EA(problem, max_generations=ITERATIONS, mp=True,algorithm='NSGAIII')
    # runs the optimization
    final_pop = ea.run()
    # optimization results
    if display:
        print(ea.dataframe())
    # save final population to file
    if filename:
        print("Saving solutions to file")
        ea.dataframe().to_csv(filename)



if __name__ == '__main__':

    N_EXP = 1

    compounds = {'TYR':'R_EX_tyr__L_e'}

    for k, v in compounds.items():
        for _ in range(N_EXP):
            millis = int(round(time() * 1000))
            ec_gecko_ou(v, 
                        display=True,
                        filename="gecko_{}_OU_{}.csv".format(k, millis))

