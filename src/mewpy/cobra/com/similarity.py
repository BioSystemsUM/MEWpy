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
##############################################################################
Models similarity measures

Author: Vitor Pereira
##############################################################################   
"""
import numpy as np
import pandas as pd
from mewpy.simulation import Simulator, get_simulator
from typing import TYPE_CHECKING, Iterable, List, Tuple, Union


if TYPE_CHECKING:
    from cobra.core import Model
    from reframed.core.cbmodel import CBModel

def get_shared_metabolites_counts(model1:Union["Model","CBModel",Simulator],
                                  model2:Union["Model","CBModel",Simulator]
                                 ) -> Tuple[int, int]:
    """Method that returns the number of unique metabolites in both models .

    :param model1: First model
    :param model2: Second model

    :return:
        int: Total number of metabolites in both models
        int: Total number of shared metabolites
    """
    sim1 = get_simulator(model1)
    sim2 = get_simulator(model2)    
    
    met1 = set([x[len(sim1._m_prefix):] for x in sim1.metabolites])
    met2 = set([x[len(sim2._m_prefix):] for x in sim2.metabolites])
    met_ids = set(met1)
    met_ids = met_ids.union(set(met2))
    common_met_ids = met1.intersection(met2)

    return len(met_ids), len(common_met_ids)


def get_shared_reactions_counts(model1:Union["Model","CBModel",Simulator],
                                model2:Union["Model","CBModel",Simulator]
                                ) -> Tuple[int, int]:
    """Computes the number of shared reactions

    
    :param model1: First model
    :param model2: Second model
    
    :return:
        int: Total number of reactions in both models
        int: Total number of shared reactions
    """
    sim1 = get_simulator(model1)
    sim2 = get_simulator(model2)    
    
    rec1 = set([x[len(sim1._r_prefix):] for x in sim1.reactions 
                if sim1.get_reaction_bounds(x)!=(0,0)])
    rec2 = set([x[len(sim2._r_prefix):] for x in sim2.reactions 
                if sim2.get_reaction_bounds(x)!=(0,0)])
    rec_ids = set(rec1)
    rec_ids = rec_ids.union(set(rec2))
    common_rec_ids = rec1.intersection(rec2)

    return len(rec_ids), len(common_rec_ids)


def jaccard_similarity(model1:Union["Model","CBModel",Simulator],
                       model2:Union["Model","CBModel",Simulator]
                       ) -> Tuple[float, float]:
    """Returns the Jacard Similarity of both models with respect to the set
    of metabolites and reactions.

   
    :param model1: First model
    :param model2: Second model

    :return:
        float: Jacard similarity of metabolite sets
        float: Jacard similarity of reaction sets
    """
    sim1 = get_simulator(model1)
    sim2 = get_simulator(model2)    
    
    total_mets, common_mets = get_shared_metabolites_counts(sim1, sim2)
    total_recs, common_recs = get_shared_reactions_counts(sim1, sim2)
    j_met = common_mets / total_mets
    j_rec = common_recs / total_recs
    return j_met, j_rec


def jaccard_similarity_matrices(
    models: Iterable[Union["Model","CBModel",Simulator]]
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The methods takes an Iterable of models and returns a dictionary
    containing all pairwise jaccard similarities for metabolites, reactions and exchange
    reactions (i.e. resource overlap).

    :param models (Iterable): List of models

    :return:
        pd.DataFrame: DataFrame of all jaccard similarities for metabolites indexed by the model ids.
        pd.DataFrame: DataFrame of all jaccard similarities for reaction indexed by the model ids.
        pd.DataFrame: DataFrame for resource overlap indexed by the model ids.
    """
    N = len(models)
    matrix_met = np.eye(N)
    matrix_rec = np.eye(N)
    matrix_ro = np.eye(N)
    index = [model.id for model in models]
    for i, model1 in enumerate(models):
        for j, model2 in zip(range(i + 1, N), models[i + 1 :]):
            ji_met, ji_rec = jaccard_similarity(model1, model2)
            ji_ro = resource_overlap(model1, model2)
            matrix_met[i, j] = ji_met
            matrix_met[j, i] = ji_met
            matrix_rec[i, j] = ji_rec
            matrix_rec[j, i] = ji_rec
            matrix_ro[i, j] = ji_ro
            matrix_ro[j, i] = ji_ro

    df1 = pd.DataFrame(matrix_met, index=index, columns=index)
    df2 = pd.DataFrame(matrix_rec, index=index, columns=index)
    df3 = pd.DataFrame(matrix_ro, index=index, columns=index)

    return df1, df2, df3


def resource_overlap(model1:Union["Model","CBModel",Simulator],
                     model2:Union["Model","CBModel",Simulator]
                    ) -> float:
    """Computes the resource overlap between two models

    
    :param model1: First model
    :param model2: Second model

    :return:
        float: Jacard index of resource overlap
    """
    sim1 = get_simulator(model1)
    sim2 = get_simulator(model2)    
    
    in_ex1 = set([x[len(sim1._r_prefix):] for x in sim1.get_uptake_reactions() 
                  if sim1.get_reaction_bounds(x)!=(0,0)])
    in_ex2 = set([x[len(sim2._r_prefix):] for x in sim2.get_uptake_reactions() 
                  if sim2.get_reaction_bounds(x)!=(0,0)])
    
    common = in_ex1.intersection(in_ex2)
    union = in_ex1.union(in_ex2)

    return len(common) / len(union)


def write_out_common_metabolites(
    models: List[Union["Model","CBModel",Simulator]], prefix: str = "common_reactions.csv"
):
    """Writes out the common reactions as excel sheet and will highligh all
    exchange reaction with yellow color

    
    :param models: List of models
    :param (str) prefix: Name of the file

    :return: DataFrame
    """

    sims = [get_simulator(model) for model in models]
    model = sims[0]
    common_metabolits = [
        model.get_metabolite(rec) for rec in model.metabolites
        if all([rec[len(model._m_prefix):] in [a[len(m._m_prefix):] for a in m.metabolites]
                for m in sims])
    ]
    # Write csv
    df_dict = {"ID": [], "NAME": [], "FORMULA": [], "COMPARTMENT": []}
    for met in common_metabolits:
        df_dict["ID"].append(met.id)
        df_dict["NAME"].append(met.name)
        df_dict["FORMULA"].append(met.formula)
        df_dict["COMPARTMENT"].append(met.compartment)
    df_common_met = pd.DataFrame(df_dict)
    df_common_met.to_csv(prefix)
    return df_common_met


def write_out_common_reactions(
    models: List[Union["Model", "CBModel", Simulator]], prefix: str = "common_metabolites"
) -> None:
    """This writes out the common reactions as excel sheet and will highligh all
    exchange reaction with yellow color

    :param models: List of models
    :param (str) prefix: Name of the file

    :return: DataFrame
    """
    model = get_simulator(models[0])
    common_reactions = [
        model.get_reaction(rec) for rec in model.reactions if all([rec in m.reactions for m in models])
    ]
    # Write csv
    df_dict = {
        "ID": [],
        "NAME": [],
        "REACTION": [],
        "LOWER_BOUND": [],
        "UPPER_BOUND": [],
    }
    for rec in common_reactions:
        df_dict["ID"].append(rec.id)
        df_dict["NAME"].append(rec.name)
        df_dict["REACTION"].append(rec.stoichiometry)
        df_dict["LOWER_BOUND"].append(rec.lb)
        df_dict["UPPER_BOUND"].append(rec.ub)


    df_common_rec = pd.DataFrame(df_dict)
    df_common_rec.to_csv(prefix)
    return df_common_rec
