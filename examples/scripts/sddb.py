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

Stain designs from SDDB (https://sddb.bio.di.uminho.pt) may be 
used to conduct computational strain optimizations, retreiving infos
on the designs through REST API.
The script illustrates how to obtain a design using a SDDB ID and 
optimize the mofification folds. This might require the addition of heterologous
genes and reactions.
"""
import requests
from cobra.io import read_sbml_model
from mewpy import get_simulator
from mewpy.util import AttrDict
from mewpy.util.crossmodel import NotationTranslator
from mewpy.problems import GOUProblem
from mewpy.optimization.evaluation import TargetFlux, BPCY
from mewpy.optimization import EA


URL = 'https://sddb.bio.di.uminho.pt/api/solution/'

ECOLI = '../models/ec/iJO1366.xml'
YEAST = '../models/yeast/iMM904.xml.gz'
MAPPING_DB = 'ids-mapping.csv'

def get_modifications_mapping(model, modifications: list, gecko: bool=False):
    
    sim = get_simulator(model)
    ids_not_mapped = []
    mds = {}
    
    if gecko:
        admissible = sim.proteins
        to_notation = "uniprot_id"
    else:
        admissible = sim.genes
        to_notation = "external_id"
        
    notation_translator = NotationTranslator(database=MAPPING_DB, 
                                             admissible=admissible,
                                             from_notation="internal_id",
                                             to_notation=to_notation,
                                             sep=",")

    for modification in modifications:
        try:
            expression = float(modification["expression"])
        except:
            expression = 1
                
        gene = modification["gene"]['name'] 
 
        if gene not in admissible:
            try:
                translation = notation_translator.translate(value=gene)
                if translation in admissible:
                    mds[translation] = expression
                else:
                    ids_not_mapped.append(gene)    
            except ValueError:
                ids_not_mapped.append(gene)    
        else:
            mds[gene] = expression

    return mds, ids_not_mapped


def optimize(solution):
    id = solution['id']
    product = AttrDict(solution['product'])
    organism = AttrDict(solution['organism'])
    modifications = solution['modifications']
    model_file = ECOLI if organism['id']==1 else YEAST
    model = read_sbml_model(model_file)
    FORMULA = product['formula']
    metabolite = None
    sim = get_simulator(model)  
    for met in sim.metabolites:
        m = sim.get_metabolite(met)
        formula = m.formula
        if FORMULA==formula and m.compartment=='e':
            metabolite=met
            break
        
    if not metabolite:
        raise RuntimeError(f"{product['name']} not found")         
    
    m_r_lookup = sim.metabolite_reaction_lookup()
    exs = [x for x in m_r_lookup[met].keys() if x in sim.get_exchange_reactions()]
    
    mds, not_mapped = get_modifications_mapping(sim, modifications)
    if len(mds)==0:
        raise RuntimeError(f'{id} No modifications')
    if len(not_mapped)>0:
        print('Unmapped modifications:',not_mapped)
    target = exs[0]
    _optimize(model,target,list(mds.keys()))
        
def _optimize(model,
              target_rxn,
              modification_targets,
              biomass= None,
              iterations = 2,
              medium = None):
    
    sim = get_simulator(model)
    BIOMASS = biomass if biomass else sim.biomass_reaction
    objectives = [
        TargetFlux(target_rxn),
        BPCY(target_rxn,BIOMASS)
    ]
    problem = GOUProblem(model, 
                         objectives, 
                         target=modification_targets,
                         candidate_min_size= len(modification_targets),
                         candidate_max_size= len(modification_targets),
                         envcond= medium
                         )
    ea = EA(problem,max_generations=iterations)
    ea.run(simplify=False)
    print(ea.dataframe())
    
def main():
    IDS = [335]
    for sid in IDS:
        try:
           response = requests.get(URL+str(sid))
           if (response.status_code == 200):
               solution = response.json()
               optimize(solution)
           else:
               print('error for sid ',sid)
        except Exception as e:
           print(e)



if __name__ == '__main__':
    main()



