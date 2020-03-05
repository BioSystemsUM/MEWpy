import os
from collections import OrderedDict
from mewpy.simulation.reframed import Simulation
from mewpy.simulation.simulation import SimulationMethod
from mewpy.utils.constants import ModelConstants
import numpy as np


def is_number(s):
    """Easy way to verify if a string is a number
    """
    return s.replace('.','').replace('-','').replace('e','').isdigit()

class Parser:
    def __init__(self, obj_labels, result_type='KO'):
        self.obj_labels = obj_labels
        self.n_obj = len(obj_labels)
        self.result_type = result_type
        self.obj_labels.append('Size')
        self.results = OrderedDict()
        self.fluxes = OrderedDict()

    def parse_line_ko(self, line, separator):
        line = line.replace(';;', ';').replace(';0;', ';')
        tokens = line.split(separator)
        tokens.remove('\n')
        fit = [abs(float(x)) for x in tokens[:self.n_obj]]
        ko = []
        for i in range(self.n_obj, len(tokens)):
            ko.append(tokens[i])
        id = len(self.results)
        if np.prod(fit) != 0:
            fit.append(len(ko))
            self.results[id] = (fit, ko)

    def parse_line_ou(self, line, separator):
        line = line.replace(';;', ';')
        tokens = line.split(separator)
        tokens.remove('\n')
        fit = [abs(float(x)) for x in tokens[:self.n_obj]]
        ko = {}
        i = self.n_obj
        while i < len(tokens)-2:
            if is_number(tokens[i+2]):
                ko[tokens[i]] = (float(tokens[i+1]), float(tokens[i+2]))
                i = i+3
            else:
                ko[tokens[i]] = (float(tokens[i+1]), float(tokens[i+1]))
                i = i+2
                
            
        id = len(self.results)
        if np.prod(fit) != 0:
            fit.append(len(ko))
            self.results[id] = (fit, ko)

    def parse_results(self, path, separator=';'):
        for r, _, fl in os.walk(path):
            for file in fl:
                if '.csv' in file:
                    # print('Parsing :',file)
                    with open(os.path.join(r, file), 'r') as f:
                        # skip the first line
                        f.readline()
                        line = f.readline()
                        while line:
                            if self.result_type == 'KO':
                                self.parse_line_ko(line, separator)
                            else:
                                self.parse_line_ou(line, separator)
                            line = f.readline()

        # TODO: ... a bit durty...
        data = []
        idx = []
        for k, v in self.results.items():
            idx.append(k)
            data.append(v[0])
        return data, idx


    
    def compute_fluxes(self, model, biomass, reactions, envcond=None):
        # TODO: reactions should be devided in product and carbon sources
        ModelConstants.RESET_SOLVER = True
        model.set_objective({biomass: 1})
        simul = Simulation(model, envcond=envcond)
        res = simul.simulate(
            objective={biomass: 1}, method=SimulationMethod.pFBA)

        self.wt_biomass = res.fluxes[biomass]
        self.wt_products ={}

        for rx in reactions:
            self.wt_products[rx] =  res.fluxes[rx]
        
        data_fba = []
        data_lmoma = []
        data_fvaMin = []
        
        for k, v in self.results.items():
            r = []
            l = []
            m = []
            if self.result_type == 'KO':
                constraints = {rxn: 0 for rxn in v[1]}
            else:
                constraints = v[1]
            res = simul.simulate(objective={biomass: 1}, method=SimulationMethod.pFBA, constraints=constraints)
            
            if res.fluxes:
                #pFBA
                biomass_value = res.fluxes[biomass]
                r.append(biomass_value)
                for rx in reactions:
                    r.append(res.fluxes[rx])
                #lMOMA
                res = simul.simulate(objective={
                    biomass: 1}, method=SimulationMethod.lMOMA, constraints=constraints, reference=res.fluxes)
                l.append(res.fluxes[biomass])
                for rx in reactions:
                    l.append(res.fluxes[rx])
                #minFVA    
                constraints[biomass] = (biomass_value*0.99, 100000.0)
                res = simul.simulate( objective={reactions[0]: 1}, constraints=constraints, maximize=False,method= SimulationMethod.pFBA)
                if res.fluxes:
                    m.append(res.fluxes[biomass])
                else:
                    m.append(0)    
                for rx in reactions:
                    if res.fluxes:
                        m.append(res.fluxes[rx])
                    else:
                        m.append(0)    
                    
            else:
                r.append(0)
                for rx in reactions:
                    r.append(0)
                l.append(0)
                for rx in reactions:
                    l.append(0)
            self.fluxes[k] = (r, l)
            data_fba.append(r)
            data_lmoma.append(l)
            data_fvaMin.append(m)

        return data_fba, data_lmoma, data_fvaMin




    def reaction_stats(self):
        stats = OrderedDict()
        for _, v in self.results.items():
            if self.result_type == 'KO':
                rxns = v[1]
            else:
                rxns = v[1].keys()
            for rx in rxns:
                if rx in stats.keys():
                    stats[rx] += 1
                else:
                    stats[rx] = 1
        self.stats = stats




if __name__ == "__main__":
    """
    import pandas as pd

    BIOMASS_ID = 'R_Ec_biomass_iJO1366_core_53p95M'
    O2 = 'R_EX_o2_LPAREN_e_RPAREN_'
    GLC = 'R_EX_glc_LPAREN_e_RPAREN_'
    PRODUCT_ID = 'R_EX_phe_DASH_L_LPAREN_e_RPAREN_'
    envcond = {GLC: (-10.0,100000.0), O2: (-9.66,100000.0)}
    parser = Parser(['BPCY','WYIELD'],result_type='OU')
    data, idx = parser.parse_results('/home/vmsapereira/Results/PHE/EColi/CBM/OU')
    df = pd.DataFrame(data,columns = parser.obj_labels,index = idx)
    df.describe()
    from reframed.io.sbml import load_cbmodel
    model = load_cbmodel('/home/vmsapereira/Models/w5/iJO1366SL/models/iJO1366SL.xml', flavor='cobra')

    data_fba, data_lmoma , data_fvaMin = parser.compute_fluxes(model,BIOMASS_ID,[PRODUCT_ID,GLC],envcond = envcond)
    print('Wild Type Biomass:', parser.wt_biomass)
    print('Wild type products:', parser.wt_products)
    """
    
    from parse_results import Parser
    import pandas as pd
    
    parser = Parser(['WYIELD','BPCY'],result_type='KO')
    data, idx = parser.parse_results('/home/vmsapereira/Results/TYR/Yeast/GECKO/KO')
    from mewpy.model.gecko import GeckoModel
    model = GeckoModel('single-pool')
    data_fba, data_lmoma, data_fva = parser.compute_fluxes(model,'r_2111',['r_1912','r_1714_REV'])
    print('Wild Type Biomass:', parser.wt_biomass)
    print('Wild type products:', parser.wt_products)