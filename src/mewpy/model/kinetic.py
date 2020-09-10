from collections import OrderedDict
from libsbml import AssignmentRule,SBMLReader
from reframed.core.model import Model
import re
import os 
import warnings



def load_sbml(filename):
    """ Loads an SBML file.
    
    :param filename: SBML file path, str.
    :returns: SBMLModel
            
    """

    if not os.path.exists(filename):
        raise IOError("Model file was not found")

    reader = SBMLReader()
    document = reader.readSBML(str(filename))
    sbml_model = document.getModel()

    if sbml_model is None:
        document.printErrors()
        raise IOError(f'Failed to load model {filename}.')

    return sbml_model


def load_ODEModel(filename, map = None):
    sbml_model = load_sbml(filename)
    from reframed.io.sbml import load_compartments,load_metabolites,load_reactions
    ode_model = ODEModel(sbml_model.getId())
    load_compartments(sbml_model, ode_model)
    load_metabolites(sbml_model, ode_model)
    ## adds constants and boundaries to metabolites
    for species in sbml_model.getListOfSpecies():
        m = ode_model.metabolites[species.getId()]
        if not hasattr(m,"constant"):
            m.constant = species.getConstant()
        if not hasattr(m,"boundary"):
            m.boundary = species.getBoundaryCondition()
    load_reactions(sbml_model, ode_model)
    _load_concentrations(sbml_model, ode_model)
    _load_global_parameters(sbml_model, ode_model)
    _load_local_parameters(sbml_model, ode_model)
    _load_ratelaws(sbml_model, ode_model)
    _load_assignment_rules(sbml_model, ode_model)
    # parse rates, rules and xdot expressions
    ode_model._set_parsed_attr()
    if isinstance(map,str):
        print("MAP =",map)
        aux = OrderedDict([(rId, re.findall(map, ratelaw)) for rId, ratelaw in ode_model.ratelaws.items()])
        #aux = OrderedDict([(rId ,  [rId+"_"+x for x in re.findall("(rmax\w*)", ratelaw)]) for rId, ratelaw in model.ratelaws.items()])#CHASSAGNOLE
        ode_model.reacParamsFactors = OrderedDict([(rId , params) for rId, params in aux.items() if len(params) > 0])
    else:
        ode_model.set_reactions_parameters_factors(map)

        
    return ode_model




def _load_concentrations(sbml_model, odemodel):
    for species in sbml_model.getListOfSpecies():
        odemodel.set_concentration(species.getId(), species.getInitialConcentration())


def _load_global_parameters(sbml_model, odemodel):
    for parameter in sbml_model.getListOfParameters():
            odemodel.set_global_parameter(parameter.getId(), parameter.getValue(), parameter.getConstant())


def _load_local_parameters(sbml_model, odemodel):
    for reaction in sbml_model.getListOfReactions():
        for parameter in reaction.getKineticLaw().getListOfParameters():
            odemodel.set_local_parameter(reaction.getId(), parameter.getId(), parameter.getValue())


def _load_ratelaws(sbml_model, odemodel):
    for reaction in sbml_model.getListOfReactions():
        odemodel.set_ratelaw(reaction.getId(), reaction.getKineticLaw().getFormula())


def _load_assignment_rules(sbml_model, odemodel):
    for rule in sbml_model.getListOfRules():
        if isinstance(rule, AssignmentRule):
            odemodel.set_assignment_rule(rule.getVariable(), rule.getFormula())





class ODEModel(Model):

    def __init__(self, model_id):
        """ ODE Model.
        
        :param model: a REFRAMED SBModel or COBRApy Model 
            
        """
        super(ODEModel,self).__init__(model_id)
        self.concentrations = OrderedDict()
        self.constant_params = OrderedDict()
        self.variable_params = OrderedDict()
        self.local_params = OrderedDict()
        self.ratelaws = OrderedDict()
        self.assignment_rules = OrderedDict()
        self._func_str = None
        self._constants = None
        self.reacParamsFactors = None  
        self.parsedRates = None
        self.parsedRules = None
        self.parsedXdot = None
        


    def get_reactions(self):
        return self.reactions

    def _clear_temp(self):
        self.update()
        self._func_str = None

    def add_reaction(self, reaction, replace=True, ratelaw=''):
        super(ODEModel,self).add_reaction(reaction,replace)
        self.ratelaws[reaction.id] = ratelaw
        self.local_params[reaction.id] = OrderedDict()

    def set_concentration(self, m_id, concentration):
        if m_id in self.metabolites:
            self.concentrations[m_id] = concentration
        else:
            warnings.warn("No such metabolite '{}'".format(m_id), RuntimeWarning)

    def set_ratelaw(self, r_id, ratelaw):
        if r_id in self.reactions:
            self.ratelaws[r_id] = ratelaw
        else:
            warnings.warn("No such reaction '{}'".format(r_id), RuntimeWarning)

    def set_assignment_rule(self, p_id, rule):
        if p_id in self.variable_params or p_id in self.metabolites:
            self.assignment_rules[p_id] = rule
        else:
            warnings.warn("No such variable parameter '{}'".format(p_id), RuntimeWarning)

    def set_global_parameter(self, key, value, constant=True):
        if constant:
            self.constant_params[key] = value
        else:
            self.variable_params[key] = value

    def set_local_parameter(self, r_id, p_id, value):
        if r_id in self.reactions:
            if r_id not in self.local_params.keys():
                self.local_params[r_id] = OrderedDict()
            self.local_params[r_id][p_id] = value
        else:
            warnings.warn("No such reaction '{}'".format(r_id), RuntimeWarning)

    def remove_reactions(self, id_list):
        self.remove_reactions(id_list)
        for r_id in id_list:
            del self.ratelaws[r_id]
            del self.local_params[r_id]
            

    def merge_constants(self):
        constants = OrderedDict()

        for c_id, comp in self.compartments.items():
            constants[c_id] = comp.size

        constants.update(self.constant_params)

        for r_id, params in self.local_params.items():
            for p_id, value in params.items():
                full_id = '{}_{}'.format(r_id, p_id)
                constants[full_id] = value

        self._constants = constants
        return constants

    def get_parameters(self, exclude_compartments=False):
        if not self._constants:
            self.merge_constants()

        parameters = self._constants.copy()

        if exclude_compartments:
            for c_id in self.compartments:
                del parameters[c_id]

        return parameters


    def print_balance(self, m_id):
        c_id = self.metabolites[m_id].compartment
        table = self.metabolite_reaction_lookup()
        terms = ["{:+g} * r['{}']".format(coeff, r_id) for r_id, coeff in table[m_id].items()]
        if len(terms)==0 or (self.metabolites[m_id].constant and self.metabolites[m_id].boundary):
            expr= "0"
        else:
            expr = "1/p['{}'] * ({})".format(c_id, ' '.join(terms))
        return expr


    def parse_rate(self, r_id, rate):

        symbols = '()+*-/,'
        rate = ' ' + rate + ' '
        for symbol in symbols:
            rate = rate.replace(symbol, ' ' + symbol + ' ')

        for i, m_id in enumerate(self.metabolites):
            rate = rate.replace(' ' + m_id + ' ', ' x[{}] '.format(i))

        for c_id in self.compartments:
            rate = rate.replace(' ' + c_id + ' ', " p['{}'] ".format(c_id))

        for p_id in self.constant_params:
            if p_id not in self.local_params[r_id]:
                rate = rate.replace(' ' + p_id + ' ', " p['{}'] ".format(p_id))

        for p_id in self.variable_params:
            if p_id not in self.local_params[r_id]:
                rate = rate.replace(' ' + p_id + ' ', " v['{}'] ".format(p_id))

        for p_id in self.local_params[r_id]:
            rate = rate.replace(' ' + p_id + ' ', " p['{}_{}']".format(r_id, p_id))

        return rate

    def parse_rule(self, rule, parsed_rates):

        symbols = '()+*-/,'
        rule = ' ' + rule + ' '
        for symbol in symbols:
            rule = rule.replace(symbol, ' ' + symbol + ' ')

        for i, m_id in enumerate(self.metabolites):
            rule = rule.replace(' ' + m_id + ' ', ' x[{}] '.format(i))

        for c_id in self.compartments:
            rule = rule.replace(' ' + c_id + ' ', " p['{}'] ".format(c_id))

        for p_id in self.constant_params:
            rule = rule.replace(' ' + p_id + ' ', " p['{}'] ".format(p_id))

        for p_id in self.variable_params:
            rule = rule.replace(' ' + p_id + ' ', " v['{}'] ".format(p_id))

        for r_id in self.reactions:
           rule = rule.replace(' ' + r_id + ' ', '({})'.format(parsed_rates[r_id]))

        return rule

    def _set_parsed_attr(self):
        self.parsedRates = {rId: self.parse_rate(rId, ratelaw)
                             for rId, ratelaw in self.ratelaws.items()}
        
        aux = {pId: self.parse_rule(rule, self.parsedRates)
               for pId, rule in self.assignment_rules.items()}
        
        trees = [_build_tree_rules(vId, aux) for vId in aux.keys()]
        order = _get_oder_rules(trees)
        
        self.parsedRules = OrderedDict([(id, aux[id]) for id in order])
        self.parsedXdot = {mId: self.print_balance(mId) for mId in self.metabolites}



    def build_ode(self, factors):
        """
        Builds de ODE model
        
            
        :param factors: (dict) The key is the parameter identifier and the value is the level of change values between 0 and 1 represent a under expression, above 1 a over expression and 0 to represent the knockouts.
        :returns: Returns  a string with the ode system.
        
        """

        # factors: ["vmax1": 0, "vmax2"=2, "ENZYME_ID":0]
        # divide vmax parameters from enzymes expression levels
        factorsEnz = OrderedDict([(k, v) for k, v in factors.items() if k in self.metabolites])
        factorsParam = OrderedDict([(k, v) for k, v in factors.items() if k not in factorsEnz.keys()])

        ruleExprs = ["    v['{}'] = {}".format(pId, self.parsedRules[pId])
                     for pId in self.parsedRules.keys()]

        rateExprs = []
        for rId in self.reactions:
            newExp = self.parsedRates[rId]
            if rId in self.reacParamsFactors.keys():
                toModify = set(factorsParam.keys()).intersection(self.reacParamsFactors[rId])
                if len(toModify) > 0:
                    for elem in toModify:
                        newExp = re.sub(r"([pv]\['" + elem + "'\])", str(factorsParam[elem]) + r" * \1", newExp)
            rateExprs.append("    r['{}'] = {}".format(rId, newExp))

        balances = []
        for m_id in self.metabolites:
            exp = self.parsedXdot[m_id]
            if m_id in factorsEnz.keys():
                if factors[m_id] == 0:
                    newExp = "0"
                else:
                    newExp = re.sub(r"\+\s*(\d)", r"+ \1 * " + str(factorsEnz[m_id]), exp)
                balances.append(' ' * 8 + newExp)
            else:
                balances.append(' ' * 8 + exp)

        func_str = 'def ode_func(t, x, r, p, v):\n\n' + \
                   '\n'.join(ruleExprs) + '\n\n' + \
                   '\n'.join(rateExprs) + '\n\n' + \
                   '    dxdt = [\n' + \
                   ',\n'.join(balances) + '\n' + \
                   '    ]\n\n' + \
                   '    return dxdt\n'
        return func_str


    def get_ode(self, r_dict=None, params=None, factors=None):
        """
        Build the ODE system.
        
        :param r_dict: (dict) This variable is used to store the reaction rates.
        :param params: (dict) Parameters and the new values used to replace the original parameters present in the SBML model.
        :param factors: (dict) The key is the parameter identifier and the value is the level of change values between 0 and 1 represent a under expression, above 1 a over expression and 0 to represent the knockouts.
        :returns: A function used to solve the ODE system.
            
        """

        p = self.merge_constants()
        v = self.variable_params.copy()

        if r_dict is not None:
            r = r_dict
        else:
            r = {}

        if params:
            p.update(params)

        exec(self.build_ode(factors), globals())
        ode_func = eval('ode_func')
        #print "Parameters"
        #print p
        #print v
        #print "-----"
        #print (self.build_ode(factors))
        f = lambda t, x: ode_func(t, x, r, p, v)
        return f

    def set_reactions_parameters_factors(self, map):
        """
        Set a new map with the parameters that can be changed for each reaction.
        
        :param map: (dict) The keys is the reaction identifier and the value a list of parameters which can be used to simulate modifications( KO, under/ over expression)
        
        """

        self.reacParamsFactors = OrderedDict(map) if map else OrderedDict()


    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)




# auxiliar functions to set the assignment rules by the correct order in the ODE system
def _build_tree_rules(parent, rules):
    regexp = "v\[\'(.*?)\'\]"
    children = re.findall(regexp, rules[parent])
    if len(children) == 0:
        return Tree(parent, None)
    else:
        childrenTrees = [_build_tree_rules(child, rules) for child in children]
        return Tree(parent, childrenTrees)


def _get_oder_rules(trees):
    res = []
    for tree in trees:
        new_elems = _get_order_nodes(tree)
        [res.append(item) for item in new_elems if item not in res]
    return res


def _get_order_nodes(tree):
    res = [tree.name]
    if len(tree.children) > 0:
        for child in tree.children:
            res = _get_order_nodes(child) + res
    return res

class Tree(object):
    "Generic tree node."
    def __init__(self, name='root', children=None):
        self.name = name
        self.children = []
        if children is not None:
            for child in children:
                self.add_child(child)

    def add_child(self, node):
       # assert isinstance(node, MyTree)
        self.children.append(node)

def get_order_nodes(tree):
    if tree.children is None:
        return [tree.name]
    else:
        res = []
        for child in tree.children:
            res = res + get_order_nodes(child)
        return res