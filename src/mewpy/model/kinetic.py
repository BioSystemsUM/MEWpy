import warnings
from collections import OrderedDict
import numpy as np
import copy
from mewpy.util.parsing import Arithmetic, build_tree

class Compartment(object):
    """ class for modeling compartments. """

    def __init__(self, comp_id, name=None, external=False, size=1.0):
        """
        Arguments:
            comp_id (str): a valid unique identifier
            name (str): compartment name (optional)
            external (bool): is external (default: false)
            size (float): compartment size (default: 1.0)
        """
        self.id = comp_id
        self.name = name if name is not None else comp_id
        self.size = size
        self.external = external
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)


class Metabolite(object):
    """ class for modeling metabolites. """

    def __init__(self, met_id, name=None, compartment=None):
        """
        Arguments:
            met_id (str): a valid unique identifier
            name (str): common metabolite name
            compartment (str): compartment containing the metabolite
        """
        self.id = met_id
        self.name = name if name is not None else met_id
        self.compartment = compartment
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)


def calculate_yprime(y, rate, substrates, products, substrate_names):
    """
    It takes the numpy array for y_prime,
    and adds or subtracts the amount in rate to all the substrates or products listed
    Returns the new y_prime
    Args:
        y: a numpy array for the substrate values, the same order as y
        rate: the rate calculated by the user made rate equation
        substrates: list of substrates for which rate should be subtracted
        products: list of products for which rate should be added
        substrate_names: the ordered list of substrate names in the model.\
            Used to get the position of each substrate or product in y_prime
    Returns:
        y_prime: following the addition or subtraction of rate to the specificed substrates
    """
    y_prime = np.zeros(len(y))

    for name in substrates:
        y_prime[substrate_names.index(name)] -= rate

    for name in products:
        y_prime[substrate_names.index(name)] += rate

    return y_prime


def check_positive(y_prime):
    """
    Check that substrate values are not negative when they shouldnt be.
    """

    for i in range(len(y_prime)):
        if y_prime[i] < 0:
            y_prime[i] = 0

    return y_prime


class Rule(object):
    """Base class for kinetic rules.
    """

    def __init__(self, r_id, formula: str, parameters: dict = {}):
        """Creates a new rule

        Args:
            r_id (str): Reaction/rule identifier
            formula (str): The rule string representation.
        """
        self.id = r_id
        self.formula = formula
        self._tree = None
        self.parameters = parameters

    @property
    def tree(self):
        """Parsing tree of the formula.

        Returns:
            Node: Root node of the parsing tree.
        """
        if not self._tree:
            self._tree = build_tree(self.formula, Arithmetic)
        return self._tree

    def parse_parameters(self):
        """Returns the list of parameters within the rule.

        Returns:
            list: parameters
        """
        return list(self.tree.get_parameters())

    def get_parameters(self):
        return self.parameters

    def replace(self, parameters=None, local=True, infix=True):
        """Replaces parameters with values taken from a dictionary.
        If no parameter are given for replacement, returns the string representation of the rule
        built from the parsing tree.

        Args:
            parameters (dict, optional): Replacement dictionary. Defaults to None.

        Returns:
            str: the kinetic rule.
        """
        param = parameters.copy() if parameters else dict()
        if local:
            param.update(self.parameters)
        t = self.tree.replace(param)
        if infix:
            return t.to_infix()
        else:
            return t

    def calculate_rate(self, substrates={}, parameters={}):
        param = dict()
        param.update(substrates)
        param.update(parameters)
        param.update(self.parameters)
        if len(param.keys()) != len(self.parse_parameters()):
            s = set(self.parse_parameters())-set(param.keys())
            raise ValueError(f"Values missing for parameters: {s}")
        t = self.replace(param)
        rate = eval(t)
        return rate

    def __str__(self):
        return self.formula

    def __repr__(self):
        return self.replace()


class KineticReaction(Rule):

    def __init__(self, r_id, formula: str, stoichiometry: dict = {}, parameters: dict = {}, modifiers: list = []):
        """[summary]

        Args:
            r_id (str): Reaction identifier
            formula (str): kinetic law
            parameters (dict, optional): local parameters. Defaults to dict().
            substrates (list, optional): substrates. Defaults to [].
            products (list, optional): products. Defaults to [].
        """
        super(KineticReaction, self).__init__(r_id, formula, parameters)
        self.stoichiometry = stoichiometry
        self.modifiers = modifiers
        self.parameter_distributions = {}

    @property
    def substrates(self):
        return [k for k, v in self.stoichiometry.items() if v < 0]

    @property
    def products(self):
        return [k for k, v in self.stoichiometry.items() if v > 0]

    def parse_law(self, map: dict, local=True):
        """Auxialiary method invoked by the model to build the ODE system.

        Args:
            map (dict): Dictionary of global paramameters replacements.

        Returns:
            str: kinetic rule.
        """
        m = {p_id: f"p['{self.id}_{p_id}']" for p_id in self.parameters.keys()}
        r_map = map.copy()
        r_map.update(m)

        return self.replace(r_map, local=local)

    def calculate_rate(self, substrates={}, parameters={}):
        param = {}
        param.update(substrates)
        param.update(parameters)
        param.update(self.parameters)
        if len(param.keys()) != len(self.parse_parameters()):
            s = set(self.parse_parameters())-set(param.keys())
            r = s - set(self.parameter_distributions.keys())
            if r:
                raise ValueError(f"Missing values or distribuitions for parameters: {r}")
            else:
                for p in s:
                    param[p] = self.parameter_distributions[p].rvs()
        t = self.replace(param)
        rate = eval(t)
        return rate

    def reaction(self, y, substrate_names, parameter_dict):

        if self.substrate_indexes == []:
            # need to move this to the model
            self.get_indexes(substrate_names)

        if self.run_model_parameters == []:
            self.run_model_parameters = self.get_parameters(parameter_dict)

        for modifier in self.modifiers:
            if modifier.substrate_indexes == []:
                modifier.get_substrate_indexes(self.reaction_substrate_names)
            if modifier.parameter_indexes == []:
                modifier.get_parameter_indexes(self.parameter_names)

        substrates = self.get_substrates(y)
        parameters = copy.copy(self.run_model_parameters)

        if len(self.modifiers) != 0:
            substrates, parameters = self.calculate_modifiers(substrates, parameters)

        rate = self.calculate_rate(substrates, parameters)

        y_prime = calculate_yprime(y, rate, self.substrates, self.products, substrate_names)
        y_prime = self.modify_product(y_prime, substrate_names)

        if self.check_positive:
            y_prime = check_positive(y_prime)

        return y_prime

    def set_parameter_defaults_to_mean(self):
        """Sets not defined parameters to the median of a distribution.
        """
        for name in self.parameter_distributions:
            if name not in self.parameters:
                if (type(self.parameter_distributions[name]) == list or
                        type(self.parameter_distributions[name]) == tuple):
                    self.parameters[name] = (self.parameter_distributions[name][0] +
                                             self.parameter_distributions[name][1]) / 2
                else:
                    self.parameters[name] = self.parameter_distributions[name].mean()


class ODEModel:
    def __init__(self, model_id):
        """ ODE Model.
        """
        self.id = model_id
        self.metabolites = OrderedDict()
        self.compartments = OrderedDict()
        # kinetic rule of each reaction
        self.ratelaws = OrderedDict()
        # initial concentration of metabolites
        self.concentrations = OrderedDict()
        # parameter defined as constantes
        self.constant_params = OrderedDict()
        # variable parameters
        self.variable_params = OrderedDict()
        self.assignment_rules = OrderedDict()

        self._func_str = None
        self._constants = None

        self._m_r_lookup = None

    def _clear_temp(self):
        self.update()
        self._func_str = None

    def add_compartment(self, compartment, replace=True):
        """ Add a compartment to the model.
        Arguments:
            compartment (Compartment): compartment to add
            replace (bool): replace previous compartment with same id (default: True)
        """
        if compartment.id in self.compartments and not replace:
            raise RuntimeError(f"Compartment {compartment.id} already exists.")
        self.compartments[compartment.id] = compartment

    def add_metabolite(self, metabolite, replace=True):
        """ Add a metabolite to the model.
        Arguments:
            metabolite (Metabolite): metabolite to add
            replace (bool): replace previous metabolite with same id (default: True)
        """

        if metabolite.id in self.metabolites and not replace:
            raise RuntimeError(f"Metabolite {metabolite.id} already exists.")

        if metabolite.compartment not in self.compartments:
            raise RuntimeError(f"Metabolite {metabolite.id} has invalid compartment {metabolite.compartment}.")

        self.metabolites[metabolite.id] = metabolite

    def set_concentration(self, m_id: str, concentration: float):
        """Sets a metabolite initial concentration

        Args:
            m_id (str): Metabolite identifier
            concentration (float): Initial concentration
        """
        if m_id in self.metabolites:
            self.concentrations[m_id] = concentration
        else:
            warnings.warn(f"No such metabolite '{m_id}'", RuntimeWarning)

    def set_ratelaw(self, r_id: str, law: KineticReaction):
        """Define the rate law for a given reaction.

        Args:
            r_id (str): Reaction Identifier
            law (KineticReaction): The reaction rate law.
        """
        self.ratelaws[r_id] = law

    def get_ratelaw(self, r_id):
        if r_id in self.ratelaws.keys():
            return self.ratelaws[r_id]
        else:
            raise ValueError('Reaction has no rate law.')

    def set_assignment_rule(self, p_id: str, rule: Rule):
        if p_id in self.variable_params or p_id in self.metabolites:
            self.assignment_rules[p_id] = rule
        else:
            warnings.warn(f"No such variable parameter '{p_id}'", RuntimeWarning)

    def set_global_parameter(self, key, value, constant=True):
        if constant:
            self.constant_params[key] = value
        else:
            self.variable_params[key] = value

    def merge_constants(self):
        constants = OrderedDict()

        for c_id, comp in self.compartments.items():
            constants[c_id] = comp.size

        constants.update(self.constant_params)

        for r_id, law in self.ratelaws.items():
            for p_id, value in law.parameters.items():
                full_id = f"{r_id}_{p_id}"
                constants[full_id] = value

        self._constants = constants
        return constants

    def metabolite_reaction_lookup(self):
        if not self._m_r_lookup:
            self._m_r_lookup = {m_id: {} for m_id in self.metabolites}

            for r_id, rule in self.ratelaws.items():
                for m_id, coeff in rule.stoichiometry.items():
                    self._m_r_lookup[m_id][r_id] = coeff
        return self._m_r_lookup

    def print_balance(self, m_id, factors=None):
        f = factors.get(m_id, 1) if factors else 1
        c_id = self.metabolites[m_id].compartment
        table = self.metabolite_reaction_lookup()

        terms = []
        for r_id, coeff in table[m_id].items():
            v = coeff * f if coeff > 0 else coeff
            terms.append(f"{v:+g} * r['{r_id}']")

        if f == 0 or len(terms) == 0 or (self.metabolites[m_id].constant and self.metabolites[m_id].boundary):
            expr = "0"
        else:
            expr = f"1/p['{c_id}'] * ({' '.join(terms)})"
        return expr

    def get_parameters(self, exclude_compartments=False):
        if not self._constants:
            self.merge_constants()
        parameters = self._constants.copy()
        if exclude_compartments:
            for c_id in self.compartments:
                del parameters[c_id]
        return parameters

    def deriv(self, y, t):
        m_y = OrderedDict(zip(self.metabolites, y))
        yprime = np.zeros(len(y))
        for _, reaction in self.ratelaws.items():
            yprime += reaction.reaction(m_y, self.get_parameters())
        return yprime

    def build_ode(self, factors=None, local=False):

        rmap = OrderedDict()
        m = {m_id: f"x[{i}]" for i, m_id in enumerate(self.metabolites)}
        c = {c_id: f"p['{c_id}']" for c_id in self.compartments}
        p = {p_id: f"p['{p_id}']" for p_id in self.constant_params}
        v = {p_id: f"v['{p_id}']" for p_id in self.variable_params}
        rmap.update(m)
        rmap.update(c)
        rmap.update(p)
        rmap.update(v)
        
        parsed_rates = {r_id: ratelaw.parse_law(rmap, local=local)
                        for r_id, ratelaw in self.ratelaws.items()}

        r = {r_id: f"({parsed_rates[r_id]})" for r_id in self.ratelaws.keys()}

        rmap.update(r)

        rate_exprs = [' '*4+"r['{}'] = {}".format(r_id, parsed_rates[r_id])
                      for r_id in self.ratelaws.keys()]

        balances = [' '*8 + self.print_balance(m_id, factors=factors) for m_id in self.metabolites]

        func_str = 'def ode_func(t, x, r, p, v):\n\n' + \
            '\n'.join(rate_exprs) + '\n\n' + \
            '    dxdt = [\n' + \
            ',\n'.join(balances) + '\n' + \
            '    ]\n\n' + \
            '    return dxdt\n'

        self._func_str = func_str
        return self._func_str

    def get_ode(self, r_dict=None, params=None, factors=None):
        """
        Args:
            r_dict: reaction identifiers to be modified
            params: modified parameters
            factors: factors to be applied to parameters
        """
        p = self.merge_constants()
        if params:
            p.update(params)

        if factors is not None:
            for k,v in factors.items():
                if k in p.keys():
                    p[k]= v * p[k] 

        v = self.variable_params.copy()
        r = r_dict if r_dict is not None else dict()
        
        np.seterr(divide='ignore', invalid='ignore')         
        exec(self.build_ode(factors), globals())
        ode_func = eval('ode_func')

        return lambda t, x: ode_func(t, x, r, p, v)
