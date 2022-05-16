import warnings
from collections import OrderedDict
import numpy as np
from reframed.core.model import Model
from mewpy.util.parsing import Arithmetic, build_tree
import copy


def calculate_yprime(y, rate, substrates, products, substrate_names):
    """
    This function is used by the rate classes the user creates.
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

    def __init__(self, r_id, formula: str, substrates: list = [],
                 products: list = [], parameters: dict = {}, modifiers: list = []):
        """[summary]

        Args:
            r_id (str): Reaction identifier
            formula (str): kinetic law
            parameters (dict, optional): local parameters. Defaults to dict().
            substrates (list, optional): substrates. Defaults to [].
            products (list, optional): products. Defaults to [].
        """
        super(KineticReaction, self).__init__(r_id, formula, parameters)
        self.substrates = substrates
        self.products = products
        self.modifiers = modifiers
        self.parameter_distributions = {}

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


class ODEModel(Model):
    # TODO: Generalize to work with any model
    def __init__(self, model_id):
        """ ODE Model.
        """
        super(ODEModel, self).__init__(model_id)
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

    def get_reactions(self):
        return self.reactions

    def _clear_temp(self):
        self.update()
        self._func_str = None

    def add_reaction(self, reaction, replace=True):
        super(ODEModel, self).add_reaction(reaction, replace)

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
        if r_id in self.reactions:
            self.ratelaws[r_id] = law
        else:
            warnings.warn(f"No such reaction '{r_id}'", RuntimeWarning)

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

    def remove_reactions(self, id_list):
        super(ODEModel, self).remove_reactions(id_list)
        for r_id in id_list:
            del self.ratelaws[r_id]

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

        if factors:
            for k, v in factors.items():
                exp = rmap[k]
                rmap[k] = f"{str(v)} * {exp}" if isinstance(exp, str) else v*exp

        parsed_rates = {r_id: ratelaw.parse_law(rmap, local=local)
                        for r_id, ratelaw in self.ratelaws.items()}

        r = {r_id: f"({parsed_rates[r_id]})" for r_id in self.reactions}

        rmap.update(r)

        rate_exprs = [' '*4+"r['{}'] = {}".format(r_id, parsed_rates[r_id])
                      for r_id in self.reactions]

        balances = [' '*8 + self.print_balance(m_id,     factors=factors) for m_id in self.metabolites]

        func_str = 'def ode_func(t, x, r, p, v):\n\n' + \
            '\n'.join(rate_exprs) + '\n\n' + \
            '    dxdt = [\n' + \
            ',\n'.join(balances) + '\n' + \
            '    ]\n\n' + \
            '    return dxdt\n'

        self._func_str = func_str

        return self._func_str

    def get_ode(self, r_dict=None, params=None, factors=None):

        p = self.merge_constants()
        v = self.variable_params.copy()

        if r_dict is not None:
            r = r_dict
        else:
            r = {}

        if params:
            p.update(params)

        exec('from math import log', globals())
        exec(self.build_ode(factors), globals())
        ode_func = eval('ode_func')

        return lambda t, x: ode_func(t, x, r, p, v)
