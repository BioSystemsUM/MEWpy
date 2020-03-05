from abc import ABC, abstractmethod
from mewpy.model.gecko import GeckoModel
from mewpy.simulation import SimulationMethod, get_simulator
from mewpy.utils.constants import EAConstants, ModelConstants
from mewpy.utils.parsing import Boolean, GeneEvaluator, build_tree
from collections import OrderedDict
from inspyred.ec.emo import Pareto
import numpy as np
from enum import Enum
import functools
import sys
import copy
import dill
import warnings



class Strategy(Enum):
    """
    The available optimization strategies
    """
    KO = 'KO'
    OU = 'OU'



class KOBounder(object):
    """
    A bounder of possible indexes in an enumeration
    """
    def __init__(self, lower_bound, upper_bound):
        
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.range = self.upper_bound - self.lower_bound +1
            
            
    # make it callable
    def __call__(self, candidate, args):
        bounded_candidate = set()
        for val in candidate: 
            if val > self.upper_bound or val < self.lower_bound:
                val =  val % self.range + self.lower_bound
            bounded_candidate.add(val) 
        return bounded_candidate






class OUBounder(object):
    """
    A bounder for (int,int) representations
    """
    def __init__(self, lower_bound, upper_bound):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.range = [ self.upper_bound[i]-self.lower_bound[i]+1 for i in range(len(self.lower_bound)) ]


    def __call__(self, candidate, args):
        bounded_candidate = set()
        for idx, lv in candidate: 
            if idx > self.upper_bound[0] or idx < self.lower_bound[0]:
                idx =  idx % self.range[0] + self.lower_bound[0]
            if lv > self.upper_bound[1] or idx < self.lower_bound[1]:
                lv =  lv % self.range[1] + self.lower_bound[1]
            bounded_candidate.add((idx,lv)) 
        return bounded_candidate



class IntTuppleBounder(object):
    """
    A bounder for (int,int,...) representations
    """
    def __init__(self, lower_bound, upper_bound):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.range = [ self.upper_bound[i]-self.lower_bound[i]+1 for i in range(len(self.lower_bound)) ]


    def __call__(self, candidate, args):
        bounded_candidate = set()
        for c in candidate:
            l = []
            for i in range(len(c)): 
                v =  c[i] % self.range[i] + self.lower_bound[i]
                l.append(v)
            bounded_candidate.add(tuple(l))
        return bounded_candidate




class Problem(ABC):
    """
    Optimization Problem base class
    
    parameters:

    model (metabolic model)
    fevaluation (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness
   
    """
    def __init__(self, model, fevaluation = None, **kwargs):
        self.model = model 
        self.fevaluation =  fevaluation
        self.number_of_objectives = len(self.fevaluation)

        # simulation context : defines the simulations environment
        self.simul_context = None
        self.environmental_conditions = kwargs.get('envcond', None)
        self.persistent_constraints = kwargs.get('constraints', None)
        self._reference = None
        # solution size    
        self.candidate_min_size = kwargs.get('candidate_min_size', EAConstants.MIN_SOLUTION_SIZE) 
        self.candidate_max_size = kwargs.get('candidate_max_size', EAConstants.MAX_SOLUTION_SIZE)
        # non target
        self.non_target = kwargs.get('non_target', None)
        # targets 
        # If not provided, targets are build in the context of the problem. 
        # Objectives are not automatically removed from the targets... this should be a user concern!?  
        self._trg_list = kwargs.get('target', None)
        # the EA representation bounds
        self._bounder = None
        # scaling factor 
        self.scalefactor = kwargs.get('scalefactor', None)
        # required simulations
        methods = []
        for f in self.fevaluation:
            methods.extend(f.required_simulations())
        self.methods = list(set(methods))
        
    
    
    def pre_process(self):
        """ Defines pre processing tasks
        """
        self.target_list
        self.reset_simulator()
        
    
    
    @property
    def simulator(self):
        if self.simul_context is None:
            self.simul_context = get_simulator(self.model, reference= self._reference)
            self.simul_context.environmental_conditions = self.environmental_conditions
            self.simul_context.constraints = self.persistent_constraints
        return self.simul_context

    def reset_simulator(self):
        self.simul_context = None

    def __str__(self):
        if self.number_of_objectives > 1:
            return '{0} ({1}  objectives)'.format(self.__class__.__name__ , self.number_of_objectives)
        else:
            return '{0} '.format(self.__class__.__name__)
        

    def __repr__(self):
        return self.__class__.__name__


    @abstractmethod
    def generator(self, random, args):
        """The generator function for the problem."""
        raise NotImplementedError

    
    @abstractmethod
    def decode(self, candidate):
        """The generator function for the problem."""
        raise NotImplementedError

    
    @property
    def target_list(self):
        "list of allowed OU"
        if self._trg_list is None:
            self._build_target_list()
        return self._trg_list

    @abstractmethod
    def _build_target_list(self):
        raise NotImplementedError
    
    def get_constraints(self, solution):
        """
        returns the constrainst enconded into an individual
        """
        return self.decode(solution.candidate)
    

    def get_environmental_conditions(self):
        return self.environmental_conditions


    def get_persistent_constraints(self):
        return self.persistent_constraints


    def evaluate(self,solution):
        """
        Evaluates a single solution
        """
        p = []
        # decoded constraints
        constraints = self.decode(solution)
        # pre simulation
        simulation_results = OrderedDict()
        try:
            for method in self.methods:
                simulation_result = self.simulator.simulate(constraints=constraints, method=method, scalefactor = self.scalefactor)
                simulation_results[method]= simulation_result      
            # apply the evaluation function(s)
            for f in self.fevaluation:
                p.append(f(simulation_results, solution, scalefactor = self.scalefactor))
        except Exception as e:
            for f in self.fevaluation:
                p.append(f.worst_fitness)
            warnings.warn(f"Solution couldn't be evaluated [{e}]\n {constraints}")
        
        # single objective
        if self.number_of_objectives == 1: 
            return p[0]
        # multi objective    
        else:
            return Pareto(p)



    def evaluator(self, candidates, args):
        """
        Evaluator 
        Note: shoudn't be dependent on args to ease multiprocessing
        returns a list of Pareto fitness values of a candidate list
        """
        fitness = []
        for candidate in candidates:
            p = self.evaluate(candidate)
            fitness.append(p)
        return fitness




    @property
    def is_maximization(self):
        return all([f.maximize for f in self.fevaluation])


    def simplify(self, solution, tolerance = 1e-6):
        """
        Simplify a solution by removing the modification that not affect the final fitness value.
        Args:
            solution : the solution to be simplified
            tolerance: max allowed objective difference values for two solutions to be considered diferent.
                       Tolerance may be defined by a single float value, or per objective by means of a list of floats
                       of size equal to the number of objectives.
        Returns: a list of simplified constraints
        """
        

        constraints = dict(copy.copy(self.get_constraints(solution)))
        constraints_to_remove = []
        if self.number_of_objectives == 1:
            fitness = [solution.fitness]
        else:
            fitness = list(solution.fitness)
        for key in constraints.keys():
            simul_constraints = copy.copy(constraints)
            del simul_constraints[key]
             # pre simulation
            simulation_results = OrderedDict()
            for method in self.methods:
                simulation_result = self.simulator.simulate(constraints=simul_constraints, method=method, scalefactor = self.scalefactor)
                simulation_results[method]= simulation_result      
            # apply the evaluation function(s)
            fit =[]
            for f in self.fevaluation:
                fit.append(f(simulation_results, solution, scalefactor = self.scalefactor))
            diff = np.abs(np.array(fit)-np.array(fitness))
            is_equal = False
            if isinstance(tolerance,float):
                is_equal = np.all(diff <= tolerance) 
            else:
                is_equal = np.all(diff <= np.array(tolerance)) 
            if is_equal:
                constraints_to_remove.append(key)
        for key in constraints_to_remove:
            del constraints[key]
        return constraints
    
    



    





class KOProblem(Problem):
    """
    Base class for Knockout optimization problems

    args:

        model (metabolic model): The constraint based metabolic model.
        fevaluation (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness

    
    **kwargs:

        envcond (OrderedDict): environmental conditions.
        constraints (OrderedDict): additional constraints to be applied to the model.
        candidate_min_size (int) : The candidates minimum size.
        candidate_min_size (int) : The candidates maximum size.
        target (list): List of target reactions.
        non_target (list): List of non target reactions. Not considered if a target list is provided.
        scalefactor (floaf): a scaling factor to be used in the LP formulation. 
    """
    def __init__(self, model, fevaluation = None, **kwargs ):
        super(KOProblem,self).__init__(model, fevaluation = fevaluation, **kwargs )
        self.strategy = Strategy.KO

        


   



    @property
    def bounder(self):
        """
        The KO list index bounder  
        """
        if self._bounder is None:
            max = len(self.target_list)-1
            self._bounder = KOBounder(0,max)
        return self._bounder



   



    def generator(self, random, args):        
        """
        Generates a solution, a random int set with length in range min_solution_size to max_solution_size
        """
        solution = set()
        solution_size= random.uniform(self.candidate_min_size,self.candidate_max_size)
        while len(solution)< solution_size:
            solution.add(random.randint(0,len(self.target_list)-1))
        return solution



    
    def decode(self,candidate):
        """
        Decodes a candidate, an integer set, into a dictionary of constraints
        """
        constraints = OrderedDict()
        for idx in candidate:
            try:
                constraints[self.target_list[idx]] = 0
            except IndexError:
                raise IndexError("Index out of range: {} from {}".format(idx,len(self.target_list[idx])))
        return constraints



    def encode(self,reactions):
        """
        Encodes a list of reaction ids into a set of indexes from the list of allowed KO reactions
        """
        code = set()
        for rxn in reactions:
            code.add(self.target_list.index(rxn))
        return code






class OUProblem(Problem):
    """ Base class for Over/Under expression optimization problems 

    arguments:

        * model* (metabolic model): the constraint metabolic model
        * fevaluation* (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness
    
    **args:

        envcond (OrderedDict): environmental conditions
        constraints (OrderedDict): additional constraints to be applied to the model 
        candidate_min_size : The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
        candidate_max_size : The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
        non_target (list): list of non target reactions
        levels (list): over/under expression levels (Default EAConstants.LEVELS)
        reference (dic): dictionary of reference flux values

    """
    def __init__(self, model, fevaluation = None, **kwargs ):
        super(OUProblem,self).__init__(model, fevaluation = fevaluation, **kwargs )
        self.levels = kwargs.get('levels', EAConstants.LEVELS)
        self._reference = kwargs.get('reference',None) 
        self.strategy = Strategy.OU


   


    @property
    def reference(self):
        if not self._reference:
            self._reference = self.simulator.reference
        return self._reference




    @property
    def bounder(self):
        """
        The list and levels index bounder  
        """
        if self._bounder is None:
            max_idx = len(self.target_list)-1
            max_lv = len(self.levels)-1
            self._bounder = OUBounder([0,0],[max_idx,max_lv])
        return self._bounder



    def generator(self, random, args):        
        """
        Generates a solution, a random (int,int) set with length in range min_solution_size to max_solution_size
        """
        solution = set()
        solution_size= random.uniform(self.candidate_min_size,self.candidate_max_size)
        while len(solution)< solution_size:
            idx = random.randint(0,len(self.target_list)-1)
            lv =  random.randint(0,len(self.levels)-1)
            solution.add((idx,lv))
        return solution


    def ou_constraint(self,level, wt):
            if level> 1:
                return (level*wt, ModelConstants.REACTION_UPPER_BOUND) if wt >= 0 else (-1*ModelConstants.REACTION_UPPER_BOUND,level*wt)
            else:
                return (0,level * wt) if wt >= 0 else (level * wt,0)


    def reaction_constraints(self,rxn,lv):
        constraints = {}
        fluxe_wt = self.reference[rxn]
        rev_rxn = self.simulator.reverse_reaction(rxn)
        if lv == 0: 
            # KO constraint
            constraints[rxn] = (0,0) 
        elif lv == 1: 
            # No contraint is applyed 
            pass
        elif rev_rxn is None or rev_rxn == rxn: 
            # if there is no reverse reaction
            constraints[rxn] = self.ou_constraint(lv,fluxe_wt)
        else: 
            # there's a reverse reaction... 
            # one of the two reactions needs to be KO, the one with no (or lesser) flux in the wt
            rev_fluxe_wt = self.reference[rev_rxn]
            if abs(fluxe_wt) >= abs(rev_fluxe_wt):
                ko_rxn, ou_rxn, fwt = rev_rxn, rxn, fluxe_wt  
            else:
                rxn , rev_rxn, rev_fluxe_wt
            constraints[ko_rxn] = (0,0)
            constraints[ou_rxn] = self.ou_constraint(lv,fwt)
        
        return constraints

    
    
    
    
    def decode(self,candidate):
        """
        Decodes a candidate, an set (idx,lv) into a dictionary of constraints
        Suposes that reverseble reactions have been treated and bounded with positive flux values
        """
        constraints = OrderedDict()
        for idx, lv_idx in candidate:
            try:
                rxn = self.target_list[idx]
                lv = self.levels[lv_idx]
                rev_rxn = self.simulator.reverse_reaction(rxn)
                # skips if the reverse reaction was already processed
                if rev_rxn and rev_rxn in constraints.keys():
                    continue
                elif lv < 0:
                    raise ValueError("All UO levels should be positive") 
                else:
                    constraints.update(self.reaction_constraints(rxn,lv))
            except IndexError:
                raise IndexError("Index out of range")
        return constraints






class RKOProblem(KOProblem):
    """
    Reaction Knockout Optimization Problem.

    args:

        model (metabolic model): The constraint based metabolic model.
        fevaluation (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness
    
    **kwargs:

        envcond (OrderedDict): environmental conditions.
        constraints (OrderedDict): additional constraints to be applied to the model.
        candidate_min_size (int) : The candidates minimum size.
        candidate_min_size (int) : The candidates maximum size.
        target (list): List of target reactions.
        non_target (list): List of non target reactions. Not considered if a target list is provided.
        scalefactor (floaf): a scaling factor to be used in the LP formulation. 
    """
    def __init__(self, model, fevaluation = None, **kwargs ):
        super(RKOProblem,self).__init__(model, fevaluation = fevaluation, **kwargs )
       
        

    def _build_target_list(self):
        
        reactions = set(self.simulator.reactions)
        essential = set(self.simulator.essential_reactions)
        drains = set(self.simulator.get_drains())
        target = reactions - essential - drains
        if self.non_target is not None:
            target = target - set(self.non_target)
        self._trg_list = list(target)
        
        
        




class ROUProblem(OUProblem):
    """
    Reaction Over/Under Expression Optimization Problem

    arguments:

        * model* (metabolic model): the constraint metabolic model
        * fevaluation* (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness
    
    **args:

        *envcond* (OrderedDict): environmental conditions
        *constraints* (OrderedDict): additional constraints to be applied to the model 
        *candidate_min_size* : The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
        *candidate_max_size* : The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
        *non_target* (list): list of non target reactions
        *levels* (list): over/under expression levels (Default EAConstants.LEVELS)
    """
    def __init__(self, model, fevaluation = None, **kwargs ):
        super(ROUProblem,self).__init__(model, fevaluation = fevaluation, **kwargs )
       
        

    def _build_target_list(self):
        reactions = set(self.simulator.reactions)
        essential = set(self.simulator.essential_reactions)
        drains = set(self.simulator.get_drains())
        target = reactions - essential - drains
        if self.non_target is not None:
            target = target - set(self.non_target)
        self._trg_list = list(target)




class GKOProblem(KOProblem):
    """
    Gene Knockout Optimization Problem

    arguments:

        * model* (metabolic model): the constraint metabolic model
        * fevaluation* (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness
    
    **args:

        *envcond* (OrderedDict): environmental conditions
        *constraints* (OrderedDict): additional constraints to be applied to the model 
        *candidate_min_size* : The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
        *candidate_max_size* : The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
        *non_target* (list): list of non target genes
    """
    def __init__(self, model, fevaluation = None, **kwargs ):
        super(GKOProblem,self).__init__(model, fevaluation = fevaluation, **kwargs )
       
        

    def _build_target_list(self):
        
        genes = set(self.simulator.genes)
        essential = set(self.simulator.essential_genes)
        target = genes - essential
        if self.non_target:
            target = target - set(self.non_target)
        self._trg_list = list(target)




    def decode(self,candidate):
        """
        Decodes a candidate, an set of genes into a dictionary of constraints
        """
        genes = [self.target_list[idx] for idx in candidate]
        active_genes = set(self.simulator.genes) - set(genes)
        active_reactions = self.simulator.evaluate_gprs(active_genes)
        inactive_reactions = set(self.simulator.reactions) - set(active_reactions)
        gr_constraints = { rxn: 0 for rxn in inactive_reactions}
        return gr_constraints




class GOUProblem(OUProblem):
    """
    Gene Over/Under expression Optimization Problem

    args:

        model (metabolic model): The constraint based metabolic model.
        fevaluation (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness
   
    **kwargs:

        envcond (OrderedDict): environmental conditions.
        constraints (OrderedDict): additional constraints to be applied to the model.
        candidate_min_size (int) : The candidates minimum size.
        candidate_min_size (int) : The candidates maximum size.
        target (list): List of target reactions.
        non_target (list): List of non target reactions. Not considered if a target list is provided.
        scalefactor (floaf): a scaling factor to be used in the LP formulation.
        reference (dic): Dictionary of flux values to be used in the over/under expression values computation.
        operators (tuple of function): (and, or) operations. Default (MIN, MAX) 

    """
    def __init__(self, model, fevaluation = None, **kwargs ):
        super(GOUProblem,self).__init__(model, fevaluation = fevaluation, **kwargs )
        # operators to replace 'and'/'or'. By default min/max     
        self._operators = kwargs.get('operators', None) 
       
        

    def _build_target_list(self):
       
        genes = set(self.simulator.genes)
        essential = set(self.simulator.essential_genes)
        target = genes - essential
        if self.non_target:
            target = target - set(self.non_target)
        self._trg_list = list(target)





    def decode(self,candidate):
        """
        Decodes a candidate, a set of (genes,lv) into a dictionary of reaction constraints
        """
        gr_constraints  = OrderedDict()
        genes = {self.target_list[idx] : self.levels[lv_idx] for idx,lv_idx in candidate}
        # get all reaction in which the genes intervene
        g_r = self.simulator.gene_reactions()
        reactions =[]
        for g in genes:
            reactions.extend(g_r[g])
        reactions = set(reactions)
        # evaluate gpr
        if not self._operators:
            self._operators = ( lambda x,y: min(x,y) , lambda x,y: max(x,y))
            
        evaluator = GeneEvaluator(genes, self._operators[0],self._operators[1])
        for rxn_id in self.simulator.reactions:
            gpr = self.simulator.get_gpr(rxn_id)
            if gpr :
                tree = build_tree(gpr, Boolean)
                # apply the operators to obtain a level for the reaction
                # if a gene as no level associated its factor is 1 (see GeneEvaluator)
                lv = tree.evaluate( evaluator.f_operand, evaluator.f_operator)
                # adds the reaction constraint
                rev_rxn = self.simulator.reverse_reaction(rxn_id)
                # skips if the reverse reaction was already processed
                if rev_rxn and rev_rxn in gr_constraints.keys():
                    continue
                elif lv < 0:
                    raise ValueError("All UO levels should be positive") 
                else:
                    gr_constraints.update(self.reaction_constraints(rxn_id,lv))
                
        return gr_constraints






class GeckoRKOProblem(KOProblem):
    """
    Gecko KnockOut Optimization Problem 

    args:

        model (metabolic model): The constraint based metabolic model.
        fevaluation (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness

    
    **kwargs:

        envcond (OrderedDict): environmental conditions.
        constraints (OrderedDict): additional constraints to be applied to the model.
        candidate_min_size (int) : The candidates minimum size.
        candidate_min_size (int) : The candidates maximum size.
        target (list): List of target reactions.
        non_target (list): List of non target reactions. Not considered if a target list is provided.
        scalefactor (floaf): a scaling factor to be used in the LP formulation. 
        prot_prefix (str): the protein draw reaction prefix. Default 'draw_prot_'  

    Note:
     Targets as well as non target proteins are defined using their prot id, ex 'P0351', and not by the associated draw reaction id, ex 'draw_prot_P0351'. 
    """
    def __init__(self, model, fevaluation = None, **kwargs ):
        #if isinstance(model,GeckoModel):
        super(GeckoRKOProblem,self).__init__(model, fevaluation = fevaluation, **kwargs )
        #else: 
        #    raise Exception("The model should be an instance of GeckoModel")
        # problem simulation context
        self.prot_prefix = kwargs.get('prot_prefix', 'draw_prot_')


    def _build_target_list(self):
        """
        If not provided, targets are all non essential proteins.
        """
        
        print('building target list')
        
        proteins = set(self.simulator.proteins)
        # as draw_prot_XXXXXX
        e = self.simulator.essential_proteins
        # remove 'draw_prot_'
        n = len(self.prot_prefix)
        essential = set([p[n:] for p in e])
        target = proteins - essential
        if self.non_target:
            target = target - set(self.non_target)
        self._trg_list = list(target)


    def decode(self,candidate):
        """
        Decodes a candidate, an integer set, into a dictionary of constraints
        """
        constraints = OrderedDict()
        for idx in candidate:
            try:
                constraints["{}{}".format(self.prot_prefix,self.target_list[idx])] = 0
            except IndexError:
                raise IndexError(f"Index out of range: {idx} from {len(self.target_list[idx])}")
        return constraints



class GeckoROUProblem(OUProblem):
    """
    Gecko Under/Over expression Optimization Problem 

    args:

        model (metabolic model): The constraint based metabolic model.
        fevaluation (list): a list of callable EvaluationFunctions. If none is given the flux value of the model objective is set as fitness
      
    **kwargs:

        envcond (OrderedDict): environmental conditions.
        constraints (OrderedDict): additional constraints to be applied to the model.
        candidate_min_size (int) : The candidates minimum size.
        candidate_min_size (int) : The candidates maximum size.
        target (list): List of target reactions.
        non_target (list): List of non target reactions. Not considered if a target list is provided.
        scalefactor (floaf): a scaling factor to be used in the LP formulation.
        reference (dic): Dictionary of flux values to be used in the over/under expression values computation. 
        prot_prefix (str): the protein draw reaction prefix. Default 'draw_prot_'  


    Note:
     Target as well as non target proteins are defined with their prot id, ex 'P0351', and with the associated reaction id, ex 'draw_prot_P0351'. 
    """
    def __init__(self, model, fevaluation = None, **kwargs ):
        #if isinstance(model,GeckoModel):
        super(GeckoROUProblem,self).__init__(model, fevaluation = fevaluation, **kwargs )
        #else: 
        #    raise Exception("The model should be an instance of GeckoModel")
        # problem simulation context
        self.prot_rev_reactions = None
        self.prot_prefix = kwargs.get('prot_prefix', 'draw_prot_')
        
    




    def _build_target_list(self):
        """
        If not provided, targets are all non essential proteins.
        """
        proteins = set(self.simulator.proteins)
        # as draw_prot_XXXXXX
        e = self.simulator.essential_proteins
        # remove 'draw_prot_'
        n = len(self.prot_prefix)
        essential = set([p[n:] for p in e])
        target = proteins - essential
        if self.non_target:
            target = target - set(self.non_target)
        self._trg_list = list(target)






    def decode(self,candidate):
        """
        Decodes a candidate, a set (idx,lv), into a dictionary of constraints
        Reverseble reactions associated to proteins with over expression are KO 
        according to the flux volume in the wild type.

        Note: Fluxes in Yeast7 gecko model are always non negative 
        """
        constraints = OrderedDict()
        
        if self.prot_rev_reactions is None:
            self.prot_rev_reactions = self.simulator.protein_rev_reactions
        
        for idx, lv_idx in candidate:
            try:
                prot = self.target_list[idx]
                rxn = self.prot_prefix+prot
                lv = self.levels[lv_idx]
                fluxe_wt = self.reference[rxn]
                if lv < 0:
                    raise ValueError("All UO levels should be positive") 
                # a level = 0 is interpreted as KO 
                elif lv == 0:
                    constraints[rxn] = 0.0
                # under expression    
                elif lv < 1: 
                    constraints[rxn] = (0.0,lv * fluxe_wt)
                # TODO: Define how a level 1 is tranlated into constraints...
                elif lv == 1:
                    continue
                else:
                    constraints[rxn] = (lv*fluxe_wt, ModelConstants.REACTION_UPPER_BOUND)
                    # Deals with reverse reactions associated with the protein. 
                    # Strategy: The reaction direction with no flux in the wild type (reference) is KO.
                    if prot in self.prot_rev_reactions.keys():
                        reactions = self.prot_rev_reactions[prot]
                        for r, r_rev in reactions:
                            # TODO: ... check if a rule can be set when both wt fluxes are 0.0
                            if self.reference[r] == 0 and self.reference[r_rev] == 0:
                                continue
                            elif self.reference[r] > 0 and self.reference[r_rev] == 0:
                                constraints[r_rev] = 0.0
                            elif self.reference[r] == 0 and self.reference[r_rev] > 0:
                                constraints[r] = 0.0
                            else:
                                warnings.warn(f"Reactions {r} and {r_rev}, associated with the protein {prot}, both have fluxes in the WT.")
            except IndexError:
                raise IndexError("Index out of range")

        return constraints




