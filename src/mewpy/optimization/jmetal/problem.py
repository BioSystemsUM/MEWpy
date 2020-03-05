""" JMetal Problems 
"""
from jmetal.core.solution import Solution, IntegerSolution
from jmetal.core.problem import Problem
from inspyred.ec.emo import Pareto 
import mewpy.optimization.problem as prb
from mewpy.simulation.simulation import SimulationMethod
from mewpy.utils.utilities import dominance_test
from collections import OrderedDict
from typing import Tuple, List
import random

#define EA representation for OU
IntTupple = Tuple[int]


class KOSolution(Solution[int]):
    """ Class representing a KO solution """

    def __init__(self, lower_bound: int, upper_bound: int, number_of_variables:int, number_of_objectives: int, number_of_constraints: int = 0):
        super(KOSolution, self).__init__(number_of_variables, number_of_objectives, number_of_constraints)
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound


    def __eq__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return self.variables.sort() == solution.variables.sort()
        return False
    
    
    # JMetal consideres all problems as minimization
    # Based on pareto dominance
     
    def __gt__(self,solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize= False) == 1
        return False
        
    
    def __lt__(self,solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize= False) == -1
        return False
    

    def __ge__(self,solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize= False) != -1
        return False
        
    
    def __le__(self,solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize= False) != 1
        return False
    


    def __copy__(self):
        new_solution = KOSolution(
            self.lower_bound,
            self.upper_bound,
            self.number_of_variables,
            self.number_of_objectives,
            self.number_of_constrains)
        new_solution.objectives = self.objectives[:]
        new_solution.variables = self.variables[:]
        new_solution.constraints = self.constraints[:]
        new_solution.attributes = self.attributes.copy()

        return new_solution

    @property
    def candidate(self):
        """
        Returns an inspyred representation candidate 
        """
        return set(self.variables)


    @property
    def fitness(self):
        """
        returns an inspyred fitness representation, ie, a single value or a Pareto object
        """
        if len(self.objectives)== 1:
            return self.objectives[0]
        else:
            return Pareto(self.objectives)

    def __str__(self):
        return " ".join((self.variables))









class OUSolution(Solution[IntTupple]):
    """
    Class representing a KO solution
    """
    def __init__(self, lower_bound: List[int], upper_bound: List[int], number_of_variables: int, number_of_objectives: int, number_of_constraints: int = 0):
        super(OUSolution, self).__init__(number_of_variables, number_of_objectives, number_of_constraints)
        self.upper_bound = upper_bound
        self.lower_bound = lower_bound


    def __eq__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return self.variables.sort() == solution.variables.sort()
        return False
   

    # JMetal consideres all problems as minimization

    def __gt__(self,solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize= False) == 1
        return False
        
    
    def __lt__(self,solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize= False) == -1
        return False
    

    def __ge__(self,solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize= False) != -1
        return False
        
    
    def __le__(self,solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize= False) != 1
        return False


    def __copy__(self):
        new_solution = OUSolution(
            self.lower_bound,
            self.upper_bound,
            self.number_of_variables,
            self.number_of_objectives,
            self.number_of_constrains
            )
        new_solution.objectives = self.objectives[:]
        new_solution.variables = self.variables[:]
        new_solution.constraints = self.constraints[:]
        new_solution.attributes = self.attributes.copy()

        return new_solution


    @property
    def candidate(self):
        return set(self.variables)


    @property
    def fitness(self):
        if len(self.objectives)== 1:
            return self.objectives[0]
        else:
            return Pareto(self.objectives)


    


class RKOProblem(prb.RKOProblem, Problem[KOSolution]):
    """ Class representing a reaction knockout problem. 
    
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
        super(RKOProblem,self).__init__( model ,fevaluation=fevaluation, **kwargs)
        self.number_of_constraints = 0
        self.obj_directions = []
        self.obj_labels = []
        for i in range(len(fevaluation)):
            self.obj_labels.append(str(fevaluation[i]))
            if fevaluation[i].maximize:
                self.obj_directions.append(self.MAXIMIZE)
            else:      
                self.obj_directions.append(self.MINIMIZE)
        
    def create_solution(self) -> KOSolution:
        # uses the super class generator
        solution = self.generator(random,None)
        
        new_solution = KOSolution(
            self.bounder.lower_bound,
            self.bounder.upper_bound,
            len(solution),
            self.number_of_objectives,
            self.number_of_constraints)

        new_solution.variables = list(solution)
        return new_solution


    def get_constraints(self, solution):
        return self.decode(set(solution.variables))
    

    def evaluate(self, solution: KOSolution) -> KOSolution:
        # solution constraints
        constraints = self.get_constraints(solution) 
        # pre simulation
        simulation_results = OrderedDict()
        for method in self.methods:
                simulation_result = self.simulator.simulate(constraints=constraints, method=method, scalefactor = self.scalefactor)
                simulation_results[method]= simulation_result

        for i in range(len(self.fevaluation)):
            f = self.fevaluation[i](simulation_results,solution.candidate, scalefactor = self.scalefactor)
            # JMetalPy only deals with minimization problems
            if self.obj_directions[i] == self.MAXIMIZE:
                solution.objectives[i] = -1 * f
            else:
                solution.objectives[i] = f
        return solution


    def get_name(self) -> str:
        return "Reaction KO Problem"









class GKOProblem(prb.GKOProblem, Problem[KOSolution]):
    """ Class representing a gene knockout problem. 
    
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
        super(GKOProblem,self).__init__( model ,fevaluation=fevaluation, **kwargs)
        self.number_of_constraints = 0
        self.obj_directions = []
        self.obj_labels = []
        for i in range(len(fevaluation)):
            self.obj_labels.append(str(fevaluation[i]))
            if fevaluation[i].maximize:
                self.obj_directions.append(self.MAXIMIZE)
            else:      
                self.obj_directions.append(self.MINIMIZE)
        
    def create_solution(self) -> KOSolution:
        # uses the super class generator
        solution = self.generator(random,None)
        
        new_solution = KOSolution(
            self.bounder.lower_bound,
            self.bounder.upper_bound,
            len(solution),
            self.number_of_objectives,
            self.number_of_constraints)

        new_solution.variables = list(solution)
        return new_solution


    def get_constraints(self, solution):
        return self.decode(set(solution.variables))
    

    def evaluate(self, solution: KOSolution) -> KOSolution:
         # solution constraints
        constraints = self.get_constraints(solution) 
        # pre simulation
        simulation_results = OrderedDict()
        for method in self.methods:
                simulation_result = self.simulator.simulate(constraints=constraints, method=method, scalefactor = self.scalefactor)
                simulation_results[method]= simulation_result

        for i in range(len(self.fevaluation)):
            f = self.fevaluation[i](simulation_results,solution.candidate, scalefactor = self.scalefactor)
            # JMetalPy only deals with minimization problems
            if self.obj_directions[i] == self.MAXIMIZE:
                solution.objectives[i] = -1 * f
            else:
                solution.objectives[i] = f
        return solution


    def get_name(self) -> str:
        return "Gene KO Problem"










class ROUProblem(prb.ROUProblem,Problem[OUSolution]):
    """ Class representing a reaction over/under expression problem. 
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

    """
    def __init__(self,model, fevaluation = None, **kwargs ):
        super(ROUProblem,self).__init__(model,fevaluation=fevaluation, **kwargs)
        self.number_of_constraints = 0
        directions = []
        labels =[]
        for i in range(len(fevaluation)):
            labels.append(str(fevaluation[i]))
            if fevaluation[i].maximize:
                directions.append(self.MAXIMIZE)
            else:      
                directions.append(self.MINIMIZE)
        self.directions = directions
        self.obj_directions = directions
        self.labels = labels
        self.obj_labels = labels


    def create_solution(self) -> OUSolution:
        # uses the super class generator
        solution = self.generator(random,None)
        
        new_solution = OUSolution(
            self.bounder.lower_bound,
            self.bounder.upper_bound,
            len(solution),
            self.number_of_objectives,
            self.number_of_constraints)

        new_solution.variables = list(solution)

        return new_solution
    


    def get_constraints(self, solution):
        return self.decode(set(solution.variables))
    

    def evaluate(self, solution: KOSolution) -> KOSolution:
         # solution constraints
        constraints = self.get_constraints(solution) 
        # pre simulation
        simulation_results = OrderedDict()
        for method in self.methods:
                simulation_result = self.simulator.simulate(constraints=constraints, method=method, scalefactor = self.scalefactor)
                simulation_results[method]= simulation_result

        for i in range(len(self.fevaluation)):
            f = self.fevaluation[i](simulation_results,solution.candidate, scalefactor = self.scalefactor)
            # JMetalPy only deals with minimization problems
            if self.obj_directions[i] == self.MAXIMIZE:
                solution.objectives[i] = -1 * f
            else:
                solution.objectives[i] = f
        return solution


    def get_name(self) -> str:
        return "Reaction OU Problem"









class GeckoRKOProblem(prb.GeckoRKOProblem, Problem[KOSolution]):
    """ Class representing a gecko reaction knockout problem.
    
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
        super(GeckoRKOProblem,self).__init__( model ,fevaluation=fevaluation, **kwargs)
        self.number_of_constraints = 0
        self.obj_directions = []
        self.obj_labels = []
        for i in range(len(fevaluation)):
            self.obj_labels.append(str(fevaluation[i]))
            if fevaluation[i].maximize:
                self.obj_directions.append(self.MAXIMIZE)
            else:      
                self.obj_directions.append(self.MINIMIZE)
        
    def create_solution(self) -> KOSolution:
        # uses the super class generator
        solution = self.generator(random,None)
        
        new_solution = KOSolution(
            self.bounder.lower_bound,
            self.bounder.upper_bound,
            len(solution),
            self.number_of_objectives,
            self.number_of_constraints)

        new_solution.variables = list(solution)
        return new_solution


    def get_constraints(self, solution):
        return self.decode(set(solution.variables))
    

    def evaluate(self, solution: KOSolution) -> KOSolution:
         # solution constraints
        constraints = self.get_constraints(solution) 
        # pre simulation
        simulation_results = OrderedDict()
        for method in self.methods:
                simulation_result = self.simulator.simulate(constraints=constraints, method=method, scalefactor = self.scalefactor)
                simulation_results[method]= simulation_result

        for i in range(len(self.fevaluation)):
            f = self.fevaluation[i](simulation_results,solution.candidate, scalefactor = self.scalefactor)
            # JMetalPy only deals with minimization problems
            if self.obj_directions[i] == self.MAXIMIZE:
                solution.objectives[i] = -1 * f
            else:
                solution.objectives[i] = f
        return solution


    def get_name(self) -> str:
        return "Gecko Protein KO Problem"





class GeckoROUProblem(prb.GeckoROUProblem, Problem[OUSolution]):
    """ Class representing a Gecko protein O/U problem. 
    
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

    """

    def __init__(self, model, fevaluation = None, **kwargs ):
        super(GeckoROUProblem,self).__init__( model ,fevaluation=fevaluation, **kwargs)
        self.number_of_constraints = 0
        self.obj_directions = []
        self.obj_labels = []
        for i in range(len(fevaluation)):
            self.obj_labels.append(str(fevaluation[i]))
            if fevaluation[i].maximize:
                self.obj_directions.append(self.MAXIMIZE)
            else:      
                self.obj_directions.append(self.MINIMIZE)
        
    def create_solution(self) -> KOSolution:
        # uses the super class generator
        solution = self.generator(random,None)
        
        new_solution = KOSolution(
            self.bounder.lower_bound,
            self.bounder.upper_bound,
            len(solution),
            self.number_of_objectives,
            self.number_of_constraints)

        new_solution.variables = list(solution)
        return new_solution


    def get_constraints(self, solution):
        return self.decode(set(solution.variables))
    

    def evaluate(self, solution: KOSolution) -> KOSolution:
         # solution constraints
        constraints = self.get_constraints(solution) 
        # pre simulation
        simulation_results = OrderedDict()
        for method in self.methods:
                simulation_result = self.simulator.simulate(constraints=constraints, method=method, scalefactor = self.scalefactor)
                simulation_results[method]= simulation_result

        for i in range(len(self.fevaluation)):
            f = self.fevaluation[i](simulation_results,solution.candidate, scalefactor = self.scalefactor)
            # JMetalPy only deals with minimization problems
            if self.obj_directions[i] == self.MAXIMIZE:
                solution.objectives[i] = -1 * f
            else:
                solution.objectives[i] = f
        return solution


    def get_name(self) -> str:
        return "Gecko Protein OU Problem"
