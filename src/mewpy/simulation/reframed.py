"""
Simulation for REFRAMED models
"""

from reframed.cobra.simulation import FBA, pFBA, MOMA , lMOMA, ROOM
from reframed.solvers import solver_instance
from reframed.core.cbmodel import CBModel
from reframed.solvers.solution import Solution
from reframed.solvers.solution import Status as s_status
from mewpy.model.gecko import GeckoModel
from mewpy.simulation import SimulationMethod, SStatus
from mewpy.simulation.simulation import Simulator, SimulationResult, ModelContainer
from mewpy.utils.constants import ModelConstants
from mewpy.utils.parsing import evaluate_expression_tree
from collections import OrderedDict





class CBModelContainer(ModelContainer):
    def __init__(self, model: CBModel):
        if not isinstance(model,CBModel):
            raise ValueError("The model is not an instance of ReFramed CBModel")
        self.model = model

    @property
    def reactions(self):
        return list(self.model.reactions.keys())

    @property
    def genes(self):
        return list(self.model.genes.keys())
    
    
    @property
    def metabolites(self):
        return list(self.model.metabolites.keys())
    
    
    def get_gpr(self, reaction_id):
        """Returns the gpr rule (str) for a given reaction ID.
        """
        if reaction_id not in self.reactions:
            raise ValueError(f"Reactions {reaction_id} does not exist")
        reaction = self.model.reactions[reaction_id]
        if reaction.gpr: 
            return str(reaction.gpr)
        else:
            return None
   
   
    def get_drains(self):
        return self.model.get_exchange_reactions()
   
        
        
        
    @property
    def medium(self):

        def is_active(rxn):
            """Determine if a boundary reaction permits flux towards creating
            metabolites
            """
            reaction = self.model.reactions[rxn]
            return ((bool(reaction.get_products()) and (reaction.ub > 0)) or
                    (bool(reaction.get_substrates) and (reaction.lb < 0)))


        def get_active_bound(rxn):
            """For an active boundary reaction, return the relevant bound"""
            reaction = self.model.reactions[rxn]
            if reaction.get_substrates():
                return -reaction.lb
            elif reaction.get_products():
                return reaction.ub

        return {rxn: get_active_bound(rxn) for rxn in self.get_drains()
                if is_active(rxn)}




class Simulation(CBModelContainer, Simulator):
    """
        Generic Simulation class for reframed CBModel.
        Defines the simulation conditions.
    """

    def __init__(self, model: CBModel , objective = None, envcond = None, constraints = None,  solver = None, reference = None):
        
        if not isinstance(model,CBModel): 
            raise ValueError("Model is None or is not an instance of REFRAMED CBModel")
        
        self.model = model
        self.objective = self.model.get_objective() if objective is None else objective
        self.environmental_conditions = OrderedDict() if envcond is None else envcond
        self.constraints = OrderedDict() if constraints is None else constraints
        self.solver = solver
        self._essential_reactions = None 
        self._essential_genes = None 
        self._reference = reference
        self._gene_to_reaction = None

        self.__status_mapping = {
            s_status.OPTIMAL: SStatus.OPTIMAL,
            s_status.UNBOUNDED: SStatus.UNBOUNDED,
            s_status.INFEASIBLE: SStatus.INFEASIBLE,
            s_status.INF_OR_UNB: SStatus.INF_OR_UNB,
            s_status.UNKNOWN: SStatus.UNKNOWN,
            s_status.SUBOPTIMAL: SStatus.SUBOPTIMAL
        }

        self.solver = solver
        self._reset_solver = ModelConstants.RESET_SOLVER
        self.reverse_sintax = [('_b','_f')]
        

    

    @property
    def reference(self):
        if self._reference is None:
            self._reference = self.simulate(method = SimulationMethod.pFBA).fluxes
        return self._reference



    @property
    def essential_reactions(self, min_growth = 0.01):
        if self._essential_reactions is not None:
            return self._essential_reactions
        wt_solution = self.simulate()
        wt_growth = wt_solution.objective_value
        reactions = self.model.reactions.keys()
        self._essential_reactions =[]
        for rxn in reactions:
            res = self.simulate(constraints={rxn : 0})
            if res: 
                if (res.status == SStatus.OPTIMAL and res.objective_value < wt_growth * min_growth) or res.status == SStatus.INFEASIBLE:
                    self._essential_reactions.append(rxn)
        return self._essential_reactions





    @property
    def essential_genes(self, min_growth = 0.01):
        if self._essential_genes is not None:
            return self._essential_genes
        self._essential_genes = []
        wt_solution = self.simulate()
        wt_growth = wt_solution.objective_value
        genes = self.model.genes
        for gene in genes:
            active_genes = set(self.model.genes) - set([gene])
            active_reactions = self.evaluate_gprs(active_genes)
            inactive_reactions = set(self.model.reactions) - set(active_reactions)
            gr_constraints = { rxn: 0 for rxn in inactive_reactions}
            res = self.simulate(constraints = gr_constraints)
            if res: 
                if (res.status == SStatus.OPTIMAL and res.objective_value < wt_growth * min_growth) or res.status == SStatus.INFEASIBLE:
                    self._essential_genes.append(gene)
        return self._essential_genes



    def evaluate_gprs(self,active_genes):
        """Returns the list of active reactions for a given list of active genes.
        """
        active_reactions =[]
        reactions = self.model.reactions
        for r_id, reaction in reactions.items():
            if reaction.gpr: 
                if evaluate_expression_tree(str(reaction.gpr),active_genes):
                    active_reactions.append(r_id)
            else:
                active_reactions.append(r_id)
        return active_reactions


        
        

    


    def get_uptake_reactions(self):
        
        drains = self.get_drains()
        reacs = [r for r in drains if self.model.reactions[r].reversible or
                    ((self.model.reactions[r].lb is None or self.model.reactions[r].lb<0 )and len(self.model.reactions[r].get_substrates())>0) or
                    ((self.model.reactions[r].ub is None or self.model.reactions[r].ub>0 )and len(self.model.reactions[r].get_products()))>0]
        return reacs





    def reverse_reaction(self,reaction_id):

        """
        Identify if a reaction is reversible and returns the 
        reverse reaction if it is the case.

        Returns 
            reaction identifier or None
        TODO: ... use regex instead  
        """
        rxn = self.model.reactions[reaction_id]
        reactions = self.model.reactions
        if rxn.lb < 0:
            rxn.reversible = True
            return reaction_id
        # The model might have been converted to irreversible by REFRAMED in which case reversible reactions
        # are decoupled into forward (reaction_id+'_f') and backward (reaction_id+'_b') reactions
        # or migth be using some other identifier which must be included in self.reverse_sufix
        else:
            for a,b in self.reverse_sintax:
                n = len(reaction_id)-len(a)
                m = len(reaction_id)-len(b)
                if reaction_id[n:] == a and reactions[reaction_id[:n]+b]:
                    return reaction_id[:n]+b
                elif reaction_id[m:] == b and reactions[reaction_id[:m]+a]:
                    return reaction_id[:m]+a
                else:
                    continue
            return None
    
    
    
    def gene_reactions(self):
        """
        returns a map of genes to reactions
        """
        if not self._gene_to_reaction:
            gr = OrderedDict()
            for rxn_id in self.reactions:
                rxn = self.model.reactions[rxn_id]
                if rxn.gpr:
                    genes = rxn.gpr.get_genes()
                    for g in genes:
                        if g in gr.keys():
                            gr[g].append(rxn_id)
                        else:
                            gr[g] = [rxn_id]
            self._gene_to_reaction = gr
        return self._gene_to_reaction
    
    



    


     


    



    def simulate(self, objective = None , method = SimulationMethod.FBA , maximize = True, constraints = None, reference = None , scalefactor = None, solver = None):
        '''
            Simulates the application of constraints using the specified method.

            arguments:
            *objective* (dic): the simulation objective. If none, the model objective is used.
            *method* (SimulationMethod):
            *maximize* (boolean) : the optimization direction
            *contraints* (dic): contraints to be applied to the model.
            *reference* 
            *scalefactor* (float) : a positive scaling factor. Default None 
        '''
        
        a_solver = solver

        if not objective:
            objective = self.objective

        simul_constraints = OrderedDict()
        if constraints:
            simul_constraints.update(constraints)
        if self.constraints: 
            simul_constraints.update(self.constraints)
        if self.environmental_conditions:
            simul_constraints.update(self.environmental_conditions) 
        
        if not a_solver and not self._reset_solver:
            if self.solver is None:
                self.solver = solver_instance(self.model)
            a_solver = self.solver

        # scales the model if a scalling factor is defined. 
        # ... for now resets the solver...
        # ... scalling should be impelemented at the solver level...
        # ... messing with the model is not a good idea, maybe using a clone instead.
        if scalefactor:
            a_solver = None
            for _, rxn in self.model.reactions.items():
                rxn.lb = rxn.lb * scalefactor
                rxn.ub = rxn.ub * scalefactor
            if simul_constraints:    
                for idx, constraint in simul_constraints.items():
                    if  isinstance(constraint, (int, float)):
                        simul_constraints[idx] = constraint * scalefactor
                    elif isinstance(constraint,tuple):
                        simul_constraints[idx] = tuple( x * scalefactor for x in constraint)
                    else:
                        raise ValueError("Could not scale the model")

        # TODO: simplifly ...
        if method in [SimulationMethod.lMOMA,SimulationMethod.MOMA,SimulationMethod.ROOM] and reference is None:
            reference = self.reference
        
        if method == SimulationMethod.FBA:
            solution = FBA(self.model,  objective= objective, minimize= not maximize, constraints=simul_constraints, solver=a_solver)
        elif method == SimulationMethod.pFBA:
            solution = pFBA(self.model, objective= objective, minimize= not maximize, constraints=simul_constraints,solver=a_solver,obj_frac=0.999)
        elif method == SimulationMethod.lMOMA:
            solution = lMOMA(self.model, constraints=simul_constraints,reference=reference,solver=a_solver)
        elif method == SimulationMethod.MOMA:
            solution = MOMA(self.model,  constraints=simul_constraints,reference=reference,solver=a_solver)
        elif method == SimulationMethod.ROOM:
            solution = ROOM(self.model,  constraints=simul_constraints,reference=reference,solver=a_solver)
        # Special case in which only the simulation context is required without any simulatin result 
        elif method == SimulationMethod.NONE: 
            solution = Solution(status=s_status.UNKNOWN, message=None, fobj=None, values=None)    
        else:
            raise Exception(
                "Unknown method to perform the simulation.")
        
        # undoes the model scalling 
        if scalefactor:
            for _, rxn in self.model.reactions.items():
                rxn.lb = rxn.lb / scalefactor
                rxn.ub = rxn.ub / scalefactor
            if solution.status in (s_status.OPTIMAL, s_status.SUBOPTIMAL):
                solution.fobj =  solution.fobj / scalefactor
                for x, y in solution.values.items():
                    solution.values[x] = y / scalefactor

        status = self.__status_mapping[solution.status]
        
        result = SimulationResult(self.model, solution.fobj , fluxes= solution.values, status= status, 
                                 envcond= self.environmental_conditions, model_constraints= self.constraints , 
                                 simul_constraints= constraints, maximize= maximize)
        return result






class GeckoSimulation(Simulation):

    def __init__(self, model : GeckoModel, objective = None, envcond = None, constraints = None,  solver = None, reference = None):
       super(GeckoSimulation,self).__init__(model,objective,envcond,constraints,solver,reference)        
       self._essential_proteins = None
    
    
    @property
    def proteins(self):
        return self.model.proteins
    
    
    @property
    def protein_rev_reactions(self):
        return self.model.protein_rev_reactions
    

    def adjust_pool_bounds(self, min_objective=0.05, inplace=False, tolerance=1e-9):
 
        """Adjust protein pool bounds minimally to make model feasible.

        Bounds from measurements can make the model non-viable or even infeasible. Adjust these minimally by minimizing
        the positive deviation from the measured values.

        Parameters
        ----------
        min_objective : float
            The minimum value of for the ojective for calling the model viable.
        inplace : bool
            Apply the adjustments to the model.
        tolerance : float
            Minimum non-zero value. Solver specific value.

        Returns
        -------
        pd.DataFrame
            Data frame with the series 'original' bounds and the new 'adjusted' bound, and the optimized 'addition'.

        """

        solver = solver_instance(self.model)
        solver.add_constraint('constraint_objective', self.model.get_objective, sense='>', rhs=min_objective)
        for pool in self.model.individual_protein_exchanges:
            solver.add_variable('pool_diff_' + pool.id,lb=0)
            solver.add_variable('measured_bound_' + pool.id, lb=pool.upper_bound, ub=pool.upper_bound)
        
        """
        with self.model as model:
            problem = model.problem
            constraint_objective = problem.Constraint(model.objective.expression, name='constraint_objective',
                                                      lb=min_objective)
            to_add = [constraint_objective]
            new_objective = S.Zero
            for pool in model.individual_protein_exchanges:
                ub_diff = problem.Variable('pool_diff_' + pool.id, lb=0, ub=None)
                current_ub = problem.Variable('measured_bound_' + pool.id, lb=pool.upper_bound, ub=pool.upper_bound)
                constraint = problem.Constraint(pool.forward_variable - current_ub - ub_diff, ub=0,
                                                name='pool_ub_' + pool.id)
                to_add.extend([ub_diff, current_ub, constraint])
                new_objective += ub_diff
                pool.bounds = 0, 1000.
            model.add_cons_vars(to_add)
            model.objective = problem.Objective(new_objective, direction='min')
            model.slim_optimize(error_value=None)
            primal_values = model.solver.primal_values
        adjustments = [(pool.id, primal_values['pool_diff_' + pool.id], pool.upper_bound)
                       for pool in model.individual_protein_exchanges
                       if primal_values['pool_diff_' + pool.id] > tolerance]
        result = pd.DataFrame(adjustments, columns=['reaction', 'addition', 'original'])
        result['adjusted'] = result['addition'] + result['original']
        if inplace:
            for adj in result.itertuples():
                model.reactions.get_by_id(adj.reaction).upper_bound = adj.adjusted
        return result
        """



    
    @property
    def essential_proteins(self, min_growth = 0.01):
        if self._essential_proteins is not None:
            return self._essential_proteins
        wt_solution = self.simulate()
        wt_growth = wt_solution.objective_value
        self._essential_proteins =[]
        proteins = self.model.proteins
        for p in proteins:
            rxn = "draw_prot_{}".format(p)    
            res = self.simulate( constraints = { rxn : 0} )
            if res: 
                if (res.status == SStatus.OPTIMAL and res.objective_value < wt_growth * min_growth) or res.status == SStatus.INFEASIBLE:
                    self._essential_proteins.append(rxn)
        return self._essential_proteins
   


    def protein_reactions(self,protein):
        """
        Returns the list of reactions associated to a protein
        """
        reactions =[]
        for r_id, rxn in self.model.reactions.items():
                    lsub = rxn.get_substrates()
                    for m in lsub:
                        if protein in m:
                            reactions.append(r_id)
        return reactions



    def reverse_reaction(self, reaction_id):
        """
        Identify if a reaction is reversible and returns the 
        reverse reaction if it is the case

        Returns 
            reaction identifier or None
        """
        f,d = zip(*self.model.protein_rev_reactions.values())
        if reaction_id in f:
            return d[f.index(reaction_id)]
        elif reaction_id in d:
            return f[d.index(reaction_id)]
        else:
            return None
        
    @staticmethod    
    def __call__(model):
        return GeckoModel(model)



if __name__ == "__main__":
    from mewpy.model.gecko import GeckoModel
    model = GeckoModel('single-pool')
    simul = GeckoSimulation(model)
    ref = simul.reference