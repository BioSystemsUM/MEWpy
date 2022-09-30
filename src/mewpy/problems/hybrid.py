from mewpy.problems import GOUProblem

class HybridGOUProblem(GOUProblem):

    def __init__(self, model, hconstraints:dict, fevaluation:list=None,**kwargs):
        """ Overrides GOUProblem by applying constraints resulting from
        sampling a kinetic model.

        :param model: The constraint metabolic model.
        :param list fevaluation: A list of callable EvaluationFunctions.
        :param dict hconstraints: The hybrid constraints definind kinetic model solution space. 
        Optional:

        :param OrderedDict envcond: Environmental conditions.
        :param OrderedDict constraints: Additional constraints to be applied to the model.
        :param int candidate_min_size: The candidate minimum size (Default EAConstants.MIN_SOLUTION_SIZE)
        :param int candidate_max_size: The candidate maximum size (Default EAConstants.MAX_SOLUTION_SIZE)
        :param list target: List of modification target genes.
        :param list non_target: List of non target genes. Not considered if a target list is provided.
        :param float scalefactor: A scaling factor to be used in the LP formulation.
        :param dic reference: Dictionary of flux values to be used in the over/under expression values computation.
        :param tuple operators: (and, or) operations. Default (MIN, MAX).
        :param list levels: Over/under expression levels (Default EAConstants.LEVELS).
        :param boolean twostep: If deletions should be applied before identifiying reference flux values.

        """
        super().__init__(model, fevaluation, **kwargs)
        self.hconstraints=hconstraints

    def solution_to_constraints(self, candidate):
        constraints = super().solution_to_constraints(candidate)
        # Apply the hybrid contraints:
        # Finds the intersection ranges of the kinetic and steady-state 
        # constraints genetic modifications space.
        # If not overlapping, the genetic modifications are not applied. 
        for r,v in  self.hconstraints.items():
            if r in constraints.keys():
                x = constraints[r]
                if isinstance(x,tuple):
                    lb,ub = x
                else:
                    lb=x
                    ub=x  
                hlb,hub = v
                l=max(lb,hlb)
                u=min(ub,hub)
                if l<=u:
                    constraints[r]=(l,u)
                else:
                    constraints[r]=(hlb,hub)
            else:
                constraints[r]=v
        return constraints    

