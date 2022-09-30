from mewpy.problems import GOUProblem

class HybridGOUProblem(GOUProblem):

    def __init__(self, model, hconstraints, fevaluation=None,**kwargs):
        """ Overrides GOUProblem by applying constraints resulting from
        sampling a kinetic model.
        """
        super().__init__(model, fevaluation, **kwargs)
        self.hconstraints=hconstraints

    def solution_to_constraints(self, candidate):
        constraints = super().solution_to_constraints(candidate)
        # apply the hybrid contraints:
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

