import copy

from reframed.core.cbmodel import CBModel


class SMomentModel(CBModel):

    def __init__(self, model, enzyme_reaction_prefix="R_ENZYME_DELIVERY_"):
        """SMomentModel an extension of REFRAMED CBModel for sMOMENT models.

        :param model: a path to a CBM xml file or an instance of REFRAMED SBModel
        :param enzyme_reaction_prefix: str enzyme prefix
                   GECKO_ANALOGON: "R_ENZYME_DELIVERY_"
                   GECKO: "R__TG_ER_"

        """
        if isinstance(model, str):
            from reframed.io.sbml import load_cbmodel
            model = load_cbmodel(model)
        elif not isinstance(model, CBModel):
            raise ValueError("The model should be a path or a CBModel")

        super(SMomentModel, self).__init__(model.id)
        self.enzyme_prefix = enzyme_reaction_prefix
        # import CBModel's data
        self.compartments = copy.deepcopy(model.compartments)
        self.metabolites = copy.deepcopy(model.metabolites)
        self.reactions = copy.deepcopy(model.reactions)
        self.genes = copy.deepcopy(model.genes)
        self.enzymes = []
        for rx in self.reactions:
            if rx.startswith(self.enzyme_prefix):
                self.enzymes.append(rx[len(self.enzyme_prefix):])

    @property
    def protein_rev_reactions(self):
        return {}

    @property
    def proteins(self):
        return self.enzymes

    def set_pool_bounds(self, lb, up):
        pass
