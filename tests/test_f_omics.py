import unittest

MODELS_PATH = 'tests/data/'
EC_CORE_MODEL = MODELS_PATH + 'e_coli_core.xml.gz'


class TestExpressionSet(unittest.TestCase):

    def setUp(self):
        import numpy as np
        self.genes = ["b0356", "b1478", "b3734", "b3731"]
        self.conditions = ["Exp#1", "Exp#2", "Exp#3"]
        self.expression = np.array(([0.17, 0.20, 0.93],
                                    [0.36, 0.83, 0.77],
                                    [0.87, 0.65, 0.07],
                                    [0.55, 0.49, 0.52]))
        from mewpy.omics import ExpressionSet
        self.expr = ExpressionSet(self.genes, self.conditions, self.expression)
        from cobra.io.sbml import read_sbml_model
        model = read_sbml_model(EC_CORE_MODEL)
        from mewpy.simulation import get_simulator
        self.sim = get_simulator(model)

    def test_GIMME(self):
        from mewpy.omics import GIMME
        GIMME(self.sim, self.expr)

    def test_eFlux(self):
        from mewpy.omics import eFlux
        eFlux(self.sim, self.expr)

    def test_iMAT(self):
        from mewpy.omics import iMAT
        iMAT(self.sim, self.expr)
