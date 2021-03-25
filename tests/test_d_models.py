import unittest

MODELS_PATH = 'tests/data/'
EC_CORE_MODEL = MODELS_PATH + 'ecoli_core_model.xml'
EC_CORE_REG_MODEL = MODELS_PATH + 'e_coli_core_trn.csv'


class TestMewModel(unittest.TestCase):
    """ Tests a MewModel
    """

    def setUp(self):
        """Set up
        Loads a model
        """

        from mewpy.io import Reader, Engines

        self.metabolic_reader = Reader(Engines.MetabolicSBML, EC_CORE_MODEL)
        self.regulatory_reader = Reader(Engines.RegulatoryCSV,
                                        EC_CORE_REG_MODEL,
                                        sep=',',
                                        id_col=1,
                                        rule_col=2,
                                        aliases_cols=[0],
                                        header=0)

    def test_algebra(self):
        """
        Tests algebra expressions for models and variables
        """
        rule = '((NOT (o2(e)>0)) AND (NOT (no3(e)>0)) AND (NOT (no2(e)>0)) AND (NOT (4tmao(e)>0)) AND (NOT (' \
               'dmso(e)>0)) AND (for(e)>0)) OR (b0001 == 1)'

        from mewpy.algebra import Expression, parse_expression, Symbolic, And, Or
        from mewpy.variables import Regulator

        # parsing
        symbolic = parse_expression(rule)

        # symbols to variables
        variables = {symbol.name: Regulator(symbol.name) for symbol in symbolic.atoms(symbols_only=True)}

        # join symbols, expression and variables
        expr = Expression(symbolic, variables)

        # regular evaluation
        state = {identifier: 1 if identifier == 'b0001' else 0
                 for identifier in expr.symbols}
        eval_result = expr.evaluate(values=state)
        self.assertEqual(eval_result, 1)

        for symbolic in expr.walk(reverse=True):
            self.assertTrue(isinstance(symbolic, Symbolic))

        for symbolic in expr.walk():
            self.assertTrue(isinstance(symbolic, Symbolic))

        active_st = expr.truth_table(active_states=True)
        full = expr.truth_table(active_states=False)

        variables.get('b0001').active_coefficient = 1

        active_st = expr.truth_table(active_states=True)
        full = expr.truth_table(active_states=False)

        rule = 'A and (B or (C and D) or F) and G'

        symbolic = parse_expression(rule)
        variables = {symbol.name: Regulator(symbol.name) for symbol in symbolic.atoms(symbols_only=True)}
        expr = Expression(symbolic, variables)

        # custom evaluation
        values = {'A': 100,
                  'B': 80,
                  'C': 90,
                  'D': 95,
                  'F': 93,
                  'G': 300}
        operators = {And: min, Or: max}

        self.assertEqual(expr.evaluate(values=values, operators=operators), 93)

    def test_model(self):
        """
        Tests model and variables object
        """
        from mewpy.variables import Gene, Metabolite, Reaction, Target, Interaction, Regulator

        # integrated model
        gene = Gene(2)
        metabolite = Metabolite(2)
        reaction = Reaction(2)
        target = Target(2)
        interaction = Interaction(2)
        regulator = Regulator(2)

        from mewpy.model import Model
        integrated_model = Model.from_types(('metabolic', 'regulatory'),
                                            identifier='IntegratedModel',
                                            genes={'1': gene},
                                            metabolites={'1': metabolite},
                                            reactions={'1': reaction},
                                            targets={'1': target},
                                            interactions={'1': interaction},
                                            regulators={'1': regulator})

        self.assertEqual(integrated_model.id, 'IntegratedModel')
        self.assertEqual(integrated_model.types, {'regulatory', 'metabolic'})
        self.assertEqual(integrated_model.genes, {2: gene})
        self.assertEqual(integrated_model.metabolites, {2: metabolite})
        self.assertEqual(integrated_model.reactions, {2: reaction})
        self.assertEqual(integrated_model.targets, {2: target})
        self.assertEqual(integrated_model.interactions, {2: interaction})

        # metabolic model
        from mewpy.model import MetabolicModel
        model = MetabolicModel('MetabolicModel')

        metabolite1 = Metabolite('o2')
        metabolite2 = Metabolite('h2o2')
        gene1 = Gene('b0001')

        from mewpy.algebra import parse_expression, Expression
        expr1 = parse_expression('b0001')

        gpr1 = Expression(expr1, {gene1.id: gene1})

        reaction1 = Reaction(identifier='R0001',
                             stoichiometry={metabolite1: -1},
                             bounds=(0.0, 999999),
                             gpr=gpr1)

        rule = 'b0002 and (b0003 or b0004)'

        reaction2 = Reaction.from_gpr_string(rule,
                                             identifier='R0002',
                                             bounds=(0.0, 999999),
                                             stoichiometry={metabolite1: -1, metabolite2: 1})

        model.add([reaction1, reaction2])

        rxns = {'R0001': reaction1, 'R0002': reaction2}
        mets = {'o2': metabolite1, 'h2o2': metabolite2}
        genes = {**{gene1.id: gene1}, **model.reactions['R0002'].genes}
        self.assertEqual(model.reactions, rxns)
        self.assertEqual(model.metabolites, mets)
        self.assertEqual(model.genes, genes)

        self.assertEqual(len(model.reactions['R0001'].metabolites), 1)
        self.assertEqual(len(model.reactions['R0002'].metabolites), 2)
        self.assertEqual(len(model.reactions['R0001'].genes), 1)
        self.assertEqual(len(model.reactions['R0002'].genes), 3)

        # regulatory model
        from mewpy.model import RegulatoryModel
        model = RegulatoryModel('RegulatoryModel')

        target = Target('b0001')
        regulator = Regulator.from_types(('regulator', 'target'), identifier='b0002')

        expr1 = parse_expression('b0002')

        reg_event1 = Expression(expr1, {'b0002': regulator})

        interaction1 = Interaction('I_b0001', target=target, regulatory_events={1.0: reg_event1})

        rule = '((NOT (o2(e)>0)) AND (NOT (no3(e)>0)) AND (NOT (no2(e)>0)) AND (NOT (4tmao(e)>0)) AND (NOT (dmso(' \
               'e)>0)) AND (for(e)>0)) '

        interaction2 = Interaction.from_string('I_b0002', stringify_rule=rule, target=regulator)

        model.add([interaction1, interaction2])

        self.assertEqual(model.interactions, {'I_b0001': interaction1, 'I_b0002': interaction2})
        self.assertEqual(model.regulators, {**{'b0002': regulator}, **model.interactions['I_b0002'].regulators})
        self.assertEqual(model.targets, {'b0001': target, 'b0002': regulator})

    def test_read(self):
        """
        Tests read model
        """

        from mewpy.io import read_model, Reader, Engines
        # from sbml
        model = read_model(self.metabolic_reader)

        self.assertEqual(len(model.reactions), 95)
        self.assertEqual(len(model.metabolites), 72)
        self.assertEqual(len(model.genes), 137)

        regulatory_reader = Reader(Engines.RegulatorySBML, MODELS_PATH + 'e_coli_lac.xml')
        model = read_model(regulatory_reader)

        self.assertEqual(len(model.interactions), 27)
        self.assertEqual(len(model.targets), 27)
        self.assertEqual(len(model.regulators), 20)

        # from csv
        model = read_model(self.regulatory_reader)

        self.assertEqual(len(model.interactions), 159)
        self.assertEqual(len(model.targets), 159)
        self.assertEqual(len(model.regulators), 45)

        # from sbml + sbml
        regulatory_reader = Reader(Engines.RegulatorySBML, MODELS_PATH + 'e_coli_lac.xml')
        model = read_model(regulatory_reader, self.metabolic_reader)

        self.assertEqual(len(model.interactions), 27)
        self.assertEqual(len(model.targets), 27)
        self.assertEqual(len(model.regulators), 20)
        self.assertEqual(len(model.reactions), 95)
        self.assertEqual(len(model.metabolites), 72)
        self.assertEqual(len(model.genes), 137)

        # from sbml + csv
        model = read_model(self.regulatory_reader, self.metabolic_reader)

        self.assertEqual(len(model.interactions), 159)
        self.assertEqual(len(model.targets), 159)
        self.assertEqual(len(model.regulators), 45)
        self.assertEqual(len(model.reactions), 95)
        self.assertEqual(len(model.metabolites), 72)
        self.assertEqual(len(model.genes), 137)

        # from cobra
        from cobra.io import read_sbml_model
        cobra_ecoli_core_model = read_sbml_model(EC_CORE_MODEL)
        metabolic_reader = Reader(Engines.CobrapyModel, cobra_ecoli_core_model)

        model = read_model(metabolic_reader)
        self.assertEqual(len(model.reactions), 95)
        self.assertEqual(len(model.metabolites), 72)
        self.assertEqual(len(model.genes), 137)

        # from reframed
        from reframed import load_cbmodel
        reframed_ecoli_core_model = load_cbmodel(EC_CORE_MODEL)
        metabolic_reader = Reader(Engines.ReframedModel, reframed_ecoli_core_model)

        model = read_model(metabolic_reader)
        self.assertEqual(len(model.reactions), 95)
        self.assertEqual(len(model.metabolites), 72)
        self.assertEqual(len(model.genes), 137)

        # from json
        model_reader = Reader(Engines.JSON, MODELS_PATH + 'iMC1010.json')

        model = read_model(model_reader)
        self.assertEqual(len(model.interactions), 1010)
        self.assertEqual(len(model.targets), 1010)
        self.assertEqual(len(model.regulators), 232)
        self.assertEqual(len(model.reactions), 1083)
        self.assertEqual(len(model.metabolites), 768)
        self.assertEqual(len(model.genes), 904)

    def test_write(self):
        """
        Tests write model
        """

        from mewpy.io import read_model, Writer, Engines, write_model
        model = read_model(self.regulatory_reader, self.metabolic_reader)

        # to sbml
        metabolic_writer = Writer(Engines.MetabolicSBML,
                                  io=MODELS_PATH + 'e_coli_core_write.xml',
                                  model=model)

        regulatory_writer = Writer(Engines.RegulatorySBML,
                                   io=MODELS_PATH + 'e_coli_lac_write.xml',
                                   model=model)

        write_model(regulatory_writer, metabolic_writer)

        # to json
        model_writer = Writer(Engines.JSON,
                              io='e_coli_core_write.json',
                              model=model)

        write_model(model_writer)

    def test_analysis(self):
        """
        Tests model analysis
        """

        from mewpy.io import read_model

        # metabolic analysis
        model = read_model(self.metabolic_reader)
        model.objective = {'Biomass_Ecoli_core': 1}

        # fba
        from mewpy.analysis import FBA
        simulator = FBA(model)
        sol = simulator.optimize()
        self.assertGreater(sol.objective_value, 0)

        # pfba
        from mewpy.analysis import pFBA
        simulator = pFBA(model)
        sol = simulator.optimize()
        self.assertGreater(sol.objective_value, 0)

        # milpfba
        from mewpy.analysis import milpFBA
        simulator = milpFBA(model)
        sol = simulator.optimize()
        self.assertGreater(sol.objective_value, 0)

        # deletions
        from mewpy.analysis import single_reaction_deletion, single_gene_deletion
        reactions_deletion = single_reaction_deletion(model=model,
                                                      method='pfba',
                                                      reactions=list(model.reactions.keys())[0:10])
        self.assertGreater(len(reactions_deletion), 0)

        genes_deletion = single_gene_deletion(model=model,
                                              method='fba',
                                              genes=list(model.genes.keys())[0:10])
        self.assertGreater(len(genes_deletion), 0)

        from mewpy.analysis import fva
        fva_sol = fva(model=model, method='pfba', fraction=0.9, reactions=list(model.reactions.keys())[0:10])
        self.assertGreater(len(fva_sol), 0)

        # regulatory analysis
        model = read_model(self.regulatory_reader)

        # milpBool
        from mewpy.analysis import milpBool
        simulator = milpBool(model)
        sol = simulator.optimize()
        self.assertEqual(sol.objective_value, 0)

        # SimBool
        from mewpy.analysis import SimBool
        simulator = SimBool(model)
        sol = simulator.optimize()
        self.assertEqual(sol.objective_value, 0)

        # deletions
        from mewpy.analysis import single_regulator_deletion
        sol = single_regulator_deletion(model)
        self.assertGreater(len(sol), 0)

        # truth table/regulatory events
        from mewpy.analysis import regulatory_events
        reg_events = regulatory_events(model=model)
        self.assertGreater(len(reg_events), 0)

        # integrated analysis
        model = read_model(self.regulatory_reader, self.metabolic_reader)

        _BIOMASS_ID = 'Biomass_Ecoli_core'
        _O2 = 'EX_o2_e'
        _GLC = 'EX_glc__D_e'
        _FUM = 'EX_fum_e'

        constraints = {_GLC: (-10.0, 100000.0), _O2: (-30, 100000.0), _FUM: (-10, 100000.0)}

        for rxn, bds in constraints.items():
            model.get(rxn).bounds = bds

        model.objective = {_BIOMASS_ID: 1}

        from mewpy.analysis import SRFBA
        simulator = SRFBA(model)
        sol = simulator.optimize()
        self.assertGreater(sol.objective_value, 0)

        from mewpy.analysis import RFBA
        simulator = RFBA(model)
        sol = simulator.optimize()
        self.assertIsNotNone(sol)

        sol = simulator.optimize(dynamic=True)
        self.assertIsNotNone(sol)

        # ifva
        from mewpy.analysis import ifva
        sol = ifva(model, method='srfba', reactions=list(model.reactions.keys())[0:10])
        self.assertGreater(len(sol), 0)

        from mewpy.analysis import isingle_reaction_deletion
        sol = isingle_reaction_deletion(model, method='srfba', reactions=list(model.reactions.keys())[0:10])
        self.assertGreater(len(sol), 0)

        from mewpy.analysis import isingle_gene_deletion
        sol = isingle_gene_deletion(model, method='srfba', genes=list(model.genes.keys())[0:10])
        self.assertGreater(len(sol), 0)

        from mewpy.analysis import isingle_regulator_deletion
        sol = isingle_regulator_deletion(model, method='srfba', regulators=list(model.regulators.keys())[0:10])
        self.assertGreater(len(sol), 0)

    def test_workflow(self):
        """
        Tests model workflow
        """

        # TODO: missing

        pass


if __name__ == '__main__':
    unittest.main()
