import unittest
from pathlib import Path
import pytest

MODELS_PATH = Path(__file__).parent.joinpath('data')
EC_CORE_MODEL = MODELS_PATH.joinpath('e_coli_core.xml')
EC_CORE_REG_MODEL = MODELS_PATH.joinpath('e_coli_core_trn.csv')
SAMPLE_MODEL = MODELS_PATH.joinpath('SampleNet.xml')
SAMPLE_REG_MODEL = MODELS_PATH.joinpath('SampleRegNet.csv')


class TestGERMModel(unittest.TestCase):
    """ Tests a GERMModel
    """

    def setUp(self):
        """Set up
        Loads a model
        """

        from mewpy.io import Reader, Engines

        self.metabolic_reader = Reader(Engines.MetabolicSBML, EC_CORE_MODEL)
        self.regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                                        EC_CORE_REG_MODEL,
                                        sep=',',
                                        id_col=0,
                                        rule_col=2,
                                        aliases_cols=[1],
                                        header=0)

    def test_algebra(self):
        """
        Tests algebra expressions for models and variables
        """
        rule = '((NOT (o2(e)>0)) AND (NOT (no3(e)>0)) AND (NOT (no2(e)>0)) AND (NOT (4tmao(e)>0)) AND (NOT (' \
               'dmso(e)>0)) AND (for(e)>0)) OR (b0001 == 1)'

        from mewpy.germ.algebra import Expression, parse_expression, Symbolic, And, Or
        from mewpy.germ.variables import Regulator

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

        expr.truth_table(strategy='max')
        expr.truth_table(strategy='all')

        variables.get('b0001').coefficients = (1, )

        expr.truth_table(strategy='max')
        expr.truth_table(strategy='all')

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
        from mewpy.germ.variables import Gene, Metabolite, Reaction, Target, Interaction, Regulator

        # integrated model
        gene = Gene(2)
        metabolite = Metabolite(2)
        reaction = Reaction(2)
        target = Target(2)
        interaction = Interaction(2)
        regulator = Regulator(2)

        from mewpy.germ.models import Model
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
        from mewpy.germ.models import MetabolicModel
        model = MetabolicModel('MetabolicModel')

        metabolite1 = Metabolite('o2')
        metabolite2 = Metabolite('h2o2')
        gene1 = Gene('b0001')

        from mewpy.germ.algebra import parse_expression, Expression
        expr1 = parse_expression('b0001')

        gpr1 = Expression(expr1, {gene1.id: gene1})

        reaction1 = Reaction(identifier='R0001',
                             stoichiometry={metabolite1: -1},
                             bounds=(0.0, 999999),
                             gpr=gpr1)

        rule = 'b0002 and (b0003 or b0004)'

        reaction2 = Reaction.from_gpr_string(identifier='R0002',
                                             rule=rule,
                                             bounds=(0.0, 999999),
                                             stoichiometry={metabolite1: -1, metabolite2: 1})

        model.add(reaction1, reaction2)

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
        from mewpy.germ.models import RegulatoryModel
        model = RegulatoryModel('RegulatoryModel')

        target = Target('b0001')
        regulator = Regulator.from_types(('regulator', 'target'), identifier='b0002')

        expr1 = parse_expression('b0002')

        reg_event1 = Expression(expr1, {'b0002': regulator})

        interaction1 = Interaction('I_b0001', target=target, regulatory_events={1.0: reg_event1})

        rule = '((NOT (o2(e)>0)) AND (NOT (no3(e)>0)) AND (NOT (no2(e)>0)) AND (NOT (4tmao(e)>0)) AND (NOT (dmso(' \
               'e)>0)) AND (for(e)>0)) '

        interaction2 = Interaction.from_string('I_b0002', rule=rule, target=regulator)

        model.add(interaction1, interaction2)

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

        regulatory_reader = Reader(Engines.RegulatorySBML, MODELS_PATH.joinpath('e_coli_lac.xml'))
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
        regulatory_reader = Reader(Engines.RegulatorySBML, MODELS_PATH.joinpath('e_coli_lac.xml'))
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
        cobra_ecoli_core_model = read_sbml_model(str(EC_CORE_MODEL))
        metabolic_reader = Reader(Engines.CobraModel, cobra_ecoli_core_model)

        model = read_model(metabolic_reader)
        self.assertEqual(len(model.reactions), 95)
        self.assertEqual(len(model.metabolites), 72)
        self.assertEqual(len(model.genes), 137)

        # from reframed
        from reframed import load_cbmodel
        reframed_ecoli_core_model = load_cbmodel(str(EC_CORE_MODEL))
        metabolic_reader = Reader(Engines.ReframedModel, reframed_ecoli_core_model)

        model = read_model(metabolic_reader)
        self.assertEqual(len(model.reactions), 95)
        self.assertEqual(len(model.metabolites), 72)
        self.assertEqual(len(model.genes), 137)

        # from json
        model_reader = Reader(Engines.JSON, MODELS_PATH.joinpath('e_coli_core.json'))

        model = read_model(model_reader)
        self.assertEqual(len(model.interactions), 159)
        self.assertEqual(len(model.targets), 159)
        self.assertEqual(len(model.regulators), 45)
        self.assertEqual(len(model.reactions), 95)
        self.assertEqual(len(model.metabolites), 72)
        self.assertEqual(len(model.genes), 137)

    #@pytest.mark.xfail
    def test_write(self):
        """
        Tests write model
        """

        from mewpy.io import read_model, Writer, Engines, write_model
        import os 
        model = read_model(self.regulatory_reader, self.metabolic_reader)

        # to sbml
        metabolic_writer = Writer(Engines.MetabolicSBML,
                                  io=MODELS_PATH.joinpath('e_coli_core_write.xml'),
                                  model=model)

        regulatory_writer = Writer(Engines.RegulatorySBML,
                                   io=MODELS_PATH.joinpath('e_coli_lac_write.xml'),
                                   model=model)

        write_model(regulatory_writer, metabolic_writer)
        os.remove(MODELS_PATH.joinpath('e_coli_core_write.xml'))
        os.remove(MODELS_PATH.joinpath('e_coli_lac_write.xml'))
                  

        # to json
        model_writer = Writer(Engines.JSON,
                              io=MODELS_PATH.joinpath('e_coli_core_write.json'),
                              model=model)

        write_model(model_writer)
        os.remove(MODELS_PATH.joinpath('e_coli_core_write.json'))

    @pytest.mark.xfail
    def test_analysis(self):
        """
        Tests model analysis
        """

        from mewpy.io import read_model

        # metabolic analysis
        model = read_model(self.metabolic_reader)
        model.objective = {'Biomass_Ecoli_core': 1}

        # fba
        from mewpy.germ.analysis import FBA, slim_fba
        simulator = FBA(model)
        sol = simulator.optimize()
        self.assertGreater(sol.objective_value, 0)
        self.assertGreater(slim_fba(model), 0)

        # pfba
        from mewpy.germ.analysis import pFBA, slim_pfba
        simulator = pFBA(model)
        sol = simulator.optimize()
        self.assertGreater(sol.objective_value, 0)
        self.assertGreater(slim_pfba(model), 0)

        # deletions
        from mewpy.germ.analysis import single_reaction_deletion, single_gene_deletion
        reactions_deletion = single_reaction_deletion(model=model, reactions=list(model.reactions.keys())[0:10])
        self.assertGreater(len(reactions_deletion), 0)

        genes_deletion = single_gene_deletion(model=model, genes=list(model.genes.keys())[0:10])
        self.assertGreater(len(genes_deletion), 0)

        from mewpy.germ.analysis import fva
        fva_sol = fva(model=model, fraction=0.9, reactions=list(model.reactions.keys())[0:10])
        self.assertGreater(len(fva_sol), 0)

        # regulatory analysis
        model = read_model(self.regulatory_reader)

        # truth table/regulatory events
        from mewpy.germ.analysis import regulatory_truth_table
        truth_table = regulatory_truth_table(model=model, initial_state={'b4401': 0, 'b1334': 0})
        self.assertGreater(len(truth_table), 0)
        self.assertEqual(truth_table.loc['b2276', 'result'], 1)

        # integrated analysis
        # model = read_model(self.regulatory_reader, self.metabolic_reader)
        # changing to sample because of CPLEX community edition
        from mewpy.io import Reader, Engines
        metabolic_reader = Reader(Engines.MetabolicSBML, SAMPLE_MODEL)
        regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                                   SAMPLE_REG_MODEL,
                                   sep=',',
                                   id_col=0,
                                   rule_col=1)

        model = read_model(regulatory_reader, metabolic_reader)

        _BIOMASS_ID = 'r11'

        model.objective = {_BIOMASS_ID: 1}
        model.get('pH').coefficients = (0, 14)

        from mewpy.germ.analysis import slim_srfba
        self.assertGreater(slim_srfba(model), 0)

        from mewpy.germ.analysis import slim_rfba
        self.assertIsNotNone(slim_rfba(model))

        sol = simulator.optimize(dynamic=True)
        self.assertIsNotNone(sol)

        # ifva
        from mewpy.germ.analysis import ifva
        sol = ifva(model, reactions=list(model.reactions.keys())[0:10])
        self.assertGreater(len(sol), 0)

        from mewpy.germ.analysis import isingle_reaction_deletion
        sol = isingle_reaction_deletion(model, reactions=list(model.reactions.keys())[0:10])
        self.assertGreater(len(sol), 0)

        from mewpy.germ.analysis import isingle_gene_deletion
        sol = isingle_gene_deletion(model, genes=list(model.genes.keys())[0:10])
        self.assertGreater(len(sol), 0)

        from mewpy.germ.analysis import isingle_regulator_deletion
        sol = isingle_regulator_deletion(model, regulators=list(model.regulators.keys())[0:10])
        self.assertGreater(len(sol), 0)

    #@pytest.mark.xfail
    def test_analysis_expression(self):
        """
        It tests model analysis with methods of expression
        """
        from mewpy.io import Reader, Engines, read_model
        metabolic_reader = Reader(Engines.MetabolicSBML, SAMPLE_MODEL)
        regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                                   SAMPLE_REG_MODEL,
                                   sep=',',
                                   id_col=0,
                                   rule_col=1)

        model = read_model(regulatory_reader, metabolic_reader)

        model.objective = {'r11': 1}
        probabilities = {
            ('g29', 'g10'): 0.1,
            ('g29', 'g11'): 0.1,
            ('g29', 'g12'): 0.1,
            ('g30', 'g10'): 0.9,
            ('g30', 'g11'): 0.9,
            ('g30', 'g12'): 0.9,
            ('g35', 'g34'): 0.1,
        }

        from mewpy.germ.analysis import PROM
        simulator = PROM(model).build()
        sol = simulator.optimize(initial_state=probabilities, regulators=['g29', 'g30', 'g35'])
        self.assertGreater(sol.solutions['ko_g35'].objective_value, 0)

        predicted_expression = {
            'g10': 2,
            'g11': 2.3,
            'g12': 2.3,
            'g34': 0.8,
        }

        from mewpy.germ.analysis import CoRegFlux
        simulator = CoRegFlux(model).build()
        sol = simulator.optimize(initial_state=predicted_expression)
        self.assertGreater(sol.objective_value, 0)

    @pytest.mark.xfail
    def test_simulation(self):
        """
        Tests model simulation
        """

        from mewpy.io import Reader, Engines, read_model

        metabolic_reader = Reader(Engines.MetabolicSBML, SAMPLE_MODEL)
        regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                                   SAMPLE_REG_MODEL,
                                   sep=',',
                                   id_col=0,
                                   rule_col=1)

        model = read_model(regulatory_reader, metabolic_reader)

        # pH (ph > 5) controls g20, which belongs to the r11 gpr
        model.get('pH').coefficients = (0, 14)

        self.assertEqual(model.id, 'COBRAModel')
        self.assertEqual(model.name, 'Model Exported from COBRA Toolbox')

        self.assertEqual(len(model.types), 2)
        self.assertEqual(len(model.simulators), 0)
        self.assertEqual(model.history.history.shape, (1, 4))
        self.assertEqual(len(model.contexts), 0)

        self.assertEqual(len(model.containers), 10)

        self.assertEqual(len(model.interactions), 36)
        self.assertEqual(len(model.targets), 36)
        self.assertEqual(len(model.regulators), 19)
        self.assertEqual(len(model.environmental_stimuli), 1)

        self.assertEqual(model.objective, {model.get('r11'): 1})
        self.assertEqual(len(model.reactions), 17)
        self.assertEqual(len(model.metabolites), 14)
        self.assertEqual(len(model.genes), 22)
        self.assertEqual(len(model.sinks), 0)
        self.assertEqual(len(model.exchanges), 2)
        self.assertEqual(len(model.demands), 0)

        self.assertEqual(len(model.compartments), 1)
        self.assertEqual(model.external_compartment, 'c')

        # fba
        from mewpy.germ.analysis import FBA
        fba = FBA(model)
        sol = fba.optimize()
        self.assertGreater(sol.x.get('r11'), 0)
        self.assertEqual(len(model.simulators), 0)

        # pfba
        from mewpy.germ.analysis import pFBA
        pfba = pFBA(model)
        sol = pfba.optimize()
        self.assertGreater(sol.x.get('r11'), 0)
        self.assertEqual(len(model.simulators), 0)

        # RFBA
        from mewpy.germ.analysis import RFBA
        initial_state = {
            'pH': 7,
            'g21': 1,
            'g22': 0,
            'g23': 1,
            'g24': 1,
            'g25': 1,
            'g26': 1,
            'g27': 1,
            'g28': 1,
            'g29': 1,
            'g30': 0,
            'g31': 0,
            'g32': 0,
            'g33': 1,
            'g36': 1,
            'r3': 100,
            'r15': 0,
            'r6': 100}

        rfba = RFBA(model)
        sol = rfba.optimize(initial_state=initial_state)
        self.assertGreater(sol.x.get('r11'), 0)
        self.assertEqual(len(model.simulators), 0)

        # SRFBA
        from mewpy.germ.analysis import SRFBA
        srfba = SRFBA(model)
        sol = srfba.optimize()
        self.assertGreater(sol.x.get('r11'), 0)
        self.assertEqual(len(model.simulators), 0)

        # multiple simulators attached
        fba = FBA(model, attach=True)
        pfba = pFBA(model, attach=True)
        srfba = SRFBA(model, attach=True)

        model.get('r16').ko()

        fba_sol = fba.optimize()
        pfba_sol = pfba.optimize()
        srfba_sol = srfba.optimize()

        # if r16 is knocked-out, only r8 can have flux, and vice-versa
        self.assertGreater(fba_sol.x.get('r8'), 333)
        self.assertGreater(pfba_sol.x.get('r8'), 333)
        self.assertGreater(srfba_sol.x.get('r8'), 333)

        model.undo()

        solver_kwargs = {'constraints': {'r16': (0, 0), 'r8': (0, 0)}}

        fba_sol = fba.optimize(solver_kwargs=solver_kwargs)
        pfba_sol = pfba.optimize(solver_kwargs=solver_kwargs)
        srfba_sol = srfba.optimize(solver_kwargs=solver_kwargs)
        self.assertEqual(fba_sol.objective_value, 0.0)
        self.assertEqual(pfba_sol.x.get('r11'), 0.0)
        self.assertEqual(srfba_sol.objective_value, 0.0)

        fba_sol = fba.optimize()
        pfba_sol = pfba.optimize()
        srfba_sol = srfba.optimize()
        self.assertGreater(fba_sol.objective_value, 0.0)
        self.assertGreater(pfba_sol.x.get('r11'), 0.0)
        self.assertGreater(srfba_sol.objective_value, 0.0)

    @pytest.mark.xfail
    def test_bounds_coefficients(self):
        """
        Tests model bounds and coefficients workflow
        """
        from mewpy.io import Reader, Engines, read_model

        metabolic_reader = Reader(Engines.MetabolicSBML, SAMPLE_MODEL)
        regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                                   SAMPLE_REG_MODEL,
                                   sep=',',
                                   id_col=0,
                                   rule_col=1)

        model = read_model(regulatory_reader, metabolic_reader)

        # pH (ph > 5) controls g20, which belongs to the r11 gpr
        model.get('pH').coefficients = (0, 14)

        # multiple simulators attached
        from mewpy.germ.analysis import FBA, pFBA, SRFBA
        fba = FBA(model, attach=True)
        pfba = pFBA(model, attach=True)
        srfba = SRFBA(model, attach=True)

        old_lb, old_ub = tuple(model.reactions.get('r16').bounds)
        model.get('r16').ko()
        self.assertEqual(model.get('r16').bounds, (0, 0))

        fba_sol = fba.optimize()
        pfba_sol = pfba.optimize()
        srfba_sol = srfba.optimize()

        # if r16 is knocked-out, only r8 can have flux, and vice-versa
        self.assertGreater(fba_sol.x.get('r8'), 333)
        self.assertGreater(pfba_sol.x.get('r8'), 333)
        self.assertGreater(srfba_sol.x.get('r8'), 333)

        # revert bound change
        model.undo()
        self.assertEqual(model.get('r16').bounds, (old_lb, old_ub))

        # using model context so all changes made to the model are reverted upon exiting the context
        with model:
            min_coef, max_coef = model.get('g14').coefficients
            lb, ub = model.get('g14').reactions.get('r8').bounds

            # a gene ko does not change the bounds of the associated reactions
            model.get('g14').ko()

            self.assertEqual(model.get('g14').coefficients, (0.0, 0.0))
            self.assertEqual(model.get('r8').bounds, (lb, ub))

            fba_sol = fba.optimize()
            pfba_sol = pfba.optimize()
            srfba_sol = srfba.optimize()

            # But the analysis methods can retrieve such change
            self.assertLess(fba_sol.x.get('r8'), 333)
            self.assertGreater(fba_sol.x.get('r16'), 333)

            self.assertLess(pfba_sol.x.get('r8'), 333)
            self.assertGreater(pfba_sol.x.get('r16'), 333)

            self.assertLess(srfba_sol.x.get('r8'), 333)
            self.assertGreater(srfba_sol.x.get('r16'), 333)
            # the gene is a variable of the srfba formulation
            self.assertEqual(srfba_sol.x.get('g14'), 0)

            # we can have as many contexts as we want
            with model:
                model.get('r16').ko()

                fba_sol = fba.optimize()
                pfba_sol = pfba.optimize()
                srfba_sol = srfba.optimize()

                # we had knocked-out all reactions that can produce I, which is reactant of r11
                self.assertLess(fba_sol.x.get('r11'), 333)
                self.assertLess(pfba_sol.x.get('r11'), 333)
                self.assertLess(srfba_sol.x.get('r11'), 333)

            # exiting the second context
            self.assertEqual(model.get('r16').bounds, (old_lb, old_ub))

            # we can use undo and redo within a context. However, this only undoes or redoes changes made to the
            # model within the given context.
            # In this case, we had revert model.get('g14').ko()
            model.undo()
            self.assertEqual(model.get('g14').coefficients, (min_coef, max_coef))

            # using the regulatory network
            # This affects the target g34 which is also the metabolic gene for r16
            model.get('g35').ko()
            # This affects the target g14 which is also the metabolic gene for r8. However, the regulatory
            # interaction for the g14 target is the following: not g31. So, we are activating the g14 and
            # respectively the reaction r8
            model.get('g31').ko()

            fba.optimize()
            pfba.optimize()
            srfba_sol = srfba.optimize()

            # fba is not affected by the regulatory share,
            # so that we cannot test whether there is flux through r8 or r16
            self.assertEqual(srfba_sol.x.get('g14'), 1)
            self.assertEqual(srfba_sol.x.get('g31'), 0)
            self.assertEqual(srfba_sol.x.get('g34'), 0)
            self.assertEqual(srfba_sol.x.get('g35'), 0)
            self.assertGreater(srfba_sol.x.get('r8'), 333)
            self.assertLess(srfba_sol.x.get('r16'), 333)

        # exiting the first context
        # redo the knock-out to the r16. All changes made to the model within a context are not recorded in the
        # model.history. Once again, these changes are also reverted upon exiting the context
        model.redo()

        fba_sol = fba.optimize()
        pfba_sol = pfba.optimize()
        srfba_sol = srfba.optimize()

        # if r16 is knocked-out, only r8 can have flux, and vice-versa
        self.assertGreater(fba_sol.x.get('r8'), 333)
        self.assertGreater(pfba_sol.x.get('r8'), 333)
        self.assertGreater(srfba_sol.x.get('r8'), 333)

    @pytest.mark.xfail
    def test_manipulation(self):
        """
        Tests model manipulation workflow
        """

        from mewpy.io import Reader, Engines, read_model

        metabolic_reader = Reader(Engines.MetabolicSBML, SAMPLE_MODEL)
        regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                                   SAMPLE_REG_MODEL,
                                   sep=',',
                                   id_col=0,
                                   rule_col=1)

        model = read_model(regulatory_reader, metabolic_reader)

        # pH (ph > 5) controls g20, which belongs to the r11 gpr
        model.get('pH').coefficients = (0, 14)

        # for rfba
        initial_state = {
            'pH': 7,
            'g21': 1,
            'g22': 0,
            'g23': 1,
            'g24': 1,
            'g25': 1,
            'g26': 1,
            'g27': 1,
            'g28': 1,
            'g29': 1,
            'g30': 0,
            'g31': 0,
            'g32': 0,
            'g33': 1,
            'g36': 1,
            'r3': 100,
            'r15': 0,
            'r6': 100}

        # multiple simulators attached
        from mewpy.germ.analysis import FBA, pFBA, RFBA, SRFBA
        fba = FBA(model, attach=True)
        pfba = pFBA(model, attach=True)
        rfba = RFBA(model, attach=True)
        srfba = SRFBA(model, attach=True)
        simulators = [fba, pfba, rfba, srfba]

        # adding new reactions, metabolites, genes, interactions, targets and regulators
        # g37: g39 and not g40
        # g38: g39 and not g40
        # r17: n <-> 2o | g37 or g38
        # r18: n <-
        # r19: o ->
        from mewpy.germ.variables import Variable, Regulator, Interaction, Metabolite, Reaction
        g39 = Regulator(identifier='g39', coefficients=(1, 1))
        g40 = Regulator(identifier='g40', coefficients=(0, 0))

        # the targets g37 and g38 will be also metabolic genes, so that regulators g39 and g40 will control the flux
        # expression of r17
        g37 = Variable.from_types(types=('target', 'gene'), identifier='g37')
        g38 = Variable.from_types(types=('target', 'gene'), identifier='g38')

        from mewpy.germ.algebra import Expression, parse_expression
        i_g37_expression = Expression(parse_expression('g39 and not g40'), {'g39': g39, 'g40': g40})
        i_g38_expression = Expression(parse_expression('g39 and not g40'), {'g39': g39, 'g40': g40})

        # it is always a good practice to build the expression of a given interaction first, and then use it in the
        # Interaction constructor. Otherwise, interaction has alternative constructors (from_expression or from_string)
        i_g37 = Interaction(identifier='Interaction_g37',
                            regulatory_events={1.0: i_g37_expression},
                            target=g37)

        i_g38 = Interaction(identifier='Interaction_g38',
                            regulatory_events={1.0: i_g38_expression},
                            target=g38)

        n = Metabolite(identifier='N',
                       charge=0,
                       compartment='c',
                       formula='C12H24O6')

        o = Metabolite(identifier='O',
                       charge=2,
                       compartment='c',
                       formula='C12H24O12')

        r17_gpr = Expression(parse_expression('g37 or g38'), {'g37': g37, 'g38': g38})

        # it is always a good practice to build the gpr of a given reaction first, and then use it in the Reaction
        # constructor. Otherwise, reaction has alternative constructors (from_gpr_expression or from_gpr_string)
        r17 = Reaction(identifier='r17',
                       bounds=(-1000, 1000),
                       gpr=r17_gpr,
                       stoichiometry={n: -1, o: 2})

        # exchanges
        r18 = Reaction(identifier='r18',
                       bounds=(-1000, 0),
                       stoichiometry={n: -1})

        r19 = Reaction(identifier='r19',
                       bounds=(0, 1000),
                       stoichiometry={o: -1})

        # note that, although the metabolic genes of r17 are linked to the interactions i_g37 and i_g38, we still had
        # to add both interactions and reactions to the model, so that the model comprehends the regulatory and
        # metabolic phenomena. Models are always modular and the main containers static, so that it is possible to
        # work with either the metabolic or regulatory share, even if the variables are linked
        model.add(i_g37, i_g38, r17, r18, r19)

        self.assertEqual(len(model.interactions), 38)
        self.assertEqual(len(model.targets), 38)
        self.assertEqual(len(model.regulators), 21)
        self.assertEqual(len(model.environmental_stimuli), 3)

        self.assertEqual(model.objective, {model.get('r11'): 1})
        self.assertEqual(len(model.reactions), 20)
        self.assertEqual(len(model.metabolites), 16)
        self.assertEqual(len(model.genes), 24)
        self.assertEqual(len(model.exchanges), 4)

        model.get('r16').ko()

        # updating for the add interactions and reactions
        for simulator in simulators:
            simulator.update()

        fba_sol = fba.optimize()
        pfba_sol = pfba.optimize()
        rfba_sol = rfba.optimize(initial_state=initial_state)
        srfba_sol = srfba.optimize()

        # if r16 is knocked-out, only r8 can have flux, and vice-versa.
        # r17, r18 and r19 do not affect the rest of the network
        self.assertGreater(fba_sol.x.get('r8'), 333)
        self.assertGreater(pfba_sol.x.get('r8'), 333)
        self.assertGreater(rfba_sol.x.get('r8'), 333)
        self.assertGreater(srfba_sol.x.get('r8'), 333)

        model.objective = {'r17': 1}
        # otherwise, the remaining network can have flux or not
        model.reactions.get('r0').bounds = (0, 0)

        # updating for the add interactions and reactions
        for simulator in simulators:
            simulator.update()

        fba_sol = fba.optimize()
        pfba_sol = pfba.optimize()
        rfba_sol = rfba.optimize(initial_state=initial_state)
        srfba_sol = srfba.optimize()

        # different objective
        self.assertLess(fba_sol.x.get('r8'), 333)
        self.assertLess(pfba_sol.x.get('r8'), 333)
        self.assertLess(rfba_sol.x.get('r8'), 333)
        self.assertLess(srfba_sol.x.get('r8'), 333)

        self.assertGreater(fba_sol.x.get('r17'), 450)
        self.assertGreater(pfba_sol.x.get('r17'), 450)
        self.assertGreater(rfba_sol.x.get('r17'), 450)
        self.assertGreater(srfba_sol.x.get('r17'), 450)

        # resetting the model to the first state. This involves removing the interactions and reactions added earlier
        model.reset()

        # This was also reset
        model.get('pH').coefficients = (0, 14)
        model.get('r16').ko()

        for simulator in simulators:
            simulator.update()

        fba_sol = fba.optimize()
        pfba_sol = pfba.optimize()
        rfba_sol = rfba.optimize(initial_state=initial_state)
        srfba_sol = srfba.optimize()

        self.assertGreater(fba_sol.x.get('r8'), 333)
        self.assertGreater(pfba_sol.x.get('r8'), 333)
        self.assertGreater(rfba_sol.x.get('r8'), 333)
        self.assertGreater(srfba_sol.x.get('r8'), 333)

        # adding a blocked by-product
        n = Metabolite(identifier='N',
                       charge=0,
                       compartment='c',
                       formula='C12H24O6')

        model.reactions.get('r8').add_metabolites(stoichiometry={n: 1})

        for simulator in simulators:
            simulator.update()

        fba_sol = fba.optimize()
        pfba_sol = pfba.optimize()
        rfba_sol = rfba.optimize(initial_state=initial_state)
        srfba_sol = srfba.optimize()

        self.assertLess(fba_sol.x.get('r11'), 1)
        self.assertLess(pfba_sol.x.get('r11'), 1)
        self.assertLess(rfba_sol.x.get('r11'), 1)
        self.assertLess(srfba_sol.x.get('r11'), 1)

        model.undo()
        for simulator in simulators:
            simulator.update()

        fba_sol = fba.optimize()
        pfba_sol = pfba.optimize()
        rfba_sol = rfba.optimize(initial_state=initial_state)
        srfba_sol = srfba.optimize()

        self.assertGreater(fba_sol.x.get('r11'), 1)
        self.assertGreater(pfba_sol.x.get('r11'), 1)
        self.assertGreater(rfba_sol.x.get('r11'), 1)
        self.assertGreater(srfba_sol.x.get('r11'), 1)

        # adding a repressor to the objective reaction
        # r11 gpr is g18 & g19 & g20
        # g18 and g19 are equally regulated by g33
        reg_g33 = model.get('g33')
        reg_g33.coefficients = (1,)

        from mewpy.germ.algebra import Not, Symbol, Expression
        symbolic = Not(variables=[Symbol(value=reg_g33.id)])
        variables = {'g33': reg_g33}
        regulatory_event = Expression(symbolic=symbolic, variables=variables)

        i_g18 = model.get('g18').interaction
        i_g18_reg_event = i_g18.regulatory_events[1]
        # this will replace the regulatory event that determines a target coefficient of 1
        i_g18.add_regulatory_event(coefficient=1, expression=regulatory_event)

        for simulator in simulators:
            simulator.update()

        fba_sol = fba.optimize()
        pfba_sol = pfba.optimize()
        rfba_sol = rfba.optimize(initial_state=initial_state)
        srfba_sol = srfba.optimize()

        self.assertGreater(fba_sol.x.get('r11'), 1)
        self.assertGreater(pfba_sol.x.get('r11'), 1)

        # only analysis methods that consider the regulatory network are affected
        self.assertLess(rfba_sol.x.get('r11'), 1)
        self.assertLess(srfba_sol.x.get('r11'), 1)

        # this will replace the regulatory event that determines a target coefficient of 1
        i_g18.add_regulatory_event(coefficient=1, expression=i_g18_reg_event)
        for simulator in simulators:
            simulator.update()

        fba_sol = fba.optimize()
        pfba_sol = pfba.optimize()
        rfba_sol = rfba.optimize(initial_state=initial_state)
        srfba_sol = srfba.optimize()

        self.assertGreater(fba_sol.x.get('r11'), 1)
        self.assertGreater(pfba_sol.x.get('r11'), 1)
        self.assertGreater(rfba_sol.x.get('r11'), 1)
        self.assertGreater(srfba_sol.x.get('r11'), 1)

    def test_serialization(self):
        """
        Tests model serialization workflow
        """

        from mewpy.io import Reader, Engines, read_model

        metabolic_reader = Reader(Engines.MetabolicSBML, SAMPLE_MODEL)
        regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                                   SAMPLE_REG_MODEL,
                                   sep=',',
                                   id_col=0,
                                   rule_col=1)

        model = read_model(regulatory_reader, metabolic_reader)

        res = model.to_dict(variables=False)
        self.assertEqual(len(res), 10)

        dict_model = model.from_dict(res)
        self.assertEqual(len(model.reactions), len(dict_model.reactions))

        shallow_copy_model = dict_model.copy()
        deep_copy_model = dict_model.deepcopy()

        self.assertIs(shallow_copy_model.get('r6'), dict_model.get('r6'))
        self.assertIsNot(deep_copy_model.get('r6'), dict_model.get('r6'))


if __name__ == '__main__':
    unittest.main()
