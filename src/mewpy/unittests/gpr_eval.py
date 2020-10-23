from mewpy.utils.parsing import build_tree, BooleanEvaluator, Boolean, GeneEvaluator, Arithmetic, ArithmeticEvaluator


def test_0():
    # aritmetic example
    t = build_tree(" 1 + 2 + 3 + ( 2 * 2 )", Arithmetic)
    res = t.evaluate(ArithmeticEvaluator.f_operand, ArithmeticEvaluator.f_operator)
    print(t, ' = ', res)


def test_1():
    # boolean example with conditions
    expression = "( (Lrp AND NOT (leu_L_e_>0)) OR NOT(((GlnG AND GlnB AND GlnD) AND RpoN) AND ((glu_L_e_>0) OR (arg_L_e_>0) OR (asp_L_e_>0) OR (his_L_e_>0) OR (pro_L_e_>0) )))"
    # expression = " (x > 0 or C or B) and not ph == 5 "
    t = build_tree(expression, Boolean)

    print(expression, ' <<==>> ', t)
    t.print_node()
    # list of true expressions
    true_list = ['GlnG']
    # dict of variables values
    v = {'leu_L_e_': 1, 'glu_L_e_': 7, "arg_L_e_": 0.5, "asp_L_e_": 2, "his_L_e_": 0, "pro_L_e_": 0}
    # propositions in the list are evaluated as True and the remaining are False
    evaluator = BooleanEvaluator(true_list, v)
    res = t.evaluate(evaluator.f_operand, evaluator.f_operator)
    print('evaluation: ', res)
    print('True list:', true_list)
    print('variables:', v)
    print('operands: ', t.get_operands())
    print('conditions: ', t.get_conditions())


def test_2():
    from urllib.request import urlretrieve
    from cobra.io import read_sbml_model
    import random

    path, _ = urlretrieve('http://bigg.ucsd.edu/static/models/RECON1.xml')
    model = read_sbml_model(path)
    ogpr = model.reactions.ATPS4m.gene_name_reaction_rule
    gpr = ogpr
    print(gpr)
    t = build_tree(gpr, Boolean)
    print(t)
    genes = list(t.get_operands())
    print("GENES:\n", genes)
    print("Evaluations:")
    evaluator = BooleanEvaluator(genes)
    res = t.evaluate(evaluator.f_operand, evaluator.f_operator)
    print(evaluator.true_list, " ==> ", res)
    for _ in range(20):
        l = []
        n = random.randint(1, len(genes))
        for _ in range(n):
            i = random.randint(0, len(genes) - 1)
            l.append(genes[i])
        evaluator.set_true_list(l)
        res = t.evaluate(evaluator.f_operand, evaluator.f_operator)
        print(evaluator.true_list, " ==> ", res)


def test_3():
    # Gene OU example
    gpr = "((G_YIL043C and G_YMR015C and G_YNL111C) or (G_YKL150W and G_YMR015C and G_YNL111C))"
    genes = {'G_YER091C': 0, 'G_YMR105C': 0.03125, 'G_YNL117W': 0.5, 'G_YNL111C': 0.125, 'G_YJR158W': 0.0625,
             'G_YLR355C': 0.5}
    operators = (lambda x, y: min(x, y), lambda x, y: max(x, y))
    evaluator = GeneEvaluator(genes, operators[0], operators[1])
    tree = build_tree(gpr, Boolean)
    tree.print_node()
    lv = tree.evaluate(evaluator.f_operand, evaluator.f_operator)
    print(f"{gpr}\n{genes}\nlevel:{lv}\n")


if __name__ == '__main__':
    test_2()
    # print()
    # test_2()
    # print()
    # test_3()
