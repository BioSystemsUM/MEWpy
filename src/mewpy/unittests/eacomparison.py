"""
EA comparison
"""
from jmetal.algorithm.multiobjective.nsgaiii import NSGAIII, UniformReferenceDirectionFactory
from jmetal.algorithm.multiobjective.nsgaii import NSGAII
from jmetal.algorithm.multiobjective.spea2 import SPEA2
from jmetal.algorithm.multiobjective.ibea import IBEA
from jmetal.core.quality_indicator import GenerationalDistance, EpsilonIndicator, HyperVolume
from jmetal.lab.experiment import Experiment, Job, generate_summary_from_experiment
from jmetal.util.termination_criterion import StoppingByEvaluations
from jmetal.operator import BinaryTournamentSelection
from jmetal.util.comparator import Comparator, MultiComparator
from jmetal.util.evaluator import MultiprocessEvaluator



candidate_max_size = 6
max_evaluations = 100
N_CPU = 2


def configure_experiment(problems: dict, n_run: int):

    from mewpy.optimization.jmetal.operators import UniformCrossoverOU, GrowMutationOU, ShrinkMutation,SingleMutationOU, MutationContainer
    crossover = UniformCrossoverOU(0.5,max_size=candidate_max_size)
    mutators = []
    mutators.append(GrowMutationOU(1.0, max_size=candidate_max_size))
    mutators.append(ShrinkMutation(1.0, min_size=candidate_max_size))
    mutators.append(SingleMutationOU(1.0))
    mutation = MutationContainer(0.3, mutators=mutators)


    jobs = []
    

    for run in range(n_run):
        for problem_tag, problem in problems.items():
            jobs.append(
                Job(
                    algorithm=NSGAII(
                        problem=problem,
                        population_evaluator=MultiprocessEvaluator(N_CPU),
                        population_size=100,
                        offspring_population_size=100,
                        mutation= mutation,
                        crossover=crossover,
                        termination_criterion=StoppingByEvaluations(max_evaluations=max_evaluations)
                        ),
                    algorithm_tag='NSGAII',
                    problem_tag=problem_tag,
                    run=run,
                )
            )
            jobs.append(
                Job(
                    algorithm=SPEA2(
                        problem=problem,
                        population_evaluator=MultiprocessEvaluator(N_CPU),
                        population_size=100,
                        offspring_population_size=100,
                        mutation= mutation,
                        crossover=crossover,
                        termination_criterion=StoppingByEvaluations(max_evaluations=max_evaluations)
                        ),
                    algorithm_tag='SPEA2',
                    problem_tag=problem_tag,
                    run=run
                    )
            )
            jobs.append(
                Job(
                    algorithm=IBEA(
                        problem=problem,
                        population_evaluator=MultiprocessEvaluator(N_CPU),
                        kappa=1.,
                        population_size = 100, 
                        offspring_population_size = 100,
                        mutation=mutation,
                        crossover=crossover,
                        termination_criterion=StoppingByEvaluations(max_evaluations=max_evaluations)
                        ),
                    algorithm_tag='IBEA',
                    problem_tag=problem_tag,
                    run=run,
                )
            )
        

    return jobs


if __name__ == '__main__':

    # Run the study
    output_directory = 'data'

    
    # Configure the experiments
    from mewpy.model.gecko import GeckoModel 
    from collections import OrderedDict 
    compound =  'r_1913'
    model = GeckoModel('single-pool', biomass_reaction_id='r_2111')
    model.set_objective({'r_2111': 1.0, 'r_4041': 0.0})
    envcond = OrderedDict()

    from mewpy.optimization.evaluation import BPCY,WYIELD
    from mewpy.simulation import SimulationMethod
    # the evaluation (objective) functions
    evaluator_1 = BPCY("r_2111", compound, method=SimulationMethod.lMOMA)
    evaluator_2 = WYIELD("r_2111", compound)

    from mewpy.problems import GeckoOUProblem
    p = GeckoOUProblem(model,
                              fevaluation=[evaluator_1, evaluator_2],
                              candidate_max_size=candidate_max_size)

    from mewpy.optimization.jmetal.problem import JMetalOUProblem
    problem = JMetalOUProblem(p)
    problem.number_of_objectives = 2
    
    problems = {'gecko':problem}
    
    jobs = configure_experiment(problems=problems, n_run=1)

    
    experiment = Experiment(output_dir=output_directory, jobs=jobs)
    experiment.run()

    import numpy
    import os
    from mewpy.optimization.ea import non_dominated_population, Solution

    for prb in problems.keys():
        # Builds a pareto front approximation from the results
        population =[]
        for r, _, fl in os.walk(output_directory):
                for file in fl:
                    if 'FUN' in file and prb in r:
                        with open(os.path.join(r, file), 'r') as f:
                            line = f.readline()
                            while line:
                                tokens = line.split()
                                fitness = [float(x) for x in tokens]
                                s = Solution(None,fitness)
                                population.append(s)
                                line = f.readline()
            
        population = non_dominated_population(population, maximize=False, filter_duplicate=False)
        # save to file
        pf_file = os.path.dirname(output_directory).join(prb+".pf")
        with open(pf_file,'w') as f:
            for s in population:
                f.write("".join([str(v)+"\t" for v in s.fitness]))
                f.write("\n")

    # Generate summary file
    generate_summary_from_experiment(
        input_dir=output_directory,
        quality_indicators=[GenerationalDistance(reference_front = output_directory), 
                            EpsilonIndicator(reference_front = output_directory),
                            HyperVolume([0,0])
                            ] 
    )