from abc import ABC, abstractmethod
from mewpy.utils.constants import EAConstants

import copy
# pathos
try:
    import pathos.multiprocesssing
    from pathos.multiprocessing import Pool
except ImportError:
    import multiprocessing
    from multiprocessing.pool import Pool
# dask
try:
    import dask
except ImportError:
    pass
# pyspark
try:
    from pyspark import SparkConf, SparkContext
except ImportError:
    pass


def cpu_count():
    """ The number of cpus
        Return EAConstants.NUM_CPUS if it is set to a positive number
        otherwise return half of the available cpus.
    """
    if EAConstants.NUM_CPUS > 0:
        return EAConstants.NUM_CPUS
    else:
        try:
            return multiprocessing.cpu_count()//2
        except (ImportError, NotImplementedError):
            return 1


class Evaluator(ABC):

    @abstractmethod
    def evaluate(self, candidates, args):
        raise NotImplementedError


class MultiProcessorEvaluator(Evaluator):

    def __init__(self, evaluator, mp_num_cpus):
        """
        Evaluate the candidates in parallel using ``multiprocessing``.
        """
        self.pool = Pool(mp_num_cpus)
        self.evaluator = evaluator

    def evaluate(self, candidates, args):
        """
        Values in args will be ignored and not passed to the evaluator to avoid unnecessary pickling in inspyred.
        """
        results = self.pool.map(self.evaluator, candidates)
        return results


class DaskEvaluator(Evaluator):

    def __init__(self, evaluator, mp_num_cpus, scheduler='processes'):
        self.evaluator = evaluator
        self.scheduler = scheduler

    def evaluate(self, candidates, args):
        with dask.config.set(scheduler=self.scheduler):
            return list(dask.compute(*[dask.delayed(self.evaluator)(c) for c in candidates]))


class SparkEvaluator(Evaluator):
    def __init__(self, evaluator, mp_num_cpus):
        self.evaluator = evaluator
        self.spark_conf = SparkConf().setAppName(
            "mewpy").setMaster(f"local[{mp_num_cpus}]")
        self.spark_context = SparkContext(conf=self.spark_conf)

    def evaluate(self, candidates, args):
        solutions_to_evaluate = self.spark_context.parallelize(candidates)

        return solutions_to_evaluate.map(lambda s: self.evaluator(s)).collect()


# ray
try:
    import ray
except ImportError:
    pass
else:

    @ray.remote
    class RayActor:
        """
        Each actor (worker) has a solver instance to overcome the need to serialize solvers which may not be pickable.
        The solver is not reset before each evaluation.
        """

        def __init__(self, problem):
            self.problem = copy.deepcopy(problem)

        def evaluate_candidates(self, candidates):
            """Evaluates a sublist of candidates
            """
            return self.problem.evaluator(candidates, None)

    class RayEvaluator(Evaluator):

        def __init__(self, problem, number_of_actors):
            ray.init(ignore_reinit_error=True)
            self.actors = [RayActor.remote(problem)
                           for _ in range(number_of_actors)]
            self.number_of_actors = len(self.actors)
            print(f"Using {self.number_of_actors} workers.")

        def evaluate(self, candidates, args):
            """
            Divides the candidates into sublists to be evaluated by each actor
            """
            size = len(candidates) // self.number_of_actors
            if len(candidates) % self.number_of_actors != 0:
                size += 1
            sub_lists = [candidates[i:i+size]
                         for i in range(0, len(candidates), size)]
            values = []
            for i in range(self.number_of_actors):
                actor = self.actors[i]
                values.append(actor.evaluate_candidates.remote(sub_lists[i]))
            r = ray.get(values)
            result = []
            for x in r:
                for y in x:
                    result.append(y)
            return result
