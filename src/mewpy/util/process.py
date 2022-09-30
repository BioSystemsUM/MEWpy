import copy
from abc import ABC, abstractmethod

from .constants import EAConstants, ModelConstants


MP_Evaluators = []

from multiprocessing import Process
from multiprocessing.pool import Pool as MPPool
MP_Evaluators.append('nodaemon')

# pathos
try:
    import pathos.multiprocessing as multiprocessing
    from pathos.multiprocessing import Pool
    MP_Evaluators.append('mp')

except ImportError as e:
    import multiprocessing
    from multiprocessing.pool import Pool
    MP_Evaluators.append('mp')
# dask
try:
    import dask
    MP_Evaluators.append('dask')
except ImportError:
    pass
# pyspark
try:
    from pyspark import SparkConf, SparkContext
    MP_Evaluators.append('spark')
except ImportError:
    pass


class NoDaemonProcess(Process):

    def _get_daemon(self):
        return False
    def _set_daemon(self,value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class NoDaemonProcessPool(MPPool):
    Process = NoDaemonProcess


def cpu_count():
    """ The number of cpus
        Return EAConstants.NUM_CPUS if it is set to a positive number
        otherwise return half of the available cpus.
    """
    if EAConstants.NUM_CPUS > 0:
        return EAConstants.NUM_CPUS
    else:
        try:
            return multiprocessing.cpu_count() // 2
        except (ImportError, NotImplementedError):
            return 1


class Evaluable(ABC):

    @abstractmethod
    def evaluator(self, candidates, *args):
        raise NotImplementedError


class Evaluator(ABC):

    """An interface for multiprocessing evaluators

    Raises:
        NotImplementedError: Requires an evaluated method to
        be implemented.
    """
    @abstractmethod
    def evaluate(self, candidates, args):
        raise NotImplementedError


class MultiProcessorEvaluator(Evaluator):

    def __init__(self, evaluator, mp_num_cpus):
        """A multiprocessing evaluator

        Args:
            evaluator(function): Evaluation function.
            mp_num_cpus(int): Number of CPUs
        """
        self.pool = Pool(mp_num_cpus)
        self.evaluator = evaluator
        self.__name__ = self.__class__.__name__

    def evaluate(self, candidates, args):
        """
        Values in args will be ignored and not passed to the evaluator to avoid unnecessary pickling in inspyred.
        """
        results = self.pool.map(self.evaluator, candidates)
        return results

    def __call__(self, candidates, args):
        return self.evaluate(candidates, args)


class NoDaemonMultiProcessorEvaluator(Evaluator):

    def __init__(self, evaluator, mp_num_cpus):
        """A multiprocessing evaluator

        Args:
            evaluator(function): Evaluation function.
            mp_num_cpus(int): Number of CPUs
        """
        self.pool = NoDaemonProcessPool(mp_num_cpus)
        self.evaluator = evaluator
        self.__name__ = self.__class__.__name__
        print('nodaemon')

    def evaluate(self, candidates, args):
        """
        Values in args will be ignored and not passed to the evaluator to avoid unnecessary pickling in inspyred.
        """
        results = self.pool.map(self.evaluator, candidates)
        return results

    def __call__(self, candidates, args):
        return self.evaluate(candidates, args)


class DaskEvaluator(Evaluator):

    def __init__(self, evaluator, mp_num_cpus, scheduler='processes'):
        """A Dask multiprocessing evaluator

        Args:
            evaluator (function): Evaluation function.
            mp_num_cpus (int): Number of CPUs.
        """
        self.evaluator = evaluator
        self.scheduler = scheduler
        self.__name__ = self.__class__.__name__

    def evaluate(self, candidates, args):
        with dask.config.set(scheduler=self.scheduler):
            return list(dask.compute(*[dask.delayed(self.evaluator)(c) for c in candidates]))

    def __call__(self, candidates, args):
        return self.evaluate(candidates, args)


class SparkEvaluator(Evaluator):
    def __init__(self, evaluator, mp_num_cpus):
        """A Spark multiprocessing evaluator

        Args:
            evaluator (function): Evaluation function.
            mp_num_cpus (int): Number of CPUs.
        """
        self.evaluator = evaluator
        self.spark_conf = SparkConf().setAppName(
            "mewpy").setMaster(f"local[{mp_num_cpus}]")
        self.spark_context = SparkContext(conf=self.spark_conf)
        self.__name__ = self.__class__.__name__

    def evaluate(self, candidates, args):
        solutions_to_evaluate = self.spark_context.parallelize(candidates)
        return solutions_to_evaluate.map(lambda s: self.evaluator(s)).collect()

    def __call__(self, candidates, args):
        return self.evaluate(candidates, args)


# ray
try:
    import ray
    MP_Evaluators.append('ray')
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
            if getattr(self.problem, "evaluator", False):
                return self.problem.evaluator(candidates, None)
            else:
                res = []
                for candidate in candidates:
                    res.append(self.problem.evaluate(candidate))
                return res

    @ray.remote
    class RayActorF:
        """
        Each actor (worker) has a solver instance to overcome the need to serialize solvers which may not be pickable.
        The solver is not reset before each evaluation.
        """

        def __init__(self, func):
            self.func = copy.deepcopy(func)

        def evaluate_candidates(self, candidates):
            """Evaluates a sublist of candidates
            """
            res = []
            for candidate in candidates:
                res.append(self.func(candidate))
            return res

    class RayEvaluator(Evaluator):
        def __init__(self, problem, number_of_actors, isfunc=False):
            """A ray actor responsible for performing evaluations.

            Args:
                problem: A class implementing an evaluator(list_of_candidates,**kwargs)
                number_of_actors (int): Number of workers
            """
            ray.init(ignore_reinit_error=True)
            if isfunc:
                self.actors = [RayActorF.remote(problem)
                               for _ in range(number_of_actors)]
            else:
                self.actors = [RayActor.remote(problem)
                               for _ in range(number_of_actors)]
            self.number_of_actors = len(self.actors)
            self.__name__ = self.__class__.__name__
            print(f"Using {self.number_of_actors} workers.")

        def evaluate(self, candidates, args):
            """
            Divides the candidates into sublists to be evaluated by each actor
            """
            size = len(candidates) // self.number_of_actors
            n = len(candidates) % self.number_of_actors
            sub_lists = []
            p = 0
            for _ in range(self.number_of_actors):
                d = size + 1 if n > 0 else size
                sub_lists.append(candidates[p:p + d])
                p = p + d
                n -= 1
            values = []
            for i in range(self.number_of_actors):
                actor = self.actors[i]
                values.append(actor.evaluate_candidates.remote(sub_lists[i]))
            r = ray.get(values)
            result = []
            for x in r:
                result += x
            return result

        def __call__(self, candidates, args):
            return self.evaluate(candidates, args)


def get_mp_evaluators():
    """"Returns the list of available multiprocessing evaluators.
    """
    return MP_Evaluators


def get_evaluator(problem, n_mp=cpu_count(), evaluator=ModelConstants.MP_EVALUATOR):
    """Retuns a multiprocessing evaluator

    Args:
        problem: a class implementing an evaluate(candidate) function
        n_mp (int, optional): The number of cpus. Defaults to cpu_count().
        evaluator (str, optional): The evaluator name: options 'ray','dask','spark'.\
            Defaults to ModelConstants.MP_EVALUATOR.

    Returns:
        [type]: [description]
    """
    if evaluator == 'ray' and 'ray' in MP_Evaluators:
        return RayEvaluator(problem, n_mp)
    elif evaluator == 'nodaemon' and 'nodaemon' in MP_Evaluators:
        return NoDaemonMultiProcessorEvaluator(problem,n_mp) 
    elif evaluator == 'dask' and 'dask' in MP_Evaluators:
        return DaskEvaluator(problem.evaluate, n_mp)
    elif evaluator == 'spark' and 'spark' in MP_Evaluators:
        return SparkEvaluator(problem.evaluate, n_mp)
    else:
        return MultiProcessorEvaluator(problem.evaluate, n_mp)


def get_fevaluator(func, n_mp=cpu_count(), evaluator=ModelConstants.MP_EVALUATOR):
    """Retuns a multiprocessing evaluator

    Args:
        problem: a class implementing an evaluate(candidate) function
        n_mp (int, optional): The number of cpus. Defaults to cpu_count().
        evaluator (str, optional): The evaluator name: options 'ray','dask','spark'.\
            Defaults to ModelConstants.MP_EVALUATOR.

    Returns:
        [type]: [description]
    """
    if evaluator == 'ray' and 'ray' in MP_Evaluators:
        return RayEvaluator(func, n_mp, isfunc=True)
    elif evaluator == 'nodaemon' and 'nodaemon' in MP_Evaluators:
        return NoDaemonMultiProcessorEvaluator(func,n_mp) 
    elif evaluator == 'dask' and 'dask' in MP_Evaluators:
        return DaskEvaluator(func, n_mp)
    elif evaluator == 'spark' and 'spark' in MP_Evaluators:
        return SparkEvaluator(func, n_mp)
    else:
        return MultiProcessorEvaluator(func, n_mp)
