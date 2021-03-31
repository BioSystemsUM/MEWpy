# TODO: this module largely depends on pandas dataframes. Should it be set as package requirement?
# noinspection PyPackageRequirements
from pandas import concat


# TODO: methods stubs and type hinting
# TODO: maybe the following solution objects should inherit from the same multisolution object!
class MultiSolution:

    def __init__(self, *solutions):

        if not solutions:
            _solutions = {}
        else:
            _solutions = {}

            for solution in solutions:
                setattr(self, f'{solution.method}', solution)
                _solutions[solution.method] = solution

        self._solutions = _solutions

    @property
    def solutions(self):
        return self._solutions

    def to_frame(self):
        frames = []
        columns = []

        for method, solution in self._solutions.items():
            frames.append(solution.to_frame().frame)
            columns.append(method)

        df = concat(frames, axis=1, join='outer', keys=columns)

        return df

    def to_summary(self):

        frames = []
        columns = []

        for method, solution in self._solutions.items():
            frames.append(solution.to_summary().frame)
            columns.append(method)

        df = concat(frames, axis=1, join='outer', keys=columns)

        return df

    def __repr__(self):
        return 'MultiSolution'

    def __str__(self):
        return 'MultiSolution:' ','.join([method for method in self._solutions])


# TODO: methods stubs and type hinting
class DynamicSolution:

    def __init__(self, *solutions, time=None):

        if not solutions:
            _solutions = {}
            time = []

        else:
            _solutions = {}

            if time is None:
                time = [i for i in range(len(solutions))]

            else:
                time = list(time)

            for t, solution in zip(time, solutions):
                setattr(self, f't_{t}', solution)
                _solutions[f't_{t}'] = solution

        self._solutions = _solutions
        self._time = time

    @property
    def solutions(self):
        return self._solutions

    def to_frame(self):
        frames = []
        columns = []

        for time, solution in self._solutions.items():
            frames.append(solution.to_frame().frame)
            columns.append(time)

        df = concat(frames, axis=1, join='outer', keys=columns)

        return df

    def to_summary(self):

        frames = []
        columns = []

        for time, solution in self._solutions.items():
            frames.append(solution.to_summary().frame)
            columns.append(time)

        df = concat(frames, axis=1, join='outer', keys=columns)

        return df

    def __repr__(self):
        return 'DynamicSolution'

    def __str__(self):
        return 'DynamicSolution:' ','.join([time for time in self._solutions])


# TODO: methods stubs and type hinting
class KOSolution:

    def __init__(self, solutions, kos=None):

        if not solutions:
            _solutions = {}
            kos = []

        else:
            _solutions = {}

            if not kos:

                if isinstance(solutions, dict):

                    kos = list(solutions.keys())
                    solutions = list(solutions.values())

                else:

                    kos = [i for i, _ in enumerate(solutions)]

            else:

                if isinstance(solutions, dict):
                    solutions = list(solutions.values())

            for ko, solution in zip(kos, solutions):
                setattr(self, f'ko_{ko}', solution)
                _solutions[f'ko_{ko}'] = solution

        self._solutions = _solutions
        self._time = kos

    @property
    def solutions(self):
        return self._solutions

    def to_frame(self):
        frames = []
        columns = []

        for ko, solution in self._solutions.items():
            frames.append(solution.to_frame().frame)
            columns.append(ko)

        df = concat(frames, axis=1, join='outer', keys=columns)

        return df

    def to_summary(self):

        frames = []
        columns = []

        for ko, solution in self._solutions.items():
            frames.append(solution.to_summary().frame)
            columns.append(ko)

        df = concat(frames, axis=1, join='outer', keys=columns)

        return df

    def __repr__(self):
        return 'KOSolution'

    def __str__(self):
        return 'KOSolution:' ','.join([time for time in self._solutions])
