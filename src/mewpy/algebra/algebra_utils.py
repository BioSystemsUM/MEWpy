def solution_decode(solution, decoder=None):
    if not decoder:
        decoder = {True: 1,
                   False: 0,
                   1: 1,
                   0: 0,
                   -1: 1,
                   -2: 0
                   }

    return decoder.get(solution, solution)
