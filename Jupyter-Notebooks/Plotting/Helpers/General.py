import numpy as np
def LogTimePoints(Start, End, N = 1e5, IncludeZero = True, IncludePowersof10 = True):
    '''
    Generates log spaced points for use as plot axes.
    
    Start, End, N = 1e5, IncludeZero = True, IncludePowersof10 = True
    '''
    # Overkill method for finding what powers of 10 are missing
    Time = np.logspace(Start,End,num=N)

    if (0.0 not in Time) and IncludeZero:
        Time = np.hstack((0.0,Time))

    if IncludePowersof10:
        Add = 10.0**np.arange(Start,End+1)[
            np.array(
                [(10.0**i not in Time) for i in range(Start,End+1)]
                )
            ]
        if len(Add) > 0:
            Time = np.sort(np.hstack([Time,Add]))
    return Time