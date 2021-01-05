import numpy as np
from numba import njit


@njit(cache=True)
def mann_kenddall_test_calculate_s(x):
    """ Calculate the sum of signs

    Args:
        x (np.array): input data

    Returns:
        int : sum of signs
    """
    n = len(x)
    # calculate S
    s = 0
    for k in range(n - 1):
        # s += np.nansum(np.sign(x[k + 1:] - x[k]))
        for j in range(k + 1, n):
            diff = x[j] - x[k]
            if diff == diff:
                if diff > 0:
                    s += 1
                elif diff < 0:
                    s -= 1
                else:
                    pass
                # s += np.sign(x[j] - x[k])
    return s


@njit(cache=True)
def in1di(x, y):
    out = np.empty(x.shape, dtype=np.int32)
    out.fill(0)
    for i in range(x.shape[0]):
        for j in range(y.shape[0]):
            if x[i] == y[j]:
                out[i] = j
                break

    return out
