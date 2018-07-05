import numpy as np

def reast(a):
    """Roll an array to the east
    Args:
        a:

    Returns:

    """
    return np.roll(a, -1, 1)

def rwest(a):
    return np.roll(a, 1, 1)

def rnorth(a):
    return np.roll(a, 1, 0)

def rsouth(a):
    return np.roll(a, -1, 0)

def rmul(a, b):
    """
    row multiply
    Args:
        a:
        b:

    Returns:

    """
    return (a.T * b).T

def gradient(a: np.ndarray, dx):
    north = np.zeros(a.shape)
    south = np.zeros(a.shape)
    east = np.zeros(a.shape)
    west = np.zeros(a.shape)

    for yy in range(-1, a.shape[0] - 1):
        for xx in range(-1, a.shape[1] - 1):
            north[yy, xx] = (a[yy][xx] - a[yy-1][xx]) / dx
            south[yy, xx] = (a[yy][xx] - a[yy+1][xx]) / dx
            east[yy, xx] = (a[yy][xx] - a[yy][xx+1]) / dx
            west[yy, xx] = (a[yy][xx] - a[yy][xx-1]) / dx

    return (north, south, east, west)




