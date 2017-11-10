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
