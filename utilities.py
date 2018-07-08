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

def edge_averages(a: np.ndarray):
    """
    Finds the average value of the array on the south and east edges of the tile
    Args:
        a:

    Returns:

    """
    south = np.zeros(a.shape)
    east = np.zeros(a.shape)
    for yy in range(-1, a.shape[0] - 1):
        for xx in range(-1, a.shape[1] - 1):
            south[yy, xx] = (a[yy][xx] + a[yy+1][xx]) / 2
            east[yy, xx] = (a[yy][xx] + a[yy][xx+1]) / 2

    return np.asarray([east, south])

def gradient_c(a: np.ndarray, dx: float=1.0):
    """
    Now finds the gradient on the south face and the east face of the tile
    Args:
        a:
        dx:
    Returns:
    """
    south = np.zeros(a.shape)
    east = np.zeros(a.shape)
    for yy in range(-1, a.shape[0] - 1):
        for xx in range(-1, a.shape[1] - 1):
            south[yy, xx] = (a[yy][xx] - a[yy+1][xx])
            east[yy, xx] = (a[yy][xx] - a[yy][xx+1])

    return np.asarray([east, south]) / dx

def gradient_a(a: np.ndarray, dx: float=1.0):
    """
    Finds the gradient for arakawa's A scheme
    Args:
        a:
        dx:
    Returns:
    """
    south = np.zeros(a.shape)
    east = np.zeros(a.shape)

    for yy in range(-1, a.shape[0] - 1):
        for xx in range(-1, a.shape[1] - 1):
            south[yy, xx] = (a[yy-1][xx] - a[yy+1][xx]) / (2 * dx)
            east[yy, xx] = (a[yy][xx-1] - a[yy][xx+1]) / (2 * dx)

    return [east, south]

gradient = gradient_c


def divergence_c(u: np.ndarray, dx: float=1.0):
    """
    Given a gradient, finds the divergence of a tile.
    Add up the gradients on all edges of the tile
    Args:
        u:
        v:
        dx:

    Returns:

    """
    div = np.zeros(u.shape[1:])
    grad_u = u[0]
    grad_v = u[1]
    for yy in range(u.shape[1]):
        for xx in range(u.shape[2]):
            north = grad_v[yy-1][xx]
            south = grad_v[yy][xx]
            east = grad_u[yy][xx]
            west = grad_u[yy][xx-1]
            #     west                     east             north              south
            # div[yy][xx] = grad_u[yy][xx-1] - grad_u[yy][xx] + grad_v[yy-1][xx] - grad_v[yy][xx]
            div[yy][xx] = west - east + north - south

    return div

def divergence_a(u: np.ndarray, dx: float=1.0):
    div = np.zeros(u.shape[1:])
    grad_u = u[0]
    grad_v = u[1]
    for yy in range(-1, u.shape[0] - 1):
        for xx in range(-1, u.shape[1] - 1):
            north = grad_v[yy-1][xx]
            south = grad_v[yy+1][xx]
            east = grad_u[yy][xx+1]
            west = grad_u[yy][xx-1]
            div[yy][xx] = west + east + north + south
    return div

divergence = divergence_c




