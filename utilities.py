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

def gradient(a: np.ndarray, dx: float=1.0):
    """
    Now finds the gradient on the south face and the east face of the tile
    Args:
        a:
        dx:
    Returns:
    """
    # north = np.zeros(a.shape)
    south = np.zeros(a.shape)
    east = np.zeros(a.shape)
    # west = np.zeros(a.shape)

    for yy in range(-1, a.shape[0] - 1):
        for xx in range(-1, a.shape[1] - 1):
            # north[yy, xx] = (a[yy][xx] - a[yy-1][xx]) / dx
            south[yy, xx] = (a[yy][xx] - a[yy+1][xx]) / dx
            east[yy, xx] = (a[yy][xx] - a[yy][xx+1]) / dx
            # west[yy, xx] = (a[yy][xx] - a[yy][xx-1]) / dx

    # return [north, south, east, west]
    return [south, east]

def divergence(u: np.ndarray, dx: float=1.0):
    """
    Add up the gradients on all edges of the tile
    Args:
        u:
        v:
        dx:

    Returns:

    """
    div = np.zeros(u.shape[1:])
    grad_u = gradient(u[1], dx)[1]
    grad_v = gradient(u[0], dx)[0]
    grad_u = u[1]
    grad_v = u[0]
    print("U:")
    print(grad_u)
    print("V:")
    print(grad_v)
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




