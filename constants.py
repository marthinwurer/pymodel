
radius = 6371000.0
#radius = 637100.0
gravity = 9.8
R_dry = 287.1
cp_dry = 1004.0
kappa = 0.286
ptop = 10.0 # the pressure at the top of the atmosphere
r_rate_earth = 7.292115e-5 # rads/sec, rotation rate of earth (sidereal day)
#r_rate_earth = 0.0
coriolis_m = 2 * r_rate_earth # will need to be multiplied by the sin of the latitude
sigma = [0.5] # proportion of pressure at that layer
dsig = [1.0] # change from this pressure layer to the one above it
# whooo globals
numx = 8
numy = 4
# numx = 36
# numy = 24
lay = 0
# numx = 72
# numy = 36
# numx = 360
# numy = 180
timestep = 900.0
# timestep = 120.0
# timestep = 20.0
#torroid = True
torroid = False


def rmul(a, b):
    """
    row multiply
    Args:
        a:
        b:

    Returns:

    """
    return (a.T * b).T