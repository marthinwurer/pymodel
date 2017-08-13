




import matplotlib.pyplot as plt
import numpy as np
from math import *

radius = 6371000.0
gravity = 9.8
R_dry = 287.1
cp_dry = 1004.0
kappa = 0.286
# whooo globals
numx = 360
numy = 180
timestep = 900.0

def dynam(u, v, p, p_next, h, dx, dy, dt):
    """
    Find the momentum gradient at the edges for the change in pressure
    """
    du = np.zeros((numy, numx))
    dv = np.zeros((numy, numx))
    for y in range(numy - 1): # don't do the last row, north-south there should be zero
        for x in range(-1, numx - 1):
            du[y][x] = (u[y-1][x] + u[y][x]) * (p[y][x] + p[y][x+1])
            dv[y][x] = (v[y][x-1] + v[y][x]) * (p[y][x] + p[y+1][x])

    # find the convergence and change the pressure
    # these look to the north and west, not south and east
    for y in range(numy): 
        for x in range(numx ):
            convergence = 0
            convergence += du[y][x-1] 
            convergence -= du[y][x] 
            convergence += dv[y-1][x] 
            convergence -= dv[y][x]

            change =  (dt * convergence /(dx[y] * dy * radius))

            if change != 0.0:
                print(change)

            p_next[y][x] = p[y][x] + change 

    # do geopotential now
    for y in range(numy - 1): # don't do the last row, north-south there should be zero
        for x in range(-1, numx - 1):
            du[y][x] = (u[y-1][x] + u[y][x]) * (p[y][x] + p[y][x+1])

            dv[y][x] = (v[y][x-1] + v[y][x]) * (p[y][x] + p[y+1][x])



            


def main():

    # max and min stuff

    # make the arrays
    dx = np.zeros(numy)
    dxc = np.zeros(numy)
    height = np.random.rand(numy, numx) * np.random.rand(numy, numx) * 200
    u = np.zeros((numy, numx))
    v = np.zeros((numy, numx))
    pressure = np.zeros((numy, numx))
    p_next = np.zeros((numy, numx))
    temp = np.zeros((numy, numx)) # potential temperature

    
    # calculate geometry
    # find the radians for each height difference: 180deg = pi rads
    # so find dx at the equator by doing 2pi / numx
    dxe = 2 * pi / numx
    # dy is 180deg so pi / numy
    dye = pi / numy

    # we want to find the differences between the centers as well
    # as the SE corners, so two different arrays are needed
    # do this by finding the starting value (near the north pole)
    # and counting down

    # do the centers first
    # this finds the latitude of the first center from the top
    lat = pi / 2 - dye / 2

    # now loop and do the calculations
    for ii in range(numy):
        lat_cur = lat - dye * ii
        dx[ii] = cos(lat_cur) * dxe
        print(lat_cur/pi * 180, dx[ii], dx[ii] * radius)

    # do the corners next
    # this finds the latitude of the first center from the top
    lat = pi / 2 - dye

    # now loop and do the calculations
    for ii in range(numy):
        lat_cur = lat - dye * ii
        dxc[ii] = cos(lat_cur) * dxe
        print(lat_cur/pi * 180, dxc[ii], dxc[ii] * radius)


    # make a random peak 
    height[numy // 3, numx // 3] = 2000

    # make an initial push
    u[numy // 2, numx // 2] = 10


    # geometry has been calculated. set up the initial pressure
    for y in range(numy):
        for x in range(numx):
            temp[y][x] = 273.15
            pressure[y][x] = 1013.25 / exp( gravity * height[y][x] / ( R_dry * temp[y][x]))

    # run dynam once 
    dynam(u, v, pressure, p_next, height, dx, dye, timestep)

    plt.imshow(p_next)
    plt.gray()
    plt.show()


    











if __name__ == "__main__":
    main()

