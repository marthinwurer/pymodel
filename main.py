




import matplotlib.pyplot as plt
import numpy as np
from math import *
from matplotlib import animation
import time

radius = 6371000.0
gravity = 9.8
R_dry = 287.1
cp_dry = 1004.0
kappa = 0.286
# whooo globals
#numx = 80
#numy = 40
numx = 360
numy = 180
timestep = 900.0



def dynam(u, v, p, p_next, h, t, dx, dxc, dy, dt):

    """
    Find the momentum gradient at the edges for the change in pressure
    """
    du = np.zeros((numy, numx))
    dv = np.zeros((numy, numx))
    spa = np.zeros((numy, numx))
    maxp = p[0][0]
    minp = maxp
    maxu = u[0][0]
    minu = maxu
    maxc = 0
    minc = 0
    maxconv = 0
    maxg = 0
    for y in range(numy - 1): # don't do the last row, north-south there should be zero
        for x in range(-1, numx - 1):
            maxu = max(maxu, u[y][x])
            minu = min(minu, u[y][x])
            du[y][x] = (u[y-1][x] + u[y][x]) * (p[y][x] + p[y][x+1])
            dv[y][x] = (v[y][x-1] + v[y][x]) * (p[y][x] + p[y+1][x])

    # find the convergence and change the pressure
    # these look to the north and west, not south and east
    for y in range(numy): 
        for x in range(numx ):
            maxp = max(maxp, p[y][x])
            minp = min(minp, p[y][x])
            convergence = 0
            convergence += du[y][x-1] 
            convergence -= du[y][x] 
            convergence += dv[y-1][x] 
            convergence -= dv[y][x]

            change =  (dt * convergence /(dx[y] * dy * radius))

            maxconv = max(maxconv, change)

            p_next[y][x] = p[y][x] + change 


    # compute spa now
    for y in range(numy): 
        for x in range(numx ):
            spa[y][x] = p[y][x]**kappa * t[y][x] * R_dry


    # do geopotential now
    for y in range(-1, numy - 1): # don't do the last row, north-south there should be zero
        for x in range(-1, numx - 1):
            p_east = (p[y][x+1] + p[y+1][x+1]) / 2
            p_west = (p[y][x] + p[y+1][x]) / 2
            g_east = (h[y][x+1] + h[y+1][x+1]) / 2
            g_west = (h[y][x] + h[y+1][x]) / 2
            #du[y][x] = (p_east + p_west) * (g_west - g_east)
            change = (g_west - g_east) / dxc[y]
            maxg = max(maxg, change)
            du[y][x] = change

            #dv[y][x] = (v[y][x-1] + v[y][x]) * (p[y][x] + p[y+1][x])


            # do the pressure gradient force
            #spa_east = (spa[y][x+1] + spa[y+1][x+1]) / 2
            #spa_west = (spa[y][x] + spa[y+1][x]) / 2
            #change = (spa_west - spa_east) / dx[y]
            change = (p_west - p_east) / dxc[y]
            maxc = max(maxc, change)
            minc = min(minc, change)
            du[y][x] += change


    # apply the change in velocity
    u += du * dt

    print(maxu, maxp, maxc, maxconv)
    print(minu, minp, minc, maxg)


            


def main():

    # max and min stuff

    # make the arrays
    dx = np.zeros(numy)
    dxc = np.zeros(numy)
    #height = np.random.rand(numy, numx) * np.random.rand(numy, numx) * 200
    height = np.zeros((numy, numx))
    u = np.zeros((numy, numx))
    v = np.zeros((numy, numx))
    pressure = np.zeros((numy, numx))
    p_next = np.zeros((numy, numx))
    temp = np.full((numy, numx), 273.15) # potential temperature

    
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
        dx[ii] = cos(lat_cur) * dxe * radius
        print(lat_cur/pi * 180, dx[ii], dx[ii])

    # do the corners next
    # this finds the latitude of the first center from the top
    lat = pi / 2 - dye

    # now loop and do the calculations
    for ii in range(numy):
        lat_cur = lat - dye * ii
        dxc[ii] = cos(lat_cur) * dxe * radius
        print(lat_cur/pi * 180, dxc[ii], dxc[ii])


    # make a random peak 
    height[numy // 3, numx // 3] = 200

    # make an initial push
    u[numy // 2, numx // 2] = 10


    # geometry has been calculated. set up the initial pressure
    for y in range(numy):
        for x in range(numx):
            temp[y][x] = 273.15
            pressure[y][x] = 1013.25 / exp( gravity * height[y][x] / ( R_dry * temp[y][x]))

        
    fig,ax = plt.subplots(1,1)
    image = ax.imshow(u, cmap='gray')
    fig.canvas.draw()
    plt.show(block=False)
    for i in range(100):
        print("frame", i)
        # run dynam once 
        dynam(u, v, pressure, p_next, height, temp, dx, dxc, dye, timestep)
        pressure = p_next


        image.set_data(pressure)
        fig.canvas.draw()
        plt.draw()
        #time.sleep(5)

        #plt.show()

    











if __name__ == "__main__":
    main()

