




import matplotlib.pyplot as plt
import numpy as np
from math import *
from matplotlib import animation
import time
import traceback

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
numx = 36
numy = 24
lay = 0
#numx = 72
#numy = 36
#numx = 360
#numy = 180
timestep = 900.0
#timestep = 120.0
#torroid = True
torroid = False



def dynam(u, v, un, vn, du, dv, dt):

    # apply the change in velocity
    np.copyto(un, u)
    np.copyto(vn, v)
    un += du * dt
    vn += dv * dt

def aflux(u, v, p, pu, pv):
    for y in range(-1,numy - 1): # don't do the last row, north-south there should be zero
        for x in range(-1, numx - 1):
            #maxu = max(maxu, u[y][x])
            #minu = min(minu, u[y][x])
            pu[y][x] = (u[y-1][x] + u[y][x]) * (p[y][x] + p[y][x+1]) / 4
            pv[y][x] = (v[y][x-1] + v[y][x]) * (p[y][x] + p[y+1][x]) / 4

            if y == -1:
                u[y][x] = 0
                v[y][x] = 0
                pv[y][x] = 0


def advecm(p, pu, pv, p_next, convergence, dx, dym, dt):
    # find the convergence and change the pressure
    # these look to the north and west, not south and east
    for y in range(numy): 
        for x in range(numx ):
            #maxp = max(maxp, p[y][x])
            #minp = min(minp, p[y][x])
            #convergence = 0
            #convergence += du[y][x-1] 
            #convergence -= du[y][x] 
            #convergence += dv[y-1][x] 
            #convergence -= dv[y][x]

            #change =  (dt * convergence /(dx[y] * dy * radius))

            ew = (pu[y][x-1] - pu[y][x]) / dx[y]
            ns = (pv[y-1][x] - pv[y][x]) / dym
            change = ew + ns
            conv =(pu[y][x-1] - pu[y][x] + pv[y-1][x] - pv[y][x])

            #maxconv = max(maxconv, change)

            convergence[y][x] = change

            p_next[y][x] = p[y][x] + change * dt
            #p_next[y][x] = p[y][x] + (conv * dt / (dx[y] * dym))


def advectracer(pu, pv, tracer, tracer_next, dx, dym, dt):
    tracer_next.fill(0.0)
    # advect the tracer
    for y in range(-1, numy - 1): 
        for x in range(-1, numx - 1):
            # compute the amount moved in each direction.
            tr_east = pu[y][x] * tracer[y][x] / dx[y] * dt
            tr_west = -pu[y][x-1] * tracer[y][x] / dx[y] * dt
            tr_north = -pv[y][x] * tracer[y][x] / dx[y] * dt
            tr_south = pv[y-1][x] * tracer[y][x] / dx[y] * dt

            # make sure that you only do the ones that are positive
            # and make sure that you only move 1/4 of the tracer in any direction
            maxmove = tracer[y][x] / 4
            tr_east = max(0, min(tr_east, maxmove))
            tr_west = max(0, min(tr_west, maxmove))
            tr_north = max(0, min(tr_north, maxmove))
            tr_south = max(0, min(tr_south, maxmove))

            # compute how much will be left
            left = tracer[y][x] - tr_east - tr_west - tr_north - tr_south

            # move the tracer
            tracer_next[y][x] += left
            tracer_next[y-1][x] += tr_north
            tracer_next[y+1][x] += tr_south
            tracer_next[y][x-1] += tr_west
            tracer_next[y][x+1] += tr_east


def pgf(du, dv, p, p_c, h, t, spa, dxc, dym):
    ## compute spa now
    for y in range(numy): 
        for x in range(numx ):
            spa[y][x] = (p[y][x] + ptop)**kappa * t[y][x] * R_dry


    # do geopotential now
    for y in range(-1, numy - 1): 
        for x in range(-1, numx - 1):
            p_east = (p[y][x+1] + p[y+1][x+1]) / 2
            p_west = (p[y][x] + p[y+1][x]) / 2
            p_south = (p[y+1][x] + p[y+1][x+1]) / 2
            p_north = (p[y][x] + p[y][x+1]) / 2
            g_east = (h[y][x+1] + h[y+1][x+1]) / 2
            g_west = (h[y][x] + h[y+1][x]) / 2
            g_south = (h[y+1][x] + h[y+1][x+1]) / 2
            g_north = (h[y][x] + h[y][x+1]) / 2
            #du[y][x] = (p_east + p_west) * (g_west - g_east)
            change = (g_west - g_east) / dxc[y]
            #maxg = max(maxg, change)
            du[y][x] = change
            dv[y][x] = (g_north - g_south) / dym

            #dv[y][x] = (v[y][x-1] + v[y][x]) * (p[y][x] + p[y+1][x])


            # do the pressure gradient force
            spa_east = (spa[y][x+1] + spa[y+1][x+1]) / 2
            spa_west = (spa[y][x] + spa[y+1][x]) / 2
            spa_north = (spa[y][x] + spa[y][x+1]) / 2
            spa_south = (spa[y+1][x] + spa[y+1][x+1]) / 2
            #change = (spa_west - spa_east) / dx[y]
            p_center = ((p_west + p_east) / 2)
            p_c[y][x] = p_center
            change = (spa_west - spa_east) / dxc[y] / p_center
            #maxc = max(maxc, change)
            #minc = min(minc, change)
            du[y][x] += change
            dv[y][x] += (spa_north - spa_south) / dym / p_center



def advecv(u, v, du, dv, p, dxc, dym):
    # do advection of momentum
    for y in range(-1, numy - 1): 
        for x in range(-1, numx - 1):

            # eastward: find the pressure and velocity to the east 
            p_east = (p[y][x+1] + p[y+1][x+1]) / 2
            u_east = (u[y][x] + u[y][x+1]) / 2
            v_east = (v[y][x] + v[y][x+1]) / 2

            flux = p_east * u_east

            fluxu = flux * u_east / dxc[y] / p_east
            du[y][x] -= fluxu
            du[y][x+1] += fluxu

            fluxv = flux * v_east / dxc[y] / p_east
            dv[y][x] -= fluxv
            dv[y][x+1] += fluxv

            # southward
            p_south = (p[y+1][x] + p[y+1][x+1]) / 2
            v_south = (v[y][x] + v[y+1][x]) / 2
            u_south = (u[y][x] + u[y+1][x]) / 2

            flux = p_south * v_south

            fluxu = flux * u_south / dym / p_south
            du[y][x] -= fluxu
            du[y+1][x] += fluxu

            fluxv = flux * v_south / dym / p_south
            dv[y][x] -= fluxv
            dv[y+1][x] += fluxv


def coriolis(u, v, du, dv, latc):

    for y in range(numy): 
        for x in range(numx ):
            f = coriolis_m * sin(latc[y])
            du[y][x] += v[y][x] * f
            dv[y][x] -= u[y][x] * f


def friction(u, v, du, dv, p, spa, dxc, dym, dt):
    vel = ((u ** 2 + v ** 2) ** (1/2))
    drag_coef_flat = 0.005
    for y in range(numy): 
        for x in range(numx):
            drag = drag_coef_flat * (vel[y][x] ** 2)# * spa[y][x]
            if drag > vel[y][x]:
                coeff = 1
            else:
                coeff = (drag / (vel[y][x] + 0.000001))
            if x == 0:
                print(coeff, u[y][x] * coeff, v[y][x] * coeff)
            du[y][x] -= (u[y][x] * coeff) / dt
            dv[y][x] -= (v[y][x] * coeff) / dt


def rad(t, tp1, lat, lon, time, dt):
    # find the center of sunlight
    # total time in seconds % day length to find day seconds
    # day seconds / seconds in day to find longitude proportion
    # multiply by 2pi to get angle around earth
    ds = time % 86400
    cx = ds / 86400.0 * 2 * pi
    
    hr = 0

    for y in range(numy): 
        for x in range(numx):
            tlon = x * lon - cx
            inc = cos(lat[y]) * cos(tlon) * hr
            tp1[y][x] = t[y][x] + inc * dt


def sdrag(u, v, p, p_center, t, dt):
    """
    Stratospheric drag (ported from GCMII)
    """
    #c1 = 0.0005
    c1 = 0.0005
    #c2 = 0.00005
    c2 = 50
    vel = ((u ** 2 + v ** 2) ** (1/2))
    maxp = 0.0
    for y in range(numy): 
        for x in range(numx):
            pk = (p[y][x] + ptop) ** kappa
            rho = ptop / (R_dry * t[y][x] * pk) # density at top of atmosphere
            cdn = c1 + c2 * vel[y][x]
            prop = dt * rho * cdn * vel[y][x] * gravity / (p_center[y][x])
            maxp = max(maxp, prop)
            prop = 1.0 - prop
            if prop < 0:
                prop = 0.0
            u[y][x] *= prop
            v[y][x] *= prop
    print("maxp", maxp)
                




    




            


def main():
    # u is eastward
    # v is southward

    # max and min stuff
    np.seterr(invalid='raise')
    np.set_printoptions(threshold=np.inf)


    # make the arrays
    # constants
    dx = np.zeros(numy)
    dxc = np.zeros(numy)
    lat = np.zeros(numy)
    latc = np.zeros(numy)
    height = np.zeros((numy, numx))
    #height = np.random.rand(numy, numx) * np.random.rand(numy, numx)# * 10

    # dynamic data: m1 = u-dt, _ = ut, p1 = u+dt
    um1 = np.zeros((numy, numx))
    vm1 = np.zeros((numy, numx))
    u = np.zeros((numy, numx))
    v = np.zeros((numy, numx))
    up1 = np.zeros((numy, numx))
    vp1 = np.zeros((numy, numx))

    pm1 = np.zeros((numy, numx))
    p = np.zeros((numy, numx))
    pp1 = np.zeros((numy, numx))

    tm1 = np.full((numy, numx), 273.15) # potential temperature
    t = np.full((numy, numx), 273.15) # potential temperature
    tp1  = np.full((numy, numx), 273.15) # potential temperature

    tracer_m1 = np.full((numy, numx), 10.0) 
    tracer = np.full((numy, numx), 10.0) 
    tracer_p1 = np.full((numy, numx), 10.0) 
    
    # working data
    p_center = np.zeros((numy, numx))
    convergence = np.zeros((numy, numx))
    pu = np.zeros((numy, numx)) # momentum
    pv = np.zeros((numy, numx))
    du = np.zeros((numy, numx)) # change in velocity
    dv = np.zeros((numy, numx)) 
    spa = np.zeros((numy, numx))


    
    # calculate geometry
    # find the radians for each height difference: 180deg = pi rads
    # so find dx at the equator by doing 2pi / numx
    dxe = 2 * pi / numx
    # dy is 180deg so pi / numy
    dye = pi / numy
    # find it in meters
    dym = dye * radius

    # we want to find the differences between the centers as well
    # as the SE corners, so two different arrays are needed
    # do this by finding the starting value (near the north pole)
    # and counting down

    # do the centers first
    # this finds the latitude of the first center from the top
    maxlat = pi / 2 - dye / 2

    if not torroid:
        # now loop and do the calculations
        for ii in range(numy):
            lat_cur = maxlat - dye * ii
            lat[ii] = lat_cur
            dx[ii] = cos(lat_cur) * dxe * radius
        #    print(lat_cur/pi * 180, dx[ii], dx[ii])
        dx[-1] = dx[1] # prevent explosions
    else:
        # torroid
        for ii in range(numy):
            dx[ii] = dxe * radius
    print(lat)
    print(dx)


    # do the corners next
    # this finds the latitude of the first center from the top
    maxlat = pi / 2 - dye

    if not torroid:
        # now loop and do the calculations
        for ii in range(numy):
            lat_cur = maxlat - dye * ii
            latc[ii] = lat_cur
            dxc[ii] = cos(lat_cur) * dxe * radius
        #    print(lat_cur/pi * 180, dxc[ii], dxc[ii])
        dxc[-1] = dxc[1] # prevent explosions
        latc[-1] = 0
    else:
        # torroid
        for ii in range(numy):
            dxc[ii] = dxe * radius
    print(latc)
    print(dxc)


    # make a random peak 
    height[numy // 3, numx // 3] = 200
    #height[numy // 4 * 3, numx // 3] = 2

    # make walls at the north and south as well as the west
    #for i in range(numx):
    #    height[0][i] = 2000
    #    height[-1][i] = 2000
    #for i in range(numy):
    #    height[i][0] = 2000

    # make an initial push
    #u[numy // 3, numx // 2] = -10
    #u[numy // 3, numx // 3] = 10
    u[numy // 3 * 2, numx // 3] = -10
    v[numy // 3, numx // 2] = 10

    #for i in range(numx):
    #    u[numy // 2][i] = 10
    #    #u[numy // 2 + 1][i] = 10
    #u[numy // 2][9] = 9
    


    # geometry has been calculated. set up the initial pressure
    for y in range(numy):
        for x in range(numx):
            p[y][x] = 1013.25 / exp( gravity * height[y][x] / ( R_dry * t[y][x])) - ptop
            pp1[y][x] = p[y][x]

    pp1[0][0] = 1100.0
    pp1[0][1] = 900.0

    #p[numy // 3 * 2, numx // 3] -= 10

    tracer_p1[0][0] = 25
    tracer_p1[0][1] = 0

    dt = timestep
    dt2 = dt * 2

    #((u ** 2 + v ** 2) ** (1/2))

    fig,ax = plt.subplots(1,1)
    #image = ax.imshow(((u ** 2 + v ** 2) ** (1/2)), cmap='gray')
    #image = ax.imshow(tracer_p1, cmap='gray')
    image = ax.imshow(pp1, cmap='gray')
    fig.canvas.draw()
    plt.show(block=False)
    i = -1
    while True:
        i+=1
        print()
        print()
        print()
        print("frame", i)
        # run dynam once 
        try:

            # leapfrog:
            # ut+Δt = ut−Δt + 2Δt · f(ut)

            # if it is an initial forward step, do one timestep
            #int_type = True
            int_type = 0

            #friction(u, v, du, dv, p, spa)
            if int_type == 0:
                # forward step
                aflux(u, v, p, pu, pv)
                advecm(p, pu, pv, pp1, convergence, dx, dym, dt)
                advectracer(pu, pv, tracer, tracer_p1, dx, dym, dt)
                #advectracer(pu, pv, t, tp1, dx, dym, dt)
                pgf(du, dv, p, p_center, height, t, spa, dxc, dym)
                advecv(u, v, du, dv, p, dxc, dym)
                coriolis(u, v, du, dv, latc)
                #friction(u, v, du, dv, p, spa, dxc, dym, dt)
                sdrag(u, v, p, p_center, t, dt)
                dynam(u, v, up1, vp1, du, dv, timestep)

                # do array swapping
                um1, u, up1 = u, up1, um1
                vm1, v, vp1 = v, vp1, vm1
                p, pp1 = pp1, p
                tracer, tracer_p1 = tracer_p1, tracer
                t, tp1 = tp1, t

            elif int_type == 1:
                # leapfrog
                if( i % 24 == 0):
                    aflux(u, v, p, pu, pv)
                    advecm(p, pu, pv, pp1, convergence, dx, dym, dt)
                    advectracer(pu, pv, tracer, tracer_p1, dx, dym, dt)
                    pgf(du, dv, p, p_center, height, t, spa, dxc, dym)
                    advecv(u, v, du, dv, p, dxc, dym)
                    #coriolis(u, v, du, dv, latc)
                    dynam(u, v, up1, vp1, du, dv, timestep)

                    # do array swapping
                    um1, u, up1 = u, up1, um1
                    vm1, v, vp1 = v, vp1, vm1
                    p, pp1 = pp1, p
                    tracer, tracer_p1 = tracer_p1, tracer
                else:
                    aflux(u, v, p, pu, pv)
                    advecm(p, pu, pv, pp1, convergence, dx, dym, dt)
                    advectracer(pu, pv, tracer, tracer_p1, dx, dym, dt)
                    pgf(du, dv, p, p_center, height, t, spa, dxc, dym)
                    advecv(u, v, du, dv, p, dxc, dym)
                    #coriolis(u, v, du, dv, latc)
                    dynam(um1, vm1, up1, vp1, du, dv, dt2)

                    # array swapping
                    um1, u, up1 = u, up1, um1
                    vm1, v, vp1 = v, vp1, vm1
                    p, pp1 = pp1, p
                    tracer, tracer_p1 = tracer_p1, tracer
            elif int_type == 2:
                # backward euler
                dtthird = dt * 2/3
                aflux(u, v, p, pu, pv)
                advecm(p, pu, pv, pm1, convergence, dx, dym, dtthird)
                advectracer(pu, pv, tracer, tracer_p1, dx, dym, dtthird)
                pgf(du, dv, p, p_center, height, t, spa, dxc, dym)
                advecv(u, v, du, dv, p, dxc, dym)
                #coriolis(u, v, du, dv, latc)
                dynam(u, v, um1, vm1, du, dv, dtthird)

                # do the second step now
                aflux(um1, vm1, pm1, pu, pv)
                advecm(p, pu, pv, pp1, convergence, dx, dym, dt)
                advectracer(pu, pv, tracer, tracer_p1, dx, dym, dt)
                pgf(du, dv, pm1, p_center, height, t, spa, dxc, dym)
                advecv(um1, vm1, du, dv, pm1, dxc, dym)
                #coriolis(um1, vm1, du, dv, latc)
                dynam(u, v, up1, vp1, du, dv, timestep)

                um1, u, up1 = u, up1, um1
                vm1, v, vp1 = v, vp1, vm1
                p, pp1 = pp1, p
                tracer, tracer_p1 = tracer_p1, tracer


            else:
                # DOUBLE backward euler
                dtthird = dt * 2/3
                aflux(u, v, p, pu, pv)
                advecm(p, pu, pv, pm1, convergence, dx, dym, dtthird)
                advectracer(pu, pv, tracer, tracer_p1, dx, dym, dtthird)
                pgf(du, dv, p, p_center, height, t, spa, dxc, dym)
                advecv(u, v, du, dv, p, dxc, dym)
                #coriolis(u, v, du, dv, latc)
                dynam(u, v, um1, vm1, du, dv, dtthird)

                # do the second step now
                aflux(um1, vm1, pm1, pu, pv)
                advecm(p, pu, pv, pp1, convergence, dx, dym, dt)
                advectracer(pu, pv, tracer, tracer_p1, dx, dym, dt)
                pgf(du, dv, pm1, p_center, height, t, spa, dxc, dym)
                advecv(um1, vm1, du, dv, pm1, dxc, dym)
                #coriolis(um1, vm1, du, dv, latc)
                dynam(u, v, up1, vp1, du, dv, timestep)

                um1, up1 = up1, um1
                vm1, vp1 = vp1, vm1
                pm1, pp1 = pp1, pm1

                # do ANOTHER backward euler step!
                aflux(um1, vm1, pm1, pu, pv)
                advecm(p, pu, pv, pp1, convergence, dx, dym, dt)
                advectracer(pu, pv, tracer, tracer_p1, dx, dym, dt)
                pgf(du, dv, pm1, p_center, height, t, spa, dxc, dym)
                advecv(um1, vm1, du, dv, pm1, dxc, dym)
                #coriolis(um1, vm1, du, dv, latc)
                dynam(u, v, up1, vp1, du, dv, timestep)

                um1, u, up1 = u, up1, um1
                vm1, v, vp1 = v, vp1, vm1
                p, pp1 = pp1, p
                tracer, tracer_p1 = tracer_p1, tracer

            # do crappy evaporation and precipitation
            for y in range(numy): 
                for x in range(numx ):
                    if (tracer_p1[y][x] > 25):
                        tracer_p1[y][x] = 25.0
                    elif tracer_p1[y][x] < 12:
                        tracer_p1[y][x] += 0.5 / dym

        except FloatingPointError:
            #print(u)
            #print(v)
            #print(pressure)
            #print(p_next)
            traceback.print_exc()
            exit()


        # diagnostics: compute total momentum
        m = p_center * ((u ** 2 + v ** 2) ** (1/2))



        print( np.amax(u), np.amax(v), np.amax(du), np.amax(dv), np.amax(convergence), np.sum(m))
        print(np.amin(p), np.amax(p), np.amax(spa))


        #image.set_data(((u ** 2 + v ** 2) ** (1/2)))
        #image.set_data(tracer)
        image.set_data(p)
        fig.canvas.draw()
        #plt.draw()
        #time.sleep(5)

        #plt.imshow(tracer, cmap='gray')
        #plt.show()

    











if __name__ == "__main__":
    main()

