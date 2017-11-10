from constants import *
import numpy as np
from utilities import reast, rwest, rnorth, rsouth, rmul


def advectemp(pu: np.ndarray, pv, t, tp1, p: np.ndarray, pp1: np.ndarray,  dx, dym, dt):
    tp1.fill(0.0)
    # p_east = np.roll(p, -1, 1)
    # p_west = np.roll(p,  1, 1)
    # p_north = np.roll(p,  1, 0)
    # p_south = np.roll(p, -1, 0)

    # adjust temp
    adjtemp = t * rmul(p, dx) * dym

    tfornext = adjtemp * dt

    # print(adjtemp)

    # compute the amount to move
    t_east = tfornext * pu
    t_west = tfornext * -pu
    t_north = tfornext * -pv
    t_south = tfornext * pv

    # make sure that you only move the max amount
    maxmove = adjtemp / 4
    t_east = np.maximum(0, np.minimum(t_east, maxmove))
    t_west = np.maximum(0, np.minimum(t_west, maxmove))
    t_north = np.maximum(0, np.minimum(t_north, maxmove))
    t_south = np.maximum(0, np.minimum(t_south, maxmove))

    left = adjtemp - t_east - t_west - t_north - t_south

    t_next = left + reast(t_east) + rwest(t_west) + rnorth(t_north) + rsouth(t_south)

    # unadjust temp
    t_next = t_next / (rmul(pp1, dx) * dym)
    # print(t_next)

    return t_next








    # advect the tracer
    # for y in range(-1, numy - 1):
    #     for x in range(-1, numx - 1):
    #         # compute the amount moved in each direction.
    #         t_east = pu[y][x] * t[y][x] / dx[y] * dt
    #         t_west = -pu[y][x-1] * t[y][x] / dx[y] * dt
    #         t_north = -pv[y][x] * t[y][x] / dx[y] * dt
    #         t_south = pv[y-1][x] * t[y][x] / dx[y] * dt
    #
    #         # make sure that you only do the ones that are positive
    #         # and make sure that you only move 1/4 of the tracer in any direction
    #         maxmove = t[y][x] / 4
    #         t_east = max(0, min(t_east, maxmove))
    #         t_west = max(0, min(t_west, maxmove))
    #         t_north = max(0, min(t_north, maxmove))
    #         t_south = max(0, min(t_south, maxmove))
    #
    #         # compute how much will be left
    #         left = t[y][x] - t_east - t_west - t_north - t_south
    #
    #         # move the tracer
    #         tp1[y][x] += left
    #         tp1[y-1][x] += t_north
    #         tp1[y+1][x] += t_south
    #         tp1[y][x-1] += t_west
    #         tp1[y][x+1] += t_east

