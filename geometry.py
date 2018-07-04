from math import pi

import numpy as np




class Geometry:
    """
    This class stores the geometrical data for a world
    """

    def __init__(self, x_cells: int, y_cells: int, radius: float, torroid: bool=False):
        self.radius = radius
        self.torroid = torroid
        self.x_cells = x_cells
        self.y_cells = y_cells
        # make the arrays
        # constants
        self.dx = np.zeros(y_cells)
        self.dxc = np.zeros(y_cells)
        self.lat = np.zeros(y_cells)
        self.latc = np.zeros(y_cells)


        # calculate geometry
        # find the radians for each height difference: 180deg = pi rads
        # so find dx at the equator by doing 2pi / numx
        dxe = 2 * pi / y_cells
        # dy is 180deg so pi / numy
        dye = pi / y_cells
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
            for ii in range(y_cells):
                lat_cur = maxlat - dye * ii
                self.lat[ii] = lat_cur
                self.dx[ii] = np.cos(lat_cur) * dxe * radius
            #    print(lat_cur/pi * 180, dx[ii], dx[ii])
        else:
            # torroid
            for ii in range(y_cells):
                self.dx[ii] = dxe * radius
        print(self.lat)
        print(self.dx)


        # do the corners next
        # this finds the latitude of the first center from the top
        maxlat = pi / 2 - dye

        if not torroid:
            # now loop and do the calculations
            for ii in range(y_cells):
                lat_cur = maxlat - dye * ii
                self.latc[ii] = lat_cur
                self.dxc[ii] = np.cos(lat_cur) * dxe * radius
            #    print(lat_cur/pi * 180, dxc[ii], dxc[ii])
            self.dxc[-1] = self.dxc[1] # prevent explosions
            self.latc[-1] = 0
        else:
            # torroid
            for ii in range(y_cells):
                self.dxc[ii] = dxe * radius
        print(self.latc)
        print(self.dxc)
