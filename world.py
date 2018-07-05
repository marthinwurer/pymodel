
import numpy as np
import constants
from geometry import Geometry
from utilities import gradient


class World:
    """
    This class defines the atmosphere data for a world.
    It stores the potential temperature, the atmospheric density, the moisture, and the wind velocities.
    It also stores the world geometry.
    The world is broken up into cells with a given height and width, as well as a depth.
    Values stored are the cellular average.
    """

    def __init__(self, y_cells: int, x_cells: int, toy:bool=False):
        self.cell_height = 33000 # the top boundary of the top cell in a column. About 1/3 of the atmosphere.
        self.layers = 1 # the number of layers in the model
        self.x_cells = x_cells # in columns
        self.y_cells = y_cells # in columns
        self.dim = (y_cells, x_cells)
        self.u = np.zeros(self.dim) # eastward velocity: m/s
        self.v = np.zeros((y_cells, x_cells)) # southward velocity: m/s
        self.p = np.full((y_cells, x_cells), constants.air_density) # density: kg/m^3 # from wikipedia for the
        self.t = np.full((y_cells, x_cells), 273.15) # potential temperature: Kelvin
        self.m = np.zeros((y_cells, x_cells)) # moisture: kg/m^3
        self.height = np.zeros((y_cells, x_cells)) # geopotential height: average land height in the column
        self.g_height = np.zeros(self.dim)

        if not toy:
            self.geometry = Geometry(x_cells, y_cells, constants.radius)
        else:
            self.geometry = Geometry(x_cells, y_cells, constants.radius, torroid=True, toy=True)




    def ideal_gas_pressure(self):
        # kg
        #---------------------------------
        #
        return self.p * constants.R_dry * self.t / 100 # / 100 turns into hPa

    def geopotential(self):
        geo_north = np.zeros(self.dim)
        geo_south = np.zeros(self.dim)
        geo_east = np.zeros(self.dim)
        geo_west = np.zeros(self.dim)
        p = self.p

        g = self.height * constants.gravity
        dx = self.geometry.dym

        g_grad = gradient(g, dx)

        out = [grad * p for grad in g_grad]

        return out

    def pgf(self):
        pass




