
import numpy as np
import constants
from geometry import Geometry


class World:
    """
    This class defines the atmosphere data for a world.
    It stores the potential temperature, the atmospheric density, the moisture, and the wind velocities.
    It also stores the world geometry.
    The world is broken up into cells with a given height and width, as well as a depth.
    Values stored are the cellular average.
    """

    def __init__(self, x_cells: int, y_cells: int):
        self.cell_height = 33000 # the top boundary of the top cell in a column. About 1/3 of the atmosphere.
        self.layers = 1 # the number of layers in the model
        self.x_cells = x_cells # in columns
        self.y_cells = y_cells # in columns
        self.u = np.zeros((y_cells, x_cells)) # eastward velocity: m/s
        self.v = np.zeros((y_cells, x_cells)) # southward velocity: m/s
        self.p = np.full((y_cells, x_cells), constants.air_density) # density: g/m^3 # from wikipedia for the
        self.t = np.full((y_cells, x_cells), 273.15) # potential temperature: Kelvin
        self.m = np.zeros((y_cells, x_cells)) # moisture: kg/m^3
        self.height = np.zeros((y_cells, x_cells)) # geopotential height: average land height in the column
        self.geometry = Geometry(x_cells, y_cells, constants.radius)




    def ideal_gas_pressure(self):
        return self.p * constants.R_dry * self.t

