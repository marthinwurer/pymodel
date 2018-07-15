
import numpy as np
import constants
from geometry import Geometry
from utilities import gradient, edge_averages, divergence, center_averages


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
        # self.p = np.full((y_cells, x_cells), constants.air_density) # density: kg/m^3 # from wikipedia for the
        self.p = np.full((y_cells, x_cells), constants.air_pressure) # pressure: Pa = kg/(m*s^2)
        self.t = np.full((y_cells, x_cells), 273.15) # potential temperature: Kelvin at 1000 hPa
        self.m = np.zeros((y_cells, x_cells)) # moisture: kg/m^3
        self.height = np.zeros((y_cells, x_cells)) # geopotential height: average land height in the column
        self.g_height = np.zeros(self.dim)
        # The values at the next timestep
        self.p_next = np.zeros(self.dim)
        self.t_next = np.zeros(self.dim)
        self.u_next = np.zeros(self.dim)
        self.v_next = np.zeros(self.dim)

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
        p = self.p

        g = self.height * constants.gravity
        dx = self.geometry.dym

        g_grad = gradient(g, dx)


        return g_grad * p

    def temperature(self):
        return self.t / ((constants.reference_pressure / self.p) ** constants.kappa)

    def density(self):
        """
        calculates the density of the air at the center of the cells from the
        potential temperature and the pressure.
        Units: K
        """
        return self.p / (constants.R_dry * self.temperature())

    def pgf(self):
        pressure_gradient = gradient(self.p)
        p_edge = edge_averages(self.p)
        d_edge = edge_averages(self.density())

        scalar = (p_edge/d_edge)

        return pressure_gradient * scalar

    def flux(self):
        u_edge = edge_averages(self.u)
        v_edge = edge_averages(self.v)
        p_edge = edge_averages(self.p)
        pu_edge = (p_edge * u_edge)[0] # just want the flux at east edge
        pv_edge = (p_edge * v_edge)[1] # just want the flux at south edge

        return np.asarray([pu_edge, pv_edge])

    
    def advect_velocity(self):
        u_edge = edge_averages(self.u)
        v_edge = edge_averages(self.v)
        p_edge = edge_averages(self.p)
        pu_edge = (p_edge * u_edge)[0] # just want the flux at east edge
        pv_edge = (p_edge * v_edge)[1] # just want the flux at south edge

        pu_center = self.p * self.u
        pv_center = self.p * self.v

        del_dot_pu = divergence(self.flux())

        u_flux_east = pu_edge * u_edge[0]
        v_flux_east = pu_edge * v_edge[0]
        u_flux_south = pv_edge * u_edge[1]
        v_flux_south = pv_edge * v_edge[1]

        u_new = del_dot_pu





        return u_new

    def update_pressure(self):
        self.p_next = self.p + divergence(self.flux()) * constants.timestep / self.geometry.dym**2
        return self.p_next


    def update_velocity(self):
        advection = self.advect_velocity()
        # advection = np.zeros(self.p.shape)
        # print(advection)
        pgf = center_averages(self.pgf())# * self.geometry.dym
        print("pgf")
        print(pgf)
        geo = center_averages(self.geopotential()) #* self.geometry.dym
        print("geo")
        print(geo)
        dU = (advection + pgf + geo) * constants.timestep / (self.p ) / ( self.geometry.dym**2)
        # dU = (advection + pgf + geo) * constants.timestep / (self.p )
        print("dU")
        print(dU)
        self.u_next = self.u + dU[0]
        self.v_next = self.v + dU[1]
        return dU

    def do_timestep(self):
        self.update_velocity()
        self.update_pressure()

        self.u = self.u_next
        self.v = self.v_next
        self.p = self.p_next

        return self.p







