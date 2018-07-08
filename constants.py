
radius = 6371000.0 # radius of the earth, in meters
#radius = 637100.0
gravity = 9.8 # in m/s^2
R_dry = 287.1 # J/(kg*K) = (kg*m^2/s^2)/(kg*K)
cp_dry = 1004.0
kappa = 0.286
ptop = 10.0 # the pressure at the top of the atmosphere
r_rate_earth = 7.292115e-5 # rads/sec, rotation rate of earth (sidereal day)
#r_rate_earth = 0.0
coriolis_m = 2 * r_rate_earth # will need to be multiplied by the sin of the latitude
sigma = [0.5] # proportion of pressure at that layer
dsig = [1.0] # change from this pressure layer to the one above it
# whooo globals
# numx = 8
# numy = 4
numx = 36
numy = 24
lay = 0
# numx = 72
# numy = 36
# numx = 360
# numy = 180
timestep = 900.0
# timestep = 120.0
# timestep = 20.0
# torroid = True
torroid = False

reference_pressure = 100000 # Pa . Reference pressure for potential temperature calculation
air_pressure = 101325 # Pa - Sea level standard atmospheric pressure in pascals
                      # pa = kg/(m*s^2)
air_density = 1.2754 # kg/m^3 : From wikipedia for dry air at STP : https://en.wikipedia.org/wiki/Density_of_air


