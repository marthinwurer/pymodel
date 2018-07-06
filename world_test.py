from world import World
import numpy as np
from constants import *
from utilities import *

w = World(5, 5, toy=True)

print(w.ideal_gas_pressure())

w.height[2][2] = 100

print(w.geopotential()[0])
print(w.geopotential()[1])
geo = np.array(w.geopotential())
print(geo.shape)
print(geo)

print(divergence(geo))
