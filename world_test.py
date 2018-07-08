from world import World
import numpy as np
from constants import *
from utilities import *

w = World(5, 5, toy=True)

print(w.ideal_gas_pressure())

w.height[2][2] = 100

geo = w.geopotential()
print(geo[0])
print(geo[1])
print(geo.shape)
print(geo)

print("divergence")
print(divergence(geo))
