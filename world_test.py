from world import World
import numpy as np
from constants import *
from utilities import *

w = World(5, 5, toy=True)

print(w.ideal_gas_pressure())

w.height[2][2] = 100
w.t[1][1] = 256

geo = w.geopotential()
print(geo[0])
print(geo[1])
print(geo.shape)
print(geo)

print("divergence")
print(divergence(geo))

print()
print("Temperature")
print(w.temperature())


print()
print("Density")
print(w.density())


print()
print("pgf")
print(w.pgf())












