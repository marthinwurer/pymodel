from world import World
import numpy as np
from constants import *
from utilities import *

w = World(5, 5, toy=True)

np.set_printoptions(threshold=np.nan, linewidth=120)

print(w.ideal_gas_pressure())

w.height[2][2] = 100
w.t[1][1] = 256
w.u[3][3] = 1
w.v[3][1] = 1

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

print()
print("momentum")
print(w.u)
temp = w.advect_velocity()
print(temp)

print()
print("averaging")
u_edge = edge_averages(w.u)[0]
v_edge = edge_averages(w.v)[1]
u_center = center_averages(np.asarray([u_edge, v_edge]))
print(u_edge)
print(v_edge)
print(u_center)

print()
print("velocity update")
du = w.update_velocity()
print(du)

print()
print("pressure update")
dp = w.update_pressure()
print(dp)

print()
print("do timestep")
new_p = w.do_timestep()
print(new_p)

for i in range(5000):
    print()
    print(w.do_timestep()/100)







