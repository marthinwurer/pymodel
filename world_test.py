from density_dynamics import World

w = World(36, 24)

print(w.ideal_gas_pressure())

w.height[10][10] = 100

print(w.geopotential()[2])