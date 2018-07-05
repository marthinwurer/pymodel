from world import World

w = World(3, 3, toy=True)

print(w.ideal_gas_pressure())

w.height[1][1] = 100

print(w.geopotential()[2])