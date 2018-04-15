import numpy as np
import rebound

sim = rebound.Simulation()
sim.units = ('AU', 'yr2pi', 'Msun')
#sim.units = ('s', 'm', 'kg')
sim.integrator = "mercurius"

solar_system_objects = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]

sim.add(solar_system_objects, date="2017-12-01 12:00")

print sim.status()

sim.N_active = 9
sim.testparticle_type = 0

sim.move_to_com()

sim.save("checkpoint_mercurius.bin")

print sim.status()
