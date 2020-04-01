# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-02-05 10:05:44
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-02-27 12:46:55

# IMPORT
# Import modules
from pyswmm_lite import environment
import numpy as np
import matplotlib.pyplot as plt

# SIMULATION
# Initialize the environment
env = environment("./tank_test_forcomparison.inp", ctrl=False)
done = False
conc = []
conc2 = []
depth = []
flow = []

env._setNodePollutant("Tank", 0, 5.0)

inflow = np.zeros(50000)
inflow[0:500] = 1.0
ggg = 1
# Run Simulation
while not done:
    env.sim._model.setNodeInflow("Tank", inflow[ggg])

    # Steps th esimulation
    done = env.step()

    ggg += 1

    conc.append(env._getNodePollutant("Tank", 0))
    conc2.append(env._getLinkPollutant("Valve", 0))



# End Simulation & Close SWMM
env.sim._model.swmm_end()
env.sim._model.swmm_close()

"""
TO DO: Create a class for graphing treatment results
"""
plt.subplot(1,4,1)
plt.plot(conc)
plt.subplot(1,4,2)
plt.plot(conc2)

#plt.xlabel("Time (s)")
plt.ylabel("Concentration")
plt.show()
