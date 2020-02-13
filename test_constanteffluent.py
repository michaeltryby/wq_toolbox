# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-02-05 10:05:44
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-02-13 11:30:32

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

env._setNodePollutant("Tank", 0, 5)

inflow = np.zeros(50000)
inflow[0:5000] = 5.0
ggg = 1
# Run Simulation
while not done:
    env.sim._model.setNodeInflow("Tank", inflow[ggg])
    env._setNodePollutant("Tank", 0, 5)

    # Steps th esimulation
    done = env.step()

    ggg += 1
    conc.append(env._getNodePollutant("Tank", 0))
    conc2.append(env._getLinkPollutant("Valve", 0))
    depth.append(env._getNodeDepth("Tank"))
    flow.append(env._getNodeInflow("Tank"))


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
plt.subplot(1,4,3)
plt.plot(depth)
plt.subplot(1,4,4)
plt.plot(flow)
#plt.xlabel("Time (s)")
plt.ylabel("Concentration")
plt.show()
