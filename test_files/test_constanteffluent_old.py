# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-02-05 10:05:44
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-02-05 11:33:09

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
depth = []
flow = []

# Run Simulation
while not done:
    # Steps the simulation
    done = env.step()

    # Set new concentration
    #env._setNodePollutant("Tank","1", 5)
    conc.append(env._getNodePollutant("Tank","P1"))
    depth.append(env._getNodeDepth("Tank"))
    flow.append(env._getNodeInflow("Outfall"))


# End Simulation & Close SWMM
env.sim._model.swmm_end()
env.sim._model.swmm_close()

"""
TO DO: Create a class for graphing treatment results
"""
plt.plot(conc)
plt.plot(depth)
plt.plot(flow)
plt.xlabel("Time (s)")
#plt.ylabel("Concentration")
plt.show()
