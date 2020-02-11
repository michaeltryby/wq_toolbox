# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-02-05 10:05:44
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-02-11 08:30:26

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

# Run Simulation
while not done:
    # Steps the simulation
    c = env._getNodePollutant("Tank", 0)
    print(c)
    # Set new concentration
    
    if 0 <= c <= 5: 
        pass
    elif c > 5:
        env._setNodePollutant("Tank", 0, 5)

    done = env.step()

    #d = min(c, 5.0)
    #print(d)
    #env._setNodePollutant("Tank", 0, d)

    #env._getNodePollutant("Tank", 0)
    #env._setNodePollutant("Tank", 0, 5)
    conc.append(env._getNodePollutant("Tank", 0))
    #conc2.append(env._getLinkPollutant("Valve", 0))
    #depth.append(env._getNodeDepth("Tank"))
    #flow.append(env._getNodeInflow("Outfall"))


# End Simulation & Close SWMM
env.sim._model.swmm_end()
env.sim._model.swmm_close()

"""
TO DO: Create a class for graphing treatment results
"""
plt.plot(conc)
#plt.plot(conc2)
plt.xlabel("Time (s)")
plt.ylabel("Concentration")
plt.show()
