# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-02-05 10:05:44
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-02-11 09:47:39

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
#depth = []
#flow = []

C_star = 10

# Run Simulation
while not done:
    # Compute the time step 
    t0 = env.sim._model.getCurrentSimulationTime()
    
    # Steps the simulation
    done = env.step()

    # call current time step
    t1 = env.sim._model.getCurrentSimulationTime()

    # calculate difference between two times (det)
    dt = t1 - t0
    dt = dt.seconds

    # Get current concentration
    C = env._getNodePollutant("Tank", 0)
    
    # Set new concentration
    C_new = C - 0.1*(C**2.0)*dt

    env._setNodePollutant("Tank", 0, C_new)

    conc.append(C)
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
plt.xlabel("Time (s)")
plt.ylabel("Concentration")
plt.show()
