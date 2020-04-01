# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-02-05 10:05:44
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-03-26 10:47:01

# IMPORT
# Import modules
from pyswmm_lite import environment
import numpy as np
import matplotlib.pyplot as plt
import faulthandler; faulthandler.enable()

# SIMULATION
# Initialize the environment
#env = environment("./tank_test_removaliszero.inp", ctrl=False)
#env = environment("./tank_test_percentremoval.inp", ctrl=False)
env = environment("./tank_test_forcomparison.inp", ctrl=False)
c = env._getNodeParam("Tank", 0)
print(c)

done = False

conc = []
conc2 = []

# Run Simulation
while not done:
    #print("step")
    # Get current concentration
    #c = env._getNodePollutant("Tank", 0)
    
    # Set new concentration
    #d = c * 0.5
    #env._setNodePollutant("Tank", 0, d)

    # Steps the simulation
    done = env.step()

    conc.append(env._getNodePollutant("Tank", 0))


# End Simulation & Close SWMM
env.sim._model.swmm_end()
env.sim._model.swmm_close()

plt.plot(conc)
plt.xlabel("Time (s)")
plt.ylabel("Concentration")
plt.show()
