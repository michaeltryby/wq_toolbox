# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-20 12:07:35
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-20 16:03:53

# IMPORT 
# Import SWMM Modules
from pyswmm_lite import environment
import numpy as np
#import wq_toolbox

#SETUP
# Build the hydraulics configuration dictionary
config1 = {
    "swmm_input": "./test_single.inp",
    "states": [("P1", "pollutantN", "1")],
    }

# SIMULATION
# Initialize the environment
env = environment("./test_single.inp", ctrl=False)
#wq = wq_toolbox(config2, env)
done = False

# Run Simulation
while not done:
	# Run water quality work
	#Qin = environment._getNodeInflow("P1")
	#Qout = environment._getLinkFlow("7")
	#Cin = environment._getNodePollutant("P1")
	#k = 0.5
	#V = environment._getNodeVolume("P1")

	# CSTR equation
	#dCdt = Qin/V*Cin - Qout/V*Cout - k*C

	#env._setNodePollutant()

    # Steps the simulation
    done = env.step()

    
# End Simulation & Close SWMM
env.sim._model.swmm_end()
env.sim._model.swmm_close()

