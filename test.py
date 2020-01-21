# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-20 12:07:35
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-21 11:16:23

# IMPORT 
# Import SWMM Modules
from pyswmm_lite import environment
import numpy as np
from scipy.integrate import odeint
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
done = False

#SETUP WQ FUNCTION
def CSTR(C, t, k, V, Qin, Qout, Cin):
	# Qin = env._getNodeInflow("P1", "1")
	# Qout = env._getLinkFlow("7")
	# Cin = env._getNodePollutant("P1")
	# V = env.sim._model.getNodeResult("P1",3)

	# CSTR equation
	dCdt = (Qin*Cin - Qout*C)/V - k*C
	return dCdt
	# Set new concentration
	env._setNodePollutant("P1", dCdt)


# Run Simulation
while not done:
	# call current time
	env.sim.swmm_model.getCurrentSimualationTime()
	# Steps the simulation
	done = env.step()

    # call current time step
    # calculate difference between two times (det)

	#Run water quality work
	#odeint(CSTR, 0.0, np.linspace(0, 40, 100), args=(0.10, env.sim._model.getNodeResult("P1",3), env._getNodeInflow("P1"), env._getLinkFlow("7"), env._getNodePollutant("P1", "1")))
    

# End Simulation & Close SWMM
env.sim._model.swmm_end()
env.sim._model.swmm_close()

