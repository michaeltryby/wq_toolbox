# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-20 12:07:35
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-21 09:18:42

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
def CSTR(Co,t,k):
	Qin = env._getNodeInflow("P1")
	Qout = env._getLinkFlow("7")
	Cin = env._getNodePollutant("P1")
	Cout =  env._getLinkPollutant("7")
	#V = environment._getNodeVolume("P1")
	V = 1

	# CSTR equation
	dCdt = Qin/V*Cin - Qout/V*Cout - k*Co

	#  Set new concentration
	env._setNodePollutant("P1", dCdt)

	return dCdt

# Run Simulation
while not done:
	# Steps the simulation
    done = env.step()

	# Run water quality work
	odeint(CSTR, 0.0, np.linspace(0,200), k=(0.5,))

    
# End Simulation & Close SWMM
env.sim._model.swmm_end()
env.sim._model.swmm_close()

