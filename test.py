# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-20 12:07:35
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-22 12:19:13

# IMPORT 
# Import modules
from pyswmm_lite import environment
import numpy as np
from scipy.integrate import odeint
import datetime
import matplotlib.pyplot as plt

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

	# CSTR equation
	if V!=0:
		dCdt = (Qin*Cin - Qout*C)/V - k*C
		if np.abs(dCdt) < 0.001:
			dCdt = 0
	else:
		dCdt = 0
	return dCdt

conc = [] 

# Run Simulation
while not done:
	# call current time
	t0 = env.sim._model.getCurrentSimulationTime()
	
	# Steps the simulation
	done = env.step()

	# call current time step
	t1 = env.sim._model.getCurrentSimulationTime()

	# calculate difference between two times (det)
	dt = t1 - t0
	dt = dt.seconds

	#Run water quality work
	sol = odeint(CSTR, 1.0, np.array([0,dt]), args=(0.10, env.sim._model.getNodeResult("P1",3), env._getNodeInflow("P1"), env._getLinkFlow("7"), env._getNodePollutant("P1", "1")))
	c = float(sol[-1])
	conc.append(c)
	print("Conc set to:", c)

	# Set new concentration
	env._setNodePollutant("P1","1", c)


# End Simulation & Close SWMM
env.sim._model.swmm_end()
env.sim._model.swmm_close()


time = np.arange(0, len(conc))
plt.plot(time, conc)
plt.xlabel("Time (s)")
plt.ylabel("Concentration")
plt.show()
