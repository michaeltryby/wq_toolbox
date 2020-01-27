# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-20 12:07:35
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-27 11:04:52

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
	"swmm_input": "./tank_test.inp",
	"states": [("P1", "pollutantN", "1")],
	}


# SIMULATION
# Initialize the environment
env = environment("./test_single.inp", ctrl=False)
done = False

# Setup CSTR equation
def CSTR(C, t, k, V, Qin, Qout, Cin):

	# CSTR equation
	if V!=0:
		dCdt = (Qin*Cin - Qout*C)/V - k*C
		if np.abs(dCdt) < 0.001:
			dCdt = 0
	else:
		dCdt = 0
	return dCdt

# Create Treatment Class
class Treatment:
	def __init__(self, environment):
		self.env = environment

	def step(self, dt):
		sol = odeint(CSTR, 1.0, np.array([0,dt]), 
			args=(0.10, self.env.sim._model.getNodeResult("P1",3), 
			self.env._getNodeInflow("P1"), self.env._getLinkFlow("7"), 
			self.env._getNodePollutant("P1", "1")))
		return sol

conc = [] 
treat = Treatment(env)
done = False

# Run Simulation
while not done:
	"""
	TO DO 1: build step class: needs to calculate dt, input it into solver,
	run treatment, calls solver, and then sets pollutant concentration.

	Needs to import the modules it needs to run. 

	Use the step() PySWMM_Lite code for inspiration
	"""
	# call current time

	# Compute the time step 
	t0 = env.sim._model.getCurrentSimulationTime()
	
	# Steps the simulation
	done = env.step()

	# call current time step
	t1 = env.sim._model.getCurrentSimulationTime()

	# calculate difference between two times (det)
	dt = t1 - t0
	dt = dt.seconds
	# Run water quality work
	"""
	TO DO 2:
	Start building cstr class, figure out how to call swmm fcns directly
	so user doesn't have to input like below, should default get Co from 
	swmm input file or user input

	Things to Consider:
	Think about the code you’re currently working on. What are the 
	properties: the contracts and invariants? Can you use 
	property-based testing framework to verify these automatically?
	"""
	#sol = odeint(CSTR, 1.0, np.array([0,dt]), 
	#	args=(0.10, env.sim._model.getNodeResult("P1",3), 
	#		env._getNodeInflow("P1"), env._getLinkFlow("7"), 
	#		env._getNodePollutant("P1", "1")))
	sol = treat.step(dt)

	c = float(sol[-1])
	conc.append(c)
	print("Conc set to:", c)

	# Set new concentration
	env._setNodePollutant("P1","1", c)


# End Simulation & Close SWMM
env.sim._model.swmm_end()
env.sim._model.swmm_close()

"""
TO DO: Create a class for graphing treatment results
"""
time = np.arange(0, len(conc))
plt.plot(time, conc)
plt.xlabel("Time (s)")
plt.ylabel("Concentration")
plt.show()
