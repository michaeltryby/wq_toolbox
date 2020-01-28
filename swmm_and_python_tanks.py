# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-28 10:31:55
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-28 10:53:17

# IMPORT 
# Import modules
from pyswmm_lite import environment
import numpy as np
from scipy.integrate import odeint
import datetime
import matplotlib.pyplot as plt


# SWMM MODEL
# SETUP
# Build the hydraulics configuration dictionary
config1 = {
    "swmm_input": "./tank_test.inp",
    "states": [("Tank", "pollutantN", "1")],
    }


# SIMULATION
# Initialize the environment
env = environment("./tank_test.inp", ctrl=False)
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
        sol = odeint(CSTR, 10.0, np.array([0,dt]), 
            args=(0.10, self.env.sim._model.getNodeResult("Tank",3), 
            self.env._getNodeInflow("Tank"), self.env._getLinkFlow("Valve"), 
            self.env._getNodePollutant("Tank", "1")))
        return sol

# Set up data containers
inflows = []
depth = []
outflows = []
conc = [] 
time = []

# Initiate treatment class
treat = Treatment(env)
done = False

# Run Simulation
while not done:
    # Compute the time step 
    t0 = env.sim._model.getCurrentSimulationTime()
    
    # Steps the simulation
    done = env.step()

    # call current time step
    t1 = env.sim._model.getCurrentSimulationTime()

    # Calculate difference between two times (det)
    dt = t1 - t0
    dt = dt.seconds
    time.append(dt)
    
    # Run water quality work
    sol = treat.step(dt)
    c = float(sol[-1])

    # Append data to data containers
    conc.append(c)
    inflows.append(env.getNodeInflow("Tank"))
    depth.append(env.getNodeDepth("Tank"))
    outflows.append(env.getLinkFlow("Valve"))

    # Set new concentration
    env._setNodePollutant("Tank","1", c)

# End Simulation & Close SWMM
env.sim._model.swmm_end()
env.sim._model.swmm_close()


# PYTHON TANK MODEL
# Setup tank volume/height equations
def tank(A, Qin, Qout, V_current):
    dvdt = Qin - Qout
    V_nextstep = V_current + dvdt
    if V_nextstep < 10**-5:
        V_nextstep = 0.0
    h = V_nextstep / A
    return V_nextstep, h

# Parameters
A = 100
V = 0
h = V / 100
Cd = 1.0
valve_area = 1.0
Qin = 5.0
Cin = 0

inflows2 = []
depth2 = []
outflows2 = []
conc2 = []

# Run simulation
for i in time:
    # Calculates Qout from tank
    if h < 10**-5:
        Qout = 0.0
    else:
        Qout = valve_area * Cd * np.sqrt(2 * 9.80 * h)
    V, h = tank(A, Qin, Qout, V)
    # Calculate CSTR concentration
    sol = odeint(CSTR, 10.0, np.array([0,i]), args=(0.10, V, Qin, Qout, Cin))
    c = float(sol[-1])
    conc2.append(c)
    inflows2.append(Qin)
    depth2.append(h)
    outflows2.append(Qout)

# Graph results
plt.subplot(4, 1, 1)
plt.plot(inflows, 'b', label="SWMM")
plt.plot(inflows2, 'r--', label="Python")
plt.ylabel("Inflows")
plt.subplots_adjust(hspace= 0.5)

plt.subplot(4, 1, 2)
plt.plot(depth, 'b', label="SWMM")
plt.plot(depth2, 'r--', label="Python")
plt.ylabel("Depth")
plt.subplots_adjust(hspace= 0.5)

plt.subplot(4, 1, 3)
plt.plot(outflows, 'b', label="SWMM")
plt.ylabel("Outflow", 'r--', label="Python")
plt.subplots_adjust(hspace= 0.5)

plt.subplot(4, 1, 4)
plt.plot(conc, 'b', label="SWMM")
plt.plot(conc2, 'r--', label="Python")
plt.xlabel("Time (s)")
plt.ylabel("Conc")

plt.legend()
plt.show()

