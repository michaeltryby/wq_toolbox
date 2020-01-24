# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-24 09:34:15
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-24 13:10:40

# IMPORT 
# Import modules
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Define simulation time
def simulation(days, dt):
    time_steps = days * 24 * 3600 / dt
    return [int(dt), int(time_steps)]

# Define tank and calculate tank dynamics
def tank(z, t, A, Qin):
    V = z[0]
    Qout = z[1]

    dVdt = Qin - Qout
    Qout = np.sqrt(2 * 9.8 * (V / A))
    return [dVdt, Qout]

"""
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

conc = []
"""

water_vol = []
Qout_all = []
dt, time_steps = simulation(1, 10)
z0 = [0,0]

# Run simulation
for i in range(0, time_steps):
    z = odeint(tank, z0, np.array([0,dt]), args=(1000, 10))
    water_vol.append(z[:,0])
    print("Tank Volume: ", z[1,0])
    Qout_all.append(z[:,1])
    print("Tank Outflow: ", z[1,1])

    #sol = CSTR()
    #conc.append(sol)
    #print("Conc set to:", sol)



# Graph results
#time = np.arange(0, len(water_vol))
#plt.plot(time, water_vol, label="Water Vol (m3)")
#plt.xlabel("Time (s)")
#plt.show()

"""
plt.plot(time, conc)
plt.xlabel("Time (s)")
plt.ylabel("Concentration")
"""