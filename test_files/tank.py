# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-24 09:34:15
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-28 10:38:19

# IMPORT 
# Import modules
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Setup tank volume/height equations
def tank(A, Qin, Qout, V_current):
    dvdt = Qin - Qout
    V_nextstep = V_current + dvdt
    if V_nextstep < 10**-5:
        V_nextstep = 0.0
    h = V_nextstep / A
    return V_nextstep, h

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

# Parameters
A = 100
V = 0
h = V / 100
Cd = 1.0
valve_area = 1.0
Qin = 5.0
Cin = 0

# Set up data containers
inflows = []
depth = []
volume = []
outflows = []
conc = []

# Run simulation
for i in range(0, 3600):
    # Calculates Qout from tank
    if h < 10**-5:
        Qout = 0.0
    else:
        Qout = valve_area * Cd * np.sqrt(2 * 9.80 * h)
    V, h = tank(A, Qin, Qout, V)
    # Calculate CSTR concentration
    sol = odeint(CSTR, 10.0, np.array([0,1.0]), args=(0.10, V, Qin, Qout, Cin))
    c = float(sol[-1])
    # Append data to data containers
    conc.append(c)
    inflows.append(Qin)
    depth.append(h)
    volume.append(V)
    outflows.append(Qout)

# Graph hydraulic results
"""
plt.subplot(5, 1, 1)
plt.plot(inflows)
plt.ylabel("Inflows")
plt.ylim(ymin= 0)
plt.xlim(xmin= 0)

plt.subplot(5, 1, 2)
plt.plot(depth)
plt.ylabel("Depth")
plt.ylim(ymin= 0)
plt.xlim(xmin= 0)

plt.subplot(5, 1, 3)
plt.plot(volume)
plt.ylabel("Volume")
plt.ylim(ymin= 0)
plt.xlim(xmin= 0)

plt.subplot(5, 1, 4)
plt.plot(outflows)
plt.ylabel("Outflow")
plt.subplots_adjust(hspace= 0.5)
plt.ylim(ymin= 0)
plt.xlim(xmin= 0)
"""

# Graph pollutant results
plt.plot(conc)
plt.ylabel("Concentration")
plt.xlabel("Time")
plt.show()
