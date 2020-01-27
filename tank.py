# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-24 09:34:15
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-27 10:55:29

# IMPORT 
# Import modules
import numpy as np
import matplotlib.pyplot as plt

def tank(A, Qin, Qout, V_current):
    dvdt = Qin - Qout
    V_nextstep = V_current + dvdt
    if V_nextstep < 10**-5:
        V_nextstep = 0.0
    h = V_nextstep / A
    return V_nextstep, h

A = 100
V = 0
h = V / 100
Cd = 1.0
valve_area = 1.0
Qin = 5.0

inflows = []
depth = []
volume = []
outflows = []

for i in range(0, 3600):
    if h < 10**-5:
        Qout = 0.0
    else:
        Qout = valve_area * Cd * np.sqrt(2 * 9.80 * h)
    V, h = tank(A, Qin, Qout, V)

    inflows.append(Qin)
    depth.append(h)
    volume.append(V)
    outflows.append(Qout)



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

    #sol = CSTR()
    #conc.append(sol)
    #print("Conc set to:", sol)
"""

# Graph results
plt.subplot(4, 1, 1)
plt.plot(inflows)
plt.ylabel("Inflows")
plt.ylim(ymin= 0)
plt.xlim(xmin= 0)

plt.subplot(4, 1, 2)
plt.plot(depth)
plt.ylabel("Depth")
plt.ylim(ymin= 0)
plt.xlim(xmin= 0)

plt.subplot(4, 1, 3)
plt.plot(volume)
plt.ylabel("Volume")
plt.ylim(ymin= 0)
plt.xlim(xmin= 0)

plt.subplot(4, 1, 4)
plt.plot(outflows)
plt.ylabel("Outflow")
plt.subplots_adjust(hspace= 0.5)
plt.ylim(ymin= 0)
plt.xlim(xmin= 0)
plt.show()

"""
plt.plot(time, conc)
plt.xlabel("Time (s)")
plt.ylabel("Concentration")
"""
