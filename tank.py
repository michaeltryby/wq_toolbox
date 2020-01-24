# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-24 09:34:15
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-24 15:28:13

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
V = 500
h = V / 100
valve = 1.0
Qin = 18.0

depth = []
outflows = []
volume = []

for i in range(0, 1000):
    if h < 10**-5:
        Qout = 0.0
    else:
        Qout = valve * np.sqrt(2 * 9.80 * h)
    V, h = tank(A, Qin, Qout, V)

    depth.append(h)
    outflows.append(Qout)
    volume.append(V)

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
plt.subplot(3, 1, 1)
plt.plot(depth)
plt.ylabel("Depth")
plt.subplot(3, 1, 2)
plt.plot(outflows)
plt.ylabel("Outflow")
plt.subplot(3, 1, 3)
plt.plot(volume)
plt.ylabel("Volume")
plt.show()

"""
plt.plot(time, conc)
plt.xlabel("Time (s)")
plt.ylabel("Concentration")
"""
