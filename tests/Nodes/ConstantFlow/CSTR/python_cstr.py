# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-20 15:01:08
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-04-29 15:35:56

"""
# METHOD 1 
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint 
import csv

Qout = np.genfromtxt('outflow.csv', delimiter=',')
Qin = np.genfromtxt('inflow.csv', delimiter=',')
V = np.genfromtxt('volume.csv', delimiter=',')
Cin = np.genfromtxt('Cin.csv', delimiter=',')

k = 2.0

def CSTR(C, t, Qin, Cin, Qout, V, k): 
    # CSTR differential equation
    dCdt = (Qin*Cin - Qout*C)/V - k*C
    return dCdt

# initial concentration
C0 = 10
# start, end, number of timesteps
t = np.linspace(0, 20, 50) 
 
# solve equation
C = odeint(CSTR, C0, t, (Qin[0], Cin[0], Qout[0], V[0], k))

# plot results
plt.plot(t,C)
plt.xlabel('time')
plt.ylabel('concentration')
plt.show()

# METHOD 2
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import csv

Qout = np.genfromtxt('outflow.csv', delimiter=',')
Qin = np.genfromtxt('inflow.csv', delimiter=',')
V = np.genfromtxt('volume.csv', delimiter=',')
Cin = np.genfromtxt('Cin.csv', delimiter=',')
k = 2.0
c0, t0 = 10.0, 0.0

def CSTR(t, C, Qin, Cin, Qout, V, k): 
    ''' Continuouslyd stirred tank reactor
    Qin     = flow into reactor
    Qout    = flow out of reactor
    Cin     = concentration flowing into reactor
    V       = volume of water in reactor
    k       = reaction rate
    '''
    # CSTR differential equation
    dCdt = (Qin*Cin - Qout*C)/V - k*C
    return dCdt

r = ode(CSTR)
r.set_f_params(Qin[0],Cin[0],Qout[0],V[0],k)
r.set_initial_value(c0, t0)
dt = 1

for i in range(0,100):
    r.integrate(r.t+dt)
    r.set_initial_value(r.y, r.t)
    r.set_f_params(Qin[i],Cin[i],Qout[i],V[i],k)
    print(r.t, r.y)
"""

# METHOD 3
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import csv

Qout = np.genfromtxt('outflow.csv', delimiter=',')
Qin = np.genfromtxt('inflow.csv', delimiter=',')
V = np.genfromtxt('volume.csv', delimiter=',')
Cin = np.genfromtxt('Cin.csv', delimiter=',')
k = 0.2

def CSTR(t, C, Qin, Cin, Qout, V, k): 
    ''' Continuouslyd stirred tank reactor
    Qin     = flow into reactor
    Qout    = flow out of reactor
    Cin     = concentration flowing into reactor
    V       = volume of water in reactor
    k       = reaction rate
    '''
    # CSTR differential equation
    dCdt = (Qin*Cin - Qout*C)/V - k*C
    return dCdt

solver = ode(CSTR)
C = np.zeros(1800)
dt=1

for i in range(0,1799):
    solver.set_f_params(Qin[i],Cin[i],Qout[i],V[i],k)
    solver.set_initial_value(C[i], i)
    C[i]=solver.integrate(i+dt)

    """
    if i == 0:
        solver.set_initial_value(0.0, 0.0)
        C[i]=solver.integrate(i)
    else: 
        solver.set_initial_value(C[i], i)
        C[i]=solver.integrate(i+dt)
    """

# Plot results
plt.subplot(311)
plt.plot(Qin, 'b', label="Qin")
plt.plot(Qout, 'r', label="Qout")
plt.ylabel("Flow")

plt.subplot(312)
plt.plot(C)
plt.ylabel("Concentration")

plt.subplot(313)
plt.plot(V)
plt.ylabel("Volume")
plt.xlabel("Time (s)")
plt.show()

