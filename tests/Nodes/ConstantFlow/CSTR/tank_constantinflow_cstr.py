# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-04-27 14:20:21
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-01 11:16:46

from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

# Single Tank, Constant Inflow, CSTR
# PySWMM Toolbox
# Create simulation
conc1 = []
flow1 = []

def tank(t, C, Qin, Cin, Qout, V, k, n): 
    # CSTR differential equation
    dCdt = (Qin*Cin - Qout*C)/V - k*C**n
    return dCdt

k = 0.2
n = 1
c0 = 0.0
solver = ode(tank)

with Simulation("./tank_constantinflow_notreatment.inp") as sim:
    Tank = Nodes(sim)["Tank"]
    start_time = sim.start_time
    last_timestep = sim.start_time

    # Step through the simulation    
    for index,step in enumerate(sim):
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - last_timestep).total_seconds()
        # Updating reference step
        last_timestep = current_step
        
        # Get newQual for Tank
        conc = Tank.pollut_quality
        conc1.append(conc['P1'])

        # Get parameters
        Qin = sim._model.getNodeResult('Tank',0)
        Cin = sim._model.getNodeCin('Tank',0)
        Qout = sim._model.getNodeResult('Tank',1)
        V = sim._model.getNodeResult('Tank',3)

        t = (sim.current_time-sim.start_time).total_seconds()
        solver.set_f_params(Qin,Cin,Qout,V,k,n)

        if index == 0:
            solver.set_initial_value(c0, 0.0)
            solver.integrate(t+dt)
            print("start",solver.y)  
        else: 
            solver.set_initial_value(solver.y, t)
            solver.integrate(t+dt)
            print(solver.y)

        # Set new concentration
        #sim._model.setNodePollutant('Tank', 0, Cnew)
        

# Plot Results
# Concentration Results
plt.plot(conc1, 'r--', label='toolbox')
plt.ylabel("Concentration")
plt.xlabel("Time (s)")
plt.legend()
plt.show()
