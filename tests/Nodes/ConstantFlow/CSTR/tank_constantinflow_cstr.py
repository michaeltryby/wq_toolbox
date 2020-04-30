# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-04-27 14:20:21
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-04-29 15:39:20

from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

# Single Tank, Constant Inflow, CSTR
# PySWMM Toolbox
# Create simulation
conc1 = []
flow1 = []

def CSTR(t, C, Qin, Cin, Qout, V, k, x): 
    ''' Continuouslyd stirred tank reactor
    Qin     = flow into reactor
    Qout    = flow out of reactor
    Cin     = concentration flowing into reactor
    V       = volume of water in reactor
    k       = reaction rate
    x       = order of reaction rate ( x = 1 if first order )
    '''
    # CSTR differential equation
    dCdt = (Qin*Cin - Qout*C)/V - k*C**x
    return dCdt

with Simulation("./tank_constantinflow_notreatment.inp") as sim:
    Tank = Nodes(sim)["Tank"]
    start_time = sim.start_time
    last_timestep = sim.start_time
    k = 0.2
    x = 1
    T = (sim.end_time-sim.start_time).total_seconds()
    C = np.zeros(int(T))
    solver = ode(CSTR)

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
        # Get flow for Tank
        flow1.append(sim._model.getNodeResult("Tank",0))

        # Get parameters
        Qin = sim._model.getNodeResult('Tank',0)
        Cin = sim._model.getNodeCin('Tank',0)
        Qout = sim._model.getNodeResult('Tank',1)
        V = sim._model.getNodeResult('Tank',3)
        solver.set_f_params(Qin,Cin,Qout,V,k,x)
        t = (sim.current_time-sim.start_time).total_seconds()
        
        #Solve ODE
        solver.set_initial_value(C[index],t)
        C[index]=solver.integrate(t+dt)

        """
        if index == 1:
            # Solve ODE
            solver.set_initial_value(0.0, 0.0)
            C[index]=solver.integrate(dt)  
        else: 
            # Solve ODE
            solver.set_initial_value(C[index],t)
            C[index]=solver.integrate(t+dt)
        """

        # Set new concentration
        sim._model.setNodePollutant('Tank', 0, C[index])
        

# Plot Results
# Concentration Results
plt.subplot(211)
plt.plot(conc1, 'r--', label='toolbox')
plt.ylabel("Concentration")
plt.legend()

# Flow Results
plt.subplot(212)
plt.plot(flow1, 'r--', label='toolbox')
plt.xlabel("Time (s)")
plt.ylabel("Flow")
plt.legend()
plt.show()
