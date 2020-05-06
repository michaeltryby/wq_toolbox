# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-04-27 14:20:21
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-06 08:53:14

from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from wq_toolbox.nodes import Node_Treatment

# Single Tank, Constant Inflow, CSTR
# PySWMM Toolbox
# Create simulation

conc1 = []

def tank(t, C, Qin, Cin, Qout, V, k, n): 
    # CSTR differential equation
    dCdt = (Qin*Cin - Qout*C)/V + k*C**n
    return dCdt

with Simulation("./tank_constantinflow_notreatment.inp") as sim:
    k = -0.2
    n = 1
    c0 = 0.0
    solver = ode(tank)
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
        solver.set_f_params(Qin,Cin,Qout,V,k,n)

        if index == 0:
            solver.set_initial_value(c0, 0.0)
            solver.integrate(solver.t+dt)
        else: 
            solver.set_initial_value(solver.y, solver.t)
            solver.integrate(solver.t+dt)
        # Set new concentration
        sim._model.setNodePollutant('Tank', 0, solver.y[0])
                

conc2 = []
dict1 = {'Tank': {0: [-0.2, 1.0, 0.0]}}

with Simulation("./tank_constantinflow_notreatment.inp") as sim:
    Tank = Nodes(sim)["Tank"]
    TR = Node_Treatment(sim, dict1)
    # Step through the simulation    
    for index,step in enumerate(sim):
        # Get newQual for Tank
        co = Tank.pollut_quality
        conc2.append(co['P1'])
        # Run treatment each time step
        TR.CSTR_solver()

# Plot Results
# Concentration Results
plt.plot(conc1, 'r--', label='no-class')
plt.plot(conc2, 'b:', label='class')
plt.ylabel("Concentration")
plt.xlabel("Time (s)")
plt.legend()
plt.show()
