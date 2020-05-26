# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-26 10:41:30

# Import required modules
from pyswmm import Simulation, Nodes
from wq_toolbox.nodes import Node_Quality
from scipy.integrate import ode 
import numpy as np
import matplotlib.pyplot as plt

# Nitrate 3 CSTRs in Series
def CSTR_tank(t, C, Qin, Cin, Qout, V, k):
    dCdt = (Qin*Cin - Qout*C)/V + k*C
    return dCdt

# Make dictionaries for each water quality method
# From Reddy et al. 1980
dict4 = {'Tank2': {0: [0.0003, 10]}}
dict5 = {'Tank2': {0: [0.0003, 10]}}
dict6 = {'Tank2': {0: [0.0003, 10]}}

# Lists to store results
tank2_flow = []
tank2_vol = []
tank2_depth = []
tank2_N = []
tank2_DO1 = []
tank2_DO2 = []
tank2_DO3 = []
Qin = []
Qout = []
Cin = []
V = []

# Setup toolbox simulation
with Simulation("./tanks_channel_NO.inp") as sim:
    # Setup toolbox methods
    #DO_Layer1 = Node_Quality(sim, dict4)
    #DO_Layer2 = Node_Quality(sim, dict5)
    #DO_Layer3 = Node_Quality(sim, dict6)
    
    # Get asset information
    Tank2 = Nodes(sim)["Tank2"]
    
    # Setup dt calculation        
    start_time = sim.start_time
    last_timestep = start_time

    # Setup CSTR solver
    solver1 = ode(CSTR_tank)
    solver1.set_integrator("dopri5")
    solver2 = ode(CSTR_tank)
    solver2.set_integrator("dopri5")
    solver3 = ode(CSTR_tank)
    solver3.set_integrator("dopri5")

    # Step through the simulation    
    for index,step in enumerate(sim):
        # Calculate dt
        current_step = sim.current_time
        dt = (current_step - last_timestep).total_seconds()
        last_timestep = current_step

        # Calculate DO concentration in tank layers
        # Get depth to calculate DO
        """
        depth = sim._model.getNodeResult("Tank2", 5)
        tank2_depth.append(depth)
        
        if depth <= 1.66:
            # Calculate DO concentration in first layer
            DO1 = DO_Layer1.DO_solver(index)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 < 1.0:
                k_ni = 0.00014  # [1/s]
            else:
                k_ni = 0.0
            tank2_DO1.append(DO1)
        
        elif 1.66 < depth <= 3.33:
            # Calculate DO concentration in first two layers
            DO1 = DO_Layer1.DO_solver(index)
            DO2 = DO_Layer2.DO_solver(index)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 < 1.0:
                k_ni1 = 0.00014
            else:
                k_ni1 = 0.0
            if DO2 < 1.0:
                k_ni2 = 0.00014
            else:
                k_ni2 = 0.0
            k_ni = k_ni1 + k_ni2
            tank2_DO1.append(DO1)
            tank2_DO2.append(DO2)
        else:
            # Calculate DO concentration in all three layers
            DO1 = DO_Layer1.DO_solver(index)
            DO2 = DO_Layer2.DO_solver(index)
            DO3 = DO_Layer3.DO_solver(index)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 < 1.0:
                k_ni1 = 0.00014
            else:
                k_ni1 = 0.0
            if DO2 < 1.0:
                k_ni2 = 0.00014
            else:
                k_ni2 = 0.0
            if DO3 < 1.0:
                k_ni3 = 0.00014
            else:
                k_nil3 = 0.0
            k_ni = k_ni1 + k_ni2 + k_ni3
            tank2_DO1.append(DO1)
            tank2_DO2.append(DO2)
            tank2_DO3.append(DO3)
        """
        # Calculate NO concnetration in tanks
        # Get parameters to calculate NO
        Qin=sim._model.getNodeResult("Tank2",0)
        Cin=sim._model.getNodeCin("Tank2",0)
        Qout=sim._model.getNodeResult("Tank2",1)
        V=sim._model.getNodeResult("Tank2",3)
        k_ni = 0.00014
        

        #Solve ODE
        if index == 0:
            solver1.set_initial_value(0.0, 0.0)
            solver1.set_f_params(Qin,Cin,Qout,V,k_ni)
            solver1.integrate(solver1.t+dt)
            solver2.set_initial_value(0.0, 0.0)
            solver2.set_f_params(Qin,Cin,Qout,V,k_ni)
            solver2.integrate(solver2.t+dt)
            solver3.set_initial_value(0.0, 0.0)
            solver3.set_f_params(Qin,Cin,Qout,V,k_ni)
            solver3.integrate(solver3.t+dt)
        else:
            solver1.set_initial_value(solver1.y, solver1.t)
            solver1.set_f_params(Qin,Cin,Qout,V,k_ni)
            solver1.integrate(solver1.t+dt)
            solver2.set_f_params(Qin,solver1.y,Qout,V,k_ni)
            solver2.set_initial_value(solver2.y, solver2.t)
            solver2.integrate(solver2.t+dt)
            solver3.set_f_params(Qin,solver2.y,Qout,V,k_ni)
            solver3.set_initial_value(solver3.y, solver3.t)
            solver3.integrate(solver3.t+dt)
        
        # Set new concentration
        sim._model.setNodePollutant("Tank2", 0, solver3.y[0])

        # Get Nitrate concentration        
        c3_N = Tank2.pollut_quality
        tank2_N.append(c3_N['NO'])

        # Get flow for each asset
        tank2_flow.append(sim._model.getNodeResult("Tank2", 0))
        tank2_vol.append(sim._model.getNodeResult("Tank2", 3))


        
#----------------------------------------------------------------------#
# Confirm CSTR matches steady state equilibrium concentration
# Found average V and Q for simulation, r same as k above 
V = sum(tank2_vol)/len(tank2_vol)
Q = sum(tank2_flow)/len(tank2_flow)
k = 0.00014
Co = 10
Cn = Co /((1+(k*V/Q))**3)
print(Cn)

#----------------------------------------------------------------------#
# Plot Results

# Calculate load
Nload3 = [a*b for a,b in zip(tank2_N,tank2_flow)]
Ncum_load3 = np.cumsum(Nload3)

# SS graphs
SS = [9.88]*len(Ncum_load3)

plt.subplot(311)
plt.plot(tank2_N, 'm', label='Tank2')
plt.plot(SS, 'k', label='SteadyState')
plt.ylabel("NO Conc")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(312)
plt.plot(tank2_flow, 'm', label='Tank2')
plt.ylabel("Flow")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(313)
plt.plot(Ncum_load3, 'm:', label='Tank2')
plt.ylabel("NO Cumm Load")
plt.xlabel("Time (s)")
plt.legend()
plt.show()

