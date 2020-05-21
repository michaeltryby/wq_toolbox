# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-21 07:48:53

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
#dict4 = {'Tank2': {0: [0.000015, 0.0100]}}
#dict5 = {'Tank2': {0: [0.000015, 0.0053]}}
#dict6 = {'Tank2': {0: [0.000015, 0.0100]}}

# Lists to store results
tank2_flow = []
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
    #solver2 = ode(CSTR_tank)
    #solver3 = ode(CSTR_tank)

    # Step through the simulation    
    for index,step in enumerate(sim):
        # Calculate dt
        current_step = sim.current_time
        dt = (current_step - last_timestep).total_seconds()
        last_timestep = current_step

        # Calculate DO concentration in tank layers
        # Get depth to calculate DO
        
        depth = sim._model.getNodeResult("Tank2", 5)
        tank2_depth.append(depth)
        """
        if depth <= 6.0:
            # Calculate DO concentration in first layer
            DO1 = DO_Layer1.DO_solver(index)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 < 1.0:
                k_ni = 0.00042  # [1/s]
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
        Qin.append(sim._model.getNodeResult("Tank2",0))
        Cin.append(sim._model.getNodeC2("Tank2",0))
        Qout.append(sim._model.getNodeResult("Tank2",1))
        V.append(sim._model.getNodeResult("Tank2",3))
        k_ni = -0.0005
        """
        solver1.set_f_params(Qin,Cin,Qout,V,k_ni)
        #Solve ODE
        if index == 0:
            solver1.set_initial_value(0.0, 0.0)
            solver1.integrate(solver1.t+dt)
            #solver2.set_f_params(Qin,Cin,Qout,V,k_ni)
            #solver2.set_initial_value(0.0, 0.0)
            #solver2.integrate(solver2.t+dt)
            #solver3.set_f_params(Qin,Cin,Qout,V,k_ni)
            #solver3.set_initial_value(0.0, 0.0)
            #solver3.integrate(solver3.t+dt)
        else:
            solver1.set_initial_value(solver1.y, solver1.t)
            solver1.integrate(solver1.t+dt)
            #solver2.set_f_params(Qin,solver1.y[0],Qout,V,k_ni)
            #solver2.set_initial_value(solver2.y, solver2.t)
            #solver2.integrate(solver2.t+dt)
            #solver3.set_f_params(Qin,solver2.y[0],Qout,V,k_ni)
            #solver3.set_initial_value(solver3.y, solver3.t)
            #solver3.integrate(solver3.t+dt)
        """
        # Set new concentration
        sim._model.setNodePollutant("Tank2", 0, 5)

        # Get Nitrate concentration        
        c3_N = Tank2.pollut_quality
        tank2_N.append(c3_N['NO'])

        # Get flow for  each asset
        tank2_flow.append(sim._model.getNodeResult("Tank2", 0))

np.savetxt('Qin.csv', Qin, delimiter=',')
np.savetxt('Qout.csv', Qout, delimiter=',')
np.savetxt('V.csv', V, delimiter=',')
np.savetxt('Cin.csv', Cin, delimiter=',')


"""        
#----------------------------------------------------------------------#
# Confirm CSTR matches steady state equilibrium concentration
# Found average V and Q for simulation, r same as k above 
V = 100
Q = 10
r = 0.01
C0 = 10
Cn = Co /((1+(k*V/Q))**3)
"""
#----------------------------------------------------------------------#
# Plot Results

# Calculate load
Nload3 = [a*b for a,b in zip(tank2_N,tank2_flow)]
Ncum_load3 = np.cumsum(Nload3)


# Plot Results
plt.subplot(511)
plt.plot(tank2_N, 'm', label='Tank2')
plt.ylabel("NO Conc")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(512)
plt.plot(tank2_DO1, 'g:', label='Tank2_DO1')
plt.plot(tank2_DO2, 'b:', label='Tank2_DO2')
plt.plot(tank2_DO3, 'm:', label='Tank2_DO3')
plt.ylabel("DO Conc")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(513)
plt.plot(tank2_depth, 'm', label='Tank2')
plt.ylabel("Depth")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(514)
plt.plot(tank2_flow, 'm', label='Tank2')
plt.ylabel("Flow")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(515)
plt.plot(Ncum_load3, 'm', label='Tank2')
plt.ylabel("NO Cumm Load")
plt.xlabel("Time (s)")
plt.legend()
plt.show()
