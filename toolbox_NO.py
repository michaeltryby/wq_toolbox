# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-28 13:33:14

# Import required modules
from pyswmm import Simulation, Nodes, Links
from wq_toolbox.nodes import Node_Quality
from scipy.integrate import ode 
import numpy as np
import matplotlib.pyplot as plt

# Nitrate 3 CSTRs in Series
def CSTR_tank(t, C, Qin, Cin, Qout, V, k):
    dCdt = (Qin*Cin - Qout*C)/V + k*C
    return dCdt

# Lists to store results
Wetland_conc = []
Wetland_flow = []
Wetland_vol = []
Wetland_depth = []
Wetland_N = []
Wetland_DO1 = []
Wetland_DO2 = []
Wetland_DO3 = []
Qin = []
Qout = []
Cin = []
V = []

# DO Variables
k_DO = 0.000015
Co_DO1 = 10.0
Co_DO2 = 5.5
Co_DO3 = 3.0
t = 0

# Tank Inflows
Qin1 = np.genfromtxt('Qin1.txt', delimiter=',')
Qin2 = np.genfromtxt('Qin2.txt', delimiter=',')
Qin3 = np.genfromtxt('Qin3.txt', delimiter=',')

# Setup toolbox simulation
with Simulation("./Ellsworth_Doyle_NO.inp") as sim:
    # Get asset information
    Ellsworth = Nodes(sim)["93-50408"]
    Doyle_Basin = Nodes(sim)["93-50404"]
    Channel = Links(sim)["95-69044"]
    Dwn_Channel = Nodes(sim)["97-50264"]
    Wetland = Nodes(sim)["93-49759"]
    
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
                
        # Set inflow
        sim._model.setNodeInflow("93-50408", Qin1[index])
        sim._model.setNodeInflow("93-50404", Qin2[index])
        sim._model.setNodeInflow("93-49759", Qin3[index])

        # Calculate dt
        current_step = sim.current_time
        dt = (current_step - last_timestep).total_seconds()
        last_timestep = current_step
        t += dt

        # Calculate DO concentration in tank layers
        # Get depth to calculate DO
        depth = sim._model.getNodeResult("93-49759", 5)
        Wetland_depth.append(depth)

        # reset DO if tank is empty
        if depth <= 0.01:
            t = 0

        if depth <= 3.00:
            # Calculate DO concentration in first layer
            DO1 = Co_DO1*np.exp(-k_DO*t)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 < 1.0:
                k_ni = 0.000006  # [1/s]
            else:
                k_ni = 0.0
            Wetland_DO1.append(DO1)
        
        elif 3.0 < depth <= 6.0:
            # Calculate DO concentration in first two layers
            DO1 = Co_DO1*np.exp(-k_DO*t)
            DO2 = Co_DO2*np.exp(-k_DO*t)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 < 1.0:
                k_ni1 = 0.000006
            else:
                k_ni1 = 0.0
            if DO2 < 1.0:
                k_ni2 = 0.000006
            else:
                k_ni2 = 0.0
            k_ni = k_ni1 + k_ni2
            Wetland_DO1.append(DO1)
            Wetland_DO2.append(DO2)
        else:
            # Calculate DO concentration in all three layers
            DO1 = Co_DO1*np.exp(-k_DO*t)
            DO2 = Co_DO2*np.exp(-k_DO*t)
            DO3 = Co_DO3*np.exp(-k_DO*t)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 < 1.0:
                k_ni1 = 0.000006
            else:
                k_ni1 = 0.0
            if DO2 < 1.0:
                k_ni2 = 0.000006
            else:
                k_ni2 = 0.0
            if DO3 < 1.0:
                k_ni3 = 0.000006
            else:
                k_ni3 = 0.0
            k_ni = k_ni1 + k_ni2 + k_ni3
            Wetland_DO1.append(DO1)
            Wetland_DO2.append(DO2)
            Wetland_DO3.append(DO3)
        
        # Calculate NO concentration in tanks
        # Get parameters to calculate NO
        Qin=sim._model.getNodeResult("93-49759",0)
        Cin=sim._model.getNodeCin("93-49759",0)
        Qout=sim._model.getNodeResult("93-49759",1)
        V=sim._model.getNodeResult("93-49759",3)

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
        sim._model.setNodePollutant("93-49759", 0, solver3.y[0])

        # Get Nitrate concentration        
        c3_N = Wetland.pollut_quality
        Wetland_N.append(c3_N['NO'])

        # Get flow for each asset
        Wetland_flow.append(sim._model.getNodeResult("93-49759", 0))
        Wetland_vol.append(sim._model.getNodeResult("93-49759", 3))

#----------------------------------------------------------------------#
# Confirm CSTR matches steady state equilibrium concentration
# Found average V and Q for simulation, r same as k above 
"""
V = sum(Wetland_vol)/len(Wetland_vol)
Q = sum(Wetland_flow)/len(Wetland_flow)
k = 0.000014
Co = 10
Cn = Co /((1+(k*V/Q))**2)
print("Cn", Cn)
"""
#----------------------------------------------------------------------#
# Plot Results

# Calculate load
Nload = [a*b for a,b in zip(Wetland_N,Wetland_flow)]
Ncum_load = np.cumsum(Nload)

# SS graphs
#SS = [Cn]*len(Ncum_load3)

plt.subplot(511)
plt.plot(Wetland_N, 'm', label='Wetland')
#plt.plot(SS, 'k', label='SteadyState')
plt.ylabel("NO Conc")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(512)
plt.plot(Wetland_DO1, 'b:', label='DO_Layer1')
plt.plot(Wetland_DO2, 'r:', label='DO_Layer2')
plt.plot(Wetland_DO3, 'g:', label='DO_Layer3')
plt.ylabel("DO Conc")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(513)
plt.plot(Wetland_depth, 'm', label='Wetland')
plt.ylabel("Depth")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(514)
plt.plot(Wetland_flow, 'm', label='Wetland')
plt.ylabel("Flow")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(515)
plt.plot(Ncum_load, 'm', label='Wetland')
plt.ylabel("NO Cumm Load")
plt.xlabel("Time (s)")
plt.legend()
plt.show()

