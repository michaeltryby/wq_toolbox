# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-06-03 14:38:04

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

# Setup toolbox simulation
with Simulation("./MBDoyle_NO.inp") as sim:
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

    #Temp count for control actions every 30 minutes
    _tempcount = 30

    # Step through the simulation    
    for index,step in enumerate(sim):

        # Calculate dt
        current_step = sim.current_time
        dt = (current_step - last_timestep).total_seconds()
        last_timestep = current_step
        t += dt

        # Calculate DO concentration in tank layers
        # Get depth to calculate DO
        depth = Wetland.depth
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
        Qin=Wetland.total_inflow
        Cin=sim._model.getNodeCin("93-49759",0)
        Qout=Wetland.total_outflow
        V=Wetland.volume

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
        Wlnd_N = Wetland.pollut_quality['NO']
        Wetland_N.append(Wlnd_N)

        # Get flow for each asset
        Wetland_flow.append(Qout)
        #Wetland_vol.append(V)

        # Control Actions (every 30 mins)
        if _tempcount == 30:
            # If depth <= 6.75 ft then close valve
            if depth <= 6.75:
                Ells_valve.target_setting = 0.0
            # After 36 hours, proportional release but no more than 33% open
            else:    
                Ells_valve.target_setting = min(0.33, Qout/(1.00*np.sqrt(2.0*32.2*depth)))
            _tempcount = 0
        _tempcount+=1
        print(Ells_valve.target_setting)
        

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
# Convert flow rate from cfs to m3/s
conv_cfs_cms = [0.02832]*len(Wetland_flow)
Wetland_flow_m = [a*b for a,b in zip(Wetland_flow,conv_cfs_cms)]

# Convert depth from ft to m
conv_ft_m = [0.3048]*len(Wetland_flow)
Wetland_depth_m = [a*b for a,b in zip(Wetland_depth,conv_ft_m)]

# Calculate load each timestep
conv_mgs_kgs = [0.000001]*len(Wetland_flow)
Nload = [a*b*c*d for a,b,c,d in zip(Wetland_N,Wetland_flow,conv_cfs_cms,conv_mgs_kgs)]
Ncum_load = np.cumsum(Nload)

# SS graphs
#SS = [Cn]*len(Ncum_load3)

#----------------------------------------------------------------------#
#Plot Results
fig, ax = plt.subplots(5, sharex=True)
ax[0].plot(Wetland_N)
#ax[0].plot(SS, 'k', label='SteadyState')
ax[0].set_ylabel("NO (mg/L)")
ax[0].set_title("Wetland")
ax[1].plot(Wetland_DO1, "m", label='DO_L1')
ax[1].plot(Wetland_DO2, label='DO_L2')
ax[1].plot(Wetland_DO3, label='DO_L3')
ax[1].legend(loc="upper right")
ax[1].set_ylabel("DO (mg/L)")
ax[2].plot(Wetland_depth_m)
ax[2].set_ylabel("Depth (m)")
ax[3].plot(Wetland_flow_m)
ax[3].set_ylabel("Flow (m³/s)")
ax[4].plot(Ncum_load)
ax[4].set_ylabel("NO Cum. Load (kg)")
ax[4].set_xlabel("Time (min)")
plt.show()
