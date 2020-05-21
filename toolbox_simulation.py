# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-19 10:06:58

# Import required modules
from pyswmm import Simulation, Nodes, Links
import numpy as np
import matplotlib.pyplot as plt

# Make dictionaries for each water quality method
# C_s from gravity Settling example in SWMM manual, v_s from Troutman et al. 2020
#dict1 = {'Tank1': {1: [0.00419, 20]}}
# From ?? 
#dict3 = {'Channel': {1: [10.0, 0.0001, 2.68, 0.45]}}
# From Reddy et al. 1980
#dict4 = {'Tank2': {0: [0.000015, 0.0030]}}
#dict5 = {'Tank2': {0: [0.000015, 0.0053]}}
#dict6 = {'Tank2': {0: [0.000015, 0.0100]}}

# Lists to store results
ch1 = []
tank1 = []
tank2 = []
ch1_flow =[]
tank1_flow = []
tank2_flow = []
ch1_N = []
tank1_N = []
tank2_N = []

# Setup toolbox simulation
with Simulation("./tanks_channel.inp") as sim:
    # Setup dt calculation
    start_time = sim.start_time
    last_timestep = start_time
    # Get asset information
    Tank1 = Nodes(sim)["Tank1"]
    Tank2 = Nodes(sim)["Tank2"]
    Ch1 = Links(sim)["Channel"]

    # Step through the simulation    
    for index,step in enumerate(sim):
        # Calculate dt
        current_step = sim.current_time
        dt = (current_step - last_timestep).total_seconds()
        last_timestep = current_step
        
        # Calculate gravity settling 
        Qin_T1 = sim._model.getNodeResult("Tank1", 0)
        k_GS = 0.001
        Cs_GS = 20.0
        Cin_T1 = sim._model.getNodeCin("Tank1", 0)
        print(Cin_T1)
        d_T1 = sim._model.getNodeResult("Tank1", 5)
        
        if d_T1 != 0.0:
            Cnew_T1 = np.heaviside((0.1-Qin_T1),0)*(Cs_GS+(Cin_T1-Cs_GS)*np.exp(-k_GS/d_T1*dt/3600))+(1-np.heaviside((0.1-Qin_T1),0))*Cin_T1
        else:
            Cnew_T1 = np.heaviside((0.1-Qin_T1),0)*Cs_GS+(Cin_T1-Cs_GS)+(1-np.heaviside((0.1-Qin_T1),0))*Cin_T1
        sim._model.setNodePollutant("Tank1", 0, Cnew_T1)
        """
        # Calculate erosion produced
        Cin_Ch = sim._model.getLinkC2("Channel", 1)
        Q_Ch = sim._model.getLinkResult("Channel", 0)
        d_Ch = sim._model.getLinkResult("Channel", 1)
        v_Ch = sim._model.getConduitVelocity("Channel")
        w = 10.0
        So = 0.0001
        Ss = 2.68
        d50 = 0.45
        g = 9.81    # m/s^2
        yw = 1000   # kg/m^3
        theta = (d_Ch*So/((Ss-1)*d50))*(1/0.001)   # unitless
        if v_Ch != 0.0:
            f = (2*g*So*d_Ch)/(v_Ch*0.305)**2     # unitless
            qs = 0.1*(1/f)*theta**(5/2)*yw*((Ss-1)*g*(d50*0.001)**3)**(1/2) # kg/m-s
            Qs = w*qs   # kg/s
        else:
            Qs = 0.0
        if Q_Ch != 0.0:
            Cnew_Ch = ((Qs/Q_Ch)*1000)  # mg/L
            Cnew_Ch = max(Cin_Ch, Cin_Ch + Cnew_Ch)
            # Set new concentration
            sim._model.setLinkPollutant("Channel", 1, Cnew_Ch)
        """
        """
        # Calculate DO concentration in tank layers
        # Get depth to calculate DO/Nitrate
        depth = sim._model.getNodeResult("Tank2", 5)
        if depth <= 3.3:
            # Calculate DO concentration in first layer
            DO1 = DO_Layer1.DO_solver(index)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 < 1.0:
                k_ni1 = 0.00014  # [1/s]
            else:
                k_ni1 = 0.0
        elif 3.3 < depth <= 6.6:
            # Calculate DO concentration in first two layers
            DO1 = DO_Layer1.DO_solver(index)
            DO2 = DO_Layer2.DO_solver(index)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 < 1.0:
                k_nil1 = 0.00014
            else:
                k_nil1 = 0.0
            if DO2 < 1.0:
                k_nil2 = 0.00014
            else:
                k_nil2 = 0.0
            k_nil = k_nil1 + k_nil2
        else:
            # Calculate DO concentration in all three layers
            DO1 = DO_Layer1.DO_solver(index)
            DO2 = DO_Layer2.DO_solver(index)
            DO3 = DO_Layer3.DO_solver(index)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 < 1.0:
                k_nil1 = 0.00014
            else:
                k_nil1 = 0.0
            if DO2 < 1.0:
                k_nil2 = 0.00014
            else:
                k_nil2 = 0.0
            if DO3 < 1.0:
                k_nil3 = 0.00014
            else:
                k_nil3 = 0.0
            k_nil = k_nil1 + k_nil2 + k_nil3
        """
        # Get TSS conc for each asset        
        c1 = Ch1.pollut_quality
        ch1.append(c1['P2'])
        c2 = Tank1.pollut_quality
        tank1.append(c2['P2'])
        c3 = Tank2.pollut_quality
        tank2.append(c3['P2'])

        # Get Nitrate conc for each asset        
        c1_N = Ch1.pollut_quality
        ch1_N.append(c1_N['P1'])
        c2_N = Tank1.pollut_quality
        tank1_N.append(c2_N['P1'])
        c3_N = Tank2.pollut_quality
        tank2_N.append(c3_N['P1'])

        # Get flow for  each asset
        ch1_flow.append(sim._model.getLinkResult("Channel", 0))
        tank1_flow.append(sim._model.getNodeResult("Tank1", 0))
        tank2_flow.append(sim._model.getNodeResult("Tank2", 0))


"""        
#----------------------------------------------------------------------#
# Confirm gravity settling matches toolbox simulation
swmm_c = []

with Simulation("./test_gravsettling.inp") as sim:
    Tank = Nodes(sim)["2"]
    # Step through the simulation    
    for step in sim:
        # Get newQual for Tank
        co = Tank.pollut_quality
        swmm_c.append(co['TSS'])
        # Get flow for Tank
        flo2.append(sim._model.getNodeResult("2", 0))

#----------------------------------------------------------------------#
# Confirm erosion mass balance
# Calculate load each timestep
load1 = [a*b for a,b in zip(ch1,ch1_flow)]
load2 = [a*b for a,b in zip(ch1ds,ch1ds_flow)]

# Calculate cumulative load
cum_load1 = np.cumsum(load1)
cum_load2 = np.cumsum(load2)

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

# Calculate loads
load1 = [a*b for a,b in zip(ch1,ch1_flow)]
load2 = [a*b for a,b in zip(tank1,tank1_flow)]
load3 = [a*b for a,b in zip(tank2,tank2_flow)]

Nload1 = [a*b for a,b in zip(ch1_N,ch1_flow)]
Nload2 = [a*b for a,b in zip(tank1_N,tank1_flow)]
Nload3 = [a*b for a,b in zip(tank2_N,tank2_flow)]

# Plot TSS load
plt.subplot(411)
plt.plot(load1, 'g--', label='Channel')
plt.plot(load2, 'b:', label='Tank1')
plt.plot(load3, 'm:', label='Tank2')
plt.ylabel("TSS Load")
plt.xlabel("Time (s)")
plt.legend()

# Plot Nitate load
plt.subplot(412)
plt.plot(Nload1, 'g--', label='Channel')
plt.plot(Nload2, 'b:', label='Tank1')
plt.plot(Nload3, 'm:', label='Tank2')
plt.ylabel("NO Load")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(413)
plt.plot(ch1, 'g--', label='Channel')
plt.plot(tank1, 'b:', label='Tank1')
plt.plot(tank2, 'm:', label='Tank2')
plt.ylabel("TSS Conc")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(414)
plt.plot(ch1_N, 'g--', label='Channel')
plt.plot(tank1_N, 'b:', label='Tank1')
plt.plot(tank2_N, 'm:', label='Tank2')
plt.ylabel("NO Conc")
plt.xlabel("Time (s)")
plt.legend()

plt.show()
