# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-07-09 14:26:57

# Import required modules
from pyswmm import Simulation, Nodes, Links
from wq_toolbox.nodes import Node_Quality
from scipy.integrate import ode 
import numpy as np
import matplotlib.pyplot as plt

# Nitrate 3 CSTRs in Series
def CSTR_tank(t, C, Qin, Cin, Qout, V, k):
    dCdt = (Qin*Cin - Qout*C)/V - k*C
    return dCdt

# DO Variables
k_DO = 0.000075 # rate/5 sec
Co_DO = 10.0    # mg/L

#----------------------------------------------------------------------#

# Uncontrolled Simulation
# Lists to store results
Ellsworth_inflow = []
Ellsworth_conc = []
Ellsworth_depth = []
#Ellsworth_flooding = []
Ellsworth_outflow = []
Ellsworth_cumload = []

DBasin_inflow = []
DBasin_conc = []
DBasin_depth = []
DBasin_outflow = []
DBasin_cumload = []

Wetland_inflow = []
Wetland_conc = []
Wetland_depth = []
#Wetland_flooding= []
Wetland_volume = []
Wetland_outflow = []
Wetland_cumload = []
Wetland_DO1 = []
Wetland_DO2 = []
Wetland_DO3 = []   
#Wtlnd_bp_inflows = []

Channel_flow = []
Channel_conc = []
Channel_depth = []
#Channel_flooding = []
Channel_cumload = []

# Setup toolbox simulation
with Simulation("./modifiedMBDoyle_NO2.inp") as sim:
    # Get asset information
    Ellsworth = Nodes(sim)["93-50408"]
    DBasin = Nodes(sim)["93-50404"]
    Wetland = Nodes(sim)["93-49759"]
    Wtlnd_bypass = Links(sim)["95-70294"]
    Channel = Links(sim)["95-70277"]
    
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

    # Tracking time for DO reaction
    t1 = 0
    t2 = 0
    t3 = 0

    # Step through the simulation    
    for index,step in enumerate(sim):

        # Calculate dt
        current_step = sim.current_time
        dt = (current_step - last_timestep).total_seconds()
        last_timestep = current_step

        # Get NO conc for each asset        
        Ell_p = Ellsworth.pollut_quality['NO']
        Ellsworth_conc.append(Ell_p)
        DB_p = DBasin.pollut_quality['NO']
        DBasin_conc.append(DB_p)
        Wt_p = Wetland.pollut_quality['NO']
        Wetland_conc.append(Wt_p)
        Ch_p = Channel.pollut_quality['NO']
        Channel_conc.append(Ch_p)

        # Get flow for  each asset
        Ell_if = Ellsworth.total_inflow
        Ellsworth_inflow.append(Ell_if)
        Ell_d = Ellsworth.depth
        Ellsworth_depth.append(Ell_d)
        #Ell_fl = Ellsworth.flooding
        #Ellsworth_flooding.append(Ell_fl)
        Ell_of = Ellsworth.total_outflow
        Ellsworth_outflow.append(Ell_of)
        DB_if = DBasin.total_inflow
        DBasin_inflow.append(DB_if)
        DB_d = DBasin.depth
        DBasin_depth.append(DB_d)
        DB_of = DBasin.total_outflow
        DBasin_outflow.append(DB_of)
        Wt_if = Wetland.total_inflow
        Wetland_inflow.append(Wt_if)
        Wt_d = Wetland.depth
        Wetland_depth.append(Wt_d)
        #Wt_fl = Wetland.flooding
        #Wetland_flooding.append(Wt_fl)
        Wt_v = Wetland.volume
        Wetland_volume.append(Wt_v)
        Wt_of = Wetland.total_outflow
        Wetland_outflow.append(DB_of)
        #Wt_bp = Wtlnd_bypass.flow
        #Wtlnd_bp_inflows.append(Wt_bp)
        Ch_f = Channel.flow
        Channel_flow.append(Ch_f)
        Ch_d = Channel.depth
        Channel_depth.append(Ch_d)

        # Calculate DO concentration in tank layers
        # reset DO if tank is empty
        if DB_d <= 0.01:
            t1 = 0
            t2 = 0
            t3 = 0
            DO1 = 0.0
            DO2 = 0.0
            DO3 = 0.0
            Wetland_DO1.append(DO1)
            Wetland_DO2.append(DO2)
            Wetland_DO3.append(DO3)
            k_ni = 0.0
        
        elif 0.01 < DB_d <= 3.00:
            # Calculate DO concentration in first layer
            t1 += dt
            t2 = 0
            t3 = 0
            DO1 = Co_DO*np.exp(-k_DO*t1)
            DO2 = 0.0
            DO3 = 0.0
            # Calculate nitate reaction rate based on DO concentration
            if DO1 <= 1.0:
                k_ni = 0.000029  # [1/5 sec]
            else:
                k_ni = 0.0
            Wetland_DO1.append(DO1)
            Wetland_DO2.append(DO2)
            Wetland_DO3.append(DO3)
        
        elif 3.0 < DB_d <= 6.0:
            # Calculate DO concentration in first two layers
            t1 += dt
            t2 += dt
            t3 = 0
            DO1 = Co_DO*np.exp(-k_DO*t1)
            DO2 = Co_DO*np.exp(-k_DO*t2)
            DO3 = 0.0
            # Calculate nitate reaction rate based on DO concentration
            if DO1 <= 1.0:
                k_ni1 = 0.000029
            else:
                k_ni1 = 0.0
            if DO2 <= 1.0:
                k_ni2 = 0.000029
            else:
                k_ni2 = 0.0
            Wetland_DO1.append(DO1)
            Wetland_DO2.append(DO2)
            Wetland_DO3.append(DO3)
        else:
            # Calculate DO concentration in all three layers
            t1 += dt
            t2 += dt
            t3 += dt
            DO1 = Co_DO*np.exp(-k_DO*t1)
            DO2 = Co_DO*np.exp(-k_DO*t2)
            DO3 = Co_DO*np.exp(-k_DO*t3)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 <= 1.0:
                k_ni1 = 0.000029
            else:
                k_ni1 = 0.0
            if DO2 <= 1.0:
                k_ni2 = 0.000029
            else:
                k_ni2 = 0.0
            if DO3 <= 1.0:
                k_ni3 = 0.000029
            else:
                k_ni3 = 0.0
            k_ni = k_ni1 + k_ni2 + k_ni3
            Wetland_DO1.append(DO1)
            Wetland_DO2.append(DO2)
            Wetland_DO3.append(DO3)
        
        # Calculate NO concentration in tanks
        # Get parameters to calculate NO
        Cin=sim._model.getNodeCin("93-49759",0)

        #Solve ODE
        if index == 0:
            solver1.set_initial_value(0.0, 0.0)
            solver1.set_f_params(Wt_if,Cin,Wt_of,Wt_v,k_ni)
            solver1.integrate(solver1.t+dt)
            solver2.set_initial_value(0.0, 0.0)
            solver2.set_f_params(Wt_if,Cin,Wt_of,Wt_v,k_ni)
            solver2.integrate(solver2.t+dt)
            solver3.set_initial_value(0.0, 0.0)
            solver3.set_f_params(Wt_if,Cin,Wt_of,Wt_v,k_ni)
            solver3.integrate(solver3.t+dt)
        else:
            solver1.set_initial_value(solver1.y, solver1.t)
            solver1.set_f_params(Wt_if,Cin,Wt_of,Wt_v,k_ni)
            solver1.integrate(solver1.t+dt)
            solver2.set_initial_value(solver2.y, solver2.t)
            solver2.set_f_params(Wt_if,solver1.y,Wt_of,Wt_v,k_ni)
            solver2.integrate(solver2.t+dt)
            solver3.set_initial_value(solver3.y, solver3.t)
            solver3.set_f_params(Wt_if,solver2.y,Wt_of,Wt_v,k_ni)
            solver3.integrate(solver3.t+dt)
        
        # Set new concentration
        sim._model.setNodePollutant("93-49759", 0, solver3.y[0])

    sim._model.swmm_end()
    print(sim.runoff_error)
    print(sim.flow_routing_error)
    print(sim.quality_error)

# Convert inflow rate from cfs to m3/s
conv_cfs_cms = [0.02832]*len(Ellsworth_inflow)
Ellsworth_inflow_m = [a*b for a,b in zip(Ellsworth_inflow,conv_cfs_cms)]
DBasin_inflow_m = [a*b for a,b in zip(DBasin_inflow,conv_cfs_cms)]
Wetland_inflow_m = [a*b for a,b in zip(Wetland_inflow,conv_cfs_cms)]
Channel_flow_m = [a*b for a,b in zip(Channel_flow,conv_cfs_cms)]

# Convert outflow rate from cfs to m3/s
conv_cfs_cms = [0.02832]*len(Ellsworth_inflow)
Ellsworth_outflow_m = [a*b for a,b in zip(Ellsworth_outflow,conv_cfs_cms)]
DBasin_outflow_m = [a*b for a,b in zip(DBasin_outflow,conv_cfs_cms)]
Wetland_outflow_m = [a*b for a,b in zip(Wetland_outflow,conv_cfs_cms)]

# Convert flooding rate from cfs to m3/s
#conv_cfs_cms = [0.02832]*len(Ellsworth_inflow)
#Ellsworth_flooding_m = [a*b for a,b in zip(Ellsworth_flooding,conv_cfs_cms)]
#Wetland_flooding_m = [a*b for a,b in zip(Wetland_flooding,conv_cfs_cms)]
#Wtlnd_bypass_m = [a*b for a,b in zip(Wtlnd_bp_inflows,conv_cfs_cms)]

# Convert depth from ft to m
conv_ft_m = [0.3048]*len(Ellsworth_inflow)
Ellsworth_depth_m = [a*b for a,b in zip(Ellsworth_depth,conv_ft_m)]
DBasin_depth_m = [a*b for a,b in zip(DBasin_depth,conv_ft_m)]
Wetland_depth_m = [a*b for a,b in zip(Wetland_depth,conv_ft_m)]
Channel_depth_m = [a*b for a,b in zip(Channel_depth,conv_ft_m)]

# Calculate load each timestep
conv_mgs_kgs = [0.000001]*len(Ellsworth_inflow)
timestep = [5]*len(Ellsworth_inflow)
Ellsworth_load = [a*b*c*d*e for a,b,c,d,e in zip(Ellsworth_conc,Ellsworth_outflow,conv_cfs_cms, conv_mgs_kgs,timestep)]
DBasin_load = [a*b*c*d*e for a,b,c,d,e in zip(DBasin_conc,DBasin_outflow,conv_cfs_cms,conv_mgs_kgs,timestep)]
Wetland_load = [a*b*c*d*e for a,b,c,d,e in zip(Wetland_conc,Wetland_outflow,conv_cfs_cms, conv_mgs_kgs,timestep)]
Channel_load = [a*b*c*d*e for a,b,c,d,e in zip(Channel_conc,Channel_flow,conv_cfs_cms,conv_mgs_kgs,timestep)]

# Calculate cumulative load (dt = 1)
Ellsworth_cumload = np.cumsum(Ellsworth_load)
DBasin_cumload = np.cumsum(DBasin_load)
Wetland_cumload = np.cumsum(Wetland_load)
Channel_cumload = np.cumsum(Channel_load)

#----------------------------------------------------------------------#

# Print final load released
print("Ellsworth:", Ellsworth_cumload[-1])
print("Doyle Basin:", DBasin_cumload[-1])
print("Wetland:", Wetland_cumload[-1])
print("Channel to Outfall:", Channel_cumload[-1]) 
print("Percent Load Reduction", ((1.6770261144428948-Channel_cumload[-1])/1.6770261144428948))

#----------------------------------------------------------------------#

# Data for flooding line
Ell_flood = [6.0960]*len(Ellsworth_inflow)
Wt_bypass = [2.8956]*len(Ellsworth_inflow)
Wt_flood  = [2.7432]*len(Ellsworth_inflow)

# Plot Result
fig, ax = plt.subplots(7, 4)
ax[0,0].plot(Ellsworth_inflow_m, 'k--')
ax[0,0].set_xticks([])
ax[0,0].set_yticks([0,3,6,9])
ax[0,0].set_yticklabels(["0","3","6","9"], size=6)
ax[0,0].set_ylim(-0.1,9.05)
ax[0,0].set_xlim(0,190080)
ax[0,0].set_ylabel("Inflow (m³/s)", size=6)

ax[1,0].plot(Ellsworth_conc, 'k--')
ax[1,0].set_xticks([])
ax[1,0].set_yticks([0,20,40])
ax[1,0].set_yticklabels(["0","20","40"], size=6)
ax[1,0].set_ylim(-0.1,50)
ax[1,0].set_xlim(0,190080)
ax[1,0].set_ylabel("NO (mg/L)", size=6)

ax[2,0].set_ylabel("DO (mg/L)", size=6)
ax[2,0].set_xticks([])
ax[2,0].set_yticks([0,5,10])
ax[2,0].set_yticklabels(["0","5","10"], size=6)
ax[2,0].set_ylim(-0.1,10.5)
ax[2,0].set_xlim(0,190080)

ax[3,0].plot(Ell_flood, color='#818282')
ax[3,0].plot(Ellsworth_depth_m, 'k--')
ax[3,0].set_xticks([])
ax[3,0].set_yticks([0,3,6])
ax[3,0].set_yticklabels(["0","3","6"], size=6)
ax[3,0].set_ylim(-0.1,6.5)
ax[3,0].set_xlim(0,190080)
ax[3,0].set_ylabel("Depth (m)", size=6)

ax[4,0].set_xticks([])
ax[4,0].set_ylabel("Valve Position", size=6)
ax[4,0].set_yticks([0,0.5,1])
ax[4,0].set_yticklabels(["0","0.5","1.0"], size=6)
ax[4,0].set_ylim(-0.1,1.05)
ax[4,0].set_xlim(0,190080)

ax[5,0].plot(Ellsworth_outflow_m, 'k--')
ax[5,0].set_xticks([])
ax[5,0].set_yticks([0,4,8])
ax[5,0].set_yticklabels(["0","4","8"], size=6)
ax[5,0].set_ylim(-0.1,8.5)
ax[5,0].set_xlim(0,190080)
ax[5,0].set_ylabel("Outflow (m³/s)", size=6)

ax[6,0].plot(Ellsworth_cumload, 'k--')
ax[6,0].set_yticks([0,2,4])
ax[6,0].set_yticklabels(["0","2","4"], size=6)
ax[6,0].set_ylim(-0.1,4)
ax[6,0].set_xlim(0,190080)
ax[6,0].set_xticks([0,34560,69120,103680,138240,172800])
ax[6,0].set_xticklabels(["0","2","4","6","8","10"], size=6)
ax[6,0].set_ylabel("Cum. Load (kg)", size=6)
ax[6,0].set_xlabel("Time (days)", size=6)

ax[0,1].plot(DBasin_inflow_m, "k--")
ax[0,1].set_xticks([])
ax[0,1].set_yticks([0,3,6,9])
ax[0,1].set_yticklabels(["0","3","6","9"], size=6)
ax[0,1].set_ylim(-0.1,9.05)
ax[0,1].set_xlim(0,190080)

ax[1,1].plot(DBasin_conc, "k--")
ax[1,1].set_xticks([])
ax[1,1].set_yticks([0,20,40])
ax[1,1].set_yticklabels(["0","20","40"], size=6)
ax[1,1].set_ylim(-0.1,50)
ax[1,1].set_xlim(0,190080)

ax[2,1].set_xticks([])
ax[2,1].set_yticks([0,5,10])
ax[2,1].set_yticklabels(["0","5","10"], size=6)
ax[2,1].set_ylim(-0.1,10.5)
ax[2,1].set_xlim(0,190080)

ax[3,1].plot(Wt_bypass, color='#818282')
ax[3,1].plot(DBasin_depth_m, "k--")
ax[3,1].set_xticks([])
ax[3,1].set_yticks([0,1.5,3])
ax[3,1].set_yticklabels(["0","1.5","3"], size=6)
ax[3,1].set_ylim(-0.1,3.25)
ax[3,1].set_xlim(0,190080)

ax[4,1].set_xticks([])
ax[4,1].set_yticks([0,0.5,1])
ax[4,1].set_yticklabels(["0","0.5","1.0"], size=6)
ax[4,1].set_ylim(-0.1,1.05)
ax[4,1].set_xlim(0,190080)

ax[5,1].plot(DBasin_outflow_m, "k--")
ax[5,1].set_xticks([])
ax[5,1].set_yticks([0,4,8])
ax[5,1].set_yticklabels(["0","4","8"], size=6)
ax[5,1].set_ylim(-0.1,8.5)
ax[5,1].set_xlim(0,190080)

ax[6,1].plot(DBasin_cumload, "k--")
ax[6,1].set_yticks([0,2,4])
ax[6,1].set_yticklabels(["0","2","4"], size=6)
ax[6,1].set_ylim(-0.1,4)
ax[6,1].set_xlim(0,190080)
ax[6,1].set_xticks([0,34560,69120,103680,138240,172800])
ax[6,1].set_xticklabels(["0","2","4","6","8","10"], size=6)
ax[6,1].set_xlabel("Time (days)", size=6)

ax[0,2].plot(Wetland_inflow_m, "k--")
ax[0,2].set_xticks([])
ax[0,2].set_yticks([0,3,6,9])
ax[0,2].set_yticklabels(["0","3","6","9"], size=6)
ax[0,2].set_ylim(-0.1,9.05)
ax[0,2].set_xlim(0,190080)

ax[1,2].plot(Wetland_conc, "k--")
ax[1,2].set_xticks([])
ax[1,2].set_yticks([0,20,40])
ax[1,2].set_yticklabels(["0","20","40"], size=6)
ax[1,2].set_ylim(-0.1,50)
ax[1,2].set_xlim(0,190080)

ax[2,2].plot(Wetland_DO1, linestyle="--", color="#B08CA1")
ax[2,2].plot(Wetland_DO2, linestyle="--", color="#99798C")
ax[2,2].plot(Wetland_DO3, linestyle="--", color="#493C44")
ax[2,2].set_xticks([])
ax[2,2].set_yticks([0,5,10])
ax[2,2].set_yticklabels(["0","5","10"], size=6)
ax[2,2].set_ylim(-0.1,10.5)
ax[2,2].set_xlim(0,190080)

ax[3,2].plot(Wt_flood, color='#818282')
ax[3,2].plot(Wetland_depth_m, "k--")
ax[3,2].set_xticks([])
ax[3,2].set_yticks([0,1.5,3])
ax[3,2].set_yticklabels(["0","1.5","3"], size=6)
ax[3,2].set_ylim(-0.1,3.25)
ax[3,2].set_xlim(0,190080)

ax[4,2].set_xticks([])
ax[4,2].set_yticks([0,0.5,1])
ax[4,2].set_yticklabels(["0","0.5","1.0"], size=6)
ax[4,2].set_ylim(-0.1,1.05)
ax[4,2].set_xlim(0,190080)

ax[5,2].plot(Wetland_outflow_m, "k--")
ax[5,2].set_xticks([])
ax[5,2].set_yticks([0,4,8])
ax[5,2].set_yticklabels(["0","4","8"], size=6)
ax[5,2].set_ylim(-0.1,8.5)
ax[5,2].set_xlim(0,190080)

ax[6,2].plot(Wetland_cumload, "k--")
ax[6,2].set_yticks([0,2,4])
ax[6,2].set_yticklabels(["0","2","4"], size=6)
ax[6,2].set_ylim(-0.1,4)
ax[6,2].set_xlim(0,190080)
ax[6,2].set_xticks([0,34560,69120,103680,138240,172800])
ax[6,2].set_xticklabels(["0","2","4","6","8","10"], size=6)
ax[6,2].set_xlabel("Time (days)", size=6)

ax[0,3].plot(Channel_flow_m, "k--")
ax[0,3].set_xticks([])
ax[0,3].set_yticks([0,3,6,9])
ax[0,3].set_yticklabels(["0","3","6","9"], size=6)
ax[0,3].set_ylim(-0.1,9.05)
ax[0,3].set_xlim(0,190080)

ax[1,3].plot(Channel_conc, "k--")
ax[1,3].set_xticks([])
ax[1,3].set_yticks([0,20,40])
ax[1,3].set_yticklabels(["0","20","40"], size=6)
ax[1,3].set_ylim(-0.1,50)
ax[1,3].set_xlim(0,190080)

ax[2,3].set_xticks([])
ax[2,3].set_yticks([0,5,10])
ax[2,3].set_yticklabels(["0","5","10"], size=6)
ax[2,3].set_ylim(-0.1,10.5)
ax[2,3].set_xlim(0,190080)

ax[3,3].plot(Channel_depth_m, "k--")
ax[3,3].set_xticks([])
ax[3,3].set_yticks([0,1.5,3])
ax[3,3].set_yticklabels(["0","1.5","3"], size=6)
ax[3,3].set_ylim(-0.1,3.25)
ax[3,3].set_xlim(0,190080)

ax[4,3].set_xticks([])
ax[4,3].set_yticks([0,0.5,1])
ax[4,3].set_yticklabels(["0","0.5","1.0"], size=6)
ax[4,3].set_ylim(-0.1,1.05)
ax[4,3].set_xlim(0,190080)

ax[5,3].plot(Channel_flow_m, "k--")
ax[5,3].set_xticks([])
ax[5,3].set_yticks([0,4,8])
ax[5,3].set_yticklabels(["0","4","8"], size=6)
ax[5,3].set_ylim(-0.1,8.5)
ax[5,3].set_xlim(0,190080)

ax[6,3].plot(Channel_cumload, "k--")
ax[6,3].set_yticks([0,2,4])
ax[6,3].set_yticklabels(["0","2","4"], size=6)
ax[6,3].set_ylim(-0.1,4)
ax[6,3].set_xlim(0,190080)
ax[6,3].set_xticks([0,34560,69120,103680,138240,172800])
ax[6,3].set_xticklabels(["0","2","4","6","8","10"], size=6)
ax[6,3].set_xlabel("Time (days)", size=6)

fig.savefig('NO.svg', format="svg", bbox_inches="tight",pad_inches=0.1)
plt.show()
