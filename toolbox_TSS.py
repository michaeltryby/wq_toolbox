# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-06-09 14:09:13

# Import required modules
from pyswmm import Simulation, Nodes, Links
from wq_toolbox.links import Link_Quality
from wq_toolbox.nodes import Node_Quality
import numpy as np
import matplotlib.pyplot as plt

# Make dictionaries for each water quality method
# Gravity Settling
# v_s from Troutman et al. 2020, C_s from Mallets Creek data
dict1 = {'93-50408': {0: [15.08, 21.0]},'93-50404': {0: [15.08, 21.0]}, \
         '93-49759': {0: [15.08, 21.0]}}
dict3 = {'95-69044': {0: [15.08, 21.0]}, '95-70180': {0: [15.08, 21.0]}, \
         '95-51776': {0: [15.08, 21.0]}, '95-51774': {0: [15.08, 21.0]}, \
         '95-51634': {0: [15.08, 21.0]}, '95-69048': {0: [15.08, 21.0]}, \
         '95-70594': {0: [15.08, 21.0]}, '95-51760': {0: [15.08, 21.0]}, \
         '95-69050': {0: [15.08, 21.0]}, '95-51758': {0: [15.08, 21.0]}, \
         '95-51757': {0: [15.08, 21.0]}, '95-70713': {0: [15.08, 21.0]}, \
         '95-70277': {0: [15.08, 21.0]}}

# Channel Erosion
# width and slope from SWMM model
# SS and d50 from Engelund and Hansel (1967) 
dict2 = {'95-69044': {0: [1, 0.000207, 2.68, 2.0]}, \
         '95-70180': {0: [1, 0.000411, 2.68, 2.0]}, \
         '95-51776': {0: [1, 0.000382, 2.68, 2.0]}, \
         '95-51774': {0: [1, 0.000133, 2.68, 2.0]}, \
         '95-51634': {0: [1, 0.000169, 2.68, 2.0]}, \
         '95-69048': {0: [1, 0.000168, 2.68, 2.0]}, \
         '95-70594': {0: [1, 0.000172, 2.68, 2.0]}, \
         '95-51760': {0: [1, 0.000374, 2.68, 2.0]}, \
         '95-69050': {0: [1, 0.000138, 2.68, 2.0]}, \
         '95-51758': {0: [1, 0.000372, 2.68, 2.0]}, \
         '95-51757': {0: [1, 0.000297, 2.68, 2.0]}, \
         '95-70713': {0: [1, 0.000180, 2.68, 2.0]}, \
         '95-70277': {0: [1, 0.000376, 2.68, 2.0]}}

# DO Variables
k_DO = 0.000015 # rate/min
Co_DO1 = 10.0   # mg/L
Co_DO2 = 5.5    # mg/L
Co_DO3 = 3.0    # mg/L

#----------------------------------------------------------------------#
# Uncontrolled Simulation
# Lists to store results
Ellsworth_inflow = []
Ellsworth_conc = []
Ellsworth_cuminload = []
Ellsworth_depth = []
Ellsworth_outflow = []
Ellsworth_cumload = []

DBasin_inflow = []
DBasin_conc = []
DBasin_cuminload = []
DBasin_depth = []
DBasin_outflow = []
DBasin_cumload = []

Wetland_inflow = []
Wetland_conc = []
Wetland_cuminload = []
Wetland_depth = []
Wetland_outflow = []
Wetland_cumload = []
Wetland_DO1 = []
Wetland_DO2 = []
Wetland_DO3 = []

Channel_flow = []
Channel_conc = []
Cannel_cuminload = []
Channel_depth = []
Channel_cumload = []

# Setup toolbox simulation
with Simulation("./modifiedMBDoyle_TSS_AArain.inp") as sim:
    # Setup toolbox methods
    GS = Node_Quality(sim, dict1)
    ER = Link_Quality(sim, dict2)
    GS2 = Link_Quality(sim, dict3)
    
    # Get asset information
    Ellsworth = Nodes(sim)["93-50408"]
    DBasin = Nodes(sim)["93-50404"]
    Wetland = Nodes(sim)["93-49759"]
    Channel = Links(sim)["95-70277"]

    # Tracking time for DO reaction
    t = 0

    # Step through the simulation    
    for index,step in enumerate(sim):

        # Calculate gravity settling 
        GS.GravitySettling()
        GS2.GravitySettling()
        # Calculate erosion produced
        ER.Erosion()

        # Get TSS conc for each asset        
        Ell_p = Ellsworth.pollut_quality['TSS']
        Ellsworth_conc.append(Ell_p)
        DB_p = DBasin.pollut_quality['TSS']
        DBasin_conc.append(DB_p)
        Ch_p = Channel.pollut_quality['TSS']
        Channel_conc.append(Ch_p)
        Wt_p = Wetland.pollut_quality['TSS']
        Wetland_conc.append(Wt_p)

        # Get data for each asset
        Ell_if = Ellsworth.total_inflow
        Ellsworth_inflow.append(Ell_if)
        Ell_d = Ellsworth.depth
        Ellsworth_depth.append(Ell_d)
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
        Wt_of = Wetland.total_outflow
        Wetland_outflow.append(DB_of)
        Ch_f = Channel.flow
        Channel_flow.append(Ch_f)
        Ch_d = Channel.depth
        Channel_depth.append(Ch_d)

        # Calculate DO concentration in tank layers
        # reset DO if tank is empty
        if DB_d <= 0.01:
            t = 0

        if DB_d <= 3.00:
            # Calculate DO concentration in first layer
            DO1 = Co_DO1*np.exp(-k_DO*t)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 < 1.0:
                k_ni = 0.000006  # [1/s]
            else:
                k_ni = 0.0
            Wetland_DO1.append(DO1)
        
        elif 3.0 < DB_d <= 6.0:
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
        t+=1

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

# Convert depth from ft to m
conv_ft_m = [0.3048]*len(Ellsworth_inflow)
Ellsworth_depth_m = [a*b for a,b in zip(Ellsworth_depth,conv_ft_m)]
DBasin_depth_m = [a*b for a,b in zip(DBasin_depth,conv_ft_m)]
Wetland_depth_m = [a*b for a,b in zip(Wetland_depth,conv_ft_m)]
Channel_depth_m = [a*b for a,b in zip(Channel_depth,conv_ft_m)]

# Calculate inflow load each timestep
#conv_mgs_kgs = [0.000001]*len(Ellsworth_inflow)
#Ellsworth_inload = [a*b*c*d for a,b,c,d in zip(Ellsworth_conc,Ellsworth_inflow,conv_cfs_cms, conv_mgs_kgs)]
#DBasin_inload = [a*b*c*d for a,b,c,d in zip(DBasin_conc,DBasin_inflow,conv_cfs_cms,conv_mgs_kgs)]
#Wetland_inload = [a*b*c*d for a,b,c,d in zip(Wetland_conc,Wetland_inflow,conv_cfs_cms, conv_mgs_kgs)]
#Channel_inload = [a*b*c*d for a,b,c,d in zip(Channel_conc,Channel_flow,conv_cfs_cms,conv_mgs_kgs)]

# Calculate cumulative inflow load (dt = 1)
#Ellsworth_cuminload = np.cumsum(Ellsworth_inload)
#DBasin_cuminload = np.cumsum(DBasin_inload)
#Wetland_cuminload = np.cumsum(Wetland_inload)
#Channel_cuminload = np.cumsum(Channel_inload)

# Calculate load each timestep
conv_mgs_kgs = [0.000001]*len(Ellsworth_inflow)
Ellsworth_load = [a*b*c*d for a,b,c,d in zip(Ellsworth_conc,Ellsworth_outflow,conv_cfs_cms, conv_mgs_kgs)]
DBasin_load = [a*b*c*d for a,b,c,d in zip(DBasin_conc,DBasin_outflow,conv_cfs_cms,conv_mgs_kgs)]
Wetland_load = [a*b*c*d for a,b,c,d in zip(Wetland_conc,Wetland_outflow,conv_cfs_cms, conv_mgs_kgs)]
Channel_load = [a*b*c*d for a,b,c,d in zip(Channel_conc,Channel_flow,conv_cfs_cms,conv_mgs_kgs)]

# Calculate cumulative load (dt = 1)
Ellsworth_cumload = np.cumsum(Ellsworth_load)
DBasin_cumload = np.cumsum(DBasin_load)
Wetland_cumload = np.cumsum(Wetland_load)
Channel_cumload = np.cumsum(Channel_load)

#----------------------------------------------------------------------#
# Controlled Simulation 
# Lists to store results
Ellsworth_inflowC = []
Ellsworth_concC = []
Ellsworth_cuminloadC = []
Ellsworth_depthC = []
Ellsworth_valveC = []
Ellsworth_outflowC = []
Ellsworth_cumloadC = []

DBasin_inflowC = []
DBasin_concC = []
DBasin_cuminloadC = []
DBasin_depthC = []
DBasin_outflowC = []
DBasin_cumloadC = []

Wetland_inflowC = []
Wetland_concC = []
Wetland_cuminloadC = []
Wetland_depthC = []
Wetland_valveC = []
Wetland_outflowC = []
Wetland_cumloadC = []
Wetland_DO1C = []
Wetland_DO2C = []
Wetland_DO3C = []

Channel_flowC = []
Channel_concC = []
Cannel_cuminloadC = []
Channel_depthC = []
Channel_cumloadC = []

# Setup toolbox simulation
with Simulation("./modifiedMBDoyle_TSS_AArain.inp") as sim:
    # Setup toolbox methods
    GS = Node_Quality(sim, dict1)
    ER = Link_Quality(sim, dict2)
    GS2 = Link_Quality(sim, dict3)
    
    # Get asset information
    Ellsworth = Nodes(sim)["93-50408"]
    Ells_valve = Links(sim)["95-70951"]
    DBasin = Nodes(sim)["93-50404"]
    Channel = Links(sim)["95-70277"]
    Wetland = Nodes(sim)["93-49759"]
    Wtlnd_valve = Links(sim)["95-70293"]

    # Tracking time for control actions every 15 minutes
    _tempcount = 15

    # Tracking time for DO reaction
    t = 0

    # Step through the simulation    
    for index,step in enumerate(sim):

        # Calculate gravity settling 
        GS.GravitySettling()
        GS2.GravitySettling()
        # Calculate erosion produced
        ER.Erosion()

        # Get TSS conc for each asset        
        Ell_p = Ellsworth.pollut_quality['TSS']
        Ellsworth_concC.append(Ell_p)
        DB_p = DBasin.pollut_quality['TSS']
        DBasin_concC.append(DB_p)
        Ch_p = Channel.pollut_quality['TSS']
        Channel_concC.append(Ch_p)
        Wt_p = Wetland.pollut_quality['TSS']
        Wetland_concC.append(Wt_p)

        # Get data each asset
        Ell_if = Ellsworth.total_inflow
        Ellsworth_inflowC.append(Ell_if)
        Ell_d = Ellsworth.depth
        Ellsworth_depthC.append(Ell_d)
        Ell_of = Ellsworth.total_outflow
        Ellsworth_outflowC.append(Ell_of)
        DB_if = DBasin.total_inflow
        DBasin_inflowC.append(DB_if)
        DB_d = DBasin.depth
        DBasin_depthC.append(DB_d)
        DB_of = DBasin.total_outflow
        DBasin_outflowC.append(DB_of)
        Wt_if = Wetland.total_inflow
        Wetland_inflowC.append(Wt_if)
        Wt_d = Wetland.depth
        Wetland_depthC.append(Wt_d)
        Wt_of = Wetland.total_outflow
        Wetland_outflowC.append(DB_of)
        Ch_f = Channel.flow
        Channel_flowC.append(Ch_f)
        Ch_d = Channel.depth
        Channel_depthC.append(Ch_d)

        # Calculate DO concentration in tank layers
        # reset DO if tank is empty
        if Wt_d <= 0.01:
            t = 0

        if Wt_d <= 3.00:
            # Calculate DO concentration in first layer
            DO1 = Co_DO1*np.exp(-k_DO*t)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 < 1.0:
                k_ni = 0.000006  # [1/s]
            else:
                k_ni = 0.0
            Wetland_DO1C.append(DO1)
        
        elif 3.0 < Wt_d <= 6.0:
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
            Wetland_DO1C.append(DO1)
            Wetland_DO2C.append(DO2)
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
            Wetland_DO1C.append(DO1)
            Wetland_DO2C.append(DO2)
            Wetland_DO3C.append(DO3)
        t+=1
        
        # Wetland Control Actions (every 15 mins)
        if _tempcount == 15:
            # If DO levels above anoxic condition:
            if (DO1 or DO2 or DO3) > 1.0:
                # And if Wetland depth <= 6.75 ft then close valve
                if Wt_d <= 6.75:
                    Wtlnd_valve.target_setting = 0.0
                # Else proportional release but no more than 33% open
                else:
                    Wtlnd_valve.target_setting = min(0.50, 0.05*np.sqrt(2.0*32.2*Wt_d))
            # Else proportional release but no more than 33% open
            else:    
                Wtlnd_valve.target_setting = min(0.50, 0.05*np.sqrt(2.0*32.2*Wt_d))
        
        #Wtlnd_valve.target_setting = min(0.25, 0.05*np.sqrt(2.0*32.2*Wt_d))
        Wetland_valveC.append(Wtlnd_valve.target_setting)

        
        # Ellsworth Control Actions (every 15 mins)
        if _tempcount == 15:
            # If Ellsworth TSS conc < 20
            if Ell_p > 20.0:
                # And if Ellsworth depth <= 15 ft & DBasin depth >= 6 ft (~75% capacity) close valve
                if Ell_d <= 15 and DB_d >= 6:
                    Ells_valve.target_setting = 0.0
                # Else proportional release but no more than 33% open
                else:
                    Ells_valve.target_setting = min(0.33, 0.05*np.sqrt(2.0*32.2*Ell_d))
            # Else proportional release but no more than 33% open
            else:  
                Ells_valve.target_setting = min(0.33, 0.05*np.sqrt(2.0*32.2*Ell_d))
            _tempcount = 0
        _tempcount += 1
        
        #Ells_valve.target_setting = min(0.25, 0.05*np.sqrt(2.0*32.2*Ell_d))
        Ellsworth_valveC.append(Ells_valve.target_setting)

    #print(sim.runoff_error)
    #print(sim.flow_routing_error)
    #print(sim.quality_error)
        
# Convert inflow rate from cfs to m3/s
conv_cfs_cms = [0.02832]*len(Ellsworth_inflowC)
Ellsworth_inflow_mC = [a*b for a,b in zip(Ellsworth_inflowC,conv_cfs_cms)]
DBasin_inflow_mC = [a*b for a,b in zip(DBasin_inflowC,conv_cfs_cms)]
Wetland_inflow_mC = [a*b for a,b in zip(Wetland_inflowC,conv_cfs_cms)]
Channel_flow_mC = [a*b for a,b in zip(Channel_flowC,conv_cfs_cms)]

# Convert outflow rate from cfs to m3/s
conv_cfs_cms = [0.02832]*len(Ellsworth_inflowC)
Ellsworth_outflow_mC = [a*b for a,b in zip(Ellsworth_outflowC,conv_cfs_cms)]
DBasin_outflow_mC = [a*b for a,b in zip(DBasin_outflowC,conv_cfs_cms)]
Wetland_outflow_mC = [a*b for a,b in zip(Wetland_outflowC,conv_cfs_cms)]

# Convert depth from ft to m
conv_ft_m = [0.3048]*len(Ellsworth_inflowC)
Ellsworth_depth_mC = [a*b for a,b in zip(Ellsworth_depthC,conv_ft_m)]
DBasin_depth_mC = [a*b for a,b in zip(DBasin_depthC,conv_ft_m)]
Wetland_depth_mC = [a*b for a,b in zip(Wetland_depthC,conv_ft_m)]
Channel_depth_mC = [a*b for a,b in zip(Channel_depthC,conv_ft_m)]

# Calculate inflow load each timestep
#conv_mgs_kgs = [0.000001]*len(Ellsworth_inflowC)
#Ellsworth_inloadC = [a*b*c*d for a,b,c,d in zip(Ellsworth_concC,Ellsworth_inflowC,conv_cfs_cms, conv_mgs_kgs)]
#DBasin_inloadC = [a*b*c*d for a,b,c,d in zip(DBasin_concC,DBasin_inflowC,conv_cfs_cms,conv_mgs_kgs)]
#Wetland_inloadC = [a*b*c*d for a,b,c,d in zip(Wetland_concC,Wetland_inflowC,conv_cfs_cms, conv_mgs_kgs)]
#Channel_inloadC = [a*b*c*d for a,b,c,d in zip(Channel_concC,Channel_flowC,conv_cfs_cms,conv_mgs_kgs)]

# Calculate cumulative load (dt = 1)
#Ellsworth_cuminloadC = np.cumsum(Ellsworth_inloadC)
#DBasin_cuminloadC = np.cumsum(DBasin_inloadC)
#Wetland_cuminloadC = np.cumsum(Wetland_inloadC)
#Channel_cuminloadC = np.cumsum(Channel_inloadC)

# Calculate outflow load each timestep
conv_mgs_kgs = [0.000001]*len(Ellsworth_inflowC)
Ellsworth_loadC = [a*b*c*d for a,b,c,d in zip(Ellsworth_concC,Ellsworth_outflowC,conv_cfs_cms, conv_mgs_kgs)]
DBasin_loadC = [a*b*c*d for a,b,c,d in zip(DBasin_concC,DBasin_outflowC,conv_cfs_cms,conv_mgs_kgs)]
Wetland_loadC = [a*b*c*d for a,b,c,d in zip(Wetland_concC,Wetland_outflowC,conv_cfs_cms, conv_mgs_kgs)]
Channel_loadC = [a*b*c*d for a,b,c,d in zip(Channel_concC,Channel_flowC,conv_cfs_cms,conv_mgs_kgs)]

# Calculate cumulative load (dt = 1)
Ellsworth_cumloadC = np.cumsum(Ellsworth_loadC)
DBasin_cumloadC = np.cumsum(DBasin_loadC)
Wetland_cumloadC = np.cumsum(Wetland_loadC)
Channel_cumloadC = np.cumsum(Channel_loadC)

#----------------------------------------------------------------------#

# Print final load released
print("Ellsworth:", Ellsworth_cumload[-1])
print("Doyle Basin:", DBasin_cumload[-1])
print("Wetland:", Wetland_cumload[-1])
print("Channel to Outfall:", Channel_cumload[-1])
print("Ellsworth Controlled:", Ellsworth_cumloadC[-1])
print("Doyle Basin Controlled:", DBasin_cumloadC[-1])
print("Wetland Controlled:", Wetland_cumloadC[-1])
print("Channel Controlled:", Channel_cumloadC[-1])

#----------------------------------------------------------------------#

# Data for flooding line
Ell_flood = [6.0960]*len(Ellsworth_inflowC)
Wt_bypass = [2.8956]*len(Ellsworth_inflowC)
Wt_flood  = [2.7432]*len(Ellsworth_inflowC)

# Plot Result
fig, ax = plt.subplots(7, 4, sharex=True)
ax[0,0].plot(Ellsworth_inflow_m, 'k--')
ax[0,0].plot(Ellsworth_inflow_mC, 'b')
ax[0,0].set_ylabel("Inflow")
ax[0,0].set_title("Ellsworth")
ax[1,0].plot(Ellsworth_conc, 'k--')
ax[1,0].plot(Ellsworth_concC, 'b')
ax[1,0].set_ylabel("TSS (mg/L)")
ax[2,0].set_ylabel("DO (mg/L")
ax[3,0].plot(Ell_flood, 'r')
ax[3,0].plot(Ellsworth_depth_m, 'k--')
ax[3,0].plot(Ellsworth_depth_mC, 'b')
ax[3,0].set_ylabel("Depth (m)")
ax[4,0].plot(Ellsworth_valveC, 'b')
ax[4,0].set_ylabel("Valve Position")
ax[5,0].plot(Ellsworth_outflow_m, 'k--')
ax[5,0].plot(Ellsworth_outflow_mC, 'b')
ax[5,0].set_ylabel("Outflow (m³/s)")
ax[6,0].plot(Ellsworth_cumload, 'k--')
ax[6,0].plot(Ellsworth_cumloadC, 'b')
ax[6,0].set_ylabel("Cum. Load (kg)")
ax[6,0].set_xlabel("Time (min)")
ax[0,1].plot(DBasin_inflow_m, "k--")
ax[0,1].plot(DBasin_inflow_mC, "b")
ax[0,1].set_title("Basin")
ax[1,1].plot(DBasin_conc, "k--")
ax[1,1].plot(DBasin_concC, "b")
ax[3,1].plot(Wt_bypass, 'r')
ax[3,1].plot(DBasin_depth_m, "k--")
ax[3,1].plot(DBasin_depth_mC, "b")
ax[5,1].plot(DBasin_outflow_m, "k--")
ax[5,1].plot(DBasin_outflow_mC, "b")
ax[6,1].plot(DBasin_cumload, "k--")
ax[6,1].plot(DBasin_cumloadC, "b")
ax[6,1].set_xlabel("Time (min)")
ax[0,2].plot(Wetland_inflow_m, "k--")
ax[0,2].plot(Wetland_inflow_mC, "b")
ax[0,2].set_title("Wetland")
ax[1,2].plot(Wetland_conc, "k--")
ax[1,2].plot(Wetland_concC, "b")
ax[2,2].plot(Wetland_DO1, "m--")
ax[2,2].plot(Wetland_DO2, "g--")
ax[2,2].plot(Wetland_DO3, "c--")
ax[2,2].plot(Wetland_DO1C, "m")
ax[2,2].plot(Wetland_DO2C, "g")
ax[2,2].plot(Wetland_DO3C, "c")
ax[3,2].plot(Wt_flood, 'r')
ax[3,2].plot(Wetland_depth_m, "k--")
ax[3,2].plot(Wetland_depth_mC, "b")
ax[4,2].plot(Wetland_valveC, "b")
ax[5,2].plot(Wetland_outflow_m, "k--")
ax[5,2].plot(Wetland_outflow_mC, "b")
ax[6,2].plot(Wetland_cumload, "k--")
ax[6,2].plot(Wetland_cumloadC, "b")
ax[6,2].set_xlabel("Time (min)")
ax[0,3].plot(Channel_flow_m, "k--")
ax[0,3].plot(Channel_flow_mC, "b")
ax[0,3].set_title("Channel")
ax[1,3].plot(Channel_conc, "k--")
ax[1,3].plot(Channel_concC, "b")
ax[3,3].plot(Channel_depth_m, "k--")
ax[3,3].plot(Channel_depth_mC, "b")
ax[5,3].plot(Channel_flow_m, "k--")
ax[5,3].plot(Channel_flow_mC, "b")
ax[6,3].plot(Channel_cumload, "k--")
ax[6,3].plot(Channel_cumloadC, "b")
ax[6,3].set_xlabel("Time (min)")
plt.show()