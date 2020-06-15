# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-06-12 08:12:24

# Import required modules
from pyswmm import Simulation, Nodes, Links
from wq_toolbox.links import Link_Quality
from wq_toolbox.nodes import Node_Quality
import numpy as np
import matplotlib.pyplot as plt

# Make dictionaries for each water quality method
# Gravity Settling in Basins
# v_s from Troutman et al. 2020, C_s from Mallets Creek data
dict1 = {'93-50408': {0: [15.08, 21.0]},'93-50404': {0: [15.08, 21.0]}, \
         '93-49759': {0: [15.08, 21.0]}}

# Channel Erosion & Gravity Settling
# width parameter changed to match Mallet's Creek load
# Slope from SWMM model
# SS and d50 found to match v_s, v_s and C_s same as bove

dict2 = {'95-69044': {0: [0.0000005, 0.0807, 2.68, 0.04, 15.08, 21.0]}, \
         '95-70180': {0: [0.0000005, 0.0411, 2.68, 0.04, 15.08, 21.0]}, \
         '95-51776': {0: [0.0000005, 0.0582, 2.68, 0.04, 15.08, 21.0]}, \
         '95-51774': {0: [0.0000005, 0.1333, 2.68, 0.04, 15.08, 21.0]}, \
         '95-51634': {0: [0.0000005, 0.1695, 2.68, 0.04, 15.08, 21.0]}, \
         '95-69048': {0: [0.0000005, 0.0368, 2.68, 0.04, 15.08, 21.0]}, \
         '95-70594': {0: [0.0000005, 1.7213, 2.68, 0.04, 15.08, 21.0]}, \
         '95-51760': {0: [0.0000005, 3.7373, 2.68, 0.04, 15.08, 21.0]}, \
         '95-69050': {0: [0.0000005, 0.1379, 2.68, 0.04, 15.08, 21.0]}, \
         '95-51758': {0: [0.0000005, 0.3719, 2.68, 0.04, 15.08, 21.0]}, \
         '95-51757': {0: [0.0000005, 0.0497, 2.68, 0.04, 15.08, 21.0]}, \
         '95-70713': {0: [0.0000005, 1.8003, 2.68, 0.04, 15.08, 21.0]}, \
         '95-70277': {0: [0.0000005, 0.3756, 2.68, 0.04, 15.08, 21.0]}}

#----------------------------------------------------------------------#
# Uncontrolled Simulation
# Lists to store results
Ellsworth_inflow = []
Ellsworth_conc = []
Ellsworth_cuminload = []
Ellsworth_depth = []
Ellsworth_flooding = []
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
Wetland_flooding = []
Wetland_outflow = []
Wetland_cumload = []
Wtlnd_bp_inflows = []

Channel_flow = []
Channel_conc = []
Cannel_cuminload = []
Channel_depth = []
Channel_cumload = []

# Setup toolbox simulation
with Simulation("./modifiedMBDoyle_TSS_V2.inp") as sim:
    # Setup toolbox methods
    GS = Node_Quality(sim, dict1)
    ER_GS = Link_Quality(sim, dict2)
    
    # Get asset information
    Ellsworth = Nodes(sim)["93-50408"]
    DBasin = Nodes(sim)["93-50404"]
    Wetland = Nodes(sim)["93-49759"]
    Wtlnd_bypass = Links(sim)["95-70294"]
    Channel = Links(sim)["95-70277"]

    # Step through the simulation    
    for index,step in enumerate(sim):

        # Calculate gravity settling in basins
        GS.GravitySettling()
        # Calculate erosion andn gravity settling in channels
        ER_GS.Erosion_and_Settling()

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
        Ell_fl = Ellsworth.flooding
        Ellsworth_flooding.append(Ell_fl)
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
        Wt_fl = Wetland.flooding
        Wetland_flooding.append(Wt_fl)
        Wt_of = Wetland.total_outflow
        Wetland_outflow.append(DB_of)
        Wt_bp = Wtlnd_bypass.flow
        Wtlnd_bp_inflows.append(Wt_bp)
        Ch_f = Channel.flow
        Channel_flow.append(Ch_f)
        Ch_d = Channel.depth
        Channel_depth.append(Ch_d)

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
conv_cfs_cms = [0.02832]*len(Ellsworth_inflow)
Ellsworth_flooding_m = [a*b for a,b in zip(Ellsworth_flooding,conv_cfs_cms)]
Wetland_flooding_m = [a*b for a,b in zip(Wetland_flooding,conv_cfs_cms)]
Wtlnd_bypass_m = [a*b for a,b in zip(Wtlnd_bp_inflows,conv_cfs_cms)]

# Convert depth from ft to m
conv_ft_m = [0.3048]*len(Ellsworth_inflow)
Ellsworth_depth_m = [a*b for a,b in zip(Ellsworth_depth,conv_ft_m)]
DBasin_depth_m = [a*b for a,b in zip(DBasin_depth,conv_ft_m)]
Wetland_depth_m = [a*b for a,b in zip(Wetland_depth,conv_ft_m)]
Channel_depth_m = [a*b for a,b in zip(Channel_depth,conv_ft_m)]

# Convert volume from ft3 to m3
#conv_cfs_cms = [0.02832]*len(Ellsworth_inflow)
#Ellsworth_volume_m = [a*b for a,b in zip(Ellsworth_volume,conv_ft_m)]
#DBasin_volume_m = [a*b for a,b in zip(DBasin_volume,conv_ft_m)]
#Wetland_volume_m = [a*b for a,b in zip(Wetland_volume,conv_ft_m)]
#Channel_volume_m = [a*b for a,b in zip(Channel_volume,conv_ft_m)]

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
timestep = [5]*len(Ellsworth_inflow)
Ellsworth_load = [a*b*c*d*e for a,b,c,d,e in zip(Ellsworth_conc,Ellsworth_outflow,conv_cfs_cms, conv_mgs_kgs,timestep)]
DBasin_load = [a*b*c*d*e for a,b,c,d,e in zip(DBasin_conc,DBasin_outflow,conv_cfs_cms,conv_mgs_kgs,timestep)]
Wetland_load = [a*b*c*d*e for a,b,c,d,e in zip(Wetland_conc,Wetland_outflow,conv_cfs_cms, conv_mgs_kgs,timestep)]
Channel_load = [a*b*c*d*e for a,b,c,d,e in zip(Channel_conc,Channel_flow,conv_cfs_cms,conv_mgs_kgs,timestep)]

# Calculate cumulative load
Ellsworth_cumload = np.cumsum(Ellsworth_load)
DBasin_cumload = np.cumsum(DBasin_load)
Wetland_cumload = np.cumsum(Wetland_load)
Channel_cumload = np.cumsum(Channel_load)

# Calculate cumulative volume
#Ellsworth_cumvol = np.cumsum(Ellsworth_volume_m)
#DBasin_cumvol = np.cumsum(DBasin_volume_m)
#Wetland_cumvol = np.cumsum(Wetland_volume_m)
#Channel_cumvol = np.cumsum(Channel_volume_m)

#----------------------------------------------------------------------#
# Controlled Simulation 
# Lists to store results
Ellsworth_inflowC = []
Ellsworth_concC = []
Ellsworth_cuminloadC = []
Ellsworth_depthC = []
Ellsworth_floodingC = []
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
Wetland_floodingC = []
Wetland_valveC = []
Wetland_outflowC = []
Wetland_cumloadC = []
Wtlnd_bp_inflowsC = []

Channel_flowC = []
Channel_concC = []
Cannel_cuminloadC = []
Channel_depthC = []
Channel_cumloadC = []

# Setup toolbox simulation
with Simulation("./modifiedMBDoyle_TSS_V2.inp") as sim:
    # Setup toolbox methods
    GS = Node_Quality(sim, dict1)
    ER_GS = Link_Quality(sim, dict2)
    
    # Get asset information
    Ellsworth = Nodes(sim)["93-50408"]
    Ells_valve = Links(sim)["95-70951"]
    DBasin = Nodes(sim)["93-50404"]
    Channel = Links(sim)["95-70277"]
    Wetland = Nodes(sim)["93-49759"]
    Wtlnd_valve = Links(sim)["95-70293"]
    Wtlnd_bypass = Links(sim)["95-70294"]

    # Tracking time for control actions every 15 minutes (5 sec time step)
    _tempcount = 180

    # Step through the simulation    
    for index,step in enumerate(sim):

        # Calculate gravity settling in basins
        GS.GravitySettling()
        # Calculate erosion andn gravity settling in channels
        ER_GS.Erosion_and_Settling()

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
        Ell_fl = Ellsworth.flooding
        Ellsworth_floodingC.append(Ell_fl)
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
        Wt_fl = Wetland.flooding
        Wetland_floodingC.append(Wt_fl)
        Wt_of = Wetland.total_outflow
        Wetland_outflowC.append(DB_of)
        Wt_bp = Wtlnd_bypass.flow
        Wtlnd_bp_inflowsC.append(Wt_bp)
        Ch_f = Channel.flow
        Channel_flowC.append(Ch_f)
        Ch_d = Channel.depth
        Channel_depthC.append(Ch_d)

        # Wetland Control Actions (every 15 mins - 5 sec timesteps)
        if _tempcount == 180:
            # If wetland has capacity, slowly open both valves
            if Wt_d <= 9.5:
                Ells_valve.target_setting = min(0.05, 70.6/(np.sqrt(2.0*32.2*Ell_d))/25)
                Wtlnd_valve.target_setting = min(0.33, 70.6/(np.sqrt(2.0*32.2*Wt_d))/12.6)
                print("Ell and Wt cap, slow open both")
            # If wetland does not have capacity
            elif Wt_d > 9.5:
                # But Ellsworth does, close Ellsworth valve and slowly open Wetland valve
                if Ell_d <= 15:
                    Ells_valve.target_setting = 0.0
                    Wtlnd_valve.target_setting = min(0.33, 70.6/(np.sqrt(2.0*32.2*Wt_d))/12.6)
                    print("Ell cap, Wt no cap, close Ell, slow open Wt")
                # If neither has capacity, slowly open both valves
                else:
                    Ells_valve.target_setting = min(0.05, 70.6/(np.sqrt(2.0*32.2*Ell_d))/25)
                    Wtlnd_valve.target_setting = min(0.33, 70.6/(np.sqrt(2.0*32.2*Wt_d))/12.6)    
                    print("Both no cap")
            _tempcount= 0
        _tempcount+= 1

        Wetland_valveC.append(Wtlnd_valve.target_setting)
        Ellsworth_valveC.append(Ells_valve.target_setting)

        # Wetland Control Actions
        #Wtlnd_valve.target_setting = 0.05
        #Wetland_valveC.append(Wtlnd_valve.target_setting)

        # Ellsworth Control Actions
        #Ells_valve.target_setting = 0.05
        #Ellsworth_valveC.append(Ells_valve.target_setting)
    
    sim._model.swmm_end()
    print(sim.runoff_error)
    print(sim.flow_routing_error)
    print(sim.quality_error)
        
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

# Convert flooding rate from cfs to m3/s
conv_cfs_cms = [0.02832]*len(Ellsworth_inflowC)
Ellsworth_flooding_mC = [a*b for a,b in zip(Ellsworth_floodingC,conv_cfs_cms)]
Wetland_flooding_mC = [a*b for a,b in zip(Wetland_floodingC,conv_cfs_cms)]
Wtlnd_bypass_mC = [a*b for a,b in zip(Wtlnd_bp_inflowsC,conv_cfs_cms)]

# Convert depth from ft to m
conv_ft_m = [0.3048]*len(Ellsworth_inflowC)
Ellsworth_depth_mC = [a*b for a,b in zip(Ellsworth_depthC,conv_ft_m)]
DBasin_depth_mC = [a*b for a,b in zip(DBasin_depthC,conv_ft_m)]
Wetland_depth_mC = [a*b for a,b in zip(Wetland_depthC,conv_ft_m)]
Channel_depth_mC = [a*b for a,b in zip(Channel_depthC,conv_ft_m)]

# Convert volume from ft3 to m3
#conv_cfs_cms = [0.02832]*len(Ellsworth_inflow)
#Ellsworth_volume_mC = [a*b for a,b in zip(Ellsworth_volumeC,conv_ft_m)]
#DBasin_volume_mC = [a*b for a,b in zip(DBasin_volumeC,conv_ft_m)]
#Wetland_volume_mC = [a*b for a,b in zip(Wetland_volumeC,conv_ft_m)]
#Channel_volume_mC = [a*b for a,b in zip(Channel_volumeC,conv_ft_m)]

# Calculate inflow load each timestep
#conv_mgs_kgs = [0.000001]*len(Ellsworth_inflowC)
#Ellsworth_inloadC = [a*b*c*d for a,b,c,d in zip(Ellsworth_concC,Ellsworth_inflowC,conv_cfs_cms, conv_mgs_kgs)]
#DBasin_inloadC = [a*b*c*d for a,b,c,d in zip(DBasin_concC,DBasin_inflowC,conv_cfs_cms,conv_mgs_kgs)]
#Wetland_inloadC = [a*b*c*d for a,b,c,d in zip(Wetland_concC,Wetland_inflowC,conv_cfs_cms, conv_mgs_kgs)]
#Channel_inloadC = [a*b*c*d for a,b,c,d in zip(Channel_concC,Channel_flowC,conv_cfs_cms,conv_mgs_kgs)]

# Calculate cumulative load
#Ellsworth_cuminloadC = np.cumsum(Ellsworth_inloadC)
#DBasin_cuminloadC = np.cumsum(DBasin_inloadC)
#Wetland_cuminloadC = np.cumsum(Wetland_inloadC)
#Channel_cuminloadC = np.cumsum(Channel_inloadC)

# Calculate outflow load each timestep
conv_mgs_kgs = [0.000001]*len(Ellsworth_inflowC)
timestep = [5]*len(Ellsworth_inflowC)
Ellsworth_loadC = [a*b*c*d*e for a,b,c,d,e in zip(Ellsworth_concC,Ellsworth_outflowC,conv_cfs_cms, conv_mgs_kgs,timestep)]
DBasin_loadC = [a*b*c*d*e for a,b,c,d,e in zip(DBasin_concC,DBasin_outflowC,conv_cfs_cms,conv_mgs_kgs,timestep)]
Wetland_loadC = [a*b*c*d*e for a,b,c,d,e in zip(Wetland_concC,Wetland_outflowC,conv_cfs_cms, conv_mgs_kgs,timestep)]
Channel_loadC = [a*b*c*d*e for a,b,c,d,e in zip(Channel_concC,Channel_flowC,conv_cfs_cms,conv_mgs_kgs,timestep)]

# Calculate cumulative load (dt = 1)
Ellsworth_cumloadC = np.cumsum(Ellsworth_loadC)
DBasin_cumloadC = np.cumsum(DBasin_loadC)
Wetland_cumloadC = np.cumsum(Wetland_loadC)
Channel_cumloadC = np.cumsum(Channel_loadC)

# Calculate cumulative volume (dt = 1)
#Ellsworth_cumvolC = np.cumsum(Ellsworth_volume_mC)
#DBasin_cumvolC = np.cumsum(DBasin_volume_mC)
#Wetland_cumvolC = np.cumsum(Wetland_volume_mC)
#Channel_cumvolC = np.cumsum(Channel_volume_mC)

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
ax[2,0].plot(Ell_flood, 'r:')
ax[2,0].plot(Ellsworth_depth_m, 'k--')
ax[2,0].plot(Ellsworth_depth_mC, 'b')
ax[2,0].set_ylabel("Depth (m)")
ax[3,0].plot(Ellsworth_valveC, 'b')
ax[3,0].set_ylabel("Valve Position")
ax[4,0].plot(Ellsworth_outflow_m, 'k--')
ax[4,0].plot(Ellsworth_outflow_mC, 'b')
ax[4,0].set_ylabel("Outflow (m³/s)")
ax[5,0].plot(Ellsworth_cumload, 'k--')
ax[5,0].plot(Ellsworth_cumloadC, 'b')
ax[5,0].set_ylabel("Cum. Load (kg)")
ax[6,0].set_xlabel("Time ")
ax[6,0].plot(Ellsworth_flooding_m, 'k--')
ax[6,0].plot(Ellsworth_flooding_mC, 'b')
ax[6,0].set_ylabel("Flooding (m³/s)")
ax[0,1].plot(DBasin_inflow_m, "k--")
ax[0,1].plot(DBasin_inflow_mC, "b")
ax[0,1].set_title("Basin")
ax[1,1].plot(DBasin_conc, "k--")
ax[1,1].plot(DBasin_concC, "b")
ax[2,1].plot(Wt_bypass, 'r:')
ax[2,1].plot(DBasin_depth_m, "k--")
ax[2,1].plot(DBasin_depth_mC, "b")
ax[4,1].plot(DBasin_outflow_m, "k--")
ax[4,1].plot(DBasin_outflow_mC, "b")
ax[5,1].plot(DBasin_cumload, "k--")
ax[5,1].plot(DBasin_cumloadC, "b")
ax[6,1].plot(Wtlnd_bypass_m, 'k--')
ax[6,1].plot(Wtlnd_bypass_mC, 'b')
ax[6,1].set_xlabel("Time ")
ax[0,2].plot(Wetland_inflow_m, "k--")
ax[0,2].plot(Wetland_inflow_mC, "b")
ax[0,2].set_title("Wetland")
ax[1,2].plot(Wetland_conc, "k--")
ax[1,2].plot(Wetland_concC, "b")
ax[2,2].plot(Wt_flood, 'r:')
ax[2,2].plot(Wetland_depth_m, "k--")
ax[2,2].plot(Wetland_depth_mC, "b")
ax[3,2].plot(Wetland_valveC, "b")
ax[4,2].plot(Wetland_outflow_m, "k--")
ax[4,2].plot(Wetland_outflow_mC, "b")
ax[5,2].plot(Wetland_cumload, "k--")
ax[5,2].plot(Wetland_cumloadC, "b")
ax[6,2].plot(Wetland_flooding_m, 'k--')
ax[6,2].plot(Wetland_flooding_mC, 'b')
ax[6,2].set_xlabel("Time ")
ax[0,3].plot(Channel_flow_m, "k--")
ax[0,3].plot(Channel_flow_mC, "b")
ax[0,3].set_title("Channel")
ax[1,3].plot(Channel_conc, "k--")
ax[1,3].plot(Channel_concC, "b")
ax[2,3].plot(Channel_depth_m, "k--")
ax[2,3].plot(Channel_depth_mC, "b")
ax[4,3].plot(Channel_flow_m, "k--")
ax[4,3].plot(Channel_flow_mC, "b")
ax[5,3].plot(Channel_cumload, "k--")
ax[5,3].plot(Channel_cumloadC, "b")
ax[6,3].set_xlabel("Time ")
plt.show()