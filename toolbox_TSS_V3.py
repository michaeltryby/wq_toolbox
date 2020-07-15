# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-07-09 11:12:22

# Import required modules
from pyswmm import Simulation, Nodes, Links
from wq_toolbox.links import Link_Quality
from wq_toolbox.nodes import Node_Quality
import numpy as np
import matplotlib.pyplot as plt
import statistics as st

# Make dictionaries for each water quality method
# Gravity Settling in Basins
# v_s based on Mallet Creek data (0.0015 m/s = 5.27 m/hr), C_s from Mallets Creek data
dict1 = {'93-50408': {0: [5.27, 21.0]},'93-50404': {0: [5.27, 21.0]}, \
         '93-49759': {0: [5.27, 21.0]}}

# Channel Erosion & Gravity Settling
# width parameter changed to match Mallet's Creek load
# Slope from SWMM model
# SS and d50 found to match v_s, v_s and C_s same as bove

dict2 = {'95-69044': {0: [0.0000005, 0.0807, 2.68, 0.04, 5.27, 21.0]}, \
         '95-70180': {0: [0.0000005, 0.0411, 2.68, 0.04, 5.27, 21.0]}, \
         '95-51776': {0: [0.0000005, 0.0582, 2.68, 0.04, 5.27, 21.0]}, \
         '95-51774': {0: [0.0000005, 0.1333, 2.68, 0.04, 5.27, 21.0]}, \
         '95-51634': {0: [0.0000005, 0.1695, 2.68, 0.04, 5.27, 21.0]}, \
         '95-69048': {0: [0.0000005, 0.0368, 2.68, 0.04, 5.27, 21.0]}, \
         '95-70594': {0: [0.0000005, 1.7213, 2.68, 0.04, 5.27, 21.0]}, \
         '95-51760': {0: [0.0000005, 3.7373, 2.68, 0.04, 5.27, 21.0]}, \
         '95-69050': {0: [0.0000005, 0.1379, 2.68, 0.04, 5.27, 21.0]}, \
         '95-51758': {0: [0.0000005, 0.3719, 2.68, 0.04, 5.27, 21.0]}, \
         '95-51757': {0: [0.0000005, 0.0497, 2.68, 0.04, 5.27, 21.0]}, \
         '95-70713': {0: [0.0000005, 1.8003, 2.68, 0.04, 5.27, 21.0]}, \
         '95-70277': {0: [0.0000005, 0.3756, 2.68, 0.04, 5.27, 21.0]}}

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

#----------------------------------------------------------------------#

# Print final load released
print("Ellsworth:", Ellsworth_cumload[-1])
print("Doyle Basin:", DBasin_cumload[-1])
print("Wetland:", Wetland_cumload[-1])
print("Channel to Outfall:", Channel_cumload[-1])

print("Ellsworth Avg", st.mean(Ellsworth_conc))
print("Doyle Basin Avg", st.mean(DBasin_conc))
print("Wetland Avg", st.mean(Wetland_conc))
print("Channel Avg", st.mean(Channel_conc))

#----------------------------------------------------------------------#
# Plot Result
fig, ax = plt.subplots(4, 4)
ax[0,0].plot(Ellsworth_inflow_m, color='#6CC6D1', linewidth=2)
ax[0,0].set_xticks([])
ax[0,0].set_yticks([0,3,6,9])
ax[0,0].set_yticklabels(["0","3","6","9"])
ax[0,0].set_ylim(0,10)
ax[0,0].set_xlim(0,86400)
ax[0,0].set_ylabel("Inflow (m³/s)")

ax[1,0].plot(Ellsworth_conc, color='#6CC6D1', linewidth=2)
ax[1,0].set_xticks([])
ax[1,0].set_yticks([0,50,100,150])
ax[1,0].set_yticklabels(["0","50","100","150"])
ax[1,0].set_ylim(0,180)
ax[1,0].set_xlim(0,86400)
ax[1,0].set_ylabel("TSS (mg/L)")

ax[2,0].plot(Ellsworth_depth_m, color='#6CC6D1', linewidth=2)
ax[2,0].set_xticks([])
ax[2,0].set_yticks([0,1,2,3])
ax[2,0].set_yticklabels(["0","1","2","3"])
ax[2,0].set_ylim(0,3)
ax[2,0].set_xlim(0,86400)
ax[2,0].set_ylabel("Depth (m)")

ax[3,0].plot(Ellsworth_outflow_m, color='#6CC6D1', linewidth=2)
ax[3,0].set_yticks([0,3,6,9])
ax[3,0].set_yticklabels(["0","3","6","9"])
ax[3,0].set_ylim(0,9)
ax[3,0].set_ylabel("Outflow (m³/s)")
ax[3,0].set_xlim(0,86400)
ax[3,0].set_xticks([0,17280,34560,51840,69120,86400])
ax[3,0].set_xticklabels(["0","1","2","3","4","5"])
ax[3,0].set_xlabel("Time (days)")

ax[0,1].plot(DBasin_inflow_m, color='#3B4D7A', linewidth=2)
ax[0,1].set_xticks([])
ax[0,1].set_yticks([0,3,6,9])
ax[0,1].set_yticklabels(["0","3","6","9"])
ax[0,1].set_ylim(0,10)
ax[0,1].set_xlim(0,86400)

ax[1,1].plot(DBasin_conc, color='#3B4D7A', linewidth=2)
ax[1,1].set_xticks([])
ax[1,1].set_yticks([0,50,100,150])
ax[1,1].set_ylim(0,180)
ax[1,1].set_xlim(0,86400)
ax[1,1].set_yticklabels(["0","50","100","150"])

ax[2,1].plot(DBasin_depth_m, color='#3B4D7A', linewidth=2)
ax[2,1].set_yticks([0,1,2,3])
ax[2,1].set_yticklabels(["0","1","2","3"])
ax[2,1].set_xticks([])
ax[2,1].set_ylim(0,3)
ax[2,1].set_xlim(0,86400)

ax[3,1].plot(DBasin_outflow_m, color='#3B4D7A', linewidth=2)
ax[3,1].set_yticks([0,3,6,9])
ax[3,1].set_yticklabels(["0","3","6","9"])
ax[3,1].set_ylim(0,9)
ax[3,1].set_xlim(0,86400)
ax[3,1].set_xticks([0,17280,34560,51840,69120,86400])
ax[3,1].set_xticklabels(["0","1","2","3","4","5"])
ax[3,1].set_xlabel("Time (days)")

ax[0,2].plot(Wetland_inflow_m, color='#B08CA1', linewidth=2)
ax[0,2].set_yticks([0,3,6,9])
ax[0,2].set_yticklabels(["0","3","6","9"])
ax[0,2].set_xticks([])
ax[0,2].set_ylim(0,10)
ax[0,2].set_xlim(0,86400)

ax[1,2].plot(Wetland_conc, color='#B08CA1', linewidth=2)
ax[1,2].set_yticks([0,50,100,150])
ax[1,2].set_yticklabels(["0","50","100","150"])
ax[1,2].set_xticks([])
ax[1,2].set_ylim(0,180)
ax[1,2].set_xlim(0,86400)

ax[2,2].plot(Wetland_depth_m, color='#B08CA1', linewidth=2)
ax[2,2].set_yticks([0,1,2,3])
ax[2,2].set_yticklabels(["0","1","2","3"])
ax[2,2].set_xticks([])
ax[2,2].set_ylim(0,3)
ax[2,2].set_xlim(0,86400)

ax[3,2].plot(Wetland_outflow_m, color='#B08CA1', linewidth=2)
ax[3,2].set_ylim(0,9)
ax[3,2].set_xlim(0,86400)
ax[3,2].set_yticks([0,3,6,9])
ax[3,2].set_yticklabels(["0","3","6","9"])
ax[3,2].set_xticks([0,17280,34560,51840,69120,86400])
ax[3,2].set_xticklabels(["0","1","2","3","4","5"])
ax[3,2].set_xlabel("Time (days)")

ax[0,3].plot(Channel_flow_m, color='#695580', linewidth=2)
ax[0,3].set_yticks([0,3,6,9])
ax[0,3].set_yticklabels(["0","3","6","9"])
ax[0,3].set_xticks([])
ax[0,3].set_ylim(0,10)
ax[0,3].set_xlim(0,86400)

ax[1,3].plot(Channel_conc, color='#695580', linewidth=2)
ax[1,3].set_yticks([0,50,100,150])
ax[1,3].set_yticklabels(["0","50","100","150"])
ax[1,3].set_ylim(0,180)
ax[1,3].set_xlim(0,86400)
ax[1,3].set_xticks([])

ax[2,3].plot(Channel_depth_m, color='#695580', linewidth=2)
ax[2,3].set_ylim(0,3)
ax[2,3].set_xlim(0,86400)
ax[2,3].set_yticks([0,1,2,3])
ax[2,3].set_yticklabels(["0","1","2","3"])
ax[2,3].set_xticks([])

ax[3,3].plot(Channel_flow_m, color='#695580', linewidth=2)
ax[3,3].set_ylim(0,9)
ax[3,3].set_xlim(0,86400)
ax[3,3].set_yticks([0,3,6,9])
ax[3,3].set_yticklabels(["0","3","6","9"])
ax[3,3].set_xticks([0,17280,34560,51840,69120,86400])
ax[3,3].set_xticklabels(["0","1","2","3","4","5"])
ax[3,3].set_xlabel("Time (days)")

fig.savefig('TSS.svg', format="svg")
plt.show()