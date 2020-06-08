# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-06-04 16:30:20
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-06-04 16:30:28

# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-06-04 16:29:33

# Import required modules
from pyswmm import Simulation, Nodes, Links
from wq_toolbox.links import Link_Quality
from wq_toolbox.nodes import Node_Quality
import numpy as np
import matplotlib.pyplot as plt

# Make dictionaries for each water quality method
# C_s and v_s from SWMM Manual
dict1 = {'93-50408': {0: [0.01, 21.0]}, '93-50404': {0: [0.01, 21.0]}} 
# width and slope from SWMM model
# SS and d50 from "Urban Hydroinformatics" (Price & Vojinovi) pg 130 
dict2 = {'95-70180': {0: [40.0, 0.00041, 2.6, 1.0]}}

# Lists to store results
#Channel_inflow = []
#Channel_conc = []
#Channel_depth = []
#Channel_outflow = []
#Channel_cumload = []

Ellsworth_inflow = []
Ellsworth_conc = []
Ellsworth_depth = []
Ellsworth_valve = []
Ellsworth_outflow = []
Ellsworth_cumload = []

Doyle_Basin_inflow = []
Doyle_Basin_conc = []
Doyle_Basin_depth = []
Doyle_Basin_outflow = []
Doyle_Basin_cumload = []

Wetland_inflow = []
Wetland_conc = []
Wetland_depth = []
Wetland_valve = []
Wetland_outflow = []
Wetland_cumload = []


# Setup toolbox simulation
with Simulation("./modifiedMBDoyle_TSS.inp") as sim:
    # Setup toolbox methods
    GS = Node_Quality(sim, dict1)
    ER = Link_Quality(sim, dict2)
    # Get asset information
    Ellsworth = Nodes(sim)["93-50408"]
    Ells_valve = Links(sim)["95-70951"]
    Doyle_Basin = Nodes(sim)["93-50404"]
    Channel = Links(sim)["95-70180"]
    #Dwn_Channel = Nodes(sim)["97-50264"]
    Wetland = Nodes(sim)["93-49759"]
    Wtlnd_valve = Links(sim)["95-70293"]

    # Temp count for control actions every 30 minutes
    _tempcount = 30
    # Tracking time so that we open the valve after 36 hours
    _time = 0
    # Step through the simulation    
    for index,step in enumerate(sim):

        # Calculate gravity settling 
        GS.GravitySettling()
        # Calculate erosion produced
        ER.Erosion()

        # Get TSS conc for each asset        
        #Ch_p = Channel.pollut_quality['TSS']
        #Channel_conc.append(Ch_p)
        Ell_p = Ellsworth.pollut_quality['TSS']
        Ellsworth_conc.append(Ell_p)
        #DwnCh_p = Dwn_Channel.pollut_quality['TSS']
        #Dwn_Channel_conc.append(DwnCh_p)
        DB_p = Doyle_Basin.pollut_quality['TSS']
        Doyle_Basin_conc.append(DB_p)
        Wt_p = Wetland.pollut_quality['TSS']
        Wetland_conc.append(Wt_p)

        # Get flow for  each asset
        #Ch_f = Channel.flow
        #Channel_flow.append(Ch_f)
        #Ch_d = Channel.depth
        #Channel_depth.append(Ch_d)
        #Dwn_Channel_flow.append(Dwn_Channel.total_inflow)
        Ell_if = Ellsworth.total_inflow
        Ellsworth_inflow.append(Ell_if)
        Ell_d = Ellsworth.depth
        Ellsworth_depth.append(Ell_d)
        Ell_of = Ellsworth.total_outflow
        Ellsworth_outflow.append(Ell_of)
        DB_if = Doyle_Basin.total_inflow
        Doyle_Basin_inflow.append(DB_if)
        DB_d = Doyle_Basin.depth
        Doyle_Basin_depth.append(DB_d)
        DB_of = Doyle_Basin.total_outflow
        Doyle_Basin_outflow.append(DB_if)
        Wt_if = Wetland.total_inflow
        Wetland_inflow.append(Wt_if)
        Wt_d = Wetland.depth
        Wetland_depth.append(Wt_d)
        Wt_of = Wetland.total_outflow
        Wetland_outflow.append(DB_of)

        # Wetland Control Actions (every 30 mins)
        if _tempcount == 30:
            # If depth <= 6.75 ft then close valve
            if depth <= 6.75:
                Wtlnd_valve.target_setting = 0.0
            # Proportional release but no more than 33% open
            else:    
                Wtlnd_valve.target_setting = min(0.33, Wt_of/(1.00*np.sqrt(2.0*32.2*Wt_d)))
        Wetland_valve.append(Wtlnd_valve.target_setting)
        
        # Ellsworth Control Actions (every 30 mins)
        # Proportional release after 36 hours
        if _time <= 2160:
            if _tempcount == 30:
                # If concentration >= 20 or depth <= 15 close valve
                if Ell_p >= 20.0 or Ell_d <= 15.0:
                    Ells_valve.target_setting = 0.0  
                else:  
                    # Else proportional release but no more than 33% open
                    Ells_valve.target_setting = min(0.33, Ell_f/np.sqrt(2.0*32.2*Ell_d))
                _tempcount = 0
            _tempcount += 1
        else:
            # After 36 hours, proportional release but no more than 33% open
            Ells_valve.target_setting = min(0.33, Ell_f/np.sqrt(2.0*32.2*Ell_d))
        _time += 1
        Ellsworth_valve.append(Ells_valve.target_setting)

        
#----------------------------------------------------------------------#
# Confirm gravity settling matches toolbox simulation
"""
swmm_c = []
swmm_f = []

with Simulation("./MBDoyle_TSS_SWMM.inp") as sim:
    Ellsworth = Nodes(sim)["93-50408"]
    Doyle_Basin = Nodes(sim)["93-50404"]
    Wetland = Nodes(sim)["93-49759"]
    # Step through the simulation    
    for step in sim:
        # Get newQual for Tank
        co = Ellsworth.pollut_quality
        swmm_c.append(co['TSS'])
        # Get flow for Tank
        swmm_f.append(sim._model.getNodeResult("93-50408", 0))
"""
#----------------------------------------------------------------------#
# Convert flow rate from cfs to m3/s
conv_cfs_cms = [0.02832]*len(Channel_flow)
Ellsworth_flow_m = [a*b for a,b in zip(Ellsworth_flow,conv_cfs_cms)]
Channel_flow_m = [a*b for a,b in zip(Channel_flow,conv_cfs_cms)]
Doyle_Basin_flow_m = [a*b for a,b in zip(Doyle_Basin_flow,conv_cfs_cms)]

# Convert depth from ft to m
conv_ft_m = [0.3048]*len(Ellsworth_flow)
Wetland_depth_m = [a*b for a,b in zip(Ellsworth_depth,conv_ft_m)]

# Calculate load each timestep
conv_mgs_kgs = [0.000001]*len(Channel_flow)
load1 = [a*b*c*d for a,b,c,d in zip(Channel_conc,Channel_flow,conv_cfs_cms, conv_mgs_kgs)]
load2 = [a*b*c*d for a,b,c,d in zip(Ellsworth_conc,Ellsworth_flow,conv_cfs_cms, conv_mgs_kgs)]
#load3 = [a*b*c*d for a,b,c,d in zip(Dwn_Channel_conc,Dwn_Channel_flow,conv_cfs_cms, conv_mgs_kgs)]
#load4 = [a*b*c*d for a,b,c,d in zip(swmm_c,swmm_f,conv_cfs_cms,conv_mgs_kgs)]
load5 = [a*b*c*d for a,b,c,d in zip(Doyle_Basin_conc,Doyle_Basin_flow,conv_cfs_cms,conv_mgs_kgs)]

# Calculate cumulative load (dt = 1)
Channel_cumload = np.cumsum(load1)
cum_load2 = np.cumsum(load2)
#cum_load3 = np.cumsum(load3)
#cum_load4 = np.cumsum(load4)
cum_load5 = np.cumsum(load5)
print("Doyle Basin:",cum_load5[-1])
print("Ellsworth:", cum_load2[-1])

# Confirm erosion mass balance in Channel
# Calculate error
#error = (cum_load1[-1]/cum_load3[-1])/cum_load1[-1]
#print(error)

#----------------------------------------------------------------------#
# Plot Result
fig, ax = plt.subplots(6, 4, sharex=True)
ax[0,0].plot(Ellsworth_inflow)
ax[0,0].set_ylabel("Inflow")
ax[0,0].set_title("Ellsworth")
ax[1,0].plot(Ellsworth_conc)
ax[1,0].set_ylabel("TSS (mg/L)")
ax[2,0].plot(Ellsworth_inflow)
ax[1,0].plot(Ellsworth_outflow_m)
ax[1,0].set_ylabel("Outflow (m³/s)")
ax[2,0].plot(Ellsworth_depth)
ax[2,0].set_ylabel("Depth (m)")
ax[3,0].plot(cum_load2)
ax[3,0].set_ylabel("Cum. Load (kg)")
ax[3,0].set_xlabel("Time (min)")
ax[0,1].plot(Channel_conc)
ax[0,1].set_title("Channel")
ax[1,1].plot(Channel_flow_m)
ax[2,1].plot(Channel_depth)
ax[3,1].plot(Channel_cumload)
ax[3,1].set_xlabel("Time (min)")
ax[0,2].plot(Doyle_Basin_conc)
ax[0,2].set_title("MB Doyle Basin")
ax[1,2].plot(Doyle_Basin_flow_m)
ax[2,2].plot(Doyle_Basin_depth)
ax[3,2].plot(cum_load5)
ax[3,2].set_xlabel("Time (min)")
plt.show()
