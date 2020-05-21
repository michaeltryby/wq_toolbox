# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-19 11:48:27

# Import required modules
from pyswmm import Simulation, Nodes, Links
from wq_toolbox.links import Link_Quality
from wq_toolbox.nodes import Node_Quality
import numpy as np
import matplotlib.pyplot as plt

# Make dictionaries for each water quality method
# C_s from gravity Settling example in SWMM manual, v_s from Troutman et al. 2020
dict1 = {'Tank1': {0: [0.01, 20]}, 'Tank2': {0: [0.01, 20]}}
# From ?? 
dict2 = {'Channel': {0: [10.0, 0.0001, 2.68, 0.68]}}


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
with Simulation("./tanks_channel_TSS.inp") as sim:
    # Setup toolbox methods
    GS = Node_Quality(sim, dict1)
    ER = Link_Quality(sim, dict2)
    # Get asset information
    Tank1 = Nodes(sim)["Tank1"]
    Tank2 = Nodes(sim)["Tank2"]
    Ch1 = Links(sim)["Channel"]

    # Step through the simulation    
    for index,step in enumerate(sim):
        # Calculate gravity settling 
        GS.GravitySettling()
        # Calculate erosion produced
        ER.Erosion()

        # Get TSS conc for each asset        
        c1 = Ch1.pollut_quality
        ch1.append(c1['TSS'])
        c2 = Tank1.pollut_quality
        tank1.append(c2['TSS'])
        c3 = Tank2.pollut_quality
        tank2.append(c3['TSS'])

        # Get flow for  each asset
        ch1_flow.append(sim._model.getLinkResult("Channel", 0))
        tank1_flow.append(sim._model.getNodeResult("Tank1", 0))
        tank2_flow.append(sim._model.getNodeResult("Tank2", 0))

       
#----------------------------------------------------------------------#
# Confirm gravity settling matches toolbox simulation
swmm_c = []
swmm_f = []

with Simulation("./tanks_channel_TSS_SWMM.inp") as sim:
    Tank1 = Nodes(sim)["Tank1"]
    # Step through the simulation    
    for step in sim:
        # Get newQual for Tank
        co = Tank1.pollut_quality
        swmm_c.append(co['TSS'])
        # Get flow for Tank
        swmm_f.append(sim._model.getNodeResult("Tank1", 0))

#----------------------------------------------------------------------#
# Confirm erosion mass balance
# Calculate load each timestep
load1 = [a*b for a,b in zip(ch1,ch1_flow)]
load2 = [a*b for a,b in zip(tank2,tank2_flow)]

# Calculate cumulative load
cum_load1 = np.cumsum(load1)
cum_load2 = np.cumsum(load2)

# Calculate error
error = (cum_load1[-1]/cum_load2[-1])/cum_load1[-1]
print(error)

#----------------------------------------------------------------------#
# Plot Results

# Calculate load
load3 = [a*b for a,b in zip(tank2,tank2_flow)]
cum_load3 = np.cumsum(load3)

load4 = [a*b for a,b in zip(swmm_c,swmm_f)]
cum_load4 = np.cumsum(load4)

plt.subplot(311)
plt.plot(ch1, 'g--', label='Channel')
plt.plot(tank1, 'b:', label='Tank1')
plt.plot(tank2, 'm:', label='Tank2')
plt.plot(swmm_c, 'r', label="SWMM_Tank1")
plt.ylabel("TSS Conc")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(312)
plt.plot(ch1_flow, 'g--', label='Channel')
plt.plot(tank1_flow, 'b:', label='Tank1')
plt.plot(tank2_flow, 'm:', label='Tank2')
plt.plot(swmm_f, 'r', label="SWMM_Tank1")
plt.ylabel("Flow")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(313)
plt.plot(cum_load1, 'g--', label='Channel')
plt.plot(cum_load2, 'b:', label='Tank1')
plt.plot(cum_load3, 'm:', label='Tank2')
plt.plot(cum_load4, 'r', label="SWMM_Tank1")
plt.ylabel("TSS Cumm Load")
plt.xlabel("Time (s)")
plt.legend()
plt.show()
