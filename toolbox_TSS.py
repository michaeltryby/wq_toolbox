# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-28 13:00:56

# Import required modules
from pyswmm import Simulation, Nodes, Links
from wq_toolbox.links import Link_Quality
from wq_toolbox.nodes import Node_Quality
import numpy as np
import matplotlib.pyplot as plt

# Make dictionaries for each water quality method
# C_s amd v_s from gravity Settling example in SWMM manual
dict1 = {'93-50408': {0: [0.01, 20]}, '93-50404': {0: [0.01, 20]}} 
dict2 = {'95-69044': {0: [68.0, 0.0037, 1.28, 0.035]}}

# Lists to store results
Channel_conc = []
Channel_flow =[]
Ellsworth_conc = []
Ellsworth_flow = []
Dwn_Channel_conc = []
Dwn_Channel_flow = []
Doyle_Basin_conc = []
Doyle_Basin_flow = []

# Tank Inflows
Qin1 = np.genfromtxt('Qin1.txt', delimiter=',')
Qin2 = np.genfromtxt('Qin2.txt', delimiter=',')
Qin3 = np.genfromtxt('Qin3.txt', delimiter=',')

# Setup toolbox simulation
with Simulation("./Ellsworth_Doyle_TSS.inp") as sim:
    # Setup toolbox methods
    GS = Node_Quality(sim, dict1)
    ER = Link_Quality(sim, dict2)
    # Get asset information
    Ellsworth = Nodes(sim)["93-50408"]
    Doyle_Basin = Nodes(sim)["93-50404"]
    Channel = Links(sim)["95-69044"]
    Dwn_Channel = Nodes(sim)["97-50264"]
    Wetland = Nodes(sim)["93-49759"]

    # Step through the simulation    
    for index,step in enumerate(sim):
        
        # Set inflow
        sim._model.setNodeInflow("93-50408", Qin1[index])
        sim._model.setNodeInflow("93-50404", Qin2[index])
        sim._model.setNodeInflow("93-49759", Qin3[index])

        # Calculate gravity settling 
        GS.GravitySettling()
        # Calculate erosion produced
        ER.Erosion()

        # Get TSS conc for each asset        
        c1 = Channel.pollut_quality
        Channel_conc.append(c1['TSS'])
        c2 = Ellsworth.pollut_quality
        Ellsworth_conc.append(c2['TSS'])
        c3 = Dwn_Channel.pollut_quality
        Dwn_Channel_conc.append(c3['TSS'])
        c4 = Doyle_Basin.pollut_quality
        Doyle_Basin_conc.append(c4['TSS'])

        # Get flow for  each asset
        Channel_flow.append(sim._model.getLinkResult("95-69044", 0))
        Ellsworth_flow.append(sim._model.getNodeResult("93-50408", 0))
        Dwn_Channel_flow.append(sim._model.getNodeResult("97-50264", 0))
        Doyle_Basin_flow.append(sim._model.getNodeResult("93-50404", 0))
        
       
#----------------------------------------------------------------------#
# Confirm gravity settling matches toolbox simulation
swmm_c = []
swmm_f = []

with Simulation("./Ellsworth_Doyle_TSS_SWMM.inp") as sim:
    Ellsworth = Nodes(sim)["93-50408"]
    Doyle_Basin = Nodes(sim)["93-50404"]
    Wetland = Nodes(sim)["93-49759"]
    # Step through the simulation    
    for step in sim:
        # Set inflow
        sim._model.setNodeInflow("93-50408", Qin1[index])
        sim._model.setNodeInflow("93-50404", Qin2[index])
        sim._model.setNodeInflow("93-49759", Qin3[index])
        # Get newQual for Tank
        co = Ellsworth.pollut_quality
        swmm_c.append(co['TSS'])
        # Get flow for Tank
        swmm_f.append(sim._model.getNodeResult("93-50408", 0))

#----------------------------------------------------------------------#
# Calculate load each timestep
load1 = [a*b for a,b in zip(Channel_conc,Channel_flow)]
load2 = [a*b for a,b in zip(Ellsworth_conc,Ellsworth_flow)]
load3 = [a*b for a,b in zip(Dwn_Channel_conc,Dwn_Channel_flow)]
load4 = [a*b for a,b in zip(swmm_c,swmm_f)]
load5 = [a*b for a,b in zip(Doyle_Basin_conc,Doyle_Basin_flow)]

# Calculate cumulative load
cum_load1 = np.cumsum(load1)
cum_load2 = np.cumsum(load2)
cum_load3 = np.cumsum(load3)
cum_load4 = np.cumsum(load4)
cum_load5 = np.cumsum(load5)


# Confirm erosion mass balance in Channel
# Calculate error
#error = (cum_load1[-1]/cum_load3[-1])/cum_load1[-1]
#print(error)

#----------------------------------------------------------------------#
# Plot Results
plt.subplot(311)
plt.plot(Channel_conc, 'g--', label='Channel')
plt.plot(Ellsworth_conc, 'b--', label='Ellsworth')
plt.plot(Doyle_Basin_conc, 'm--', label='Doyle Basin')
#plt.plot(swmm_c, 'r--', label="SWMM_Ellsworth")
plt.ylabel("TSS Conc")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(312)
plt.plot(Channel_flow, 'g--', label='Channel')
plt.plot(Ellsworth_flow, 'b--', label='Ellsworth')
plt.plot(Doyle_Basin_flow, 'm--', label='Doyle Basin')
#plt.plot(swmm_f, 'r--', label="SWMM_Ellsworth")
plt.ylabel("Flow")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(313)
plt.plot(cum_load1, 'g--', label='Channel')
plt.plot(cum_load2, 'b--', label='Ellsworth')
plt.plot(cum_load5, 'm--', label='Doyle Basin')
#plt.plot(cum_load4, 'r--', label="SWMM_Ellsworth")
plt.ylabel("TSS Cumm Load")
plt.xlabel("Time (s)")
plt.legend()
plt.show()
