# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-02 13:36:05

from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt

# Multiple Tanks, Variable Inflow, Sedimentation & Resuspension 
# PySWMM Toolbox

class Node_Treatment:
    
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 

    def SedimentationResuspension(self):
        """
        SEDIMENTATION & RESUSPENSION (Troutman et al. 2020)
        Model considers both settling, as a function of depth, and resuspension,
        as a function of flow.

        Dictionary format: 
        dict = {'SWMM_Node_ID1': {pindex1: [v_s, a, b], pindex2: [v_s, a, b]},
                'SWMM_Node_ID2': {pindex1: [v_s, a, b], pindex2: [v_s, a, b]}}
        
        v_s = settling velcity (SI: m/s, US: ft/s)
        a   = ratio between velcity and TSS resuspension to result in 100% 
              resuspension for the maximum velociy through storage pipe
        b   = linear approximation of the ratio bewteen flow and velocity
              computed for each upstream inline storage asset
        """
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step

        # Read from user dictionary
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                Qin = sim._model.getNodeResult(node,0)
                Cin = sim._model.getNodeCin(node,pollutant)
                C = sim._model.getNodeC2(node,pollutant)
                d = sim._model.getNodeResult(node,5)
                v_s = self.node_dict[node][pollutant][0]
                a = self.node_dict[node][pollutant][1]
                b = self.node_dict[node][pollutant][2]
                # Calculate removal
                if d != 0.0 and Qin != 0.0:
                    R = 1-np.exp(-v_s*dt/d)-np.exp(-a*b/Qin)
                else:
                    R = 0
                # Calculate new concentration
                Cnew = (1-R)*Cin
                #Cnew = min(Cnew,C)
                # Set new concentration
                sim._model.setNodePollutant(node, pollutant, Cnew)

# Create simulation
co_1 = []
co_2 = []
co_3 = []
co_4 = []

dict1 = {'004': {0: [0.00419, 2.18, 27.263]}, '006': {0: [0.00419, 7.29, 10.528]}, 
        '011': {0: [0.00419, 1.61, 17.815]}, '022': {0: [0.00419, 1.38, 23.852]}}

with Simulation("./epsilon_notreatment.inp") as sim:
    NT = Node_Treatment(sim, dict1)
    Tank4 = Nodes(sim)["004"]
    Tank6 = Nodes(sim)["006"]
    Tank11 = Nodes(sim)["011"]
    Tank22 = Nodes(sim)["022"]
    # Step through the simulation    
    for step in sim:
        # Run treatment each time step
        NT.SedimentationResuspension()
        # Get newQual for Tank
        c1 = Tank4.pollut_quality
        co_1.append(c1['TSS'])
        c2 = Tank6.pollut_quality
        co_2.append(c2['TSS'])
        c3 = Tank11.pollut_quality
        co_3.append(c3['TSS'])
        c4 = Tank22.pollut_quality
        co_4.append(c4['TSS'])

# SWMM Comparision
# Create simulation
con_1 = []
con_2 = []
con_3 = []
con_4 = [] 

with Simulation("./epsilon.inp") as sim:
    Tank4 = Nodes(sim)["004"]
    Tank6 = Nodes(sim)["006"]
    Tank11 = Nodes(sim)["011"]
    Tank22 = Nodes(sim)["022"]

    # Step through the simulation    
    for step in sim:
        # Get newQual for Tank
        co1 = Tank4.pollut_quality
        con_1.append(co1['TSS'])
        co2 = Tank6.pollut_quality
        con_2.append(co2['TSS'])
        co3 = Tank11.pollut_quality
        con_3.append(co3['TSS'])
        co4 = Tank22.pollut_quality
        con_4.append(co4['TSS'])

# Plot Results
# Concentration Results
plt.plot(co_1, 'g--', label='tb-T4')
plt.plot(co_2, 'm--', label='tb-T6')
plt.plot(co_3, 'b--', label='tb-T11')
plt.plot(co_4, 'c--', label='tb-T22')
plt.plot(con_1, 'g:', label='swmm-T4')
plt.plot(con_2, 'm:', label='swmm-T6')
plt.plot(con_3, 'b:', label='swmm-T11')
plt.plot(con_4, 'c:', label='swmm-T22')
plt.ylabel("Concentration")
plt.xlabel("Time (s)")
plt.legend()
plt.show()
