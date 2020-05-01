# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-01 10:24:10

from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt

# Multiple Tanks, Variable Inflow, Co-Removal 

# PySWMM Toolbox
# Create simulation
conc2_P1 = []
conc2_P2 = []
conc5_P1 = []
conc5_P2 = []
flow2 = [] 
flow5 = [] 

class CoRemoval:
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict

    def treatment(self):
        # Read from user dictionary
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                #Get parameters
                Cin = sim._model.getNodeCin(node, pollutant)
                R1 = self.node_dict[node][pollutant][0]
                R2 = self.node_dict[node][pollutant][1]
                # Calculate new concentration
                Cnew = (1-R1*R2)*Cin
                # Set new concentration
                sim._model.setNodePollutant(node, pollutant, Cnew)

class ConstantRemoval:
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict

    def treatment(self):
        # Read from user dictionary
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                #Get parameters
                Cin = sim._model.getNodeCin(node, pollutant)
                R = self.node_dict[node][pollutant]
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration 
                sim._model.setNodePollutant(node, pollutant, Cnew)

dict1 = {'2': {0: [0.8, 0.3]}, '5': {0: [0.5, 0.3]}}
dict2 = {'2': {1: 0.3}, '5': {1: 0.3}}


with Simulation("./gamma_notreatment.inp") as sim:
    CoR = CoRemoval(sim, dict1)
    Co  = ConstantRemoval(sim, dict2)
    Tank2 = Nodes(sim)["2"]
    Tank5 = Nodes(sim)["5"]
    # Step through the simulation    
    for step in sim:
        # Run treatment each time step
        Co.treatment()
        CoR.treatment()
        # Get newQual for Tank
        c2 = Tank2.pollut_quality
        conc2_P1.append(c2['P1'])
        conc2_P2.append(c2['P2'])
        c5 = Tank5.pollut_quality
        conc5_P1.append(c5['P1'])
        conc5_P2.append(c5['P2'])
        # Get flow for Tank
        flow2.append(sim._model.getNodeResult("2", 0))
        flow5.append(sim._model.getNodeResult("5", 0))

# SWMM Comparision
# Create simulation
con2_P1 = []
con2_P2 = []
con5_P1 = []
con5_P2 = []
flo2 = [] 
flo5 = [] 

with Simulation("./gamma_coremoval.inp") as sim:
    Tank2 = Nodes(sim)["2"]
    Tank5 = Nodes(sim)["5"]
    # Step through the simulation    
    for step in sim:
        # Get newQual for Tank
        co2 = Tank2.pollut_quality
        con2_P1.append(co2['P1'])
        con2_P2.append(co2['P2'])
        co5 = Tank5.pollut_quality
        con5_P1.append(co5['P1'])
        con5_P2.append(co5['P2'])
        # Get flow for Tank
        flo2.append(sim._model.getNodeResult("2", 0))
        flo5.append(sim._model.getNodeResult("5", 0))

# Plot Results
# Concentration Results
plt.subplot(211)
plt.plot(conc2_P1, 'g--', label='tb-P1-T2')
plt.plot(conc2_P2, 'm--', label='tb-P2-T2')
plt.plot(conc5_P1, 'b--', label='tb-P1-T5')
plt.plot(conc5_P2, 'c--', label='tb-P2-T5')
plt.plot(con2_P1, 'g:', label='sw-P1-T2')
plt.plot(con2_P2, 'm:', label='sw-P2-T2')
plt.plot(con5_P1, 'b:', label='sw-P1-T5')
plt.plot(con5_P2, 'c:', label='sw-P2-T5')
plt.ylabel("Concentration")
plt.xlabel("Time (s)")
plt.legend()

# Flow Results
plt.subplot(212)
plt.plot(flow2, 'g--', label='tb-T2')
plt.plot(flow5, 'b--', label='tb-T5')
plt.plot(flo2, 'm:', label='sw-T2')
plt.plot(flo5, 'y:', label='sw-T5')
plt.xlabel("Time (s)")
plt.ylabel("Flow")
plt.legend()
plt.show()

