# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-02 11:58:54

from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt

# Multiple Tanks, Variable Inflow, Gravity Settling

# PySWMM Toolbox
# Create simulation
conc2_P1 = []
conc5_P1 = []
flow2 = [] 
flow5 = [] 

class GravitySettling:
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 

    def treatment(self):
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step
        # Read from user dictionary
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                Qin = sim._model.getNodeResult(node, 0)
                k = self.node_dict[node][pollutant][0]
                C_star = self.node_dict[node][pollutant][1]
                C = sim._model.getNodeC2(node,pollutant)
                depth = sim._model.getNodeResult(node, 5)
                if depth != 0:
                    # Calculate new concentration
                    Cnew = np.heaviside((0.1-Qin),0)*(C_star+(C-C_star)*np.exp(-k/depth*dt/3600))+(1-np.heaviside((0.1-Qin),0))*C
                else:
                    Cnew = np.heaviside((0.1-Qin),0)*C_star+(C-C_star)+(1-np.heaviside((0.1-Qin),0))*C
                # Set new concentration
                sim._model.setNodePollutant(node, pollutant, Cnew)

dict1 = {'2': {0: [0.01, 10]}, '5': {0: [0.10, 15]}}

with Simulation("./gamma_notreatment.inp") as sim:
    GS = GravitySettling(sim, dict1)
    Tank2 = Nodes(sim)["2"]
    Tank5 = Nodes(sim)["5"]
    # Step through the simulation    
    for step in sim:
        # Run treatment each time step
        GS.treatment()
        # Get newQual for Tank
        c2 = Tank2.pollut_quality
        conc2_P1.append(c2['P1'])
        c5 = Tank5.pollut_quality
        conc5_P1.append(c5['P1'])
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

with Simulation("./gamma_gravsettling.inp") as sim:
    Tank2 = Nodes(sim)["2"]
    Tank5 = Nodes(sim)["5"]
    # Step through the simulation    
    for step in sim:
        # Get newQual for Tank
        co2 = Tank2.pollut_quality
        con2_P1.append(co2['P1'])
        co5 = Tank5.pollut_quality
        con5_P1.append(co5['P1'])
        # Get flow for Tank
        flo2.append(sim._model.getNodeResult("2", 0))
        flo5.append(sim._model.getNodeResult("5", 0))

# Plot Results
# Concentration Results
plt.subplot(211)
plt.plot(conc2_P1, 'g--', label='tb-P1-T2')
plt.plot(conc5_P1, 'b--', label='tb-P1-T5')
plt.plot(con2_P1, 'g:', label='sw-P1-T2')
plt.plot(con5_P1, 'b:', label='sw-P1-T5')
plt.ylabel("Concentration")
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

