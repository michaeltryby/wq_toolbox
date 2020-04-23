# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-04-23 15:45:21

from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt

# Multiple Tanks, Variable Inflow, Nth Order Reaction

# PySWMM Toolbox
# Create simulation
conc2 = []
conc5 = []
flow2 = [] 
flow5 = []

class NthOrderReaction:
    def __init__(self, sim):
        self.sim = sim
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 

    def treatment(self, node_dict):
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step

        for node in node_dict:
            for pollutant in node_dict[node]:
                # Get current concentration
                C = sim._model.getNodePollutant(node, pollutant)
                # Calculate treatment
                Cnew = C - (node_dict[node][pollutant][0]*(C**node_dict[node][pollutant][1])*dt)
                # Set concentration each time step
                sim._model.setNodePollutant(node, pollutant, Cnew)

dict1 = {'2': {0: [0.01, 0.5]}, '5': {0: [0.01, 0.5]}}

with Simulation("./gamma_notreatment.inp") as sim:
    Tank2 = Nodes(sim)["2"]
    Tank5 = Nodes(sim)["5"]
    #last_timestep = NOR.time_tracking()
    NOR = NthOrderReaction(sim)

    # Step through the simulation    
    for step in sim:
        # Run treatment each time step
        NOR.treatment(dict1)
        # Get influent concentration
        c2 = Tank2.pollut_quality
        conc2.append(c2['P1'])
        c5 = Tank5.pollut_quality
        conc5.append(c5['P1'])
        # Get flow into Tank
        flow2.append(sim._model.getNodeResult("2", 0))
        flow5.append(sim._model.getNodeResult("5", 0))

# SWMM Comparision
# Create simulation
con2 = []
con5 = []
flo2 = [] 
flo5 = [] 

with Simulation("./gamma_nthorderreaction.inp") as sim:
    Tank2 = Nodes(sim)["2"]
    Tank5 = Nodes(sim)["5"]
    # Step through the simulation    
    for step in sim:
        # Get newQual for Tank
        co2 = Tank2.pollut_quality
        con2.append(co2['P1'])
        co5 = Tank5.pollut_quality
        con5.append(co5['P1'])
        # Get flow into Tank
        flo2.append(sim._model.getNodeResult("2", 0))
        flo5.append(sim._model.getNodeResult("5", 0))

# Plot Results
# Concentration Results
plt.subplot(211)
plt.plot(conc2, 'g--', label='tb-T2')
plt.plot(conc5, 'b--', label='tb-T5')
plt.plot(con2, 'm:', label='sw-T2')
plt.plot(con5, 'y:', label='sw-T5')
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
