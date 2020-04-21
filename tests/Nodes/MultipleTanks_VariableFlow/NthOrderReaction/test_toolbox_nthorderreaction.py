# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-04-21 15:50:54

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

# Call before simulation
def setup_time():
    # Get simulation starting time
    start_time = sim.start_time
    # Get simulation ending time
    end_time = sim.end_time
    # Initial Value for Last Step
    last_timestep = start_time
    return last_timestep

# Call during simulation
def NthOrderReaction(Node_dict):
    # Get current time
    current_step = sim.current_time
    # Calculate model dt in seconds
    dt = (current_step - last_timestep).total_seconds()
    # Updating reference step
    last_timestep = current_step

    for node in Node_dict:
        for pollutant in Node_dict[node]:
            # Get current concentration
            C = sim._model.getNodePollutant(node, pollutant)
            # Calculate treatment
            Cnew = C - (Node_dict[node][pollutant][0]*(C**Node_dict[node][pollutant][0])*dt)
            # Set concentration each time step
            sim._model.setNodePollutant(node, pollutant, Cnew)

dict1 = {'2': {0: [0.01, 0.5]}, '5': {0: [0.01, 0.5]}}

with Simulation("./gamma_notreatment.inp") as sim:
    Tank2 = Nodes(sim)["2"]
    Tank5 = Nodes(sim)["5"]
    setup_time()

    # Step through the simulation    
    for step in sim:
        # Run treatment each time step
        NthOrderReaction(dict1)
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
