# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-04-22 13:37:42

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

with Simulation("./gamma_notreatment.inp") as sim:
    Tank2 = Nodes(sim)["2"]
    Tank5 = Nodes(sim)["5"]
    # Get simulation starting time
    start_time = sim.start_time
    # Initial Value for Last Step
    last_timestep = start_time

    # Step through the simulation    
    for step in sim:
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - last_timestep).total_seconds()
        # Updating reference step
        last_timestep = current_step
        # Get current concentration
        C2 = sim._model.getNodePollutant("2", 0)
        C5 = sim._model.getNodePollutant("5", 0)
        # Calculate treatment
        Cnew2 = C2 - (0.01*(C2**0.50)*dt)
        Cnew5 = C5 - (0.01*(C5**0.50)*dt)
        # Set concentration each time step
        sim._model.setNodePollutant("2", 0, Cnew2)
        sim._model.setNodePollutant("5", 0, Cnew5)
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