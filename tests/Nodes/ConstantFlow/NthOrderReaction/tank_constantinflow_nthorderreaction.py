from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt
import faulthandler; faulthandler.enable()

# Single Tank, Constant Inflow, Nth Order Reaction

# PySWMM Toolbox
# Create simulation
conc1 = []
flow1 = [] 

with Simulation("./tank_constantinflow_notreatment.inp") as sim:
    Tank = Nodes(sim)["Tank"]
    # Get simulation starting time
    start_time = sim.start_time
    # Get simulation ending time
    end_time = sim.end_time
    # Initial Value for Last Step
    last_timestep = start_time

    # Step through the simulation    
    for index, step in enumerate(sim):
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - last_timestep).total_seconds()
        # Updating reference step
        last_timestep = current_step
        # Get current concentration
        C = sim._model.getNodePollutant("Tank", 0)
        # Calculate treatment
        Cnew = C - (0.01*(C**2.0)*dt)
        # Set concentration each time step
        sim._model.setNodePollutant("Tank", 0, Cnew)
        # Get influent concentration
        conc = Tank.pollut_quality
        c1= conc['P1']
        conc1.append(c1)
        # Get flow into Tank
        flow1.append(sim._model.getNodeResult("Tank", 0))

# SWMM Comparision
# Create simulation
conc2 = []
flow2 = []

with Simulation("./tank_constantinflow_nthorderreaction.inp") as sim:
    Tank = Nodes(sim)["Tank"]

    # Step through the simulation    
    for step in sim:
        # Get newQual for Tank
        c0 = Tank.pollut_quality
        conc2.append(c0['P1'])
        # Get flow into Tank
        flow2.append(sim._model.getNodeResult("Tank", 0))

# Plot Results
# Concentration Results
plt.subplot(211)
plt.plot(conc1, 'r--', label='toolbox')
plt.plot(conc2, 'b:', label='swmm')
plt.ylabel("Concentration")
plt.xlabel("Time (s)")
plt.legend()

# Flow Results
plt.subplot(212)
plt.plot(flow1, 'r--', label='toolbox')
plt.plot(flow2, 'b:', label='swmm')
plt.xlabel("Time (s)")
plt.ylabel("Flow")
plt.legend()
plt.show()
