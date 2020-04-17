from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt
import faulthandler; faulthandler.enable()

# Single Tank, Constant Inflow, Constant Effluent

# PySWMM Toolbox
# Create simulation
conc1 = []
flow1 = [] 

with Simulation("./tank_constantinflow_notreatment.inp") as sim:
    Tank = Nodes(sim)["Tank"]
    # Set initial concentration
    sim._model.setNodePollutant("Tank", 0, 2)

    # Step through the simulation    
    for step in sim:
        # Get newQual for Tank
        conc = Tank.pollut_quality
        # Set constant effluent concentration each time step
        sim._model.setNodePollutant("Tank", 0, 2)
        conc1.append(conc['P1'])
        # Get flow for Tank
        flow1.append(sim._model.getNodeResult("Tank", 0))

# SWMM Comparision
# Create simulation
conc2 = []
flow2 = [] 

with Simulation("./tank_constantinflow_constanteffluent.inp") as sim:
    Tank = Nodes(sim)["Tank"]

    # Step through the simulation    
    for step in sim:
        # Get newQual for Tank
        c0 = Tank.pollut_quality
        print(c0)
        conc2.append(c0['P1'])
        # Get flow for Tank
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
