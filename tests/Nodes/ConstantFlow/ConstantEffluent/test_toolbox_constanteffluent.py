from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt

# Single Tank, Constant Inflow, Constant Effluent

# PySWMM Toolbox
class Node_Treatment:
    
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 


    def EventMeanConc(self):
        """
        Event Mean Concentration Treatment (SWMM Water Quality Manual, 2016)
        Treatment results in a constant concentration.

        Dictionary format: 
        dict = {'SWMM_Node_ID1': {pindex1: C, pindex2: C},
                'SWMM_Node_ID2': {pindex1: C, pindex2: C}}
        
        C = constant treatment concentration for each pollutant (SI or US: mg/L)
        """
        # Read from user dictionary
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                # Set concentration
                sim._model.setNodePollutant(node, pollutant, self.node_dict[node][pollutant])

# Create simulation
conc1 = []
flow1 = [] 
dict1 = {'Tank': {0: 2}}

with Simulation("./tank_constantinflow_notreatment.inp") as sim:
    NT = NodeTreatment(sim, dict1)
    Tank = Nodes(sim)["Tank"]

    # Step through the simulation    
    for step in sim:
        # Set constant effluent concentration each time step
        NT.EventMeanConc()
        # Get newQual for Tank
        conc = Tank.pollut_quality
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
