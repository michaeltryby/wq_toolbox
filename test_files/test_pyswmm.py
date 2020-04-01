from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt
import faulthandler; faulthandler.enable()

conc1 = []
flow1 = []

# Create simulation 
with Simulation("./tank_test.inp") as sim:
    # Turn on externalTreatment flag
    sim._model.setNodeParam("Tank", 5, 1)
    Tank = Nodes(sim)["Tank"]

    # Step through the simulation    
    for step in sim:
        conc = Tank.pollut_quality
        print(conc)
        #c1 = conc['P1']*0.5
        sim._model.setNodePollutant("Tank", 0, 2)
        conc1.append(conc['P1'])
        flow1.append(sim._model.getNodeResult("Tank", 0))


conc2 = []
flow2 = []

# Create simulation 
with Simulation("./tank_test_percremoval.inp") as sim:
    
    # Step through the simulation    
    for step in sim:
        conc = Tank.pollut_quality
        print(conc)
        conc2.append(conc['P1'])
        flow2.append(sim._model.getNodeResult("Tank", 0))


#Plot results
plt.subplot(211)
plt.plot(conc1, 'r--', label='toolbox')
plt.plot(conc2, 'b:', label='swmm')
plt.ylabel("Concentration")
plt.xlabel("Time (s)")
plt.legend()

plt.subplot(212)
plt.plot(flow1, 'r--', label='toolbox')
plt.plot(flow2, 'b:', label='swmm')
plt.xlabel("Time (s)")
plt.ylabel("Flow")
plt.legend()
plt.show()
