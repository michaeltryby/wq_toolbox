from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt
import faulthandler; faulthandler.enable()

conc1 = []
flow1 = []

# Create simulation 
with Simulation("./tank_test.inp") as sim:
    Tank = Nodes(sim)["Tank"]

    # Step through the simulation    
    for index, step in enumerate(sim):
        #if index == 5:
        #    sim._model.setNodePollutant("Tank", 0, 20)
        conc = Tank.pollut_quality
        print(conc)
        print(sim._model.getNodePollutant("Tank", 0))
        sim._model.setNodePollutant("Tank", 0, 5)
        conc1.append(conc['P1'])
        flow1.append(sim._model.getNodeResult("Tank", 0))

conc2 = []
flow2 = []

# Create simulation 
with Simulation("./tank_test_percremoval.inp") as sim:
    Tank = Nodes(sim)["Tank"]

    # Step through the simulation    
    for step in sim:
        c0 = Tank.pollut_quality
        print(c0)
        conc2.append(c0['P1'])
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
