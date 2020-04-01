from pyswmm import Simulation
import numpy as np
import matplotlib.pyplot as plt
import faulthandler; faulthandler.enable()


conc1 = []
flow1 = []

# Create simulation 
with Simulation("./tank_test.inp") as sim:
    
    # Turn on externalTreatment flag
    #sim._model.setNodeParam("Tank", 5, 1)
    #c = sim._model.getNodeParam("Tank", 5)
    #print(c)
    
    # Step through the simulation    
    for step in sim:
        # run percent removal treatment
        d = sim._model.getNodePollutant("Tank", 0)
        
        # set new concentration in node
        #e = d * 0.5
        #print(e)
        #sim._model.setNodePollutant("Tank", 0, e)
        
        conc1.append(sim._model.getNodePollutant("Tank", 0))
        flow1.append(sim._model.getNodeResult("Tank", 0))


conc2 = []
flow2 = []


# Create simulation 
with Simulation("./tank_test_percremoval.inp") as sim:
    
    # Step through the simulation    
    for step in sim:

        conc2.append(sim._model.getNodePollutant("Tank", 0))
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
