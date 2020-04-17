from pyswmm import Simulation, Nodes, Links
import numpy as np
import matplotlib.pyplot as plt

# Testing Link Setters, Constant Inflow
# Add pollutant upstream in link and treatment in node downstream of link

# PySWMM Toolbox
# Create simulation
conc1 = []
conc2 = []
conc3 = []
flow1 = [] 
flow2 = []
flow3 = []

with Simulation("./ModifiedLinkTest_variableinflow.inp") as sim:
    inlet = Nodes(sim)["Inlet"]
    culvert = Links(sim)["Culvert"]
    outlet = Nodes(sim)["Outlet"]

    # Step through the simulation    
    for step in sim:
        # Get current concentrations
        C2 = sim._model.getLinkPollutant("Culvert", 0)
        C3 = sim._model.getNodePollutant("Outlet", 0)
        # Calculate treatment in node and added load in link
        Cnew2 = C2*(1.01)
        Cnew3 = C3*(1-0.01)
        # Set constant effluent concentration each time step
        sim._model.setLinkPollutant("Culvert", 0, Cnew2)
        sim._model.setNodePollutant("Outlet", 0, Cnew3)
        # Get newQual for Tank
        c1 = inlet.pollut_quality
        c2 = culvert.pollut_quality
        c3 = outlet.pollut_quality
        conc1.append(c1['P1'])
        conc2.append(c2['P1'])
        conc3.append(c3['P1'])
        # Get flow for Tank
        flow1.append(sim._model.getNodeResult("Inlet", 0))
        flow2.append(sim._model.getLinkResult("Culvert", 0))
        flow3.append(sim._model.getNodeResult("Outlet", 0))

load1 = [a*b for a,b in zip(conc1,flow1)]
load2 = [a*b for a,b in zip(conc2,flow2)]
load3 = [a*b for a,b in zip(conc3,flow3)]

cu_load1 = np.cumsum(load1)
cu_load2 = np.cumsum(load2)
cu_load3 = np.cumsum(load3)


# Plot Results
# Concentration Results
plt.subplot(311)
plt.plot(conc1, 'g--', label='inlet')
plt.plot(conc2, 'b', label='culvert')
plt.plot(conc3, 'm--', label='outlet')
plt.ylabel("Concen")
plt.legend()

# Flow Results
plt.subplot(312)
plt.plot(flow1, 'g--', label='inlet')
plt.plot(flow2, 'b', label='culvert')
plt.plot(flow3, 'm--', label='outlet')
plt.ylabel("Flow")
plt.legend()

# Flow Results
plt.subplot(313)
plt.plot(cu_load1, 'g--', label='inlet')
plt.plot(cu_load2, 'b', label='culvert')
plt.plot(cu_load3, 'm--', label='outlet')
plt.xlabel("Time (s)")
plt.ylabel("Cumulative Load")
plt.legend()
plt.show()
