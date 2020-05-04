from pyswmm import Simulation, Nodes, Links
import numpy as np
import matplotlib.pyplot as plt

# Testing Link Getters & Setters, Variable Inflow
# PySWMM Toolbox

class Link_Treatment:
    
    def __init__(self, sim, link_dict):
        self.sim = sim
        self.link_dict = link_dict
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 

    def Erosion(self): 
        """
        ENGELUND-HANSEN EROSION (1967)
        Engelund and Hansen (1967) developed a procedure for predicting stage-
        discharge relationships and sediment transport in alluvial streams.

        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [w, So, Ss, d50], pindex2: [w, So, Ss, d50]},
                'SWMM_Link_ID2': {pindex1: [w, So, Ss, d50], pindex2: [w, So, Ss, d50]}}
        
        w   = channel width (SI: m, US: ft)
        So  = bottom slope (SI: m/m, US: ft/ft)
        d50 = mean sediment particle diameter (SI or US: mm)
        Ss  = specific gravity of sediment (for soil usually between 2.65-2.80)
        d   = depth (SI: m, US: ft)
        qs = sediment discharge per unit width (SI: kg/m-s, US: lb/ft-s)
        Qs = sediment discharge (SI: kg/s, US: lb/s)
        """
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step

        # Read from user dictionary
        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                Cin = sim._model.getLinkC2(link,pollutant)
                Qin = sim._model.getLinkResult(link,0)
                A = sim._model.getLinkResult(link,3)
                d = sim._model.getLinkResult(link,1)
                v = Qin/A
                w = self.link_dict[link][pollutant][0]
                So = self.link_dict[link][pollutant][1]
                Ss = self.link_dict[link][pollutant][2]
                d50 = self.link_dict[link][pollutant][3]
                
                #Calculate erosion
                if sim._model.getSimUnit(0) == "US":
                    g = 32.2    # ft/s^2
                    yw = 62.4   # lb/ft^3
                    theta = (d*So/((Ss-1)*d50))*(1/0.00328) # unitless
                    if v != 0.0:
                        f = (2*g*So*d)/v**2     # unitless
                        qs = 0.1*(1/f)*theta**(5/2)*yw*((Ss-1)*g*(d50*0.00328)**3)**(1/2) # lb/ft-s
                        Qs = w*qs   # lb/s
                    else:
                        Qs = 0.0
                    if Qin !=0.0:
                        Cnew = (Qs/Qin)*(453592/28.3168) + Cin   # mg/L
                        # Set new concentration
                        sim._model.setLinkPollutant(link, pollutant, Cnew)

                else:
                    g = 9.81    # m/s^2
                    yw = 1000   # kg/m^3
                    theta = (d*So/((Ss-1)*d50))*(1/0.001)   # unitless
                    if v != 0.0:
                        f = (2*g*So*d)/v**2     # unitless
                        qs = 0.1*(1/f)*theta**(5/2)*yw*((Ss-1)*g*(d50*0.001)**3)**(1/2) # kg/m-s
                        Qs = w*qs   # kg/s
                    else:
                        Qs = 0.0
                    if Qin != 0.0:
                        Cnew = ((Qs/Qin)*1000) + Cin  # mg/L
                        # Set new concentration
                        sim._model.setLinkPollutant(link, pollutant, Cnew)

# Create simulation
conc1 = []
conc2 = []
conc3 = []
flow1 = [] 
flow2 = []
flow3 = []
dict1 = {'Channel': {0: [10.0, 2.0, 2.68, 0.19]}}

with Simulation("./LinkTest_variableinflow.inp") as sim:
    inlet = Nodes(sim)["Inlet"]
    channel = Links(sim)["Channel"]
    tailwater = Nodes(sim)["TailWater"]
    LT = Link_Treatment(sim, dict1)

    # Step through the simulation    
    for step in sim:
        # Calculate erosion in node and add load to link
        LT.Erosion()
        # Get newQual for Tank
        c1 = inlet.pollut_quality
        c2 = channel.pollut_quality
        c3 = tailwater.pollut_quality
        conc1.append(c1['P1'])
        conc2.append(c2['P1'])
        conc3.append(c3['P1'])
        # Get flow for Tank
        flow1.append(sim._model.getNodeResult("Inlet", 0))
        flow2.append(sim._model.getLinkResult("Channel", 0))
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
plt.plot(conc2, 'b', label='channel')
plt.plot(conc3, 'm--', label='tailwater')
plt.ylabel("Concen")
plt.legend()

# Flow Results
plt.subplot(312)
plt.plot(flow1, 'g--', label='inlet')
plt.plot(flow2, 'b', label='channel')
plt.plot(flow3, 'm--', label='tailwater')
plt.ylabel("Flow")
plt.legend()

# Flow Results
plt.subplot(313)
plt.plot(cu_load1, 'g--', label='inlet')
plt.plot(cu_load2, 'b', label='channel')
plt.plot(cu_load3, 'm--', label='tailwater')
plt.xlabel("Time (s)")
plt.ylabel("Cumulative Load")
plt.legend()
plt.show()
