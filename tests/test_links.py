from wq_toolbox.links import Link_Quality
from pyswmm import Simulation, Nodes, Links
import numpy as np
from sklearn.metrics import mean_squared_error as mse
import pytest
import matplotlib.pyplot as plt

"""
SWMM Water Quality Methods:
For each method, check the cummulative load in the link where the 
pollutant transformation is occurring is equivalent to the cummulative 
load directly downstream.

Additional Water Quality Methods:
For each method, check the cummulative load in the link where the 
pollutant transformation is occurring is equivalent to the cummulative 
load downstream.
"""

# SWMM WATER QUALITY METHODS
# Event Mean Concentration
def test_EventMeanConc_load():
    dict1 = {'Culvert': {0: 5.0}}
    conc = []
    conc1 = []
    flow = []
    flow1 = []
    with Simulation("./inps/LinkTest_variableinflow.inp") as sim:
        EMC = Link_Quality(sim, dict1)
        culvert = Links(sim)["Culvert"]
        outlet = Nodes(sim)["Outlet"]
        for step in sim:
            EMC.EventMeanConc()
            c = culvert.pollut_quality
            conc.append(c['P1'])
            c1 = outlet.pollut_quality
            conc1.append(c1['P1'])
            flow.append(sim._model.getLinkResult("Culvert",0))
            flow1.append(sim._model.getNodeResult("Outlet",0))
        load = [a*b for a,b in zip(conc,flow)]
        cum_load = np.cumsum(load)
        load1 = [a*b for a,b in zip(conc1,flow1)]
        cum_load1 = np.cumsum(load1)    
    error = (cum_load1[-1]/cum_load[-1])/cum_load1[-1]
    print(error)
    assert error <= 0.03


# Constant Removal
def test_ConstantRemoval_load():
    dict1 = {'Culvert': {0: 0.5}}
    conc = []
    conc1 = []
    flow = []
    flow1 = []
    with Simulation("./inps/LinkTest_variableinflow.inp") as sim:
        CR = Link_Quality(sim, dict1)
        culvert = Links(sim)["Culvert"]
        outlet = Nodes(sim)["Outlet"]
        for step in sim:
            CR.ConstantRemoval()
            c = culvert.pollut_quality
            conc.append(c['P1'])
            c1 = outlet.pollut_quality
            conc1.append(c1['P1'])
            flow.append(sim._model.getLinkResult("Culvert",0))
            flow1.append(sim._model.getNodeResult("Outlet",0))
        load = [a*b for a,b in zip(conc,flow)]
        cum_load = np.cumsum(load)
        load1 = [a*b for a,b in zip(conc1,flow1)]
        cum_load1 = np.cumsum(load1)    
    error = (cum_load1[-1]/cum_load[-1])/cum_load1[-1]
    print(error)
    assert error <= 0.03


# CoRemoval
def test_CoRemoval_load():
    dict1 = {'Culvert': {0: [0.75, 0.15]}}
    dict2 = {'Culvert': {1: 0.15}}
    conc = []
    conc1 = []
    flow = []
    flow1 = []
    with Simulation("./inps/LinkTest_variableinflow2.inp") as sim:
        CR = Link_Quality(sim, dict1)
        culvert = Links(sim)["Culvert"]
        outlet = Nodes(sim)["Outlet"]
        for step in sim:
            CR.CoRemoval()
            c = culvert.pollut_quality
            conc.append(c['P1'])
            c1 = outlet.pollut_quality
            conc1.append(c1['P1'])
            flow.append(sim._model.getLinkResult("Culvert",0))
            flow1.append(sim._model.getNodeResult("Outlet",0))
        load = [a*b for a,b in zip(conc,flow)]
        cum_load = np.cumsum(load)
        load1 = [a*b for a,b in zip(conc1,flow1)]
        cum_load1 = np.cumsum(load1)    
    error = (cum_load1[-1]/cum_load[-1])/cum_load1[-1]
    print(error)
    assert error <= 0.03
 

# ConcDependRemoval
def test_ConcDependRemoval_load():
    dict1 = {'Culvert': {0: [0.50, 10.0, 0.75]}}
    conc = []
    conc1 = []
    flow = []
    flow1 = []
    with Simulation("./inps/LinkTest_variableinflow.inp") as sim:
        CDR = Link_Quality(sim, dict1)
        culvert = Links(sim)["Culvert"]
        outlet = Nodes(sim)["Outlet"]
        for step in sim:
            CDR.ConcDependRemoval()
            c = culvert.pollut_quality
            conc.append(c['P1'])
            c1 = outlet.pollut_quality
            conc1.append(c1['P1'])
            flow.append(sim._model.getLinkResult("Culvert",0))
            flow1.append(sim._model.getNodeResult("Outlet",0))
        load = [a*b for a,b in zip(conc,flow)]
        cum_load = np.cumsum(load)
        load1 = [a*b for a,b in zip(conc1,flow1)]
        cum_load1 = np.cumsum(load1)    
    error = (cum_load1[-1]/cum_load[-1])/cum_load1[-1]
    print(error)
    assert error <= 0.03


# NthOrderReaction
def test_NthOrderReaction_load():
    dict1 = {'Culvert': {0: [0.01, 2.0]}}
    conc = []
    conc1 = []
    flow = []
    flow1 = []
    with Simulation("./inps/LinkTest_variableinflow.inp") as sim:
        NOR = Link_Quality(sim, dict1)
        culvert = Links(sim)["Culvert"]
        outlet = Nodes(sim)["Outlet"]
        for step in sim:
            NOR.NthOrderReaction()
            c = culvert.pollut_quality
            conc.append(c['P1'])
            c1 = outlet.pollut_quality
            conc1.append(c1['P1'])
            flow.append(sim._model.getLinkResult("Culvert",0))
            flow1.append(sim._model.getNodeResult("Outlet",0))
        load = [a*b for a,b in zip(conc,flow)]
        cum_load = np.cumsum(load)
        load1 = [a*b for a,b in zip(conc1,flow1)]
        cum_load1 = np.cumsum(load1)    
    error = (cum_load1[-1]/cum_load[-1])/cum_load1[-1]
    print(error)
    assert error <= 0.03

  
# GravitySettling
def test_GravitySettling_load():
    dict1 = {'Culvert': {0: [0.01, 10.0]}}
    conc = []
    conc1 = []
    flow = []
    flow1 = []
    with Simulation("./inps/LinkTest_variableinflow.inp") as sim:
        GS = Link_Quality(sim, dict1)
        culvert = Links(sim)["Culvert"]
        outlet = Nodes(sim)["Outlet"]
        for step in sim:
            GS.GravitySettling()
            c = culvert.pollut_quality
            conc.append(c['P1'])
            c1 = outlet.pollut_quality
            conc1.append(c1['P1'])
            flow.append(sim._model.getLinkResult("Culvert",0))
            flow1.append(sim._model.getNodeResult("Outlet",0))
        load = [a*b for a,b in zip(conc,flow)]
        cum_load = np.cumsum(load)
        load1 = [a*b for a,b in zip(conc1,flow1)]
        cum_load1 = np.cumsum(load1)    
    error = (cum_load1[-1]/cum_load[-1])/cum_load1[-1]
    print(error)
    assert error <= 0.03


# Erosion
def test_Erosion_load():
    dict1 = {'Channel': {0: [10.0, 0.001, 2.68, 0.7]}}
    conc = []
    conc1 = []
    flow = []
    flow1 = []
    with Simulation("./inps/LinkTest_variableinflow.inp") as sim:
        ER = Link_Quality(sim, dict1)
        channel = Links(sim)["Channel"]
        tailwater = Nodes(sim)["TailWater"]
        for step in sim:
            ER.Erosion()
            c = channel.pollut_quality
            conc.append(c['P1'])
            c1 = tailwater.pollut_quality
            conc1.append(c1['P1'])
            flow.append(sim._model.getLinkResult("Channel",0))
            flow1.append(sim._model.getNodeResult("Tailwater",0))
        load = [a*b for a,b in zip(conc,flow)]
        cum_load = np.cumsum(load)
        load1 = [a*b for a,b in zip(conc1,flow1)]
        cum_load1 = np.cumsum(load1)    
    error = (cum_load1[-1]/cum_load[-1])/cum_load1[-1]
    print(error)
    assert error <= 0.03

