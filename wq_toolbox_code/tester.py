# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-04-17 14:57:45

from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: conc1, pindex2: conc2},
        'SWMM_Node_ID2': {pindex1: conc1, pindex2: conc2}}
"""
# Example Constant Effluent Dictionary:
dict = {"'Node1'": {0: 5, 1: 10, 2: 15}, "'Node2'": {0: 2, 1: 4, 2: 6}}

def NodeConstantEffluent(dict):
    for i in dict:
        for j in dict[i]:
            # Set constant effluent concentration each time step
            sim._model.setNodePollutant(i, j, dict[i][j])

"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: %R1, pindex2: %R2},
        'SWMM_Node_ID2': {pindex1: %R1, pindex2: %R2}}
"""
# Example Constant Effluent Dictionary:
dict = {"'Node1'": {0: 5, 1: 10, 2: 15}, "'Node2'": {0: 2, 1: 4, 2: 6}}

def NodePercentRemoval(dict):
    Cin = []
    j = 0
    for i in dict:
        # Create list for Cin values
        Cin.append("Cin"+str(j))
        j += 1
        for k,l in zip(Cin, dict[i]):
            # Get influent concentration
            k = sim._model.getNodeCin(i, l)

            # Calculate treatment:

            sim._model.setNodePollutant(i, j, dict[i][j])



        # Get influent concentration
        Cin2 = sim._model.getNodeCin("2", 0)
        Cin5 = sim._model.getNodeCin("5", 0)
        # Calculate treatment
        Cnew2 = Cin2*(1-0.3)
        Cnew5 = Cin5*(1-0.1)
        # Set constant effluent concentration each time step
        sim._model.setNodePollutant("2", 0, Cnew2)
        sim._model.setNodePollutant("5", 0, Cnew5)
