# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-04-21 15:50:40

from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt

# CONSTANT EFFLUENT TREATMENT
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: conc1, pindex2: conc2},
        'SWMM_Node_ID2': {pindex1: conc1, pindex2: conc2}}
"""
# Example Constant Effluent Dictionary:
dict1 = {'Node1': {0: 5, 1: 10, 2: 15}, 'Node2': {0: 2, 1: 4, 2: 6}}

#Call during simulation
def NodeConstantEffluent(node_dict):
    # Read from user dictionary
    for node in node_dict:
        for pollutant in node_dict[node]:
            # Set constant effluent concentration each time step
            sim._model.setNodePollutant(str(node), pollutant, node_dict[node][pollutant])


########################################################################
# PERCENT REMOVAL TREAMTENT
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: %R1, pindex2: %R2},
        'SWMM_Node_ID2': {pindex1: %R1, pindex2: %R2}}
"""
# Example Percent Removal Dictionary
dict2 = {'Node1': {0: 0.5, 1: 0.10, 2: 0.15}, 'Node2': {0: 0.2, 1: 0.4, 2: 0.6}}

# Call during simulation
def NodePercentRemoval(Node_dict): 
    for node in Node_dict:
        for pollutant in Node_dict[node]:
            # Get Cin for each pollutant/node and append to Cin
            Cin = sim._model.getNodeCin(node, pollutant)
            # Calculate new concentration from percent removal treatment
            Cnew = (1 - Node_dict[node][pollutant])*Cin
            # Set new concentration each time step
            sim._model.setNodePollutant(node, pollutant, Cnew)

########################################################################
# NTH ORDER REMOVAL TREAMTENT
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: [k, n], pindex2: [k, n]},
        'SWMM_Node_ID2': {pindex1: [k, n], pindex2: [k, n]}}
"""
# Example Nth Order Reaction Dictionary
dict2 = {'Node1': {0: [0.01, 0.5], 1: [0.01, 0.5], 2: [0.01, 0.5]}, 
        'Node2': {0: [0.01, 0.5], 1: [0.01, 0.5], 2: [0.01, 0.5]}}

# Call before simulation
def setup_time():
    # Get simulation starting time
    start_time = sim.start_time
    # Get simulation ending time
    end_time = sim.end_time
    # Initial Value for Last Step
    last_timestep = start_time
    return last_timestep

# Call during simulation
def NthOrderReaction(Node_dict):
    # Get current time
    current_step = sim.current_time
    # Calculate model dt in seconds
    dt = (current_step - last_timestep).total_seconds()
    # Updating reference step
    last_timestep = current_step

    for node in Node_dict:
        for pollutant in Node_dict[node]:
            # Get current concentration
            C = sim._model.getNodePollutant(node, pollutant)
            # Calculate treatment
            Cnew = C - (Node_dict[node][pollutant][0]*(C**Node_dict[node][pollutant][0])*dt)
            # Set concentration each time step
            sim._model.setNodePollutant(node, pollutant, Cnew)












