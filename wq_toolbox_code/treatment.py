# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-04-23 15:51:07

from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt

########################################################################
# Event Mean Concentration Treatment
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: conc1, pindex2: conc2},
        'SWMM_Node_ID2': {pindex1: conc1, pindex2: conc2}}
"""
# Example EMC Dictionary:
dict1 = {'Node1': {0: 5, 1: 10, 2: 15}, 'Node2': {0: 2, 1: 4, 2: 6}}

class EventMeanConc:
    def __init__(self, sim):
        self.sim = sim
        
    def treatment(self, node_dict):
        # Read from user dictionary
        for node in node_dict:
            for pollutant in node_dict[node]:
                # Set concentration
                sim._model.setNodePollutant(node, pollutant, node_dict[node][pollutant])


########################################################################
# CONSTANT REMOVAL TREATMENT
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: %R1, pindex2: %R2},
        'SWMM_Node_ID2': {pindex1: %R1, pindex2: %R2}}
"""
# Example Constant Removal Dictionary
dict2 = {'Node1': {0: 0.5, 1: 0.10, 2: 0.15}, 'Node2': {0: 0.2, 1: 0.4, 2: 0.6}}

class ConstantRemoval:
    def __init__(self, sim):
        self.sim = sim

    def treatment(self, node_dict):
        # Read from user dictionary
        for node in node_dict:
            for pollutant in node_dict[node]:
                # Calculate new concentration
                Cnew = (1 - Node_dict[node][pollutant])*sim._model.getNodeCin(node, pollutant)
                # Set new concentration 
                sim._model.setNodePollutant(node, pollutant, Cnew)


########################################################################
# CO-REMOVAL TREATMENT
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: [%R1, %R_otherpollutant], pindex2: [%R1, %R_otherpollutant]},
        'SWMM_Node_ID2': {pindex1: [%R1, %R_otherpollutant], pindex2: [%R1, %R_otherpollutant]}}
"""
# Example Constant Removal Dictionary
dict2 = {'Node1': {0: [0.75, 0.50], 1: [0.10, 0.40]}, 
        'Node2': {0: [0.75, 0.50], 1: [0.10, 0.40]}}

class CoRemoval:
    def __init__(self, sim):
        self.sim = sim

    def treatment(self, node_dict):
        # Read from user dictionary
        for node in node_dict:
            for pollutant in node_dict[node]:
                # Calculate new concentration
                Cnew = (1 - Node_dict[node][pollutant][0]Node_dict[node][pollutant][1])*sim._model.getNodeCin(node, pollutant)
                # Set new concentration
                sim._model.setNodePollutant(node, pollutant, Cnew)


########################################################################
# CONCENTRATION-DEPENDENT REMOVAL
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: [%R_lower, boundary_conc, %R_upper], pindex2: [%R_lower, boundary_conc, %R_upper]},
        'SWMM_Node_ID2': {pindex1: [%R_lower, boundary_conc, %R_upper], pindex2: [%R_lower, boundary_conc, %R_upper]}}
"""
# Example Concentration-Dependent Removal Dictionary
dict2 = {'Node1': {0: [0.50, 50, 0.75], 1: [0.60, 20, 0.80]}, 
        'Node2': {0: [0.50, 50, 0.75], 1: [0.60, 20, 0.80]}}

class ConcDependRemoval:
    def __init__(self, sim):
        self.sim = sim

    def treatment(self, node_dict):
        # Read from user dictionary
        for node in node_dict:
            for pollutant in node_dict[node]:
                # Get Cin for each pollutant/node
                Cin = sim._model.getNodeCin(node, pollutant)
                # Calculate removal
                R = (1 - np.heaviside((Cin-Node_dict[node][pollutant][1])))*Node_dict[node][pollutant][0]+np.heaviside((Cin-Node_dict[node][pollutant][1]))*Node_dict[node][pollutant][2]
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration
                sim._model.setNodePollutant(node, pollutant, Cnew)


########################################################################
# NTH ORDER REACTION KINETICS
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: [k, n], pindex2: [k, n]},
        'SWMM_Node_ID2': {pindex1: [k, n], pindex2: [k, n]}}
"""
# Example Nth Order Reaction Kinetics Dictionary
dict2 = {'Node1': {0: [0.01, 0.5], 1: [0.01, 0.5], 2: [0.01, 0.5]}, 
        'Node2': {0: [0.01, 0.5], 1: [0.01, 0.5], 2: [0.01, 0.5]}}

class NthOrderReaction:
    def __init__(self, sim):
        self.sim = sim
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 

    def treatment(self, node_dict):
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step

        for node in node_dict:
            for pollutant in node_dict[node]:
                # Get current concentration
                C = sim._model.getNodePollutant(node, pollutant)
                # Calculate treatment
                Cnew = C - (node_dict[node][pollutant][0]*(C**node_dict[node][pollutant][1])*dt)
                # Set concentration each time step
                sim._model.setNodePollutant(node, pollutant, Cnew)


########################################################################
# K-C_STAR MODEL
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: [k, C_star], pindex2: [k, C_star]},
        'SWMM_Node_ID2': {pindex1: [k, C_star], pindex2: [k, C_star]}}
"""
# Example k-C_star Dictionary
dict2 = {'Node1': {0: [0.01, 20], 1: [0.01, 20], 2: [0.01, 20]}, 
        'Node2': {0: [0.01, 20], 1: [0.01, 20], 2: [0.01, 20]}}

class kCModel:
    def __init__(self, sim):
        self.sim = sim

    def treatment(self, node_dict):
        # Read from user dictionary
        for node in node_dict:
            for pollutant in node_dict[node]:
                # Get Cin for each pollutant/node
                Cin = sim._model.getNodeCin(node, pollutant)
                # Calculate removal
                R = np.heaviside((Cin-Node_dict[node][pollutant][1]),0)*((1 - np.exp(Node_dict[node][pollutant][0]*sim._model.getNodeHRT(node)/sim._model.getNodeResult(node, 5)))*(1-Node_dict[node][pollutant][0]/Cin))
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration
                sim._model.setNodePollutant(node, pollutant, Cnew)


########################################################################
# GRAVITY SETTLING
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: [k, C_star], pindex2: [k, C_star]},
        'SWMM_Node_ID2': {pindex1: [k, C_star], pindex2: [k, C_star]}}
"""
# Example Gravity Settling Dictionary
dict2 = {'Node1': {0: [0.01, 20], 1: [0.01, 20], 2: [0.01, 20]}, 
        'Node2': {0: [0.01, 20], 1: [0.01, 20], 2: [0.01, 20]}}

class GravitySettling:
    def __init__(self, sim):
        self.sim = sim
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 

    def treatment(self, node_dict, last_timestep):
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - last_timestep).total_seconds()
        # Updating reference step
        last_timestep = current_step

        # Read from user dictionary
        for node in node_dict:
            for pollutant in node_dict[node]:
                # Calculate new concentration
                Cnew = np.heaviside((0.1-sim._model.getNodeResult(node, 0)))*(Node_dict[node][pollutant][1]+(sim._model.getNodeCin(node, pollutant)-Node_dict[node][pollutant][1])*np.exp(-Node_dict[node][pollutant][0]/sim._model.getNodeResult(node, 5)*dt/3600))+(1-np.heaviside((0.1-sim._model.getNodeResult(node, 0))))*sim._model.getNodeCin(node, pollutant)
                # Set new concentration
                sim._model.setNodePollutant(node, pollutant, Cnew)


########################################################################
# CSTR
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: [k, C_star], pindex2: [k, C_star]},
        'SWMM_Node_ID2': {pindex1: [k, C_star], pindex2: [k, C_star]}}
"""
# Example Gravity Settling Dictionary
dict2 = {'Node1': {0: [0.01, 20], 1: [0.01, 20], 2: [0.01, 20]}, 
        'Node2': {0: [0.01, 20], 1: [0.01, 20], 2: [0.01, 20]}}

class CSTR:
    def __init__(self, sim):
        self.sim = sim
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 

    def treatment(self, node_dict):
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - last_timestep).total_seconds()
        # Updating reference step
        last_timestep = current_step

        # Read from user dictionary
        for node in node_dict:
            for pollutant in node_dict[node]:
                Cin = sim._model.getNodeCin(node, pollutant)
                Qin = sim._model.getNodeResult(node, 0)
                Qout = sim._model.getNodeResult(node, 1)
                V = sim._model.getNodeResult(node, 3)
                k = Node_dict[node][pollutant][0]
                dCdt = (Qin*Cin - Qout*C)/V - k*C


########################################################################
