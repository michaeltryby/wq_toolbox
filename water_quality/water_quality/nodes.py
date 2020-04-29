# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-04-27 14:50:28

from pyswmm.simulation import Simulation
import numpy as np
from scipy.integrate import ode 

########################################################################
# Event Mean Concentration Treatment
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: conc1, pindex2: conc2},
        'SWMM_Node_ID2': {pindex1: conc1, pindex2: conc2}}
"""
class EventMeanConc:
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict

    def treatment(self):
        # Read from user dictionary
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                # Set concentration
                sim._model.setNodePollutant(node, pollutant, self.node_dict[node][pollutant])

########################################################################
# CONSTANT REMOVAL TREATMENT
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: %R1, pindex2: %R2},
        'SWMM_Node_ID2': {pindex1: %R1, pindex2: %R2}}
"""
class ConstantRemoval:
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict

    def treatment(self):
        # Read from user dictionary
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                #Get parameters
                Cin = sim._model.getNodeCin(node, pollutant)
                R = self.node_dict[node][pollutant]
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration 
                sim._model.setNodePollutant(node, pollutant, Cnew)

########################################################################
# CO-REMOVAL TREATMENT
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: [%R1, %R_otherpollutant], pindex2: [%R1, %R_otherpollutant]},
        'SWMM_Node_ID2': {pindex1: [%R1, %R_otherpollutant], pindex2: [%R1, %R_otherpollutant]}}
"""
class CoRemoval:
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict

    def treatment(self):
        # Read from user dictionary
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                #Get parameters
                Cin = sim._model.getNodeCin(node, pollutant)
                C = sim._model.getNodePollutant(node, pollutant)
                R1 = self.node_dict[node][pollutant][0]
                R2 = self.node_dict[node][pollutant][1]
                # Calculate new concentration
                Cnew = (1-R1*R2)*Cin
                Cnew = min(Cin, C)
                # Set new concentration
                sim._model.setNodePollutant(node, pollutant, Cnew)

########################################################################
# CONCENTRATION-DEPENDENT REMOVAL
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: [%R_lower, boundary_conc, %R_upper], pindex2: [%R_lower, boundary_conc, %R_upper]},
        'SWMM_Node_ID2': {pindex1: [%R_lower, boundary_conc, %R_upper], pindex2: [%R_lower, boundary_conc, %R_upper]}}
"""
class ConcDependRemoval:
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict

    def treatment(self):
        # Read from user dictionary
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                # Get Cin for each pollutant/node
                Cin = sim._model.getNodeCin(node, pollutant)
                R_lower = self.node_dict[node][pollutant][0]
                bound_C = self.node_dict[node][pollutant][1]
                R_upper = self.node_dict[node][pollutant][2]
                # Calculate removal
                R = (1-np.heaviside((Cin-bound_C),0))*R_lower+np.heaviside((Cin-bound_C),0)*R_upper
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
class NthOrderReaction:
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 

    def treatment(self):
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step

        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                # Get parameters
                k = self.node_dict[node][pollutant][0]
                n = self.node_dict[node][pollutant][1]
                C = sim._model.getNodePollutant(node, pollutant)
                # Calculate treatment
                Cnew = C - (k*(C**n)*dt)
                # Set concentration each time step
                sim._model.setNodePollutant(node, pollutant, Cnew)

########################################################################
# K-C_STAR MODEL
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: [k, C_star], pindex2: [k, C_star]},
        'SWMM_Node_ID2': {pindex1: [k, C_star], pindex2: [k, C_star]}}
"""
class kCModel:
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict

    def treatment(self):
        # Read from user dictionary
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                # Get Cin for each pollutant/node
                Cin = sim._model.getNodeCin(node, pollutant)
                depth = sim._model.getNodeResult(node,5)
                k = self.node_dict[node][pollutant][0]
                C_star = self.node_dict[node][pollutant][1]
                hrt = sim._model.getNodeHRT(node)
                # Calculate removal
                if depth != 0.0 and Cin != 0.0:
                    R = np.heaviside((Cin-C_star),0)*((1-np.exp(-k*hrt/depth))*(1-C_star/Cin))
                else:
                    R = 0
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
class GravitySettling:
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 

    def treatment(self):
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step
        # Read from user dictionary
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                Qin = sim._model.getNodeResult(node,0)
                k = self.node_dict[node][pollutant][0]
                C_star = self.node_dict[node][pollutant][1]
                C = sim._model.getNodePollutant(node,pollutant)
                depth = sim._model.getNodeResult(node,5)
                if depth != 0.0:
                    # Calculate new concentration
                    Cnew = np.heaviside((0.1-Qin),0)*(C_star+(C-C_star)*np.exp(-k/depth*dt/3600))+(1-np.heaviside((0.1-Qin),0))*C
                else:
                    Cnew = np.heaviside((0.1-Qin),0)*C_star+(C-C_star)+(1-np.heaviside((0.1-Qin),0))*C
                # Set new concentration
                sim._model.setNodePollutant(node, pollutant, Cnew)

########################################################################
# CSTR
"""
Dictionary format: 
dict = {'SWMM_Node_ID1': {pindex1: [k, C_star], pindex2: [k, C_star]},
        'SWMM_Node_ID2': {pindex1: [k, C_star], pindex2: [k, C_star]}}
"""
class CSTR:
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 

    def tank(self, t, y):
        # Read from user dictionary
        for node in node_dict:
            for pollutant in node_dict[node]:
                Cin = sim._model.getNodeCin(node, pollutant)
                Qin = sim._model.getNodeResult(node, 0)
                Qout = sim._model.getNodeResult(node, 1)
                V = sim._model.getNodeResult(node, 3)
                k = Node_dict[node][pollutant][0]
                # Setup variables for ODE solver
                C = y[0]
                n = len(y)
                dCdt = np.zeros((n,1))
                dCdt[0] = (Qin*Cin - Qout*C)/V - k*C
                return dCdt

    def treatment(self):
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - last_timestep).total_seconds()
        # Updating reference step
        last_timestep = current_step

        # Read from user dictionary
        for node in node_dict:
            for pollutant in node_dict[node]:
                r = integrate.ode(tank).set_integrator('vode', method='bdf')
                t_start = 0.0
                t_final = dt
                num_steps = np.floor((t_final - t_start)/dt)+1

                C_t_zero = Cin
                r.set_initial_value([C_t_zero], t_start)

                t = np.zeros((num_steps, 1))
                C = np.zeros((num_steps, 1))
                t[0] = t_start
                C[0] = C_t_zero

                k = 1
                while r.successful() and k < num_steps:
                    r.integrate(r.t + dt)
                    t[k] = r.t
                    C[k] = r.y[0]
                    k += 1

            sim._model.setNodePollutant(node, pollutant, Cnew)

########################################################################