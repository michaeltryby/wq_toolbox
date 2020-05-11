# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-11 10:46:46

from pyswmm.simulation import Simulation
import numpy as np
from scipy.integrate import ode 


class Node_Quality:
    
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict
        self.start_time = self.sim.start_time
        self.last_timestep = self.start_time
        self.solver = ode(self.CSTR_tank)


    def EventMeanConc(self):
        """
        Event Mean Concentration Treatment (SWMM Water Quality Manual, 2016)
        Treatment results in a constant concentration.

        Dictionary format: 
        dict = {'SWMM_Node_ID1': {pindex1: C, pindex2: C},
                'SWMM_Node_ID2': {pindex1: C, pindex2: C}}
        
        C = constant treatment concentration for each pollutant (SI or US: mg/L)
        """
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                # Set concentration
                self.sim._model.setNodePollutant(node, pollutant, self.node_dict[node][pollutant])


    def ConstantRemoval(self):
        """
        CONSTANT REMOVAL TREATMENT (SWMM Water Quality Manual, 2016)
        Treatment results in a constant percent removal.

        Dictionary format: 
        dict = {'SWMM_Node_ID1': {pindex1: %R1, pindex2: %R2},
                'SWMM_Node_ID2': {pindex1: %R1, pindex2: %R2}}
        
        R = pollutant removal fraction (unitless)
        """
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                #Get parameters
                Cin = self.sim._model.getNodeCin(node, pollutant)
                R = self.node_dict[node][pollutant]
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration 
                self.sim._model.setNodePollutant(node, pollutant, Cnew)


    def CoRemoval(self):
        """
        CO-REMOVAL TREATMENT (SWMM Water Quality Manual, 2016)
        Removal of some pollutant is proportional to the removal of
        some other pollutant.

        Dictionary format: 
        dict = {'SWMM_Node_ID1': {pindex1: [R1, R2], pindex2: [R1, R2]},
                'SWMM_Node_ID2': {pindex1: [R1, R2], pindex2: [R1, R2]}}
        
        R1 = pollutant removal fraction (unitless) 
        R2 = pollutant removal fraction for other pollutant (unitless)
        """
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                #Get parameters
                Cin = self.sim._model.getNodeCin(node, pollutant)
                R1 = self.node_dict[node][pollutant][0]
                R2 = self.node_dict[node][pollutant][1]
                # Calculate new concentration
                Cnew = (1-R1*R2)*Cin
                # Set new concentration
                self.sim._model.setNodePollutant(node, pollutant, Cnew)


    def ConcDependRemoval(self):
        """
        CONCENTRATION-DEPENDENT REMOVAL (SWMM Water Quality Manual, 2016)
        When higher pollutant removal efficiencies occur with higher 
        influent concentrations.

        Dictionary format: 
        dict = {'SWMM_Node_ID1': {pindex1: [R_l, BC, R_u], pindex2: [R_l, BC, R_u]},
                'SWMM_Node_ID2': {pindex1: [R_l, BC, R_u], pindex2: [R_l, BC, R_u]}}
        
        R_l = lower removal rate (unitless)
        BC  = boundary concentration that determines removal rate (SI or US: mg/L)
        R_u = upper removal rate (unitless)
        """
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                # Get Cin for each pollutant/node
                Cin = self.sim._model.getNodeCin(node, pollutant)
                R_l = self.node_dict[node][pollutant][0]
                BC = self.node_dict[node][pollutant][1]
                R_u = self.node_dict[node][pollutant][2]
                # Calculate removal
                R = (1-np.heaviside((Cin-BC),0))*R_l+np.heaviside((Cin-BC),0)*R_u
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration
                self.sim._model.setNodePollutant(node, pollutant, Cnew)


    def NthOrderReaction(self):
        """
        NTH ORDER REACTION KINETICS (SWMM Water Quality Manual, 2016)
        When treatment of pollutant X exhibits n-th order reaciton kinetics
        where the instantaneous reaction rate is kC^n.

        Dictionary format: 
        dict = {'SWMM_Node_ID1': {pindex1: [k, n], pindex2: [k, n]},
                'SWMM_Node_ID2': {pindex1: [k, n], pindex2: [k, n]}}
        
        k   = reaction rate constant (SI: m/hr, US: ft/hr)
        n   = reaction order (first order, second order, etc.) (unitless)
        """
        # Get current time
        current_step = self.sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step

        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                # Get parameters
                k = self.node_dict[node][pollutant][0]
                n = self.node_dict[node][pollutant][1]
                C = self.sim._model.getNodeC2(node, pollutant)
                # Calculate treatment
                Cnew = C - (k*(C**n)*dt)
                # Set concentration each time step
                self.sim._model.setNodePollutant(node, pollutant, Cnew)


    def kCModel(self):
        """
        K-C_STAR MODEL (SWMM Water Quality Manual, 2016)
        The first-order model with bachground concnetration made popular by 
        Kadlec and Knight (1996) for long-term treatment performance of wetlands.

        Dictionary format: 
        dict = {'SWMM_Node_ID1': {pindex1: [k, C_s], pindex2: [k, C_s]},
                'SWMM_Node_ID2': {pindex1: [k, C_s], pindex2: [k, C_s]}}
        
        k   = reaction rate constant (SI: m/hr, US: ft/hr)
        C_s = constant residual concentration that always remains (SI or US: mg/L)
        """
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                # Get Cin for each pollutant/node
                Cin = self.sim._model.getNodeCin(node, pollutant)
                d = self.sim._model.getNodeResult(node,5)
                k = self.node_dict[node][pollutant][0]
                C_s = self.node_dict[node][pollutant][1]
                hrt = self.sim._model.getNodeHRT(node)
                # Calculate removal
                if d != 0.0 and Cin != 0.0:
                    R = np.heaviside((Cin-C_s),0)*((1-np.exp(-k*hrt/d))*(1-C_s/Cin))
                else:
                    R = 0
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration
                self.sim._model.setNodePollutant(node, pollutant, Cnew) 


    def GravitySettling(self):
        """
        GRAVITY SETTLING (SWMM Water Quality Manual, 2016)
        During a quiescent period of time within a storage volume, a fraction
        of suspended particles will settle out.

        Dictionary format: 
        dict = {'SWMM_Node_ID1': {pindex1: [k, C_s], pindex2: [k, C_s]},
                'SWMM_Node_ID2': {pindex1: [k, C_s], pindex2: [k, C_s]}}
        
        k   = reaction rate constant (SI: m/hr, US: ft/hr)
        C_s = constant residual concentration that always remains (SI or US: mg/L)
        """
        # Get current time
        current_step = self.sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step
        
        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                Qin = self.sim._model.getNodeResult(node,0)
                k = self.node_dict[node][pollutant][0]
                C_s = self.node_dict[node][pollutant][1]
                Cin = self.sim._model.getNodeCin(node,pollutant)
                d = self.sim._model.getNodeResult(node,5)
                if d != 0.0:
                    # Calculate new concentration
                    Cnew = np.heaviside((0.1-Qin),0)*(C_s+(Cin-C_s)*np.exp(-k/d*dt/3600))+(1-np.heaviside((0.1-Qin),0))*Cin
                else:
                    Cnew = np.heaviside((0.1-Qin),0)*C_s+(Cin-C_s)+(1-np.heaviside((0.1-Qin),0))*Cin
                # Set new concentration
                self.sim._model.setNodePollutant(node, pollutant, Cnew)


    def CSTR_tank(self, t, C, Qin, Cin, Qout, V, k, n):
        """
        UNSTEADY CONTINUOUSLY STIRRED TANK REACTOR (CSTR)
        CSTR is a common model for a chemical reactor. The behavior of a CSTR
        is modeled assuming it is not in steady state. This is because
        outflow, inflow, volume, and concentration are constantly changing.

        NOTE: You do not need to call this class, only the CSTR_solver. 
        CSTR_tank is intitalized in __init__ in Node_Treatment.  
        """
        dCdt = (Qin*Cin - Qout*C)/V + k*C**n
        return dCdt


    def CSTR_solver(self, index):
        """
        UNSTEADY CONTINUOUSLY STIRRED TANK REACTOR (CSTR) SOLVER
        CSTR is a common model for a chemical reactor. The behavior of a CSTR
        is modeled assuming it is not in steady state. This is because
        outflow, inflow, volume, and concentration are constantly changing.
        Therefore, Scipy.Integrate.ode solver is used to solve for concentration.
        
        NOTE: You only need to call this method, not CSTR_tank. CSTR_tank is
        intitalized in __init__ in Node_Treatment.  

        Dictionary format: 
        dict = {'SWMM_Node_ID1': {pindex1: [k, n, c0], pindex2: [k, n, c0]},
                'SWMM_Node_ID2': {pindex1: [k, n, c0], pindex2: [k, n, c0]}}
        
        k   = reaction rate constant (SI or US: 1/s)
        n   = reaction order (first order, second order, etc.) (unitless)
        c0  = intital concentration inside reactor (SI or US: mg/L)
        """
        # Get current time
        current_step = self.sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step

        for node in self.node_dict:
            for pollutant in self.node_dict[node]:
                # Get parameters
                Qin = self.sim._model.getNodeResult(node,0)
                Cin = self.sim._model.getNodeCin(node,pollutant)
                Qout = self.sim._model.getNodeResult(node,1)
                V = self.sim._model.getNodeResult(node,3)
                k = self.node_dict[node][pollutant][0]
                n = self.node_dict[node][pollutant][1]
                c0 = self.node_dict[node][pollutant][2]
                # Parameterize solver
                self.solver.set_f_params(Qin,Cin,Qout,V,k,n)
                # Solve ODE
                if index == 0:
                    self.solver.set_initial_value(c0, 0.0)
                    self.solver.integrate(self.solver.t+dt)
                else:
                    self.solver.set_initial_value(self.solver.y, self.solver.t)
                    self.solver.integrate(self.solver.t+dt)
                # Set new concentration
                self.sim._model.setNodePollutant(node, pollutant, self.solver.y[0])
        
