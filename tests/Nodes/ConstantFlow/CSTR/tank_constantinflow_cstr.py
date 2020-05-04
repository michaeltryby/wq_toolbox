# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-04-27 14:20:21
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-05-04 14:44:36

from pyswmm import Simulation, Nodes
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

# Single Tank, Constant Inflow, CSTR
# PySWMM Toolbox
# Create simulation
conc1 = []

def tank(t, C, Qin, Cin, Qout, V, k, n): 
    # CSTR differential equation
    dCdt = (Qin*Cin - Qout*C)/V + k*C**n
    print(dCdt)
    return dCdt

with Simulation("./tank_constantinflow_notreatment.inp") as sim:
    k = - 0.2
    n = 1
    c0 = 0.0
    solver = ode(tank)
    Tank = Nodes(sim)["Tank"]
    start_time = sim.start_time
    last_timestep = sim.start_time

    # Step through the simulation    
    for index,step in enumerate(sim):
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - last_timestep).total_seconds()
        # Updating reference step
        last_timestep = current_step
        
        # Get newQual for Tank
        conc = Tank.pollut_quality
        conc1.append(conc['P1'])

        # Get parameters
        Qin = sim._model.getNodeResult('Tank',0)
        Cin = sim._model.getNodeCin('Tank',0)
        Qout = sim._model.getNodeResult('Tank',1)
        V = sim._model.getNodeResult('Tank',3)
        solver.set_f_params(Qin,Cin,Qout,V,k,n)

        if index == 0:
            solver.set_initial_value(c0, 0.0)
            solver.integrate(solver.t+dt)
        else: 
            solver.set_initial_value(solver.y, solver.t)
            solver.integrate(solver.t+dt)
        # Set new concentration
        sim._model.setNodePollutant('Tank', 0, solver.y[0])


class Node_Treatment:
    
    def __init__(self, sim, node_dict):
        self.sim = sim
        self.node_dict = node_dict
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 

    def CSTR(self):
        """
        UNSTEADY CONTINUOUSLY STIRRED TANK REACTOR (CSTR)
        CSTR is a common model for a chemical reactor. The behavior of a CSTR
        is modeled assuming it is not in steady state. This is because
        outflow, inflow, volume, and concentration are constantly changing.
        Therefore, an ODE solver is used to solve for cocentration.s
        
        Dictionary format: 
        dict = {'SWMM_Node_ID1': {pindex1: [k, n, c0], pindex2: [k, n, c0]},
                'SWMM_Node_ID2': {pindex1: [k, n, c0], pindex2: [k, n, c0]}}
        
        k   = reaction rate constant (SI or US: 1/s)
        n   = reaction order (first order, second order, etc.) (unitless)
        c0  = intital concentration inside reactor (SI or US: mg/L)
        """
        def tank(t, C, Qin, Cin, Qout, V, k, n):
            dCdt = (Qin*Cin - Qout*C)/V + k*C**n
            return dCdt

        def treatment(node_dict, sim, last_timestep):
            # Get current time
            current_step = sim.current_time
            # Calculate model dt in seconds
            dt = (current_step - last_timestep).total_seconds()
            # Updating reference step
            last_timestep = current_step

            # Setup solver
            solver = ode(tank)

            for node in node_dict:
                for pollutant in node_dict[node]:
                    # Get parameters
                    Qin = sim._model.getNodeResult(node,0)
                    Cin = sim._model.getNodeCin(node,pollutant)
                    Qout = sim._model.getNodeResult(node,1)
                    V = sim._model.getNodeResult(node,3)
                    k = node_dict[node][pollutant][0]
                    n = node_dict[node][pollutant][1]
                    c0 = node_dict[node][pollutant][2]
                    # Parameterize solver
                    solver.set_f_params(Qin,Cin,Qout,V,k,n)
                    if index == 0:
                        solver.set_initial_value(c0, 0.0)
                        solver.integrate(solver.t+dt)
                    else:
                        solver.set_initial_value(solver.y, solver.t)
                        solver.integrate(solver.t+dt)
                    # Set new concentration
                    sim._model.setNodePollutant(node, pollutant, solver.y[0])
        treatment(self.node_dict, self.sim, self.last_timestep)

conc2 = []
dict1 = {'Tank': {0: [-0.2, 1.0, 0.0]}}

with Simulation("./tank_constantinflow_notreatment.inp") as sim:
    Tank = Nodes(sim)["Tank"]
    TR = Node_Treatment(sim, dict1)
    # Step through the simulation    
    for index,step in enumerate(sim):
        # Get newQual for Tank
        co = Tank.pollut_quality
        conc2.append(co['P1'])
        # Run treatment each time step
        TR.CSTR()

# Plot Results
# Concentration Results
plt.plot(conc1, 'r--', label='no-class')
plt.plot(conc2, 'b:', label='class')
plt.ylabel("Concentration")
plt.xlabel("Time (s)")
plt.legend()
plt.show()
