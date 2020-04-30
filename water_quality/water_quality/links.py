# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-04-30 16:53:20

from pyswmm.simulation import Simulation
import numpy as np
from scipy.integrate import ode 


class Treatment:
    
    def __init__(self, sim, link_dict):
        self.sim = sim
        self.asset_dict = link_dict
        self.start_time = sim.start_time
        self.last_timestep = self.start_time 


    def EventMeanConc(self):
        # Event Mean Concentration Treatment
        """
        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: conc1, pindex2: conc2},
                'SWMM_Link_ID2': {pindex1: conc1, pindex2: conc2}}
        """
        # Read from user dictionary
        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                # Set concentration
                sim._model.setLinkPollutant(link, pollutant, self.link_dict[link][pollutant])


    def ConstantRemoval(self):
        # CONSTANT REMOVAL TREATMENT
        """
        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: %R1, pindex2: %R2},
                'SWMM_Link_ID2': {pindex1: %R1, pindex2: %R2}}
        """
        # Read from user dictionary
        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                #Get parameters
                Cin = sim._model.getLinkPollutant(link, pollutant)
                R = self.link_dict[link][pollutant]
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration 
                sim._model.setLinkPollutant(link, pollutant, Cnew)


    def CoRemoval(self):
        # CO-REMOVAL TREATMENT
        """
        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [%R1, %R_otherpollutant], pindex2: [%R1, %R_otherpollutant]},
                'SWMM_Link_ID2': {pindex1: [%R1, %R_otherpollutant], pindex2: [%R1, %R_otherpollutant]}}
        """
        # Read from user dictionary
        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                #Get parameters
                Cin = sim._model.getLinkPollutant(link, pollutant)
                R1 = self.link_dict[link][pollutant][0]
                R2 = self.link_dict[link][pollutant][1]
                # Calculate new concentration
                Cnew = (1-R1*R2)*Cin
                # Set new concentration
                sim._model.setLinkPollutant(link, pollutant, Cnew)


    def ConcDependRemoval(self):
        # CONCENTRATION-DEPENDENT REMOVAL
        """
        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [%R_lower, boundary_conc, %R_upper], pindex2: [%R_lower, boundary_conc, %R_upper]},
                'SWMM_Link_ID2': {pindex1: [%R_lower, boundary_conc, %R_upper], pindex2: [%R_lower, boundary_conc, %R_upper]}}
        """
        # Read from user dictionary
        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                # Get Cin for each pollutant/link
                Cin = sim._model.getLinkPollutant(link, pollutant)
                R_lower = self.link_dict[link][pollutant][0]
                bound_C = self.link_dict[link][pollutant][1]
                R_upper = self.link_dict[link][pollutant][2]
                # Calculate removal
                R = (1-np.heaviside((Cin-bound_C),0))*R_lower+np.heaviside((Cin-bound_C),0)*R_upper
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration
                sim._model.setLinkPollutant(link, pollutant, Cnew)


    def NthOrderReaction(self):
        # NTH ORDER REACTION KINETICS
        """
        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [k, n], pindex2: [k, n]},
                'SWMM_Link_ID2': {pindex1: [k, n], pindex2: [k, n]}}
        """
        # Get current time
        current_step = sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step

        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                # Get parameters
                k = self.link_dict[link][pollutant][0]
                n = self.link_dict[link][pollutant][1]
                C = sim._model.getLinkPollutant(link, pollutant)
                # Calculate treatment
                Cnew = C - (k*(C**n)*dt)
                # Set concentration each time step
                sim._model.setLinkPollutant(link, pollutant, Cnew)


    def kCModel(self):
        # K-C_STAR MODEL
        """
        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [k, C_star], pindex2: [k, C_star]},
                'SWMM_Link_ID2': {pindex1: [k, C_star], pindex2: [k, C_star]}}
        """
        # Read from user dictionary
        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                # Get Cin for each pollutant/link
                Cin = sim._model.getLinkPollutant(link, pollutant)
                depth = sim._model.getLinkResult(link,1)
                k = self.link_dict[link][pollutant][0]
                C_star = self.link_dict[link][pollutant][1]
                hrt = sim._model.getLinkHRT(link)
                # Calculate removal
                if depth != 0.0 and Cin != 0.0:
                    R = np.heaviside((Cin-C_star),0)*((1-np.exp(-k*hrt/depth))*(1-C_star/Cin))
                else:
                    R = 0
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration
                sim._model.setLinkPollutant(link, pollutant, Cnew) 


    def GravitySettling(self):
        # GRAVITY SETTLING
        """
        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [k, C_star], pindex2: [k, C_star]},
                'SWMM_Link_ID2': {pindex1: [k, C_star], pindex2: [k, C_star]}}
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
                Qin = sim._model.getLinkResult(link,0)
                k = self.link_dict[link][pollutant][0]
                C_star = self.link_dict[link][pollutant][1]
                C = sim._model.getLinkPollutant(link,pollutant)
                depth = sim._model.getLinkResult(link,1)
                if depth != 0.0:
                    # Calculate new concentration
                    Cnew = np.heaviside((0.1-Qin),0)*(C_star+(C-C_star)*np.exp(-k/depth*dt/3600))+(1-np.heaviside((0.1-Qin),0))*C
                else:
                    Cnew = np.heaviside((0.1-Qin),0)*C_star+(C-C_star)+(1-np.heaviside((0.1-Qin),0))*C
                # Set new concentration
                sim._model.setLinkPollutant(link, pollutant, Cnew)

    def CSTR(self):
        # CSTR
        """
        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [k, C_star, reaction_order], pindex2: [k, C_star, reaction_order]},
                'SWMM_Link_ID2': {pindex1: [k, C_star, reaction_order], pindex2: [k, C_star, reaction_order]}}
        """
        def tank(self, t, y):
            # Read from user dictionary
            for link in link_dict:
                for pollutant in link_dict[link]:
                    Cin = sim._model.getLinkCin(link, pollutant)
                    Qin = sim._model.getLinkResult(link,0)
                    # Do not have a link flow out-- is it the same as is?
                    #Qout = sim._model.getLinkResult(link,0)
                    V = sim._model.getLinkResult(link,2)
                    k = Link_dict[link][pollutant][0]
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
            for link in link_dict:
                for pollutant in link_dict[link]:
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

                sim._model.setLinkPollutant(link, pollutant, Cnew)
    
    def SedimentationResuspension(self):
        # SEDIMENTATION & RESUSPENSION 
        """
        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [v_s, a, b], pindex2: [v_s, a, b]},
                'SWMM_Link_ID2': {pindex1: [v_s, a, b], pindex2: [v_s, a, b]}}
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
                Qin = sim._model.getLinkResult(link,0)
                Cin = sim._model.getLinkPollutant(link,pollutant)
                depth = sim._model.getLinkResult(link,1)
                v_s = self.link_dict[link][pollutant][0]
                a = self.link_dict[link][pollutant][1]
                b = self.link_dict[link][pollutant][2]
                # Calculate removal
                if depth != 0.0 and Qin != 0.0:
                    R = 1 - np.exp(-v_s*dt/depth)-np.exp(-a*b/Qin)
                else:
                    R = 0
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration
                sim._model.setLinkPollutant(link, pollutant, Cnew)

    def Erosion(self): 
        # ENGELUND-HANSEN EROSION (1967)
        """
        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [, a, b], pindex2: [v_s, a, b]},
                'SWMM_Link_ID2': {pindex1: [v_s, a, b], pindex2: [v_s, a, b]}}
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
                Cin = sim._model.getLinkPollutant(link,pollutant)
                Qin = sim._model.getLinkResult(link,0)
                depth = sim._model.getLinkResult(link,1)
                v_s = self.link_dict[link][pollutant][0]
                a = self.link_dict[link][pollutant][1]
                b = self.link_dict[link][pollutant][2]
                Ss = 
                # Set new concentration
                sim._model.setLinkPollutant(link, pollutant, Cnew)

