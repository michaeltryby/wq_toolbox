from pyswmm.simulation import Simulation
import numpy as np


class Link_Quality:
    
    def __init__(self, sim, link_dict):
        self.sim = sim
        self.link_dict = link_dict
        self.start_time = self.sim.start_time
        self.last_timestep = self.start_time


    def EventMeanConc(self):
        """
        Event Mean Concentration Treatment (SWMM Water Quality Manual, 2016)
        Treatment results in a constant concentration.

        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: C, pindex2: C},
                'SWMM_Link_ID2': {pindex1: C, pindex2: C}}
        
        C = constant treatment concentration for each pollutant (SI or US: mg/L)
        """
        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                # Set concentration
                self.sim._model.setLinkPollutant(link, pollutant, self.link_dict[link][pollutant])


    def ConstantRemoval(self):
        """
        CONSTANT REMOVAL TREATMENT (SWMM Water Quality Manual, 2016)
        Treatment results in a constant percent removal.

        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: R, pindex2: R},
                'SWMM_Link_ID2': {pindex1: R, pindex2: R}}
        
        R = pollutant removal fraction (unitless)
        """
        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                #Get parameters
                Cin = self.sim._model.getLinkC2(link, pollutant)
                R = self.link_dict[link][pollutant]
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration 
                self.sim._model.setLinkPollutant(link, pollutant, Cnew)


    def CoRemoval(self):
        """
        CO-REMOVAL TREATMENT (SWMM Water Quality Manual, 2016)
        Removal of some pollutant is proportional to the removal of
        some other pollutant.

        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [R1, R2], pindex2: [R1, R2t]},
                'SWMM_Link_ID2': {pindex1: [R1, R2], pindex2: [R1, R2]}}
        
        R1 = pollutant removal fraction (unitless) 
        R2 = pollutant removal fraction for other pollutant (unitless)
        """
        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                #Get parameters
                Cin = self.sim._model.getLinkC2(link, pollutant)
                R1 = self.link_dict[link][pollutant][0]
                R2 = self.link_dict[link][pollutant][1]
                # Calculate new concentration
                Cnew = (1-R1*R2)*Cin
                # Set new concentration
                self.sim._model.setLinkPollutant(link, pollutant, Cnew)


    def ConcDependRemoval(self): 
        """
        CONCENTRATION-DEPENDENT REMOVAL (SWMM Water Quality Manual, 2016)
        When higher pollutant removal efficiencies occur with higher 
        influent concentrations.

        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [R_l, BC, R_u], pindex2: [R_l, BC, R_u]},
                'SWMM_Link_ID2': {pindex1: [R_l, BC, R_u], pindex2: [R_l, BC, R_u]}}
        
        R_l = lower removal rate (unitless)
        BC  = boundary concentration that determines removal rate (SI or US: mg/L)
        R_u = upper removal rate (unitless)
        """
        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                # Get Cin for each pollutant/link
                Cin = self.sim._model.getLinkC2(link, pollutant)
                R_l = self.link_dict[link][pollutant][0]
                BC = self.link_dict[link][pollutant][1]
                R_u = self.link_dict[link][pollutant][2]
                # Calculate removal
                R = (1-np.heaviside((Cin-BC),0))*R_l+np.heaviside((Cin-BC),0)*R_u
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration
                self.sim._model.setLinkPollutant(link, pollutant, Cnew)


    def NthOrderReaction(self):
        """
        NTH ORDER REACTION KINETICS (SWMM Water Quality Manual, 2016)
        When treatment of pollutant X exhibits n-th order reaciton kinetics
        where the instantaneous reaction rate is kC^n.

        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [k, n], pindex2: [k, n]},
                'SWMM_Link_ID2': {pindex1: [k, n], pindex2: [k, n]}}
        
        k   = reaction rate constant (SI: m/hr, US: ft/hr)
        n   = reaction order (first order, second order, etc.) (unitless)
        """
        # Get current time
        current_step = self.sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step

        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                # Get parameters
                k = self.link_dict[link][pollutant][0]
                n = self.link_dict[link][pollutant][1]
                C = self.sim._model.getLinkC2(link, pollutant)
                # Calculate treatment
                Cnew = C - (k*(C**n)*dt)
                # Set concentration each time step
                self.sim._model.setLinkPollutant(link, pollutant, Cnew)


    def GravitySettling(self):
        """
        GRAVITY SETTLING (SWMM Water Quality Manual, 2016)
        During a quiescent period of time within a storage volume, a fraction
        of suspended particles will settle out.

        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [k, C_s], pindex2: [k, C_s]},
                'SWMM_Link_ID2': {pindex1: [k, C_s], pindex2: [k, C_s]}}
        
        k   = reaction rate constant (SI: m/hr, US: ft/hr)
        C_s = constant residual concentration that always remains (SI or US: mg/L)
        """
        # Get current time
        current_step = self.sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step
        
        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                Q = self.sim._model.getLinkResult(link,0)
                k = self.link_dict[link][pollutant][0]
                C_s = self.link_dict[link][pollutant][1]
                C = self.sim._model.getLinkC2(link,pollutant)
                d = self.sim._model.getLinkResult(link,1)
                if d != 0.0:
                    # Calculate new concentration
                    Cnew = np.heaviside((0.1-Q),0)*(C_s+(C-C_s)*np.exp(-k/d*dt/3600))+(1-np.heaviside((0.1-Q),0))*C
                else:
                    Cnew = np.heaviside((0.1-Q),0)*C_s+(C-C_s)+(1-np.heaviside((0.1-Q),0))*C
                # Set new concentration
                self.sim._model.setLinkPollutant(link, pollutant, Cnew)


    def SedimentationResuspension(self):
        """
        SEDIMENTATION & RESUSPENSION
        Model considers both settling, as a function of depth, and resuspension,
        as a function of velocity.

        Dictionary format: 
        dict = {'SWMM_Link_ID1': {pindex1: [v_s, a], pindex2: [v_s, a]},
                'SWMM_Link_ID2': {pindex1: [v_s, a], pindex2: [v_s, a]}}
        
        v_s = settling velcity (SI: m/s, US: ft/s)
        a   = ratio between velcity and pollutant resuspension to result in 
              75% resuspension for the maximum velociy through storage pipe
        """
        # Get current time
        current_step = self.sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step

        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                Q = self.sim._model.getLinkResult(link,0)
                Cin = self.sim._model.getLinkC2(link,pollutant)
                v = self.sim._model.getConduitVelocity(link)
                d = self.sim._model.getLinkResult(link,1)
                v_s = self.link_dict[link][pollutant][0]
                a = self.link_dict[link][pollutant][1]
                # Calculate removal
                if d != 0.0 and v != 0.0:
                    R = 1 - np.exp(-v_s*dt/d) - np.exp(-a/(v*0.305))
                else:
                    R = 0
                # Calculate new concentration
                Cnew = (1-R)*Cin
                # Set new concentration
                self.sim._model.setLinkPollutant(link, pollutant, Cnew)


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
        Ss  = specific gravity of sediment (for soil usually between 2.65-2.80)
        d50 = mean sediment particle diameter (SI or US: mm)
        d   = depth (SI: m, US: ft)
        qs  = sediment discharge per unit width (SI: kg/m-s, US: lb/ft-s)
        Qs  = sediment discharge (SI: kg/s, US: lb/s)
        """
        # Get current time
        current_step = self.sim.current_time
        # Calculate model dt in seconds
        dt = (current_step - self.last_timestep).total_seconds()
        # Updating reference step
        self.last_timestep = current_step

        for link in self.link_dict:
            for pollutant in self.link_dict[link]:
                Cin = self.sim._model.getLinkC2(link,pollutant)
                Q = self.sim._model.getLinkResult(link,0)
                d = self.sim._model.getLinkResult(link,1)
                v = self.sim._model.getConduitVelocity(link)
                w = self.link_dict[link][pollutant][0]
                So = self.link_dict[link][pollutant][1]
                Ss = self.link_dict[link][pollutant][2]
                d50 = self.link_dict[link][pollutant][3]
                
                #Calculate erosion
                if self.sim._model.getSimUnit(0) == "US":
                    g = 32.2    # ft/s^2
                    yw = 62.4   # lb/ft^3
                    theta = (d*So/((Ss-1)*d50))*(1/0.00328) # unitless
                    if v != 0.0:
                        f = (2*g*So*d)/v**2     # unitless
                        qs = 0.1*(1/f)*theta**(5/2)*yw*((Ss-1)*g*(d50*0.00328)**3)**(1/2) # lb/ft-s
                        Qs = w*qs   # lb/s
                    else:
                        Qs = 0.0
                    if Q !=0.0:
                        Cnew = (Qs/Q)*(453592/28.3168)   # mg/L
                        Cnew = max(Cin, Cin+Cnew)
                        # Set new concentration
                        self.sim._model.setLinkPollutant(link, pollutant, Cnew)

                else:
                    g = 9.81    # m/s^2
                    yw = 1000   # kg/m^3
                    theta = (d*So/((Ss-1)*d50))*(1/0.001)   # unitless
                    if v != 0.0:
                        f = (2*g*So*d)/(v*0.305)**2     # unitless
                        qs = 0.1*(1/f)*theta**(5/2)*yw*((Ss-1)*g*(d50*0.001)**3)**(1/2) # kg/m-s
                        Qs = w*qs   # kg/s
                    else:
                        Qs = 0.0
                    if Q != 0.0:
                        Cnew = ((Qs/Q)*1000)  # mg/L
                        Cnew = max(Cin, Cin+Cnew)
                        # Set new concentration
                        self.sim._model.setLinkPollutant(link, pollutant, Cnew)

