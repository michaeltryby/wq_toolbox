# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-20 15:01:08
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-21 08:56:15

import numpy as np
from scipy.integrate import odeint

#SETUP WQ FUNCTION
def CSTR(Co,t,k):
	#Qin = env._getNodeInflow("P1")
	#Qout = env._getLinkFlow("7")
	#Cin = env._getNodePollutant("P1")
	#V = environment._getNodeVolume("P1")
	Qin = 50
	Qout =  25
	Cin = 10
	Cout  = 5
	V = 100

	# CSTR equation
	dCdt = Qin/V*Cin - Qout/V*Cout - k*Co
	return print(dCdt)

	#env._setNodePollutant("P1", dCdt)

k = 0.5
t =  np.linspace(0,10,11)
odeint(CSTR, 0.0, t, args=(k,))