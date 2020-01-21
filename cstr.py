# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-20 15:01:08
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-21 11:26:25

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#SETUP WQ FUNCTION
def CSTR(C, t, k, V, Qin, Qout, Cin):

	# CSTR equation
	dCdt = (Qin*Cin - Qout*C)/V - k*C
	return dCdt


t = np.linspace(0, 40, 100)
sol = odeint(CSTR, 0.0, t, args=(0.10, 10.0, 1.0, 1.0, 10.0))

plt.plot(t, sol)
plt.ylim(0, 6.0)
plt.show()