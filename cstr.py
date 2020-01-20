# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-20 15:01:08
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-20 15:03:01

def CSTR(self, reaction_rate, C_initial, solver_selected):
	#Dynamic CSTR equation
	Qin = env._getNodeInflow
	Qout = env._getLinkFlow
	Cin = env._getNodePollutant
	k = reaction_rate
	V = env._getNodeVolume

	# CSTR equation
	dCdt = Qin/V*Cin - Qout/V*Cout - k*C

	Solve.solver_selected()

	env._setNodePollutant()