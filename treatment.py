# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-20 15:42:16

# IMPORT MODULES
from pyswmm_lite import environment as env
from wq_toolbox.solver import Solve
import math


# MAIN CLASS FOR WQ_TOOLBOX
class Treatment:
	"""
	Attributes
    ----------
    config1 : dict
        dictionary for pyswmm_lite with swmm_input and, action and state space `(ID, attribute)`
	"""

    # TREATMENT OPTIONS
	# provides various treatment options
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

	"""
	def PFR(self, solver_selected):
		# PFR equation
		pfr = add_equation
		Solve.solver_selected()
	
	def KC_Model(self, solver_selected, C_star, k):
		# KC Model equation

		(C_out - C_star) = math.exp(-k/environment.get_NodeInflow)
		Solve.solver_selected()	

	# Maybe user just puts equation straight into dictionary instead?
	def User_Defined(self, equation, solver_selected):
		# Allows user to define their own treatment equation
		Solve.solver_selected()

	# METHODS TO ADD LATER
	# log water quality data
	def _logger(self,):

	#  estimates the performance
	def performance(self):

	# track a pollutant through the system
	def tracker(self):
	"""

