# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-15 11:12:53

# IMPORT MODULES
from pyswmm_lite import environment


# MAIN CLASS FOR WQ_TOOLBOX
class Treatment:
	"""
	Attributes
    ----------
    config1 : dict
        dictionary for pyswmm_lite with swmm_input and, action and state space `(ID, attribute)`
	"""

	def getPollutant(self, config1, log=True):
		# get the pollutant concentration at current time step in each object
		if i==pollutantN in config1:
		    node_conc = environment._getNodePollutant
		elif i==pollutantL in config1:
		    link_conc = environment._getLinkPollutant
		elif i==pollutantS in config1:  
		    sub_conc = environment._getSubcatchPollutant

    def setPollutant(self, config1, log=True):
		# set the new pollutant concentration at current time step in each object
		if i==pollutantN in config1:
		    node_conc = environment._setNodePollutant
		elif i==pollutantL in config1:
		    link_conc = environment._setLinkPollutant
		elif i==pollutantS in config1:  
		    sub_conc = environment._setSubcatchPollutant
    

    # TREATMENT OPTIONS
	# provides various treatment options
	def PFR(self, reaction_rate, solver):
		# PFR equation
		pfr = add_equation

	
	def CSTR(self, reaction_rate, solver):
		# CSTR equation
		cstr = add_equation


	def KC_Model(self, reaction_rate, solver):
		# KC Model equation
		kc_model = add_equation

	# Maybe user just puts equation straight into dictionary instead?
	def User_Defined(self, equation, solver):
		# Allows user to define their own treatment equation



	# METHODS TO ADD LATER
	# log water quality data
	def _logger(self,):

	#  estimates the performance
	def performance(self):

	# track a pollutant through the system
	def tracker(self):


