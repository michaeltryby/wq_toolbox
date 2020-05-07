import wq_toolbox
import numpy as np
from sklearn.metrics import mean_squared_error as mse
import pytest

"""
SWMM Water Quality Methods:
For each method check the mean square deviation between the toolbox 
computed concentration value and the SWMM computed concentration value 
for the entire simulation.

For each method, check the cummulative load in the node where the 
pollutant transformation is occurring is equivalent to the cummulative
load downstream.

Additional Water Quality Methods:
For each method, check the cummulative load in the node where the 
pollutant transformation is occurring is equivalent to the cummulative 
load downstream.

Additionally, for CSTR, check the toolbox's calculated steady state 
concentration is equal to the closed form steady state CSTR equation. 
"""

# Event Mean Concentration



# Constant Removal



# CoRemoval



# ConcDependRemoval



# NthOrderReaction



# kCModel



# GravitySettling



# CSTR



# SedimentationResuspension


