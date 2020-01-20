# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-20 12:07:35
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-20 12:18:12

# IMPORT 
# Import SWMM Modules
from pyswmm_lite import environment
import wq_toolbox

# SETUP
# Generate the state tuple for hydraulics for pyswmm_lite
states = ["pollutantN", "pollutantL"]

# Build the hydraulics configuration dictionary
config1 = {
    "swmm_input": "./test_single.inp",
    "states": states,
    }

# Build the treatment configuration dictionary
config2 = {
    "pollutantN": "cstr", 1.5
    "pollutantL": "cstr", 1.5
    }

# SIMULATION
# Initialize the environment
env = environment(config1, ctrl=True)
wq = wq_toolbox(config2, env)
done = False
state = env.initial_state()

# Run Simulation
while not done:

    new_treatment, done = wq.step()

    """
    All this in the step: 
    
    # Steps water quality simulation:
    # 1st gets pollutant at current timestep
    wq.Treatment.getPollutant() 

    # 2nd runs treatment using specified solver
    wq.Treatment.CSTR(solver1)

    #  3rd sets new concentration at current time step
    wq.Treatment.setPollutant()"""

    # Steps the simulation
    new_state, done = env.step(np.ones(2))

    state = new_state
    treatment = new_treatment
    
# End Simulation & Close SWMM
env.sim._model.swmm_end()
env.sim._model.swmm_close()

