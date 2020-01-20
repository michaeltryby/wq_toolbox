# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:56:47
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-01-15 12:42:47

# IMPORT 
# Import SWMM Modules
from pyswmm_lite import environment
import wq_toolbox

#SETUP
# Generate the state tuple for hydraulics for pyswmm_lite
states = []
subcatchment = np.linspace(1, 2, 2, dtype=int)
nodes = np.linspace(1, 2, 2, dtype=int)
links = np.linspace(1, 2, 2, dtype=int)
pollutants = np.linspace(1, 2, 2, dtype=int)

#States for the subcatchment
for i in subcatchment:
    for j in pollutants:
        states.append(("S"+str(i), "pollutantS", int(j))) # 2 pollutants/subcatchment

# States for the nodes
for i in nodes:
    for j in pollutants:
        states.append(("N"+str(i), "pollutantN", int(j))) # 2 pollutants/node

#States for the links
for i in links:
    for j in pollutants:
        states.append(("L"+str(i), "pollutantL", int(j))) # 2 pollutants/link

# Build the hydraulics configuration dictionary
config1 = {
    "swmm_input": "./test.inp",
    "states": states,
    }

#Generate treatment tuple for water quality
treatments = []

# Treatments for the subcatchment
for i in subcatchment:
    for j in pollutants:
        treatments.append(("S"+str(i), "pfr", int(j)), 1.5, "flow") 
        # 2 pollutants/subcatchment, same "pfr" treatment, reaction rate, SWMM data needed

# Treatments for the nodes
for i in nodes:
    for j in pollutants:
        treatments.append(("N"+str(i), "cstr", int(j)), 0.4, "depthN") 
        # 2 pollutants/node, same "cstr" treatment, reaction rate, SWMM data needed

# Treatments for the links
for i in links:
    for j in pollutants:
        treatments.append(("L"+str(i), "kcmodel", int(j), 2.0, 5.0, "flow") 
        # 2 pollutants/link, same "kcmodel" treatment, reaction rate, background conc., SWMM data needed

# Build the treatment configuration dictionary
config2 = {
    "treatments": treatment,
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

