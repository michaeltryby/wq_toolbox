# Working from this code for wq toolbox (1/7/2020)

# IMPORT 
# Import SWMM Modules
from pyswmm_lite import environment
from pyswmm_lite import water_quality

#SETUP
# Generate the state tuple for hydraulics
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
        treatments.append(("S"+str(i), "pfr", int(j)), 1.5) 
        # 2 pollutants/subcatchment, same "pfr" treatment, reaction rate

# Treatments for the nodes
for i in nodes:
    for j in pollutants:
        treatments.append(("N"+str(i), "cstr", int(j)), 0.4) 
        # 2 pollutants/node, same "cstr" treatment, reaction rate

# Treatments for the links
for i in links:
    for j in pollutants:
        treatments.append(("L"+str(i), "kcmodel", int(j), 2.0, 5.0) 
        # 2 pollutants/link, same "kcmodel" treatment, reaction rate, background conc.

# Build the treatment configuration dictionary
config2 = {
    "treatments": treatment,
    }


# SIMULATION
# Initialize the environment
env = environment(config1, ctrl=True)
wq = water_quality(config2, ctrl=True)
done = False
state = env.initial_state()
treatment = wq.initial_state()

# Run Simulation
while not done:
    
    # Steps water quality simulation: 1st uses getters, 
    # 2nd computes treatment, 3rd uses setters
    new_treatment, done = wq.step()

    # Steps the simulation
    new_state, done = env.step(np.ones(2))

    state = new_state
    treatment = new_treatment
    
# End Simulation & Close SWMM
env.sim._model.swmm_end()
env.sim._model.swmm_close()

