#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 16:11:45 2026

@author: syedtariqishtiaq
"""


# Tutorial 4


with open("/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/docs/DefiningInitialCloneSizes.md", "r") as f:
    
    print(f.read())



import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

import numpy as np
import matplotlib.pyplot as plt
from clone_competition_simulation import Parameters, TimeParameters, PopulationParameters

# Two clones: top half clone 0, bottom half clone 1

initial_grid = np.repeat([0, 1], [50]).reshape(10, 10)

# If we wanted to instead have a 70/30 split, we'd write:
# initial_grid = np.repeat([0, 1], [70, 30]).reshape(10, 10)

p = Parameters(
    algorithm="Moran2D",
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_grid=initial_grid, cell_in_own_neighbourhood=False)
)
# "cell_in_own_neighbourhood=False" asks whether the cell can replace 
# itself: i.e are there 6 or 7 candidates for the cell's replacement?
s = p.get_simulator()
s.run_sim()

s.muller_plot(figsize=(5, 5))
plt.xlabel('Time')
plt.ylabel('Number of cells')
plt.title('Two clones competing on a 2D grid')
plt.show()




# We can also add fitness parameters


import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

import numpy as np
import matplotlib.pyplot as plt
from clone_competition_simulation import Parameters, TimeParameters, FitnessParameters, PopulationParameters

# Two clones: top half clone 0, bottom half clone 1

initial_grid = np.repeat([0, 1], [50]).reshape(10, 10)

# If we wanted to instead have a 70/30 split, we'd write:
# initial_grid = np.repeat([0, 1], [70, 30]).reshape(10, 10)

p = Parameters(
    algorithm="Moran2D",
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_grid=initial_grid, cell_in_own_neighbourhood=False),
    fitness=FitnessParameters(initial_fitness_array=[1, 1.05])  
)
# # clone 0 = normal, clone 1 = fitter


# "cell_in_own_neighbourhood=False" asks whether the cell can  
# replace itself: i.e are there 6 or are there 7 candidates for 
# the cell's replacement?
s = p.get_simulator()
s.run_sim()

s.muller_plot(figsize=(5, 5))
plt.xlabel('Time')
plt.ylabel('Number of cells')
plt.title('Two clones competing on a 2D grid (with fitness parameters)')
plt.show()

# As the fitness of the second clone is 2.05 it drastically overtakes
# the clone with fitness 1

# We could also keep the fitness of 1 and 1.05, and run the 
# simulation 100 times instead of 10, and see the effects






























