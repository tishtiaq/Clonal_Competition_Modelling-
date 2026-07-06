#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 13:03:01 2026

@author: syedtariqishtiaq
"""

# Tutorial 5

with open("/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/docs/DefiningInitialCloneFitness.md", "r") as f:
    
    print(f.read())


import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

import numpy as np
import matplotlib.pyplot as plt
from clone_competition_simulation import Parameters, TimeParameters, PopulationParameters, FitnessParameters

p = Parameters(
    algorithm='Moran2D',
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_grid=np.arange(16).reshape(4, 4), cell_in_own_neighbourhood=False),
    fitness=FitnessParameters(initial_fitness_array=np.random.uniform(0.5, 1.5, size=16))
)

# This generates 16 random fitness values between 0.5 and 1.5, one
# for each of the 16 clones 

s = p.get_simulator()
s.run_sim()

s.muller_plot(figsize=(5, 5))
plt.xlabel('Time')
plt.ylabel('Number of cells')
plt.title('16 clones with random fitness values')
plt.show()

# The plot shows how 16 different colours (each a different clone) 
# start off, with each having a different, randomly assigned fitness. 
# The colours most prominent at the end of the time will tend to be 
# those that started with the highest fitness, thereby showing how 
# clones with higher fitness tend to compete better. 






















