#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 11:26:53 2026

@author: syedtariqishtiaq
"""

# Tutorial 1


import os
docs_path = "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/docs"
for root, dirs, files in os.walk(docs_path):
    for file in files:
        print(os.path.join(root, file))
        
        
with open("/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/docs/Introduction.md", "r") as f:
    print(f.read())




# Tells spyder where to find package (need this every time)
import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

# Loading the tools we need for the package
from clone_competition_simulation import Parameters, PopulationParameters, TimeParameters, FitnessParameters
import matplotlib.pyplot as plt


p = Parameters(
    algorithm='Moran',
# Moran algorithm used 
    times=TimeParameters(max_time=25, division_rate=1.4),
# Run simulation 25 times, cells dividing 1.4 times each time unit
    population=PopulationParameters(initial_size_array=[100, 100, 100]),
    fitness=FitnessParameters(initial_fitness_array=[1, 1.02, 1.04])
)
# Three clones, each with 100 cells, different fitness values


sim = p.get_simulator()
sim.run_sim()
# Create the simulator, run it, then plot it

sim.muller_plot(figsize=(5, 5))
plt.xlabel('Time')
plt.ylabel('Number of Cells')
plt.title('Clone Competition Simulation')
plt.show()


# So we clearly see that the cells with the higher fitness have a greater 
# chance of taking over. In biological terms: "the fittest clone underwent 
# clonal expansion and came to dominate the tissue through positive selection"











