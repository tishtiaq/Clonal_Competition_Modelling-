#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 13:36:16 2026

@author: syedtariqishtiaq
"""

# Tutorial 6


with open("/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/docs/LineageTracingSimulations.md", "r") as f:
    print(f.read())
    
    
import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

import numpy as np
import matplotlib.pyplot as plt
from clone_competition_simulation import Parameters, TimeParameters, PopulationParameters

p = Parameters(
    algorithm='Moran',
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_size_array=np.ones(10000))
)
# This sets up 10,000 single-cell clones, division rate 1, and 
# it runs for 10 time steps in total

# We use the non-spatial Moran because the SPM isn't set up to 
# account for physical neighbours

s = p.get_simulator()
s.run_sim()

# The above actually runs the simulation. By the end, s contains
# full history of 10,000 clones sizes' over time


fig, ax = plt.subplots(figsize=(5, 5))
s.plot_mean_clone_size_graph_for_non_mutation(show_spm_fit=True, legend_label='Moran', legend_label_fit='SPM trend',
                                              fit_plot_kwargs={'c': 'r', 'linestyle': '--'}, ax=ax)

# This says to plot the mean clone size from the Moran simulation
# and to overlay the theoretical SPM prediction line on the same
# graph

plt.legend()
plt.xlabel('Time')
plt.ylabel('Mean clone size')
plt.show()


# The fact that the plots are always very close shows us that our
# simulation code is closely modelling the actual biology.




import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/Python Projects/UCL Coding Project 2026/clone-competition-simulation-master/src")

import numpy as np
import matplotlib.pyplot as plt
from clone_competition_simulation import (
    Parameters, TimeParameters, PopulationParameters,
    FitnessParameters, Gene, FitnessCalculator, NormalDist,
    LabelParameters, PlottingParameters, PLOT_COLOURS_EXAMPLE1
)

fitness_calculator = FitnessCalculator(
    genes=[Gene(name="Gene1", mutation_distribution=NormalDist(mean=1.1, var=0.1),
                synonymous_proportion=0.5)]
)

p = Parameters(
    algorithm='Moran',
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_cells=10000),
    fitness=FitnessParameters(
        mutation_rates=0.002,
        fitness_calculator=fitness_calculator
    ),
    labels=LabelParameters(
        label_times=3,
        label_frequencies=0.01,
        label_values=1,
        
# This tells us when to input the labels on the cells (as they're
# being added partway through the simulation)

    ),
    plotting=PlottingParameters(
        plot_colour_maps=PLOT_COLOURS_EXAMPLE1
    )
)
moran_label = p.get_simulator()
moran_label.run_sim()

moran_label.muller_plot(figsize=(7, 7))
plt.show()

'''
The code above is for the second graph which is plotted. The black 
crosses show that at time=3, the labels have been introduced and the
faint green lines show how these have been followed. 

'''


























