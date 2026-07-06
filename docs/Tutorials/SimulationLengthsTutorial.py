#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:34:23 2026

@author: syedtariqishtiaq
"""

# Tutorial 3


with open("/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/docs/SimulationLength.md", "r") as f:
    
    print(f.read())


import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

import numpy as np
import matplotlib.pyplot as plt
from clone_competition_simulation import Parameters, TimeParameters, PopulationParameters

p = Parameters(
    algorithm='Moran',
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_size_array=np.ones(1000))
)

# The above looks at a clone with division rate 1


s = p.get_simulator()
s.run_sim()

p2 = Parameters(
    algorithm='Moran',
    times=TimeParameters(max_time=10, division_rate=1.7),
    population=PopulationParameters(initial_size_array=np.ones(1000))
)

# The above looks at a clone with division rate 1.7


s2 = p2.get_simulator()
s2.run_sim()

s.plot_mean_clone_size_graph_for_non_mutation(show_spm_fit=False, legend_label=1)
s2.plot_mean_clone_size_graph_for_non_mutation(ax=plt.gca(), show_spm_fit=False, legend_label=1.7)
plt.legend(title='Division rate')
plt.xlabel('Time')
plt.ylabel('Mean clone size')
plt.title('Effect of division rate on clone dynamics')
plt.show()

# We plots both lines onto the same graph and see how they compare

# We can see that the tissue with a higher division rate shows faster
# clonal expansion and larger mean clone sizes over time


















