#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 16:20:00 2026

@author: syedtariqishtiaq
"""

# Tutorial 8

'''
Here we're adding in the idea of mutant clones having
further mutations. If clone 1 mutates further it will 
become clone 3, with a different ID. The population_array
only shows the final ID of the clone. If we wanted to see
all the mutations that a clone carried, we wouldn't be
able to in this format. That is why here, we can backtrack
mutations to see everything that was hit. 

This is very useful because it allows us to see what
proportion of mutations with P53 also had NOTCH1. 
'''

with open("/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/docs/Mutations2.md", "r") as f:
    print(f.read())



import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/Python Projects/UCL Coding Project 2026/clone-competition-simulation-master/src")

import matplotlib.pyplot as plt
import numpy as np
from clone_competition_simulation import Parameters, TimeParameters, PopulationParameters, FitnessParameters, Gene, FitnessCalculator, FixedValue

gene = Gene(name='Gene1', mutation_distribution=FixedValue(2), synonymous_proportion=0)
fit_calc = FitnessCalculator(genes=[gene])

np.random.seed(0)
p = Parameters(
    algorithm='Moran',
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_cells=50),
    fitness=FitnessParameters(
        mutation_rates=0.01,
        fitness_calculator=fit_calc,
    )
)

# This is a 1% probability per cell division that a 
# mutation happens

# We didn't specifically tell the mutation to occur at a
# certain time, but we DID say that when it did occur it 
# would have fitness 2

s = p.get_simulator()
s.run_sim()
s.muller_plot(figsize=(5, 5))

plt.xlabel('Time')
plt.ylabel('Number of cells')
plt.show()

print(s.view_clone_info())
print("Mutant clone sizes:", s.get_mutant_clone_sizes())
print("Ancestors of clone 5:", s.get_clone_ancestors(5))
print("Descendants of clone 1:", s.get_clone_descendants(1))


'''
Now we produce 4 plots showing the same simulation but 
displayed with increasing clarity 
'''

import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/Python Projects/UCL Coding Project 2026/clone-competition-simulation-master/src")

import matplotlib.pyplot as plt
import numpy as np
from clone_competition_simulation import Parameters, TimeParameters, PopulationParameters, FitnessParameters, Gene, FitnessCalculator, FixedValue

# A busier simulation with more cells and mutations
gene = Gene(name='Gene1', mutation_distribution=FixedValue(1.5), synonymous_proportion=0.5)
fit_calc = FitnessCalculator(genes=[gene])

np.random.seed(0)
p = Parameters(
    algorithm='Moran',
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_cells=10000),
    fitness=FitnessParameters(
        mutation_rates=0.01,
        fitness_calculator=fit_calc,
    )
)
s = p.get_simulator()
s.run_sim()

# Busy default plot
s.muller_plot(figsize=(5, 5))
plt.title("Default (busy)")
plt.show()

# Cleaner - no X markers
s.muller_plot(figsize=(5, 5), show_mutations_with_x=False)
plt.title("No X markers")
plt.show()

# Cleaner still - only clones that reached size 30+
s.muller_plot(figsize=(5, 5), show_mutations_with_x=False, min_size=30)
plt.title("Min size 30 only")
plt.show()

# Incomplete moment plot
fig, ax = plt.subplots(figsize=(5, 5))
s.plot_incomplete_moment(ax=ax)
plt.xlabel('Time')
plt.show()

'''
The incomplete moment plot is a statistical summary that
shows the average clone size, ONLY counting clones that 
are above a certain size threshold. The decreasing pattern
reflects clones progressively dying out or growing, which
causes the distribution to shift. 

At time x=0, the first incomplete moment is 1.0, as this
means 100% of cells belong to clones of size 0 or bigger. 
Then as x increases (we raise the size threshold) the 
y-axis tells us what fraction of cells are in clone sizes
above THAT size. 

This is a key summary to use in my model, as this 
quantitatively compares clone sizes, rather than relying 
on visual comparisons of Muller plots. 
'''























