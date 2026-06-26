#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 14:50:12 2026

@author: syedtariqishtiaq
"""


# Tutorial 10


with open("/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/docs/SimulatingBiopsies.md", "r") as f:
    print(f.read())
    
    

import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

import numpy as np
from clone_competition_simulation import Parameters, TimeParameters, PopulationParameters, biopsy_sample, Biopsy

np.random.seed(0)
p = Parameters(
    algorithm='Moran2D',
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_grid=np.arange(100).reshape(10, 10),
                                    cell_in_own_neighbourhood=False)
)
s = p.get_simulator()
s.run_sim()

# See the full grid of clone ids
print("Full grid:")
print(s.grid_results[-1])

# s.grid_results is a list containing the grid's state at each
# individual time point. We want the final time point so use -1
# to start from the end


# Count cells per clone across the WHOLE grid
print("\nWhole grid clone counts:")
print(biopsy_sample(s.grid_results[-1], s, biopsy=None, remove_initial_clones=False))

# biopsy=None means to count the entire grid rather than restrict


# Take a smaller rectangular biopsy
biopsy = Biopsy(origin=(3, 4), shape=(5, 3))
print("\nBiopsy region:")
print(s.grid_results[-1][3:3+5, 4:4+3])

# This cuts out the specific grid that we want

print("\nBiopsy clone counts:")
print(biopsy_sample(s.grid_results[-1], s, biopsy=biopsy, remove_initial_clones=False))





'''
Handling subclones
'''


import numpy as np
from clone_competition_simulation import Gene, FitnessCalculator, NormalDist, FitnessParameters

fit_calc = FitnessCalculator(
    genes=[Gene(name="Gene1", mutation_distribution=NormalDist(0.1), synonymous_proportion=0.5)],
)

np.random.seed(0)
p = Parameters(
    algorithm='Moran2D',
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_cells=100,
                                    cell_in_own_neighbourhood=False),
    fitness=FitnessParameters(
        mutation_rates=0.1,
        fitness_calculator=fit_calc
    )
)
s = p.get_simulator()
s.run_sim()

print("Final grid:")
print(s.grid_results[-1])

biopsy = Biopsy(origin=(3, 4), shape=(5, 3))
print("\nBiopsy region:")
print(s.grid_results[-1][3:3+5, 4:4+3])

print("\nBiopsy counts WITH initial clones included:")
print(biopsy_sample(s.grid_results[-1], s, biopsy, remove_initial_clones=False))

print("\nBiopsy counts WITHOUT initial clones (default):")
print(biopsy_sample(s.grid_results[-1], s, biopsy))

# The above removes clone ID=0 because this is the only wild-
# type clone (clone without mutation). That is why 0 shows up 
# in the WITH count but not in the WITHOUT count. 


''' 
We had to fix a bug since the code couldn't convert NaN to an integer. Each 
clone started with clone ID=0.
'''


with open("/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src/clone_competition_simulation/tissue_sampling/sim_sampling.py", "r") as f:
    content = f.read()

import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

import numpy as np
import matplotlib.pyplot as plt
from clone_competition_simulation import (
    Parameters, TimeParameters, PopulationParameters, FitnessParameters,
    Gene, FitnessCalculator, FixedValue, biopsy_sample, get_vafs_for_all_biopsies
)

fit_calc = FitnessCalculator(
    genes=[
        Gene(name='Gene1', mutation_distribution=FixedValue(1.2), synonymous_proportion=0.5),
        Gene(name='Gene2', mutation_distribution=FixedValue(1.4), synonymous_proportion=0.5),
        Gene(name='Gene3', mutation_distribution=FixedValue(0.7), synonymous_proportion=0.5)
    ],
    combine_mutations='add'
)

np.random.seed(0)
p = Parameters(
    algorithm='WF2D',
    times=TimeParameters(max_time=200, division_rate=1),
    population=PopulationParameters(initial_cells=10000,
                                    cell_in_own_neighbourhood=False),
    fitness=FitnessParameters(
        mutation_rates=0.01,
        fitness_calculator=fit_calc
    )
)
s = p.get_simulator()
s.run_sim()

vafs = get_vafs_for_all_biopsies(s, biopsies=None, heterozygous=False)
plt.hist(vafs['vaf'], bins=np.linspace(0, 0.5, 50))
plt.xlabel('Variant allele fraction')
plt.ylabel('Frequency')
plt.title('Exact VAFs (no sequencing noise)')
plt.show()

print(f"Number of clones detected (exact VAFs): {len(vafs)}")

# Now simulate realistic sequencing
sequenced_vafs = get_vafs_for_all_biopsies(s, biopsies=None, detection_limit=5, coverage=100)
print(f"Number of clones detected with realistic sequencing: {len(sequenced_vafs)}")
print(sequenced_vafs.head(10))



'''
Final part of tutorial 10 below. This involves everything above, but also
includes the random coverage and dN/ds calculations from sequenced data
'''

from clone_competition_simulation import get_sample_dnds
import inspect
print(inspect.signature(get_sample_dnds))

for i in range(1, 4):
    gene = f"Gene{i}"
    print(f'dN/dS for {gene}:', get_sample_dnds(vafs, s, gene=gene))


'''
The above shows the dN/ds values for Gene1, Gene2 and Gene3. Each gene type
had a fixed fitness rate, and this code indirectlycalculates the dN/ds from 
that, by using the 
'''





