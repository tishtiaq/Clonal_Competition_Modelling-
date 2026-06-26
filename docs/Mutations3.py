#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:47:07 2026

@author: syedtariqishtiaq
"""

# Tutorial 9


with open("/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/docs/Mutations3.md", "r") as f:
    print(f.read())


'''
We're going to run the Epistatic effects, modelling NOTCH1 and 
TP53 with a special combined effect when they are mutated
together.
'''

import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

import numpy as np
import matplotlib.pyplot as plt
from clone_competition_simulation import (
    Parameters, TimeParameters, PopulationParameters, FitnessParameters,
    Gene, FitnessCalculator, FixedValue, EpistaticEffect
)

# We import EpistaticEffect this time


fit_calc = FitnessCalculator(
    genes=[
        Gene(name='NOTCH1', mutation_distribution=FixedValue(1.1), synonymous_proportion=0),
        Gene(name='TP53', mutation_distribution=FixedValue(1.05), synonymous_proportion=0)
    ],
# Gene1 = NOTCH1-style gene, Gene2 = TP53-style gene 

    epistatics=[
        EpistaticEffect(
            name='NOTCH1_TP53_combo',
            gene_names=['NOTCH1', 'TP53'],
            fitness_distribution=FixedValue(3)
        )
    ],
    multi_gene_array=True,
    combine_mutations='replace',
)
# If BOTH genes are mutated, override with a special combined fitness of 3

np.random.seed(0)
p = Parameters(
    algorithm='Moran',
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_cells=6),
    fitness=FitnessParameters(
        mutation_rates=0.15,
        fitness_calculator=fit_calc,
    )
)
s = p.get_simulator()
s.run_sim()

print(s.view_clone_info(include_raw_fitness=True))

'''
This sets up a simulation where, alone, NOTCH1 has fitness 1.1, 
and TP53 has fitness 1.05. Yet when both genes are mutated, the 
fitness has a special value of 3 (greater than what you'd get
from multiplying or adding the values)

Note: NaN stands for Not a Number and is a common placeholder
'''

# Now we look at BoundedLogisticFitness, the way of capping how
# high fitness can climb even as more mutations occur


import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

import numpy as np
import matplotlib.pyplot as plt
from clone_competition_simulation import (
    Parameters, TimeParameters, PopulationParameters, FitnessParameters,
    Gene, FitnessCalculator, FixedValue, BoundedLogisticFitness
    )

b = BoundedLogisticFitness(3)
plt.plot(np.linspace(0, 10, 100), b.fitness(np.linspace(0, 10, 100)))
plt.xlabel('Raw fitness')
plt.ylabel('Transformed fitness')
plt.title('Diminishing returns: cap at 3')
plt.show()

fit_calc = FitnessCalculator(
    genes=[Gene(name='Gene1', mutation_distribution=FixedValue(1.4), synonymous_proportion=0)],
    combine_mutations='multiply',
    mutation_combination_class=BoundedLogisticFitness(3)
)

# Bounding this gives the shape of a logistic curve, and maps
# every value from the raw fitness to somewhere on this curve

np.random.seed(0)
p = Parameters(
    algorithm='Moran',
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_cells=1000),
    fitness=FitnessParameters(
        mutation_rates=0.1,
        fitness_calculator=fit_calc
    )
)
s = p.get_simulator()
s.run_sim()

print(s.view_clone_info(include_raw_fitness=True)[-10:])


'''
Slicing with the [-10:] gives us the last 10 rows. 


The bounded fitness formula comes from the formula:
    
    fitness = a / (1 + c × b^(-x))
where a is the maximum fitness ceiling we choose, b controls the
gradient of the curve, c is calculated as (a - 1) × b, so that
fitness(1)=1. 
'''












