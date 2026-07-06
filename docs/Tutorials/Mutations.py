#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 14:58:18 2026

@author: syedtariqishtiaq
"""

# Tutorial 7 


with open("/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/docs/Mutations.md", "r") as f:
    print(f.read())
    
    
import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

import matplotlib.pyplot as plt
from clone_competition_simulation import Parameters, TimeParameters, PopulationParameters, FitnessParameters, Gene, FitnessCalculator, UniformDist

gene1 = Gene(name='Gene1', mutation_distribution=UniformDist(0.5, 1.1), synonymous_proportion=0.4, weight=3)
gene2 = Gene(name='Gene2', mutation_distribution=UniformDist(1.1, 1.5), synonymous_proportion=0.5, weight=1)


# The weight=3 means a gene is 3 times more likely to be on
# the site of a mutation than that of weight=1. 


fit_calc = FitnessCalculator(genes=[gene1, gene2])
p = Parameters(
    algorithm='Moran2D',
    times=TimeParameters(max_time=10, division_rate=1),
    population=PopulationParameters(initial_cells=62500, cell_in_own_neighbourhood=False),
    fitness=FitnessParameters(
        mutation_rates=0.01,
        fitness_calculator=fit_calc
    )
)
s = p.get_simulator()
s.run_sim()

print("Overall dN/dS:", s.get_dnds())
print("Gene1 dN/dS:", s.get_dnds(gene='Gene1'))
print("Gene2 dN/dS:", s.get_dnds(gene='Gene2'))

# Above, the code separates out the dN/ds values for the 
# individual genes. This is important as it stops strongly 
# selected genes (such as NOTCH1) being diluted by other 
# genes that are being sequenced at the same time.

s.plot_dnds()
plt.xlabel('Time')
plt.ylabel('dN/ds')
plt.show()
    
'''
The early time points here are not reliable due to low 
sample size. At the very start of the simulation, only a 
handful of mutations have occurred, so the dN/dS ratio is 
essentially being calculated from a tiny, noisy sample. A 
single mutation event can swing the ratio dramatically.
'''


# The graph shows dN/ds being plotted over time. 


# dN is the ratio of nonsynonymous substitutions: mutations that 
# change the amino acid in a protein which can affect protein
# function
# ds is the ratio of synonymous substitutions: mutations that 
# do NOT change the amino acid in a protein which are invisible to 
# selection

# If dN/ds<1 this means negative selection as protein 
# function overall conserved (more invisibility)

# If dN/ds>1 then we see positive selection as protein is 
# evolving to adapt and we see more changes in protein 
# function 
    
# dN/ds=1 means neutral selection (combination of positive and 
# negative)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    