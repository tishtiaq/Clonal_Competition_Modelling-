#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:43:19 2026

@author: syedtariqishtiaq
"""

'''
import urllib.request
import zipfile
import os

url = "https://github.com/michaelhall28/clone-competition-simulation/archive/refs/heads/master.zip"

urllib.request.urlretrieve(url, "repo.zip")

with zipfile.ZipFile("repo.zip", "r") as z:
    z.extractall(".")

print("Done:", os.listdir("."))
'''
'''
# Checks what is in the clone-competition-simulation-master 
# folder
import os
print(os.listdir("clone-competition-simulation-master"))


# Install to import into spyder
import subprocess
import sys

subprocess.run([sys.executable, "-m", "pip", "install", "-e", "clone-competition-simulation-master"])


# What the README says
with open("clone-competition-simulation-master/README.md", "r") as f:
    print(f.read())
'''


import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/Python Projects/UCL Coding Project 2026/clone-competition-simulation-master/src")

from clone_competition_simulation.parameters import Parameters, TimeParameters, PopulationParameters, FitnessParameters
from clone_competition_simulation.fitness.fitness_classes import Gene, UniformDist, FitnessCalculator

print("Import successful!")

import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/Python Projects/UCL Coding Project 2026/clone-competition-simulation-master/src")

from clone_competition_simulation.parameters import Parameters, TimeParameters, PopulationParameters, FitnessParameters
from clone_competition_simulation.fitness.fitness_classes import Gene, UniformDist, FitnessCalculator

# Define fitness effects of mutations
fitness_calculator = FitnessCalculator(
    genes=[Gene(name='example_gene', mutation_distribution=UniformDist(1, 2), synonymous_proportion=0.5)],
    combine_mutations='multiply'
)

# Set up parameters
p = Parameters(
    algorithm='WF2D',
    population=PopulationParameters(grid_shape=(100, 100), cell_in_own_neighbourhood=True),
    times=TimeParameters(max_time=20, division_rate=1),
    fitness=FitnessParameters(mutation_rates=0.01, fitness_calculator=fitness_calculator)
)

# Run simulation
s = p.get_simulator()
s.run_sim()
s.muller_plot()


# Now have the plot of the code
















