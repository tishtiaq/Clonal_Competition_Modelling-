#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 12:50:10 2026

@author: syedtariqishtiaq
"""

# Beginning of new code


import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

import numpy as np
import matplotlib.pyplot as plt
from clone_competition_simulation import Parameters, TimeParameters, PopulationParameters, FitnessParameters, Moran2D

class P53ExponentialDecay(Moran2D):
    """
    This is going to be a custom version of Moran2D. TP53-mutant 
    cells (identified by clone id) have fitness that decays 
    exponentially over time, from initial_fitness towards 
    baseline 1.0. Decay is based on time since each individual 
    clone's own birth. Fitness is recalculated once per sample 
    point (not continuously).
    """

    def __init__(self, parameters, tp53_clone_id, initial_fitness, decay_rate=0.1):
        super().__init__(parameters)
# Run all of Moran2D's normal setup first, unchanged
        self.tp53_clone_id = tp53_clone_id
        self.initial_fitness = initial_fitness
        self.decay_rate = decay_rate
        self.baseline_fitness = 1.0

        # Set the TP53 clone's starting fitness explicitly
        self.clones_array[self.tp53_clone_id, self.fitness_idx] = self.initial_fitness

    def _record_results(self, i, current_data, progress, task):
        super()._record_results(i, current_data, progress, task)

        # Only recalculate if we just took a sample (plot_idx was incremented inside _take_sample)
        if i == self.sample_points[self.plot_idx - 1]:
            birth_sample_idx = int(self.clones_array[self.tp53_clone_id, self.generation_born_idx])
            birth_time = self.times[birth_sample_idx]
            current_time = self.times[self.plot_idx - 1]
            age = current_time - birth_time

            new_fitness = self.baseline_fitness + (self.initial_fitness - self.baseline_fitness) * np.exp(-self.decay_rate * age)
            self.clones_array[self.tp53_clone_id, self.fitness_idx] = new_fitness


# Simple test setup: wild-type background + one TP53 clone, no NOTCH1 involved yet
grid = np.zeros((20, 20), dtype=int)
grid[8:12, 8:12] = 1  # TP53 patch in the middle

fitness_array = [1, 1.4]  # clone 0 = wild-type, clone 1 = TP53 starting at fitness 1.4

params = Parameters(
    algorithm='Moran2D',
    times=TimeParameters(max_time=50, division_rate=1, samples=50),
    population=PopulationParameters(initial_grid=grid.copy(), cell_in_own_neighbourhood=False),
    fitness=FitnessParameters(initial_fitness_array=fitness_array)
)

np.random.seed(0)
sim = P53ExponentialDecay(parameters=params, tp53_clone_id=1, initial_fitness=1.4, decay_rate=0.1)
sim.run_sim()

tp53_sizes = sim.population_array.toarray()[1]

plt.plot(sim.times, tp53_sizes, label='TP53 clone size')
plt.xlabel('Time')
plt.ylabel('Number of cells')
plt.title('TP53 with exponentially decaying fitness')
plt.legend()
plt.show()

print("TP53 fitness at end of simulation:", sim.clones_array[1, sim.fitness_idx])




import sys
sys.path.insert(0, "/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/clone-competition-simulation-master/src")

import numpy as np
import matplotlib.pyplot as plt
from clone_competition_simulation import Parameters, TimeParameters, PopulationParameters, FitnessParameters, Moran2D

class P53ExponentialDecay(Moran2D):
    def __init__(self, parameters, tp53_clone_id, initial_fitness, decay_rate=0.1):
        super().__init__(parameters)
        self.tp53_clone_id = tp53_clone_id
        self.initial_fitness = initial_fitness
        self.decay_rate = decay_rate
        self.baseline_fitness = 1.0
        self.clones_array[self.tp53_clone_id, self.fitness_idx] = self.initial_fitness

    def _record_results(self, i, current_data, progress, task):
        super()._record_results(i, current_data, progress, task)
        if i == self.sample_points[self.plot_idx - 1]:
            birth_sample_idx = int(self.clones_array[self.tp53_clone_id, self.generation_born_idx])
            birth_time = self.times[birth_sample_idx]
            current_time = self.times[self.plot_idx - 1]
            age = current_time - birth_time
            new_fitness = self.baseline_fitness + (self.initial_fitness - self.baseline_fitness) * np.exp(-self.decay_rate * age)
# The above reflects the exponential nature of the code
            self.clones_array[self.tp53_clone_id, self.fitness_idx] = new_fitness


def build_params():
    grid = np.zeros((20, 20), dtype=int)
    grid[8:12, 8:12] = 1
    fitness_array = [1, 1.4]
    return Parameters(
        algorithm='Moran2D',
        times=TimeParameters(max_time=50, division_rate=1, samples=50),
        population=PopulationParameters(initial_grid=grid.copy(), cell_in_own_neighbourhood=False),
        fitness=FitnessParameters(initial_fitness_array=fitness_array)
    )


n_repeats = 10
all_runs = []

for seed in range(n_repeats):
    np.random.seed(seed)
    params = build_params()
    sim = P53ExponentialDecay(parameters=params, tp53_clone_id=1, initial_fitness=1.4, decay_rate=0.1)
    sim.run_sim()
    tp53_sizes = sim.population_array.toarray()[1]
    all_runs.append(tp53_sizes)
    plt.plot(sim.times, tp53_sizes, color='lightblue', alpha=0.5)

all_runs = np.array(all_runs)
mean_trajectory = all_runs.mean(axis=0)
plt.plot(sim.times, mean_trajectory, color='darkblue', linewidth=2.5, label='Mean across 10 runs')

plt.xlabel('Time')
plt.ylabel('Number of cells')
plt.title('TP53 clone size: individual runs (light) vs mean (dark)')
plt.legend()
plt.show()

