# OSCILLATING fitness function implementation for the clonal competition simulation.

import math
import matplotlib.pyplot as plt
import numpy as np
from clone_competition_simulation import (WF, Moran, Parameters, FitnessParameters, TimeParameters,
                                          PopulationParameters, TreatmentParameters)

from clone_competition_simulation.simulation_algorithms.current_data import NonSpatialCurrentData
import clone_competition_simulation

class OscillatingFitness(WF):
# Defining a class based on the sin function f(t)=a*cos(bt - c) + d

# NOTE: I am changing the function to be f(t)=a*cos(t-b) for simplicity.
    def __init__(self, parameters, a_amplitude, b_phase):
        self.amplitude = a_amplitude
        self.phase = b_phase
    
        # Defines the scaling of the function, phase shift, and vertical offset
        super().__init__(parameters)

    def get_next_generation(self, current_data:NonSpatialCurrentData) ->  np.ndarray[tuple[int], np.dtype[np.int_]]:
      # This function returns cell counts for the next generation. 
      current_time = self.i / self.division_rate
      new_fitness = self.amplitude * np.cos(current_time - self.phase)
      self.clones_array[1:, self.fitness_idx] = new_fitness
      # self.clones_array is the table listing id, label, fitness, generation, parent id. This line updates the fitness of the 
      # second clone (only, as the other is wild-type) based on the exponential function.
      return super().get_next_generation(current_data) 
      # Returns the next generation of cells based on the current data and the updated fitness values.
params = Parameters(
algorithm="WF", 
times=TimeParameters(max_time=10, division_rate=1), 
population=PopulationParameters(initial_size_array=np.concatenate([[9000], np.ones(1000, dtype=int)])),
fitness=FitnessParameters(initial_fitness_array=np.concatenate([[1.0], np.ones(1000)*1.5])),
)


import os
docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')

# Pass the parameters to the custom class

sim_oscillating= OscillatingFitness(params, a_coefficient=0.5, b_scale=1.0, c_phase=-1.0, d_offset=1.1)
sim_oscillating.run_sim()
# These exact coefficients range between roughly 0.5 and 1.5

# Run the simulation in a population of 10,000 with 1,000 wild-type
sim_constant = WF(params)
sim_constant.run_sim()

# Muller Plots Comparison (side by side)
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
sim_oscillating.muller_plot(ax=ax1) # Plot the linear simulation on one side
ax1.set_title("Oscillating Fitness")
sim_constant.muller_plot(ax=ax2) # Plot the constant simulation on the other side
ax2.set_title("Constant Fitness")
fig1.suptitle("Muller Plots")
plt.savefig(os.path.join(docs_dir, 'muller_oscillating_comparison_1000_mutants.png'), dpi=150, bbox_inches='tight')


# Create plot for survival rate
fig, ax = plt.subplots(figsize=(8, 5))  # creates 1 figure with two graphs
sim_constant.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Constant Fitness')
sim_oscillating.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Oscillating Fitness') 
ax.legend()
plt.title("Survival Comparison: Constant vs Oscillating Fitness")
plt.tight_layout
plt.savefig(os.path.join(docs_dir, 'constant_survival_plot_1000_mutants.png'), dpi=150, bbox_inches='tight')


# Calculate ratio
times = sim_constant.times
constant_surviving = np.array([np.sum(sim_constant.population_array.toarray()[:, i] > 0) for i in range(len(times))])
oscillating_surviving = np.array([np.sum(sim_oscillating.population_array.toarray()[:, i] > 0) for i in range(len(times))])

ratio = oscillating_surviving / constant_surviving
plt.figure(figsize=(8, 5))
plt.plot(times, ratio)
plt.axhline(1, color='grey', linestyle='--')
plt.xlabel('Time')
plt.ylabel('Oscillating / Constant survival')
plt.title('Survival ratio')
plt.tight_layout()
plt.savefig(os.path.join(docs_dir, 'survival_ratio.png'), dpi=150, bbox_inches='tight')


plt.show()
