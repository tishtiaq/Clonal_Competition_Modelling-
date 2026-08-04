

# Creating code for a WF Multiple Clones simulation with fitness remaining constant throughout. 
# This is to compare with my linearly decreasing model, to plot side by side.

import matplotlib.pyplot as plt 
import numpy as np

from clone_competition_simulation import (WF, Parameters, FitnessParameters, TimeParameters,
                                          PopulationParameters, TreatmentParameters)

from clone_competition_simulation.simulation_algorithms.current_data import NonSpatialCurrentData

np.random.seed(28) # could start with anything
n_clones = 10000     # number of clones in the initial population
mean_fitness = 1.3 # mean of the initial fitness values for the clones
std_fitness = 0.1 # standard deviation of the fitness values (up for changing?)

mutant_fitnesses = np.random.normal(loc=mean_fitness, scale=std_fitness, size=n_clones - 1) # Draws n-1 samples from a normal distribution 
mutant_fitnesses = np.clip(mutant_fitnesses, 1.01, None) # Rounds any fitness values outside the bounds to the interval edges.
                                                         # This specifically ensures that all fitness values are above 1.01, and None gives no ceiling
fitness_array = np.concatenate([[1], mutant_fitnesses]) # Joins the wild-type fitness value of 1 with the 99 mutant fitness values into a single array
size_array = np.ones(n_clones, dtype=int) # Means that each clone starts with 1 cell



params = Parameters(
    algorithm="WF", 
    times=TimeParameters(max_time=100, division_rate=1), 
    population=PopulationParameters(initial_size_array=size_array), # Sets up the size based on definition above
    fitness=FitnessParameters(initial_fitness_array=fitness_array), # Sets the fitness values for the clones in the initial population. The first clone is wild-type, and the rest are mutants with different fitness values.
)
    # Pass the parameters to the custom class
sim = WF(params)

sim.run_sim()

import os
docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')


sim.muller_plot(figsize=(5, 5))
plt.title("WF Competition Between One Wild-Type Cell and 9999 Fitter Cells of Constant Fitness")
plt.xlabel("Time")
plt.ylabel("Clone Size")
plt.savefig(os.path.join(docs_dir, 'muller_wf_constant10000clones.png'), dpi=150, bbox_inches='tight')

sim.plot_surviving_clones_for_non_mutation()
plt.title("Constant Fitness: Surviving Clones for Non-Mutation Case")
plt.savefig(os.path.join(docs_dir, 'survival_wf_constant10000clones.png'), dpi=150, bbox_inches='tight')
# Adding in a plot of the surviving clones for the non-mutation case.


sim.plot_mean_clone_size_graph_for_non_mutation()
plt.title("Constant Fitness: Mean Clone Size for Non-Mutation Case")
plt.savefig(os.path.join(docs_dir, 'meanclonesize_wf_constant10000clones.png'), dpi=150, bbox_inches='tight')
# Adding in a plot of the mean clone size for the non-mutation case.


plt.show()
# print("Final fitness of the second clone:", sim.clones_array[1, sim.fitness_idx])


'''
# This draws the muller plot and survival plot on the same graph:
fig, ax1 = plt.subplots(figsize=(10, 5))  # creates 1 figure with two graphs

sim.muller_plot(ax=ax1)  # draw muller plot on the first subplot
ax1.set_ylabel("Clone Size", color='black')  

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
sim.plot_surviving_clones_for_non_mutation(ax=ax2) # draw survival plot on the second subplot
ax2.set_ylabel("Surving Clones", color='blue')

plt.title("WF Constant Fitness - Clone Size and Survival")
plt.tight_layout()
plt.savefig(os.path.join(docs_dir, 'combined_wf_constant10000clones.png'), dpi=150, bbox_inches='tight')
'''







