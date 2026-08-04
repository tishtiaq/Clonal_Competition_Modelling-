


# We're looking to plot the survival plots on the same graph here


import matplotlib.pyplot as plt 
import numpy as np

from clone_competition_simulation import (WF, Parameters, FitnessParameters, TimeParameters,
                                          PopulationParameters, TreatmentParameters)

from WFLinearFitness import WFLinearFitness

import os


np.random.seed(28) # could start with anything
n_clones = 10000     # number of clones in the initial population
mean_fitness = 1.3 # mean of the initial fitness values for the clones
std_fitness = 0.1 # standard deviation of the fitness values (up for changing?)

mutant_fitnesses = np.random.normal(loc=mean_fitness, scale=std_fitness, size=int(n_clones * 0.50)) # Draws 50% of the total population from a normal distribution: leaves 50% wild-type
mutant_fitnesses = np.clip(mutant_fitnesses, 1.01, None) # Rounds any fitness values outside the bounds to the interval edges.
                                                         # This specifically ensures that all fitness values are above 1.01, and None gives no ceiling
fitness_array = np.concatenate([[1], mutant_fitnesses]) # Joins the wild-type fitness value of 1 with the 99 mutant fitness values into a single array
n_wildtype_cells = int(n_clones * 0.5) # Define how many wild-type cells we want to start with (as a percentage of total cells)
size_array = np.concatenate([[n_wildtype_cells], np.ones(int(n_clones * 0.50), dtype=int)]) # Means that each clone starts with 1 cell


params = Parameters(
    algorithm="WF", 
    times=TimeParameters(max_time=1000, division_rate=1), 
    population=PopulationParameters(initial_size_array=size_array), # Sets up the size based on definition above
    fitness=FitnessParameters(initial_fitness_array=fitness_array), # Sets the fitness values for the clones in the initial population. The first clone is wild-type, and the rest are mutants with different fitness values.
)
    # Pass the parameters to the custom class

np.random.seed(28)
sim_constant = WF(params)
sim_constant.run_sim()

np.random.seed(28)
sim_linear = WFLinearFitness(params, a_intercept=1.5)
sim_linear.run_sim()


# This draws the survival plots on the same graph:
fig, ax = plt.subplots(figsize=(8, 5))  # creates 1 figure with two graphs

sim_constant.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Constant Fitness')
sim_linear.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Linearly Decreasing Fitness') 

ax.legend()
plt.title("Survival Comparison: Constant vs Linearly Decreasing Fitness\n"
          f"({n_clones} clones, mean fitness {mean_fitness}, SD {std_fitness})")
plt.tight_layout


docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')
plt.savefig(os.path.join(docs_dir, '50percent_wild_survival_comparison.png'), dpi=150, bbox_inches='tight')
plt.show()

