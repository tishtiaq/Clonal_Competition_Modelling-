

# Creating a simple LINEAR time function to replace a float as the fitness. 

# This function always finishes with a fitness of 1, yet takes different starting values. The slope is calculated 
# based on the starting value. 

# I have updated the bottom of the code so that it now plots comparison plots between simulations which uses the linear decay model and simulations which use the old model.
# The plots show the difference in muller plots, mean clone size and survival rate of clones. You can easily adapt population sizes, ratio of wild-type : higher fitness clones, 
# and starting fitness values by changing the stored floats in the code where labelled. 

import matplotlib.pyplot as plt
import numpy as np
from clone_competition_simulation import (WF, Moran, Parameters, FitnessParameters, TimeParameters,
                                          PopulationParameters, TreatmentParameters)

from clone_competition_simulation.simulation_algorithms.current_data import NonSpatialCurrentData
import clone_competition_simulation

class WFLinearFitness(WF):
    # Defining a class based on the linear function f(t)=a-bt
    
    def __init__(self, parameters, a_intercept):
        b_slope = (a_intercept - 1) / 10
        self.slope = b_slope
        self.intercept = a_intercept
        # Defines the gradient and intercept
        super().__init__(parameters)

# This works when our max_time is 10. If we change the max_time, we will have to change the slope calculation accordingly.
# Only need to change denominator in slope calculation if we change max_time.

    def get_next_generation(self, current_data:NonSpatialCurrentData) ->  np.ndarray[tuple[int], np.dtype[np.int_]]:
      # This function returns cell counts for the next generation. 
      current_time = self.i / self.division_rate
      new_fitness = self.intercept - self.slope * current_time   # Here we define the operation 
      self.clones_array[1, self.fitness_idx] = new_fitness
      # self.clones_array is the table listing id, label, fitness, generation, parent id. This line updates the fitness of the 
      # second clone (only, as the other is wild-type) based on the linear function.
      return super().get_next_generation(current_data) 
      # Returns the next generation of cells based on the current data and the updated fitness values.
params = Parameters(
    algorithm="WF", 
    times=TimeParameters(max_time=10, division_rate=1), 
    population=PopulationParameters(initial_size_array=np.array([1000, 9000])), # Change ratio of wild-type : higher fitness HERE
    fitness=FitnessParameters(initial_fitness_array=np.array([1, 1.05])) # Change fitness of fitter clone HERE
)

import os
docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')


sim_linear = WFLinearFitness(params, a_intercept=1.3)
sim_linear.run_sim()

sim_constant = WF(params)
sim_constant.run_sim()


# Below is the Muller Plots side by side
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
sim_linear.muller_plot(ax=ax1) # Plot the linear simulation on one side
ax1.set_title("Linearly Decreasing Fitness")
sim_constant.muller_plot(ax=ax2) # Plot the constant simulation on the other side
ax2.set_title("Constant Fitness")
fig1.suptitle("Muller Plots")
plt.savefig(os.path.join(docs_dir, '10percent_1.05_muller_comparison_wf.png'), dpi=150, bbox_inches='tight')


# Plot Mean Clone Sizes
fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 5))
sim_linear.plot_mean_clone_size_graph_for_non_mutation(ax=ax3)
ax3.set_title("Linearly Decreasing Fitness")
sim_constant.plot_mean_clone_size_graph_for_non_mutation(ax=ax4)
ax4.set_title("Constant Fitness")
fig2.suptitle("Mean Clone Size")
plt.savefig(os.path.join(docs_dir, '10percent_1.05_size_comparison_wf.png'), dpi=150, bbox_inches='tight')


# Plot Survival Rate on the same graph 
fig, ax = plt.subplots(figsize=(8, 5))  # creates 1 figure with two graphs
sim_constant.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Constant Fitness')
sim_linear.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Linearly Decreasing Fitness') 
ax.legend()
plt.title("Survival Comparison: Constant vs Linearly Decreasing Fitness")
plt.tight_layout
plt.savefig(os.path.join(docs_dir, '10percent_1.05_survival_comparison_wf.png'), dpi=150, bbox_inches='tight')


plt.show()