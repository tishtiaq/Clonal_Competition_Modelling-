

# Making custom Moran (Non-Spatial)

# I have updated the bottom of the code so that it now plots comparison plots between simulations which uses the linear decay model and simulations which use the old model.
# The plots show the difference in muller plots, mean clone size and survival rate of clones. You can easily adapt population sizes, ratio of wild-type : higher fitness clones, 
# and starting fitness values by changing the stored floats in the code where labelled. 


import numpy as np 
import matplotlib.pyplot as plt

from clone_competition_simulation import (Moran, NonSpatialCurrentData, Parameters, PopulationParameters, 
                                          TimeParameters, FitnessParameters)


class MoranLinearFitness(Moran):

    def __init__(self, parameters, a_intercept):
        b_slope = (a_intercept - 1) / 10
        self.slope = b_slope
        self.intercept = a_intercept
        # Defines the gradient and intercept
        super().__init__(parameters)
    

    def get_dividing_cell(self, current_data: NonSpatialCurrentData) -> int:     # This function determines which cell will divide
# Returns the index of the cell that has highest fitness in the current_data.current_population array.

        current_time = self.i / self.division_rate
        new_fitness = max(self.intercept - self.slope * current_time, 0.01)   
        # Here we define the operation, and ensure fitness doesn't go below 0.01 to avoid negative fitness values
        self.clones_array[1, self.fitness_idx] = new_fitness

        return super().get_dividing_cell(current_data)

'''
    def get_differentiating_cell(self, current_data: NonSpatialCurrentData) -> int:  # This function determines which cell will die
        return 0
'''        
# We are commenting this out because we want to keep the parent class' code for this: we want cell death to be random. In the Moran model, usually 
# you need to overwrite 2 functions. 

params = Parameters(
    algorithm="Moran", times = TimeParameters(max_time=10, division_rate=1, samples=6),
    population=PopulationParameters(initial_size_array=np.array([5000,5000])),  # Change ratio of wild-type : higher fitness HERE
    fitness = FitnessParameters(initial_fitness_array=np.array([1, 1.3]))       # Change fitness of fitter clone HERE
    )


import os
docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')


sim_linear = MoranLinearFitness(params, a_intercept=1.3)
sim_linear.run_sim()

sim_constant = Moran(params)
sim_constant.run_sim()


# Below plots the Muller Plots side by side
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
sim_linear.muller_plot(ax=ax1) # Plot the linear simulation on one side
ax1.set_title("Linearly Decreasing Fitness")
sim_constant.muller_plot(ax=ax2) # Plot the constant simulation on the other side
ax2.set_title("Constant Fitness")
fig1.suptitle("Muller Plots")
plt.savefig(os.path.join(docs_dir, '50percent_1.3_muller_comparison_moran.png'), dpi=150, bbox_inches='tight')


# Plot Mean Clone Sizes
fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 5))
sim_linear.plot_mean_clone_size_graph_for_non_mutation(ax=ax3)
ax3.set_title("Linearly Decreasing Fitness")
sim_constant.plot_mean_clone_size_graph_for_non_mutation(ax=ax4)
ax4.set_title("Constant Fitness")
fig2.suptitle("Mean Clone Size")
plt.savefig(os.path.join(docs_dir, '50percent_1.3_size_comparison_moran.png'), dpi=150, bbox_inches='tight')


# Plot Survival Rate on the same graph 
fig, ax = plt.subplots(figsize=(8, 5))  # creates 1 figure with two graphs
sim_constant.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Constant Fitness')
sim_linear.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Linearly Decreasing Fitness') 
ax.legend()
plt.title("Survival Comparison: Constant vs Linearly Decreasing Fitness")
plt.tight_layout
plt.savefig(os.path.join(docs_dir, '50percent_1.3_survival_comparison_moran.png'), dpi=150, bbox_inches='tight')


plt.show()