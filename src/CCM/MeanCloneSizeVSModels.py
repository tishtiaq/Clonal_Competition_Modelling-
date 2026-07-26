
# Creating code to take each individual fitness class and give us the mean clone size and survival rates to compare.

import matplotlib.pyplot as plt
import numpy as np
from clone_competition_simulation import (WF, Parameters, FitnessParameters, TimeParameters,
                                          PopulationParameters)

from clone_competition_simulation.simulation_algorithms.current_data import NonSpatialCurrentData
import clone_competition_simulation

from WFStepFitness import WFStepFitness
from WFExponentialFitnessNoSpace import ExponentialFitness
from WFLinearFitnessNoSpace import WFLinearFitness
from OscillatingFitnessFunction import OscillatingFitness

params = Parameters(
    algorithm="WF", 
    times=TimeParameters(max_time=100, division_rate=1), 
    population=PopulationParameters(initial_size_array=np.concatenate([[9000], np.ones(1000, dtype=int)])), # Change ratio of wild-type : higher fitness HERE
    fitness=FitnessParameters(initial_fitness_array=np.concatenate([[1.0], np.ones(1000)*1.05])) # Change fitness of fitter clone HERE
)

import os
docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')


# 1: Constant simulation run for control variable
sim_constant = WF(params)
sim_constant.run_sim()

# 2: Run Step model 
sim_step = WFStepFitness(params, a_intercept=1.05, t_time=10) 
sim_step.run_sim()

# 3: Run Linear model
sim_linear = WFLinearFitness(params, a_intercept=1.05)
sim_linear.run_sim()

# 4: Run Exponential model
sim_exponential = ExponentialFitness(params, a_coefficient=1.05)
sim_exponential.run_sim()


# Muller Plots side by side
fig1, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(12, 5))

sim_constant.muller_plot(ax=ax1) # Plot the constant simulation 
ax1.set_title("Constant Fitness")

sim_step.muller_plot(ax=ax2) # Plot the step simulation 
ax2.set_title("Step Fitness")

sim_linear.muller_plot(ax=ax3) # ...etc
ax3.set_title("Linear Fitness")

sim_exponential.muller_plot(ax=ax4)
ax4.set_title("Exponential Fitness")

fig1.suptitle("Muller Plots")
plt.savefig(os.path.join(docs_dir, 'models_comparison_muller_plot.png'), dpi=150, bbox_inches='tight')


'''
# Mean Clone Sizes
fig2, (ax5, ax6, ax7, ax8) = plt.subplots(1, 4, figsize=(12, 5))

sim_constant.plot_mean_clone_size_graph_for_non_mutation(ax=ax5)
ax5.set_title("Constant Fitness")

sim_step.plot_mean_clone_size_graph_for_non_mutation(ax=ax6)
ax6.set_title("Step Fitness")

sim_linear.plot_mean_clone_size_graph_for_non_mutation(ax=ax7)
ax7.set_title("Linear Fitness")

sim_exponential.plot_mean_clone_size_graph_for_non_mutation(ax=ax8)
ax8.set_title("Exponential Fitness")

fig2.suptitle("Mean Clone Size")
plt.savefig(os.path.join(docs_dir, 'models_comparison_mean_clone_size.png'), dpi=150, bbox_inches='tight')
'''

# Mean Clone Sizes on the same plot:
fig3, ax_combined = plt.subplots(figsize=(8, 5))
sim_constant.plot_mean_clone_size_graph_for_non_mutation(ax=ax_combined, legend_label='Constant Fitness')
sim_step.plot_mean_clone_size_graph_for_non_mutation(ax=ax_combined, legend_label='Step Fitness')
sim_linear.plot_mean_clone_size_graph_for_non_mutation(ax=ax_combined, legend_label='Linear Fitness')
sim_exponential.plot_mean_clone_size_graph_for_non_mutation(ax=ax_combined, legend_label='Exponential Fitness')
ax_combined.legend()
plt.title("Mean Clone Size of Different Models")
plt.tight_layout()
plt.savefig(os.path.join(docs_dir, 'models_comparison_mean_clone_size_combined.png'), dpi=150, bbox_inches='tight')


# Plot Survival Rate on the same graph 
fig, ax = plt.subplots(figsize=(8, 5))  # creates 1 figure with two graphs
sim_constant.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Constant Fitness')
sim_step.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Step Fitness')
sim_linear.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Linear Fitness')
sim_exponential.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Exponential Fitness')
ax.legend()
plt.title("Clone Survival Rate")
plt.tight_layout()
plt.savefig(os.path.join(docs_dir, 'models_comparison_survival.png'), dpi=150, bbox_inches='tight')


plt.show()
