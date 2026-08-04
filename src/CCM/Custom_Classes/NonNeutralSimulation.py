

# Creating a non-neutral simulation with X wild-type cells and Y single-cell
# clones with a higher fitness (these ratios can be easily changed). This is taken from the 
# SimpleModels section of the PracticalTeacher notebook, but I am going to run one simulation 
# with constant fitness for the clones with higher fitness, and another where there fitness is 
# linearly decreasing. By comparing muller plots, mean clone size and survival rates, we will be
# able to see what effect a time-dependent fitness function has on the takeover in the tissue. 


from clone_competition_simulation import(WF, Parameters, FitnessParameters, TimeParameters, 
                                        PopulationParameters, LabelParameters)
import matplotlib.pyplot as plt
import numpy as np 

from WFLinearFitness import WFLinearFitness


MAX_TIME=25 
DIVISION_RATE=1.4 
SELECTIVE_FITNESS = 1.15   # Changing the initial fitness here so that there is a very small difference between the wild-type and non-wild-types
NEUTRAL_FITNESS = 1


# Run a simulation with 9900 wild type cells and 100 single-cell clones with higher fitness. 
params = Parameters(
    algorithm='WF',
    times=TimeParameters(max_time=MAX_TIME, division_rate=DIVISION_RATE), 
    population=PopulationParameters(
        initial_size_array=np.concatenate([
            np.array([9900]),   # The wild type cells tracked as one big clone
            np.ones(100)]   # 100 single-cell clones
    )),
    fitness=FitnessParameters(
        initial_fitness_array=np.concatenate([[1],   # Wild type cells with fitness 1
np.full(100, SELECTIVE_FITNESS)  # The next 100 clones given fitness 1.15
            ])
    ),
    labels=LabelParameters(
        initial_label_array=np.concatenate([
            [0],  # Wild type clone labelled with 0
            np.ones(100)  # Mutant clones labelled with 1
            ])
    )
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
# plt.savefig(os.path.join(docs_dir, '9900v100_1.15_muller_comparison_wf.png'), dpi=150, bbox_inches='tight')


# Plot Mean Clone Sizes
fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 5))
sim_linear.plot_mean_clone_size_graph_for_non_mutation(ax=ax3)
ax3.set_title("Linearly Decreasing Fitness")
sim_constant.plot_mean_clone_size_graph_for_non_mutation(ax=ax4)
ax4.set_title("Constant Fitness")
fig2.suptitle("Mean Clone Size")
# plt.savefig(os.path.join(docs_dir, '9900v100_1.15_size_comparison_wf.png'), dpi=150, bbox_inches='tight')


# Plot Survival Rate on the same graph 
fig, ax = plt.subplots(figsize=(8, 5))  # creates 1 figure with two graphs
sim_constant.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Constant Fitness')
sim_linear.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Linearly Decreasing Fitness') 
ax.legend()
plt.title("Survival Comparison: Constant vs Linearly Decreasing Fitness")
plt.tight_layout
# plt.savefig(os.path.join(docs_dir, '9900v100_1.15_survival_comparison_wf.png'), dpi=150, bbox_inches='tight')


plt.show()