

# Using Moran2D
 
  
# Creating a TWO DIMENSIONAL non-neutral simulation with X wild-type cells and Y single-cell
# clones with a higher fitness (these ratios can be easily changed). This is taken from the 
# SimpleModels section of the PracticalTeacher notebook, but I am going to run one simulation 
# with constant fitness for the clones with higher fitness, and another where there fitness is 
# linearly decreasing. By comparing SPATIAL PLOTS
# we will be able to see what effect a time-dependent fitness function has on the takeover in the tissue. 


from clone_competition_simulation  import (Moran2D, Parameters, SpatialCurrentData, TimeParameters,
                                                             PopulationParameters, FitnessParameters)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  
from Moran2DLinearFitness import MoranLinearFitness2D


initial_grid = np.zeros((100,100), dtype=int) 
initial_grid[:50, :] = 1      # 2D slicing syntax: first colon is rows, second is columns. This assigns 50/50 split
                              # between mutants and wild-types. 

params = Parameters(
    algorithm="Moran2D", times= TimeParameters(max_time=100, division_rate=1, samples=4),
    population= PopulationParameters(initial_grid=initial_grid, cell_in_own_neighbourhood=False),
    fitness= FitnessParameters(initial_fitness_array=np.array([1, 1.15])),        # Change mutant fitness HERE and below
)


import os
docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')


sim_linear = MoranLinearFitness2D(params, a_intercept=1.15)          # Change mutant fitness HERE and above
sim_linear.run_sim()

sim_constant = MoranLinearFitness2D(params)
sim_constant.run_sim()


# Fig1 shows the Muller plots
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
sim_linear.muller_plot(ax=ax1) # Plot the linear simulation on one side
ax1.set_title("Linearly Decreasing Fitness")
sim_constant.muller_plot(ax=ax2) # Plot the constant simulation on the other side
ax2.set_title("Constant Fitness")
fig1.suptitle("Muller Plots")
plt.tight_layout()
# plt.savefig(os.path.join(docs_dir, 'Moran2D_1.15_muller_comparison.png'), dpi=150, bbox_inches='tight')


# Fig2 shows the Spatial grid which displays the physical layout at end of simulation
fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 5))

sim_linear.plot_grid(ax=ax3)
ax3.set_aspect('auto')
ax3.axis('off')
ax3.set_title("Linear Spatial Layout")

sim_constant.plot_grid(ax=ax4)
ax4.set_aspect('auto')
ax4.axis('off')
ax4.set_title("Constant Spatial Layout")

fig2.suptitle("Spatial Grids")
plt.tight_layout()
plt.savefig(os.path.join(docs_dir, 'Moran2D_1.15_spatial_comparison.png'), dpi=150, bbox_inches='tight')


plt.show()

