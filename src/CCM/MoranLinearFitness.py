

# Making custom Moran (Non-Spatial)

import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

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
# We are commenting this out because we want to keep the parent class' code for this: we want cell death to be random

params = Parameters(
    algorithm="Moran", times = TimeParameters(max_time=10, division_rate=1, samples=6),
    population=PopulationParameters(initial_size_array=np.array([500,500])),
    fitness = FitnessParameters(initial_fitness_array=np.array([1, 1.5]))
    )

sim = MoranLinearFitness(params, a_intercept=1.5)

sim.run_sim()
sim.muller_plot(figsize=(5, 5))
plt.title("Moran Simulation Between One Wild-Type Clone and One Fitter Clone of Linearly Decreasing Fitness")
plt.xlabel("Time")
plt.ylabel("Clone Size")

print(sim.population_array.toarray())


import os
docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')
plt.savefig(os.path.join(docs_dir, 'moran_muller_linear_500_500.png'), dpi=150, bbox_inches='tight')
# Saves the file under name moran (make sure to get rid of the MullerPlotsMoran.md folder and these aren't actual Moran simulations)

plt.show()