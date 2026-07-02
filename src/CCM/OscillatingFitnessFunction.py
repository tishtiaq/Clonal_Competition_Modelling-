# OSCILLATING fitness function implementation for the clonal competition simulation.

import math
import matplotlib.pyplot as plt
import numpy as np
from clone_competition_simulation import (WF, Parameters, FitnessParameters, TimeParameters,
                                          PopulationParameters, TreatmentParameters)

from clone_competition_simulation.simulation_algorithms.current_data import NonSpatialCurrentData
import clone_competition_simulation

class OscillatingFitness(WF):
# Defining a class based on the sin function f(t)=a*cos(bt - c) + d
    def __init__(self, parameters, a_coefficient, b_scale, c_phase, d_offset):
        self.coefficient = a_coefficient
        self.scale = b_scale
        self.phase = c_phase
        self.offset = d_offset
        # Defines the scaling of the function, phase shift, and vertical offset
        super().__init__(parameters)
# TODO: # This works when our max_time is 10. If we change the max_time, we will have to change the slope calculation accordingly.
# Only need to change denominator in slope calculation if we change max_time.
    def get_next_generation(self, current_data:NonSpatialCurrentData) ->  np.ndarray[tuple[int], np.dtype[np.int_]]:
      # This function returns cell counts for the next generation. 
      current_time = self.i / self.division_rate
      new_fitness = self.coefficient * np.cos(self.scale * current_time - self.phase) + self.offset
      self.clones_array[1, self.fitness_idx] = new_fitness
      # self.clones_array is the table listing id, label, fitness, generation, parent id. This line updates the fitness of the 
      # second clone (only, as the other is wild-type) based on the exponential function.
      return super().get_next_generation(current_data) 
      # Returns the next generation of cells based on the current data and the updated fitness values.
params = Parameters(
algorithm="WF", 
times=TimeParameters(max_time=10, division_rate=1), 
population=PopulationParameters(initial_size_array=np.array((500, 500))),
fitness=FitnessParameters(initial_fitness_array=np.array([1, 1.5])),
)
    # Pass the parameters to the custom class
sim = OscillatingFitness(params, a_coefficient=0.2, b_scale=1.0, c_phase=-1.0, d_offset=1.0)
# TODO: Fix numbers to ensure we finish with a fitness of 1 for this starting fitness of 1.5

sim.run_sim()
sim.muller_plot(figsize=(5, 5))
plt.title("Competition Between One Wild-Type Clone and One Fitter Clone of Oscillating Fitness")
plt.xlabel("Time")
plt.ylabel("Clone Size")


import os
docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')
plt.savefig(os.path.join(docs_dir, 'muller_oscillating_500_500.png'), dpi=150, bbox_inches='tight')
# Saves the plot to the docs folder 

plt.show()
print("Final fitness of the second clone:", sim.clones_array[1, sim.fitness_idx])

