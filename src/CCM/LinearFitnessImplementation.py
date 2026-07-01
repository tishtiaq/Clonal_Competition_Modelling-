

# Creating a simple LINEAR time function to replace a float as the fitness. 

# This function always finishes with a fitness of 1, yet takes different starting values. The slope is calculated 
# based on the starting value. 
import os
print(os.getcwd())

import matplotlib.pyplot as plt
import numpy as np
from clone_competition_simulation import (WF, Parameters, FitnessParameters, TimeParameters,
                                          PopulationParameters, TreatmentParameters)

from clone_competition_simulation.simulation_algorithms.current_data import NonSpatialCurrentData
import clone_competition_simulation

class LinearFitness(WF):
    # Defining a class based on the linear function f(t)=a-bt
    '''
    def __init__(self, parameters, b_slope, a_intercept):
        self.slope = b_slope
        self.intercept = a_intercept
        #Defines the gradient and intercept
        super().__init__(parameters)
      # This calls the superclass (WF) to begin the simulation with their provided parameters.

      # Redefined so that it only takes the starting value, and calculates the slope based off the starting value
      # and fact that we have to finish with a fitness of 1.
    '''
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
      new_fitness = self.intercept - self.slope * current_time
      self.clones_array[1, self.fitness_idx] = new_fitness
      # self.clones_array is the table listing id, label, fitness, generation, parent id. This line updates the fitness of the 
      # second clone (only, as the other is wild-type) based on the linear function.
      return super().get_next_generation(current_data) 
      # Returns the next generation of cells based on the current data and the updated fitness values.
params = Parameters(
    algorithm="WF", 
    times=TimeParameters(max_time=10, division_rate=1), 
    population=PopulationParameters(initial_size_array=np.array([250, 750])),
    fitness=FitnessParameters(initial_fitness_array=np.array([1, 1.5])),
)
    # Pass the parameters to the custom class
sim = LinearFitness(params, a_intercept=1.5)
# These exact numbers ensure we finish with a fitness of 1 for this starting fitness of 1.5

sim.run_sim()
sim.muller_plot(figsize=(5, 5))
plt.title("Competition Between One Wild-Type Clone and One Fitter Clone of Linearly Decreasing Fitness")
plt.xlabel("Time")
plt.ylabel("Clone Size")
plt.show()
# print("Final fitness of the second clone:", sim.clones_array[1, sim.fitness_idx])

# TODO: Add calculation of total mutant takeover/ VAF



'''
    p = Parameters(
        algorithm='WF', 
        times=TimeParameters(max_time=10, division_rate=1), 
        population=PopulationParameters(initial_size_array=[500, 500]),  # Start with equal clone sizes
        fitness=FitnessParameters(initial_fitness_array=[1, 1.5]),  # The second clone is fitter than the first (at the start)
        treatment=TreatmentParameters(
            # Define the treatment. 
            treatment_timings=[4],  # Start the treatment at time 4
            treatment_effects=[
                [1, 0.5]  # One value per initial clone. 
            ],  
            treatment_replace_fitness=False   # The `treatment_effects` will multiply the fitness, rather than replace it 
            # (if we were to write True instead)
        )
    ) This is old fitness code. 
    '''










