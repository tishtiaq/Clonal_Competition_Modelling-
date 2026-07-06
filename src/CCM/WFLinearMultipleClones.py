

# In this code we generate 100 clones with different starting fitness values (selected from a normal distribution) and run a 
# WF simulation with linearly decreasing fitness for each clone. The fitness of each clone decreases linearly over time, with the 
# slope determined by its initial fitness value (as each clone must end with fitness of 1). 

# We draw muller plots, survival of clone plots and mean clone size plots for the simulation.

import matplotlib.pyplot as plt
import numpy as np
from clone_competition_simulation import (WF, Parameters, FitnessParameters, TimeParameters,
                                          PopulationParameters, TreatmentParameters)

from clone_competition_simulation.simulation_algorithms.current_data import NonSpatialCurrentData
import clone_competition_simulation

class WFLinearFitness(WF):
    # Defining a class based on the linear function f(t)=a-bt
    
    def __init__(self, parameters):
        self.initial_fitness_array = parameters.fitness.initial_fitness_array # Defines the initial fitness' of clones
        
        self.slopes = (self.initial_fitness_array - 1) / parameters.times.max_time # As the starting fitness' of the clones are different, we need to calculate the slope for each clone individually.
        super().__init__(parameters)

        # This is now defined for multiple clones with different starting fitness values. The slope is calculated for each clone based on its starting fitness value.

    def get_next_generation(self, current_data:NonSpatialCurrentData) ->  np.ndarray[tuple[int], np.dtype[np.int_]]:
      # This function returns cell counts for the next generation. 
      current_time = self.i / self.division_rate
      for clone_id in range(1, len(self.clones_array)):
         starting_fitness = self.initial_fitness_array[clone_id]
         self.clones_array[clone_id, self.fitness_idx] = max(starting_fitness - self.slopes[clone_id] * current_time, 0.01)
      
      # self.clones_array is the table listing id, label, fitness, generation, parent id. This line updates the fitness of the 
      # second clone (only, as the other is wild-type) based on the linear function.
      return super().get_next_generation(current_data) 
      # Returns the next generation of cells based on the current data and the updated fitness values.

np.random.seed(28)
n_clones = 100     # number of clones in the initial population
mean_fitness = 1.3 # mean of the initial fitness values for the clones
std_fitness = 0.1 # standard deviation of the fitness values 

mutant_fitnesses = np.random.normal(loc=mean_fitness, scale=std_fitness, size=n_clones - 1) # Draws samples from a normal distribution 
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
sim = WFLinearFitness(params)
# These exact numbers ensure we finish with a fitness of 1 for this starting fitness of 1.5

sim.run_sim()
sim.muller_plot(figsize=(5, 5))
plt.title("WF Competition Between One Wild-Type Clone (1 cell) and 99 Fitter Clones of Linearly Decreasing Fitness")
# Changed the above to include more clones in the initial population.
plt.xlabel("Time")
plt.ylabel("Clone Size")


import os
docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')
plt.savefig(os.path.join(docs_dir, 'muller_wf_linear100clones.png'), dpi=150, bbox_inches='tight')


sim.plot_surviving_clones_for_non_mutation()
plt.savefig(os.path.join(docs_dir, 'survivial_wf_linear100clones.png'), dpi=150, bbox_inches='tight')
# Adding in a plot of the surviving clones for the non-mutation case.
sim.plot_mean_clone_size_graph_for_non_mutation()
plt.savefig(os.path.join(docs_dir, 'meanclonesize_wf_linear100clones.png'), dpi=150, bbox_inches='tight')
# Adding in a plot of the mean clone size for the non-mutation case.
# To make this worthwhile, we will increase the number of clones in the initial population.


plt.show()
# print("Final fitness of the second clone:", sim.clones_array[1, sim.fitness_idx])
