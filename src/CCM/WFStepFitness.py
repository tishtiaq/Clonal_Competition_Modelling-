
# Creating custom class using a step function for fitness.
# The aim is that fitness will be at X until time t where it will drop down to neutral. 

# We want to study how well this fits the data, compared with linear and exponential fit.
# A theory is that maybe fitness does not decline per se, but instead the loss of TP53 causes
# a temporary increase in fitness that then drops to neutrality.


import matplotlib.pyplot as plt
import numpy as np
from clone_competition_simulation import (WF, Parameters, FitnessParameters, TimeParameters,
                                          PopulationParameters)

from clone_competition_simulation.simulation_algorithms.current_data import NonSpatialCurrentData
import clone_competition_simulation

class WFStepFitness(WF):
    # Defining a class based on the step function which is one of two values, dependent on the time t
    
    # a_intercept gives us the starting fitness, t_time tell us the time at which the fitness drops to neutrality.
    def __init__(self, parameters, a_intercept, t_time):
       self.intercept = a_intercept
       self.time = t_time
       super().__init__(parameters)
       for i, fitness in enumerate(parameters.fitness.initial_fitness_array):
           self.clones_array[i, self.fitness_idx] = fitness

    def get_next_generation(self, current_data:NonSpatialCurrentData) ->  np.ndarray[tuple[int], np.dtype[np.int_]]:
      # This function returns cell counts for the next generation. This is the only function we need to overwrite in the WF algorithm.
      current_time = self.i / self.division_rate
      if current_time < self.time:
         new_fitness = self.intercept
      else: 
         new_fitness = 1.0
      self.clones_array[1, self.fitness_idx] = new_fitness
      # self.clones_array is the table listing id, label, fitness, generation, parent id. This line updates the fitness of the 
      # second clone (only, as the other is wild-type) based on the linear function.
      return super().get_next_generation(current_data) 
      # Returns the next generation of cells based on the current data and the updated fitness values.

if __name__ == "__main__":

   params = Parameters(
      algorithm="WF", 
      times=TimeParameters(max_time=100, division_rate=1), 
      population=PopulationParameters(initial_size_array=np.concatenate([[9000], np.ones(1000, dtype=int)])), # Change ratio of wild-type : higher fitness HERE
      fitness=FitnessParameters(initial_fitness_array=np.concatenate([[1.0], np.ones(1000)*1.05])) # Change fitness of fitter clone HERE
   )

   import os
   docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')


   sim_step = WFStepFitness(params, a_intercept=1.05, t_time=25) # change time step and starting fitness here
   sim_step.run_sim()

   sim_constant = WF(params)
   sim_constant.run_sim()


   # Below is the Muller Plots side by side
   fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
   sim_step.muller_plot(ax=ax1) # Plot the linear simulation on one side
   ax1.set_title("Step Fitness")
   sim_constant.muller_plot(ax=ax2) # Plot the constant simulation on the other side
   ax2.set_title("Constant Fitness")
   fig1.suptitle("Muller Plots")
   plt.savefig(os.path.join(docs_dir, 't25_1.05_step_muller.png'), dpi=150, bbox_inches='tight')


   # Plot Mean Clone Sizes
   fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 5))
   sim_step.plot_mean_clone_size_graph_for_non_mutation(ax=ax3)
   ax3.set_title("Step Fitness")
   sim_constant.plot_mean_clone_size_graph_for_non_mutation(ax=ax4)
   ax4.set_title("Constant Fitness")
   fig2.suptitle("Mean Clone Size")
   plt.savefig(os.path.join(docs_dir, 't25_1.05_step_size.png'), dpi=150, bbox_inches='tight')


   # Plot Survival Rate on the same graph 
   fig, ax = plt.subplots(figsize=(8, 5))  # creates 1 figure with two graphs
   sim_constant.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Constant Fitness')
   sim_step.plot_surviving_clones_for_non_mutation(ax=ax, legend_label='Step Fitness') 
   ax.legend()
   plt.title("Survival Comparison: Constant vs Step Fitness")
   plt.tight_layout()
   plt.savefig(os.path.join(docs_dir, 't25_1.05_step_survival.png'), dpi=150, bbox_inches='tight')


   plt.show()


