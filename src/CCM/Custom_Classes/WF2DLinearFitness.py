
# Building custom 2D algorithms

# This file contains a custom implementation of a WF2D algorithm with a linearly decreasing fitness function. 

from clone_competition_simulation  import (WF2D, Parameters, SpatialCurrentData, TimeParameters,
                                                             PopulationParameters, FitnessParameters)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  

class WFLinearFitness2D(WF2D):
    def __init__(self, parameters, a_intercept):
        b_slope = (a_intercept - 1) / parameters.times.max_time
        self.slope = b_slope
        self.intercept = a_intercept
        # Defines the gradient and intercept
        super().__init__(parameters)
        for i, fitness in enumerate(parameters.fitness.initial_fitness_array):
            self.clones_array[i, self.fitness_idx] = fitness
    # Same linear function defined above. 

    def get_next_generation(self, current_data: SpatialCurrentData) -> np.ndarray[tuple[int], np.dtype[np.int_]]:
        current_time = self.i / self.division_rate
        new_fitness = self.intercept - self.slope * current_time   # Here we define the operation
        self.clones_array[1, self.fitness_idx] = new_fitness
        return super().get_next_generation(current_data) 
if __name__ == "__main__":    
    initial_grid = np.zeros((50,50), dtype=int) 
    initial_grid[:37, :] = 1      # 2D slicing syntax: first colon is rows, second is columns

    params = Parameters(
        algorithm="WF2D", times= TimeParameters(max_time=100, division_rate=1, samples=4),
        population= PopulationParameters(initial_grid=initial_grid, cell_in_own_neighbourhood=False),
        fitness= FitnessParameters(initial_fitness_array=np.array([1, 1.5])),
    )
            

    sim = WFLinearFitness2D(params, a_intercept=1.5)
    sim.run_sim()
    sim.muller_plot(figsize=(5, 5))

    plt.title("Competition Between One Wild-Type Clone and One Fitter Clone of Linearly Decreasing Fitness")
    plt.xlabel("Time")
    plt.ylabel("Clone Size")

    import os
    docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')
    plt.savefig(os.path.join(docs_dir, 'f_wf_2D_linear_75_25.png'), dpi=150, bbox_inches='tight')

    plt.show()








