

# Moran2DLinearFitness

from clone_competition_simulation import (Moran2D, Parameters, SpatialCurrentData, TimeParameters,
                                                             PopulationParameters, FitnessParameters)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

class MoranLinearFitness2D(Moran2D):
    def __init__ (self, parameters, a_intercept):
        b_slope = (a_intercept - 1) / 100
        self.slope = b_slope
        self.intercept = a_intercept
        # Defines the gradient and intercept
        super().__init__(parameters)


    def get_differentiating_cell(self, i: int, current_data: SpatialCurrentData) -> int:
        return 
    # TODO: Write in custom rules for these functions

    def get_dividing_cell(self, coord: int, current_data: SpatialCurrentData) -> int:
        return 
    # TODO: Write in custom rules for these functions
    
initial_grid = np.zeros((50,50), dtype=int) 
initial_grid[:37, :] = 1                    # 2D slicing syntax: first colon is rows, second is columns. Defines what proportion
                                            # of the grid is wild-type and what proportion is mutant. 

params = Parameters(
    algorithm="Moran2D", times= TimeParameters(max_time=100, division_rate=1, samples=4),
    population= PopulationParameters(initial_grid=initial_grid, cell_in_own_neighbourhood=False),
    fitness= FitnessParameters(initial_fitness_array=np.array([1, 1.5])),
)

sim = MoranLinearFitness2D(params, a_intercept=1.5)
sim.run_sim()
sim.muller_plot(figsize=(5, 5))




