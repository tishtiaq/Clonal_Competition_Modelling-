# Creating spatial Step class 

from clone_competition_simulation  import (WF2D, Parameters, SpatialCurrentData, TimeParameters,
                                                             PopulationParameters, FitnessParameters)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  

class WFStepFitness2D(WF2D):
    def __init__(self, parameters, a_intercept, t_time):
        self.intercept = a_intercept
        self.time = t_time
        super().__init__(parameters)
        for i, fitness in enumerate(parameters.fitness.initial_fitness_array):
            self.clones_array[i, self.fitness_idx] = fitness

    def get_next_generation(self, current_data):
        current_time = self.i / self.division_rate
        new_fitness = self.intercept if current_time < self.time else 1.0
        self.clones_array[1, self.fitness_idx] = new_fitness
        return super().get_next_generation(current_data)