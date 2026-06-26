#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:21:44 2026

@author: syedtariqishtiaq
"""

# Introducing simple feedbacks


from clone_competition_simulation import (WF, NonSpatialCurrentData, Parameters, PopulationParameters, 
                                          TimeParameters, FitnessParameters)

import numpy as np
import pandas as pd

p = Parameters(
    algorithm='Moran',  # We will run a non-spatial Moran simulation. 
    times=TimeParameters(
        max_time=25,  # Run for 25 time units
        division_rate=1.4  # Set average division rate for all cells to 1.7 per time unit
    ), 
    population=PopulationParameters(  # Define the cell population
        initial_size_array=[100, 100, 100]  # There are three initial clones, with 100 cells in each
    ),
    fitness=FitnessParameters(  # Define the cell fitness
        initial_fitness_array=[1, 1.02, 1.04]    # Each clone has a different fitness value
    )
)

sim = p.get_simulator()
sim.run_sim()


# Taken from Introduction.md
# The WF model has just one function to overwrite: get_next_generation

class DecreasingFitness(WF):
    def get_next_generation(self, current_data: NonSpatialCurrentData) -> np.ndarray[tuple[int], np.dtype[np.int_]]:
     # -> indicates what kind of output you expect from the function  
        weights = current_data.current_population * self.clones_array[current_data.non_zero_clones, self.fitness_idx]
        relative_weights = weights / weights.sum()




# Below is the get_next_generation function that is in the actual WF definition file. 

    def get_next_generation(self, current_data: NonSpatialCurrentData) -> np.ndarray[tuple[int], np.dtype[np.int_]]: 
        """
        Sample from the current cells to output the next generation of cells

        Draws from each clone in proportion to the clone size and the clone fitness

        Args:
            current_data (CurrentData): contains the current clone cell populations and the indices of the living clones

        Returns:
            np.ndarray[tuple[int], np.dtype[np.int_]]: _description_
        """
        # Draw the new generation of cells from the old generation.
        # First, calculate the relative weight of each clone (population size multiplied by the fitness)
        weights = current_data.current_population * self.clones_array[current_data.non_zero_clones, self.fitness_idx]
        relative_weights = weights / weights.sum()

        # Then draw the new population.
        new_population = np.random.multinomial(self.total_pop, relative_weights)
        return new_population
        































