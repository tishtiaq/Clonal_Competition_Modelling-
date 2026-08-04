import numpy as np
from collections import OrderedDict

from clone_competition_simulation.parameters import Parameters, TimeParameters, PopulationParameters, LabelParameters, FitnessParameters

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from WFStepFitness import WFStepFitness
# Added import of Step Fitness custom class

# Silence the progress bar — must use the module reference directly
import clone_competition_simulation.simulation_algorithms.simulation_loop as sl
from rich.console import Console
sl.console = Console(quiet=True)

def rangeModifier(valueRange,steps):
    low, high = valueRange
    gradient = (high-low)/steps
    intercept = low+gradient/2
    return(intercept,gradient)

fitnessRange = (1, 2)
inductionRange = (0.001, 0.03)
decayRange = (50, 350)

steps_3d = 10


result_3d = np.load("./step_sim_result_3d_exponential.npy")
print("3D sweep completed and saved")

best_i, best_j, best_k = np.unravel_index(np.argmax(result_3d), result_3d.shape) # np.argmax finds the highest value in the entire cube. unravel converts it back to 3D coordinates
fi3, fg3 = rangeModifier(fitnessRange, steps_3d)
ii3, ig3 = rangeModifier(inductionRange, steps_3d)
di3, dg3 = rangeModifier(decayRange, steps_3d)
best_fitness = fg3*best_i + fi3
best_induction = ig3*best_j + ii3
best_decay = dg3*best_k + di3      # converts grid indices back into actual paramters
print(f"Best: fitness={best_fitness:.4f}, induction={best_induction:.6f}, decay={best_decay:.6f}")

