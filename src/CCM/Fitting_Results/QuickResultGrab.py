
import numpy as np
import matplotlib.pyplot as plt
import os 

docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'src', 'CCM')

fitnessRange = (1,2)
inductionRange = (0.001, 0.03)
decayRange = (0,50)
steps_3d = 10

def rangeModifier(valueRange,steps):
    low, high = valueRange
    gradient = (high-low)/steps
    intercept = low+gradient/2
    return(intercept,gradient)

result_3d = np.load(os.path.join(docs_dir, '0-50_step_sim_result_3d.npy'))
print("3D sweep completed and saved")


best_i, best_j, best_k = np.unravel_index(np.argmax(result_3d), result_3d.shape) # np.argmax finds the highest value in the entire cube. unravel converts it back to 3D coordinates
fi3, fg3 = rangeModifier(fitnessRange, steps_3d)
ii3, ig3 = rangeModifier(inductionRange, steps_3d)
di3, dg3 = rangeModifier(decayRange, steps_3d)
best_fitness = fg3*best_i + fi3
best_induction = ig3*best_j + ii3
best_decay = dg3*best_k + di3      # converts grid indices back into actual paramters
print(f"Best: fitness={best_fitness:.4f}, induction={best_induction:.6f}, decay={best_decay:.6f}")