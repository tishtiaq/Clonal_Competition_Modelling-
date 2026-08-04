
# Creating file in order to show the relation between fitness and time in the different models used


import numpy as np
import matplotlib.pyplot as plt
import os 
docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')

# Step Fitness
'''
Best values: fitness= 1.05, induction= 0.002450, time = 2.5 (days)
'''
def step_function(t):
    if t < 3:
        return 1.05
    else:
        return 1

# Linear Fitness
'''
Best values: fitness= 1.05, induction= , decay = 0.00135
'''
def linear_function(t):
    return max(1.05 - 0.00135*t, 1.0)

# Exponential Fitness
'''
Best values: fitness=1.05 , induction= , decay = 0.002250
'''
def exp_function(t):
    #return 1 + (1.05 - 1)*np.exp(-0.002250*t)
    return 1.05 * np.exp(-0.002250*t)


t = np.linspace(0, 364, 1000)

fitness_step = np.array([step_function(ti) for ti in t])
fitness_linear = np.array([linear_function(ti) for ti in t])
fitness_exponential = np.array([exp_function(ti) for ti in t])

plt.figure(figsize=(8, 5))
plt.plot(t, fitness_step, label='Step', c='blue')
plt.plot(t, fitness_linear, label='Linear', c='red')
plt.plot(t, fitness_exponential, label='Exponential', c='green')
plt.axhline(1, color='grey', linestyle='--', linewidth=0.5, label='Neutral (fitness=1)')
plt.xlabel('Time (days)')
plt.ylabel('Fitness')
plt.title('Fitness over time for different decay models')
plt.legend()
plt.tight_layout()

plt.savefig(os.path.join(docs_dir, '2.5_step_fitness_time_plots.png'), dpi=150, bbox_inches='tight')

plt.show()
