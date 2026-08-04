
# Calculating residuals code
import numpy as np
import matplotlib.pyplot as plt 
import os           
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  

docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')

from TimeVaryingProbabilisticFitting import (get_sim_means, load_data_tp53, run_sim, get_grid, ERROR_OBJECT, 
                                             GRID_SHAPE, CELLS, DIVISION_RATE, LOOP_LIMITS)
from clone_competition_simulation.parameters import Parameters, TimeParameters, PopulationParameters, FitnessParameters, LabelParameters
from clone_competition_simulation import WF2D
from WFStepFitness2D import WFStepFitness2D # New spatial custom class
from WF2DLinearFitness import WFLinearFitness2D
from WF2DExponentialFitness import WFExponentialFitness2D



DATA_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', '41467_2022_33945_MOESM5_ESM.xlsx')
TP53 = load_data_tp53(DATA_FILE)
docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')


observed = np.array([np.mean(TP53[t]) for t in TP53])
times_data = np.array([t for t in TP53])


models = {
    'Constant': {'fitness': 1.065, 'induction': 0.0155},
    'Step':     {'fitness': 1.05, 'induction': 0.00245, 'decay': 95.0, 'model': 'step'},
    'Linear':   {'fitness': 1.05, 'induction': 0.01405, 'decay': 0.00135, 'model': 'linear'},
    'Exponential': {'fitness': 1.05, 'induction': 0.01405, 'decay': 0.00225, 'model': 'exponential'},
}

colours = {'Constant': 'grey', 'Linear': 'red', 'Step': 'blue', 'Exponential': 'green'}

for name, params in models.items():
    means, _, _ = get_sim_means(params, 10, TP53)
    predicted = np.array(means[1:])
    residuals = observed - predicted
    rmse = np.sqrt(np.mean(residuals**2))
    print(f"RMSE {name}: {rmse:.4f}")
    plt.scatter(times_data, residuals, label=f'{name} (RMSE={rmse:.4f})', c=colours[name])


# Plot a scatter plot with days on the x-axis and residual (observed - predicted) on the y-axis
# Points above 0 mean the model UNDERpredicts, so the simulation says fewer mutants than the data actually shows
# Points below 0 mean the model OVERpredicts, so the simulation says more mutants than the data actually shows
# Hence points close to 0 are good fits

# The RMSE gives a single-value to compare models directly 

plt.axhline(0, color='grey', linestyle='--', linewidth=0.5)
plt.xlabel('Time (days)')
plt.ylabel('Residual (observed - predicted)')
plt.title('Residuals at each timepoint')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(docs_dir, 'residuals_all_models.png'), dpi=150, bbox_inches='tight')
plt.show()
