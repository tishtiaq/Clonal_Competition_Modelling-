


# This is Michael Hall's code on fitting the early takeover data
# Uses the first 3 data points


from pyabc import (ABCSMC, RV, Distribution, PNormDistance)
from pyabc.sampler import MulticoreEvalParallelSampler, SingleCoreSampler
from functools import partial, update_wrapper
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
from collections import OrderedDict
from clone_competition_simulation.parameters import Parameters
import pyabc.visualization.credible as credible
import sys

DATA_FILE = "09-03-21 Final clonal counting dataset.xlsx"

# Functions and fixed parameters for simulations
GRID_SHAPE = (500, 500)
DIVISION_RATE = 0.27
CELLS = GRID_SHAPE[0]*GRID_SHAPE[1]

ERROR_OBJECT = {'distance': 99999}
#NUM_CLONES = 100

# Each 500x500 grid is 0.25 million cells.
# If 1million cells in the oesophagus, then 4 grids is a full oesophagus.
# Asking for 100 clones in total for each time point.
# Run simulations repeatedly until enough clones exist for each time point.
# Don't want to get stuck in a loop if the induction/fitness is too low and no clones survive
# Limit the number of loops. If 4 grids per oeosophagus, then 50 grids is 12.5 mice worth, more than sampled number.
LOOP_LIMITS = 50


def get_grid(fitness, induction, grid_shape, cells):
    # Make the initial grid with randomly placed induced cells
    initial_grid = np.zeros(grid_shape, dtype=int)
    total_mutants = int(induction*cells)
    mutant_locs = np.random.choice(grid_shape[0]*grid_shape[1], total_mutants, replace=False)
    mutant_locs = [(m // grid_shape[1], m % grid_shape[1]) for m in mutant_locs]  # Convert the random draws into array indices

    count = 0
    for i in range(total_mutants):
        initial_grid[mutant_locs[count]] = i + 1
        count += 1

    fitness_array = [1] + [fitness]*total_mutants
    label_array = [0] + [1]*total_mutants
    return initial_grid, fitness_array, label_array

def get_mutant_takeover(sim):
    mutants = sim.population_array.toarray()[1:]
    mutant_pop = mutants.sum(axis=0)
    if mutant_pop.max() == sim.total_pop:
        mutant_pop[np.argmax(mutant_pop):] = sim.total_pop
    return mutant_pop/CELLS

# replace with r^2

def distance_ks(target, sim_results):
    return ks_2samp(target, sim_results).statistic

def distance_ir2(target, sim_results):
    #consider testing for shape- passing a 2d np array gives weird results
    residuals = target - sim_results
    sumSquareResiduals = np.sum(residuals ** 2)
    diffMean = target - np.mean(target)
    sumSquareDiff = np.sum(diffMean ** 2)
    return((sumSquareResiduals/sumSquareDiff))
# Defines two different distance functions 

def transformBySampling(takeover,samples):
    return([t for (t,individualSample) in zip(takeover,samples) for i in range(individualSample) ])

def run_sim(parameters, times, samplesPerTimepoint, target_data, return_takeover=False):
    #print(parameters)
    fitness, induction = parameters['fitness'], parameters['induction']
    #times = [t for t in target_data]
    #samplesPerTimepoint = [4,4,3,4,4,3]

    try:
        #setup simulation parameters
        initial_grid, fitness_array, label_array = get_grid(fitness, induction, GRID_SHAPE, CELLS)
        if len(fitness_array) == 1:  # Induction rate too low. No mutants on grid.
            return ERROR_OBJECT
        p = Parameters(algorithm='WF2D', initial_grid=initial_grid, times=times, fitness_array=fitness_array,
                       label_array=label_array,
                       print_warnings=False, division_rate=DIVISION_RATE,
                       cell_in_own_neighbourhood=True)
        
        #clone_counts = np.zeros(len(times), dtype=int)
        takeover = np.zeros(len(times), dtype=int)
        sim_distributions = [[] for i in range(len(times))]
        if p:
            for loop in range(LOOP_LIMITS):
                s = p.get_simulator()
                s.run_sim()
                takeover = get_mutant_takeover(s)

            if return_takeover:
                return transformBySampling(takeover,samplesPerTimepoint) #takeover
            
            #transform takeover to experimental timepoints
            result = transformBySampling(takeover,samplesPerTimepoint)
            
            #total_distance = sum([distance_ks(target_data[t], res[t]) for t in target_data])
            total_distance = distance_ir2(target_data, result)
            return {'distance': total_distance}
        else:
            return ERROR_OBJECT
    except (Exception, SystemExit) as e:
        print('Error')
        print(e)
        return ERROR_OBJECT
    
    # Set up ABC run
priors = Distribution(
        fitness=RV("uniform", 0, 50),
        induction=RV("uniform", 0, 0.1)
    )

# load the data

DATA_FILE = "41467_2022_33945_MOESM5_ESM.xlsx"
times = np.array([7*i for i in [1.5, 3, 6]])#, 12, 24, 52]])
dataset = pd.read_excel(DATA_FILE, sheet_name="Supplementary Data 5", skiprows=5, skipfooter=12,
                       usecols="E", header=None, engine='openpyxl').to_numpy()[:,0]

samplesPerTimepoint = [4,4,3]#,4,4,3]


distance = PNormDistance()
sampler = MulticoreEvalParallelSampler()  # Can change the number of cores used here

f = partial(run_sim, times=times, samplesPerTimepoint=samplesPerTimepoint, target_data=dataset)
update_wrapper(f, run_sim)  # Adds the __name__ attribute that pyabc uses

# Run the fitting using PyABC
abc = ABCSMC(f, priors, distance, population_size=100, sampler=sampler)
db_path = ("sqlite:///" + "TP53Early"+'_pyabc.db')

r = abc.new(db_path, {'distance': 0})
print("RunID:", r.id)
history = abc.run(minimum_epsilon=0.1, max_nr_populations=15)

# The database is stored and can be reloaded if needed - see pyabc documentation
# Here just print the median and 95% credible intervals.
def get_estimate_and_ci_for_param(param, df, w, confidence=0.95):
    vals = np.array(df[param])
    lb, ub = credible.compute_credible_interval(vals, w, confidence)
    median = credible.compute_quantile(vals, w, 0.5)
    return {'median': median, 'CI_lower_bound': lb, 'CI_upper_bound':ub}


df, w = history.get_distribution()  # Get the accepted parameters from the last generation

for p in ['fitness', 'induction']:
    print(p, get_estimate_and_ci_for_param(p, df, w))
