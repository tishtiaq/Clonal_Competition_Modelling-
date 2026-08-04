
# 
# Previously, we have seen that the simulation in the file ProbabilisticFittingResultsTP53 fits the early timepoints well, 
# but not the later timepoints.
# In this file, I have taken the exact code from ProbabilisticFittingResultsTP53 and added in a search for an extra parameter: decay.
# By adding this parameter we are now searching 3-dimensional space instead of 2D, and we hope that this leads to a better fit. 

# This code currently uses the custom class with Exponentially Decreasing Fitness. 


import numpy as np
import pandas as pd
from collections import OrderedDict
from scipy.stats import ks_2samp
from scipy.stats import anderson_ksamp
from scipy.stats import multivariate_normal

from clone_competition_simulation.parameters import Parameters, TimeParameters, PopulationParameters, LabelParameters, FitnessParameters

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from WF2DExponentialFitness import WFExponentialFitness2D

# Silence the progress bar — must use the module reference directly
import clone_competition_simulation.simulation_algorithms.simulation_loop as sl
from rich.console import Console
sl.console = Console(quiet=True)


def _silence_progress():
    import clone_competition_simulation.simulation_algorithms.simulation_loop as sl
    from rich.console import Console
    from loguru import logger

    sl.console = Console(quiet=True)
    logger.disable("clone_competition_simulation")

# Loading data for all 6 timepoints 
docs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'docs')

DATA_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', '41467_2022_33945_MOESM5_ESM.xlsx')
# Changed to data file reaching TP53 mouse data
def load_data_tp53(data_file):
    df = pd.read_excel(data_file, sheet_name="Supplementary Data 5", skiprows=5, skipfooter=1,
                           usecols="A,E", header=None, engine='openpyxl') 
    df.columns = ['week', 'percent_gfp']
    df['week'] = df['week'].ffill()
    result = OrderedDict()
    for week in sorted(df['week'].unique()): # unique removes duplicate values
        days = int(week*7) 
        vals = df[df['week']==week]['percent_gfp'].values   # values returns data inside the dataframe without columns, headers or row indices
        result[days] = vals           
    return result                    # Hands back a dictionary with days passed as keys and a list of the takeover for each mouse in the sample as values.

def mean_data_tp53(data_file):
    data = load_data_tp53(data_file)
    return OrderedDict([(t, np.mean(v)) for t, v in data.items()])

TP53_DATA = load_data_tp53(DATA_FILE)    # Creates the dictionary defined in the function above (load_data_tp53)
TP53_MEAN = mean_data_tp53(DATA_FILE)    # Createes the same dictionary but for the values it has a single mean of the gfp% instead of a list
print(TP53_MEAN)                         # Return the ordered dictionary containing mean %gfp for each timepoint


# Functions and fixed parameters for simulations
GRID_SHAPE = (500, 500)
DIVISION_RATE = 0.27
CELLS = GRID_SHAPE[0]*GRID_SHAPE[1]

ERROR_OBJECT = {'distance': -np.inf}
NUM_CLONES = 100    

LOOP_LIMITS = 10
SIMS_PER_PARAMETER = 100
N_BOOTSTRAP = 1000


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


def get_sim_mean_pl(parameters, target_data):
    fitness = parameters['fitness']
    induction = parameters['induction']
    decay = parameters.get('decay', 0) # defaults to 0 if not provided
    times = [t for t in target_data]
    # Added decay to the above here

    try:
        initial_grid, fitness_array, label_array = get_grid(fitness, induction, GRID_SHAPE, CELLS)
        if len(fitness_array) == 1:  # Induction rate too low. No mutants on grid.
            return ERROR_OBJECT
        
    # The above function computes mean simulated takeover which gets used in the calculation of likelihood
        p = Parameters(algorithm='WF2D', 
                       population=PopulationParameters(initial_grid=initial_grid,cell_in_own_neighbourhood=True),
                       times=TimeParameters(times=times, division_rate=DIVISION_RATE), 
                       labels=LabelParameters(initial_label_array=label_array),
                       fitness=FitnessParameters(initial_fitness_array=fitness_array)
                       )

        # Just redefining get_mutant_takeover to use later
        def get_mutant_takeover(sim):
            mutants = sim.population_array.toarray()[1:]
            mutant_pop = mutants.sum(axis=0)
            if mutant_pop.max() == sim.total_pop:
                mutant_pop[np.argmax(mutant_pop):] = sim.total_pop
            return mutant_pop/CELLS
        # Above stays the same
        '''
        sim_results = []
        for loop in range(LOOP_LIMITS):
            s = p.get_simulator()
            s.run_sim()
            takeover = get_mutant_takeover(s)
            sim_results.append(takeover)

            Replaced this code with use of LinearFitness2D custom class
        '''
        sim_results = []
        for loop in range(LOOP_LIMITS):
            if decay>0 :
                s = WFExponentialFitness2D(p, a_coefficient=fitness) 
                s.rate = decay
            else:
                s = p.get_simulator() # Use normal WF2D for when decay=0
            s.run_sim()
            sim_results.append(get_mutant_takeover(s))

        sim_array = np.array(sim_results)
        mean_vec = np.mean(sim_array, axis=0)
        cov_mat = np.cov(sim_array, rowvar=False)    # calculates the covariance matrix for the variables. rowvar=False tells you it is the columns 
                                                     # (not rows) that you want to analyse
        cov_mat += np.eye(len(times)) * 1e-10        # creates an identity matrix with tiny values across the diagonal and we add it to the 
                                                     # covariance matrix to stop it from becoming singular


        observations = np.array([np.mean(target_data[t]) for t in times])
        residuals = observations - mean_vec
        sum_square_residuals = np.sum(residuals ** 2)
        diff_mean = observations - np.mean(observations) # computes how far each observation differs from the mean of all observations
        sum_square_diff = np.sum(diff_mean ** 2)
        total_distance = -(sum_square_residuals / sum_square_diff)
        # Above is what Claude wrote to try and get around the problem of the zero-width confidence intervals: to try and make 
        # it into a smooth bell curve shape instead 
    
        return {'distance': total_distance}
        # Above all stays the same 

    except (Exception, SystemExit) as e:
        print('Error:', e)
        return ERROR_OBJECT         
    # Above all stays the same

'''
total_distance = multivariate_normal.logpdf(observations, mean=mean_vec, cov=cov_mat)  # creates the probability density function
# Put this instead of total_distance if you want to revert back to the old formula
'''   

test_p = {'fitness':1.5,'induction':0.05}
# target_data = mean_data('hom', DATA_FILE)
get_sim_mean_pl(test_p, TP53_MEAN)



from joblib import Parallel, delayed

def rangeModifier(valueRange,steps):
    low, high = valueRange
    gradient = (high-low)/steps
    intercept = low+gradient/2
    return(intercept,gradient)


def pseudoLikelihoodSweep(target_data, n=10,fitnessRange=(0,10),inductionRange=(0.005,0.105)):
    #values = [[0 for i in range(n)] for j in range(n)]
    #result = np.zeros((n,n))
    fi, fg = rangeModifier(fitnessRange,n)
    ii, ig = rangeModifier(inductionRange,n)
    # fi=fitness intercept, fg=fitness gradient, ii=induction intercept, ig=induction gradient 
    # These let you convert any grid index i into the actual parameter value it represents

    # This function divides fitnessRange and inductionRange into a grid of nxn points, runs a systematic grid search, and calculates 
    # pseudo-likelihood score for each combination.

    grid = [(i, j, fg*i+fi, ig*j+ii) for i in range(n) for j in range(n)]

    def _worker(i, j, fit, ind):
        p = {'fitness': fit, 'induction': ind}
        r = get_sim_mean_pl(p, target_data)
        if r == ERROR_OBJECT:
            return i, j, np.nan
        return i, j, r['distance']  # already a log-likelihood, don't log again

    results = Parallel(n_jobs=-1, initializer=_silence_progress)(delayed(_worker)(i, j, fit, ind) for i, j, fit, ind in grid)

    result = np.zeros((n, n))
    values = [[None]*n for _ in range(n)]
    for i, j, val in results:
        result[i][j] = val
        values[i][j] = (fg*i+fi, ig*j+ii)

    return(result,values)


# Creating new 3D grid searching function
def pseudoLikelihoodSweep3D(target_data, n=10, fitnessRange=(1, 2), inductionRange=(0.001, 0.03), decayRange=(0, 0.005)):
    fi, fg = rangeModifier(fitnessRange, n)
    ii, ig = rangeModifier(inductionRange, n)
    di, dg = rangeModifier(decayRange, n)

    grid = [(i, j, k, fg*i + fi, ig*j + ii, dg*k + di)
        for i in range(n) for j in range(n) for k in range(n)] # builds a list of combination of grid positions
    
    def _worker(i, j, k, fit, ind, dec): # This function takes one specific parameter combination and runs the simulation with it
        p = {'fitness': fit, 'induction': ind, 'decay': dec}
        r = get_sim_mean_pl(p, target_data)
        if r == ERROR_OBJECT:
            return i, j, k, np.nan
        return i, j, k, r['distance']
    # Either returns (i, j, k, distance) if sim was successful or (i, j, k, nan) if unsuccessful

    results = Parallel(n_jobs=-1, initializer=_silence_progress)(  # n_jobs=-1 means use all CPU cores simulataneously. _silence_progress supresses the simulation progress bars so the terminal is not flooded with them
        delayed(_worker)(i, j, k, fit, ind, dec) # delayed(_worker) packages up jobs and sends them to workers on cores which are free
        for i, j, k, fit, ind, dec in grid) # loops through tuples in grid, creating 1 packaged job per combination
    # The above distributes parameter combinations and distributes them across CPU cores to run simulataneously

    result = np.zeros((n, n, n)) # creates an empty 3D array for where likelihood values will be stored
    for i, j, k, val in results: # loops through the 1000 returned results and places each distance value into the right position in the cube
        result[i][j][k] = val # result[i][j][k] is the cell at fitness i, induction j, decay k. After the loop the cube is fully populated

    return result   # hands back the likelihood cube


fitnessRange = (1, 2)
inductionRange = (0.001, 0.03)
steps = 20

result, values = pseudoLikelihoodSweep(TP53_MEAN,fitnessRange=fitnessRange,inductionRange=inductionRange,n=steps)
# Now the ranges are specified for the ranges of TP53


import matplotlib.pyplot as plt

def fmt(x, sf=2):
    return f"{x:.{sf}g}"

def likelihood_heatmap(data,fitnessRange,inductionRange,steps=10,expTransform=False, filename=None):
    #if plotting likelihood, convert from log 
    if expTransform:
        data = np.exp(data- np.nanmax(data))
        data = np.where(np.isnan(data), 0, data)
        data /= data.sum()

    fig, ax = plt.subplots()
    im = ax.imshow((data), cmap='plasma')

    fi, fg = rangeModifier(fitnessRange,steps)
    ii, ig = rangeModifier(inductionRange,steps)

    #float(i)+0.5,'induction':float(1+j)/100
    ax.set_xticks(range(steps), [fmt(ig*float(i)+ii) for i in range(steps)])
    ax.set_yticks(range(steps), [fmt(fg*float(i)+fi) for i in range(steps)])

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    cbar = fig.colorbar(im, ax=ax)
    if filename:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.show()

'''
likelihood_heatmap(result,fitnessRange,inductionRange,steps=steps, filename=os.path.join(docs_dir, 'tp53_heatmap_loglik.png')) 
# Plots log-likelihoods


likelihood_heatmap(result,fitnessRange,inductionRange,steps=steps,expTransform=True, filename=os.path.join(docs_dir, 'tp53_heatmap_prob.png')) 
# Exponentiates the log-likelihoods to give probabilities which allows us to calculate confidence intervals
'''

def getBestEdges(result,fitnessRange=(0,10),inductionRange=(0.005,0.105),n=10):
    best_x, best_y = np.unravel_index(np.argmax(result), result.shape)
    fi, fg = rangeModifier(fitnessRange,n)
    lowFit = fitnessRange[0]
    lowInd = inductionRange[0]
    ii, ig = rangeModifier(inductionRange,n)
    edgeX = (lowFit+fg*best_x, lowFit+fg*(best_x+1))
    edgeY = (lowInd+ig*best_y, lowInd+ig*(best_y+1))
    print(best_x,best_y)
    return(edgeX,edgeY)

# Finds the peak of the likelihood map (best parameter combination of fitness and induction)

'''
top_coarse_fitness, top_coarse_induction = getBestEdges(result,fitnessRange=fitnessRange, 
                                                        inductionRange=inductionRange,
                                                        n=steps)
# Changed the inputs to this function to include TP53 fitness and induction ranges


data = np.exp(result- np.nanmax(result))
data = np.where(np.isnan(data), 0, data)
data /= data.sum()

prob_surface = np.exp(grid - np.nanmax(grid))
prob_surface = np.where(np.isnan(prob_surface), 0, prob_surface)
prob_surface /= prob_surface.sum()
'''


def get_95ci(data,fitnessRange,inductionRange,steps):
    def credible_interval(vals, probs, ci=0.95):
        """Central credible interval from a discrete 1D distribution."""
        cdf   = np.cumsum(probs)
        lower = vals[np.searchsorted(cdf, (1 - ci) / 2)]
        upper = vals[np.searchsorted(cdf, (1 + ci) / 2)]
        return lower, upper

   # Computes credible intervals by treating exponentiated likelihood as a probability distribution and integrating it

    data = np.exp(data- np.nanmax(data))
    data = np.where(np.isnan(data), 0, data)
    data /= data.sum()

    fi, fg = rangeModifier(fitnessRange,steps)
    ii, ig = rangeModifier(inductionRange,steps)

    x_vals = [(ig*float(i)+ii) for i in range(steps)]
    y_vals = [(fg*float(i)+fi) for i in range(steps)]

    marginal_x = data.sum(axis=0)   # sum over rows → distribution over x
    marginal_y = data.sum(axis=1)   # sum over cols → distribution over y

    x_lo, x_hi = credible_interval(x_vals, marginal_x)
    y_lo, y_hi = credible_interval(y_vals, marginal_y)

    print(f"x 95% CI: [{x_lo:.2f}, {x_hi:.2f}]")
    print(f"y 95% CI: [{y_lo:.2f}, {y_hi:.2f}]")

    #bounds is a dict of tuples with a new inductionRange 
    bounds = {"inductionRange":(x_lo-ig,x_hi+ig),"fitnessRange":(y_lo-fg,y_hi+fg)}

    return({"data":data,
            "x":{"low":x_lo,"high":x_hi},
            "y":{"low":y_lo,"high":y_hi},
            "bounds":bounds,
            "x_vals":x_vals,
            "y_vals":y_vals,
            "marginal_x":marginal_x,
            "marginal_y":marginal_y})


ci_tp53 = get_95ci(result,fitnessRange,inductionRange,steps)   # changed name from ci_hom to ci_tp53


import matplotlib.pyplot as plt

def plot_marginals(ci_result, filename=None):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    x_vals = np.array(ci_result["x_vals"])
    y_vals = np.array(ci_result["y_vals"])

    marginal_x = ci_result["marginal_x"]
    marginal_y = ci_result["marginal_y"]

    x_lo = ci_result["x"]["low"]
    x_hi = ci_result["x"]["high"]
    y_lo = ci_result["y"]["low"]
    y_hi = ci_result["y"]["high"]

    # heatmap
    axes[0].imshow(ci_result["data"], origin='lower', extent=[min(x_vals), max(x_vals), min(y_vals), max(y_vals)], cmap='inferno', aspect='auto')
    axes[0].axvline(x_lo, color='w', linestyle='--', label='x CI')
    axes[0].axvline(x_hi, color='w', linestyle='--')
    axes[0].axhline(y_lo, color='c', linestyle='--', label='y CI')
    axes[0].axhline(y_hi, color='c', linestyle='--')
    axes[0].set_title('Joint distribution')
    axes[0].legend(fontsize=8)

    # marginal x
    axes[1].plot(x_vals, marginal_x)
    axes[1].axvline(x_lo, color='r', linestyle='--', label=f'95% CI')
    axes[1].axvline(x_hi, color='r', linestyle='--')
    axes[1].set_title('Marginal x')
    axes[1].legend()

    # marginal y
    axes[2].plot(y_vals, marginal_y)
    axes[2].axvline(y_lo, color='r', linestyle='--', label='95% CI')
    axes[2].axvline(y_hi, color='r', linestyle='--')
    axes[2].set_title('Marginal y')
    axes[2].legend()

    plt.tight_layout()
    if filename:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.show()


plot_marginals(ci_tp53, filename=os.path.join(docs_dir, 'tp53_marginals_coarse.png'))   # Saving plot to docs


# # Look in the area of the peak specifically for a better distribution

result_tp53_tight, values_tp53_tight = pseudoLikelihoodSweep(TP53_MEAN,fitnessRange=ci_tp53["bounds"]["fitnessRange"],inductionRange=ci_tp53["bounds"]["inductionRange"],n=steps)


likelihood_heatmap(result_tp53_tight,ci_tp53["bounds"]["fitnessRange"],ci_tp53["bounds"]["inductionRange"],steps=steps,
                   filename=os.path.join(docs_dir, 'tp53_heatmap_tight_loglik.png'))
likelihood_heatmap(result_tp53_tight,ci_tp53["bounds"]["fitnessRange"],ci_tp53["bounds"]["inductionRange"],steps=steps,expTransform=True,
                   filename=os.path.join(docs_dir, 'tp53_heatmap_tight_prob.png'))
# Also saving these plots to docs, these are the heatmaps where we have searched a finer area

ci_tp53_tight = get_95ci(result_tp53_tight,ci_tp53["bounds"]["fitnessRange"],ci_tp53["bounds"]["inductionRange"],steps)

plot_marginals(ci_tp53_tight, filename=os.path.join(docs_dir, 'tp53_marginals_tight.png'))

import pickle 
np.save('result_tp53_tight.npy', result_tp53_tight)
with open('ci_tp53_tight.pkl', 'wb') as f:
    pickle.dump(ci_tp53_tight, f)

# # Make some plots (do the best fits fit?)
from scipy.stats import sem

def run_sim(parameters, target_data, return_clone_sizes=False):
    fitness, induction = parameters['fitness'], parameters['induction']
    decay = parameters.get('decay', 0)
    times = [t for t in target_data]

    try:
        initial_grid, fitness_array, label_array = get_grid(fitness, induction, GRID_SHAPE, CELLS)
        if len(fitness_array) == 1:  # Induction rate too low. No mutants on grid.
            return ERROR_OBJECT

        p = Parameters(algorithm='WF2D', 
               population=PopulationParameters(initial_grid=initial_grid,cell_in_own_neighbourhood=True),
               times=TimeParameters(times=times, division_rate=DIVISION_RATE), 
               labels=LabelParameters(initial_label_array=label_array),
               fitness=FitnessParameters(initial_fitness_array=fitness_array)
               )

        def get_mutant_takeover(sim):         # Rewriting the same get_mutant_takeover function that we used earlier
            mutants = sim.population_array.toarray()[1:]
            mutant_pop = mutants.sum(axis=0)
            if mutant_pop.max() == sim.total_pop:
                mutant_pop[np.argmax(mutant_pop):] = sim.total_pop
            return mutant_pop / CELLS
        
        sim_results = []
        for loop in range(LOOP_LIMITS):
            if decay > 0:
                s = WFExponentialFitness2D(p, a_coefficient=fitness)
                s.rate = parameters['decay']
            else:
                s = p.get_simulator()
            s.run_sim()
            takeover = get_mutant_takeover(s)  # creates an array of takeover fractions, one per timepoint
            sim_results.append(takeover)
        # Everything that is above is what we have done previously in the simulations

        mean_takeover = np.mean(sim_results, axis=0)

        if return_clone_sizes:
            return OrderedDict([(t, np.array([mean_takeover[i]])) for i, t in enumerate(times)])
        
        total_distance = sum([abs(np.mean(target_data[t]) - mean_takeover[i])
                                  for i, t in enumerate(times)])
        return {'distance': total_distance}
    
    except (Exception, SystemExit) as e:
        print('Error:', e)
        return ERROR_OBJECT 


def get_sim_means(parameters, num, data, max_attempts_per_sample=10):
    res = [[] for t in data]

    i = 0
    attempts = 0
    successes = 0
    while successes < num:
        print(successes, end=' ')
        np.random.seed(i)
        clone_sizes = run_sim(parameters=parameters, target_data=data, return_clone_sizes=True)
        i += 1
        attempts += 1

        # Detect failure: ERROR_OBJECT won't have one entry per timepoint in `data`
        if clone_sizes is ERROR_OBJECT or set(clone_sizes.keys()) != set(data.keys()):
            if attempts - successes > max_attempts_per_sample * num:
                raise RuntimeError("Too many failed simulations — check parameters.")
            continue

        for j, (t, clones) in enumerate(clone_sizes.items()):
            res[j].append(clones.mean())
        successes += 1
    ind = parameters['induction']
    intervals_high = [ind]
    intervals_low = [ind]
    means = [ind]         # Changed all 3 of these as they were 1 before: made sense for hom, het and wt but now we are 
                          # talking about induction fraction, so need to reduce it. I have changed it so it starts at whatever
                          # is coded into TP53 data set
    for r in res:
        intervals_high.append(np.quantile(r, 0.975))
        intervals_low.append(np.quantile(r, 0.025))
        means.append(np.mean(r))

    return means, intervals_high, intervals_low

def plot_data(results, colour=None, label=None, err_stat=sem, elinewidth=2, capsize=4, markersize=3):
    x = results.keys()
    y = [results[k].mean() for k in results]
    yerr = [err_stat(results[k]) for k in results]
    plt.errorbar(x, y, yerr=yerr, label=label, c=colour, 
                 elinewidth=elinewidth, 
                 capsize=capsize, capthick=elinewidth,
                 markersize=markersize,
                 fmt='o')

def ci_to_params(ci_result):
    #returns the value half way between CIs

    #use bounds as names are unambigous
    bounds = ci_result['bounds']

    fit = np.mean(bounds["fitnessRange"])
    ind = np.mean(bounds["inductionRange"])
    return {'fitness': fit, 'induction': ind}


# Load the data
TP53 = load_data_tp53(DATA_FILE)
TP53_PARAMS = ci_to_params(ci_tp53_tight)

tp53_means, tp53_interval_high, tp53_interval_low = get_sim_means(TP53_PARAMS, 10, TP53) # changed here to 10 rather than 100 for cost

plt.figure(figsize=(6.5, 5))
data = TP53
colour_sim = 'k'
colour_data = 'k'
times = [0] + [t for t in data]
plt.plot(times, tp53_means, label='+/+ fit', c=colour_sim)
plt.fill_between(times, tp53_interval_high, tp53_interval_low, alpha=0.3, color=colour_sim)
plot_data(data, label='+/+ data', colour=colour_data)

plt.xlim(left=0)
plt.legend(bbox_to_anchor=(1.03, 1))
plt.yscale('linear')
plt.ylabel('Mutant Takeover Fraction')
plt.xlabel('Time (days)')
plt.ylim(bottom=0, top=0.3)
plt.tight_layout()
plt.savefig(os.path.join(docs_dir, 'tp53_best_fit.png'), dpi=150, bbox_inches='tight')
plt.show()

# Now we add in the 3 dimensional grid search:

decayRange = (0, 0.003)
steps_3d = 10

result_3d = pseudoLikelihoodSweep3D(TP53_MEAN,
                                    fitnessRange=fitnessRange,
                                    inductionRange=inductionRange,
                                    decayRange=decayRange,
                                    n=steps_3d)
# Calls the new pseudoLikelihoodSweep3D function

np.save('result_3d_exponential.npy', result_3d) # saves the result as soon as the sweep finishes
print("3D sweep completed and saved")

best_i, best_j, best_k = np.unravel_index(np.argmax(result_3d), result_3d.shape) # np.argmax finds the highest value in the entire cube. unravel converts it back to 3D coordinates
fi3, fg3 = rangeModifier(fitnessRange, steps_3d)
ii3, ig3 = rangeModifier(inductionRange, steps_3d)
di3, dg3 = rangeModifier(decayRange, steps_3d)
best_fitness = fg3*best_i + fi3
best_induction = ig3*best_j + ii3
best_decay = dg3*best_k + di3      # converts grid indices back into actual paramters
print(f"Best: fitness={best_fitness:.4f}, induction={best_induction:.6f}, decay={best_decay:.6f}")

# Then we need to take slices of 3 combinations of 2 parameters each time. So we take slices of fitnessXinduction, fitnessXdecay,
# inductionXdecay to find the value of each at the best of the other. It is like in 2D when we sum up inductions across fitness values
# to work out maximum induction. 



# Fitness X Induction:

from mpl_toolkits.mplot3d import Axes3D

fitness_vals = np.array([fg3*i + fi3 for i in range(steps_3d)])
induction_vals = np.array([ig3*j + ii3 for j in range(steps_3d)])
decay_vals = np.array([dg3*k + di3 for k in range(steps_3d)])

slice_fit_ind = result_3d[:, :, best_k]
X, Y = np.meshgrid(induction_vals, fitness_vals)   # meshgrid takes 1d coordinate arrays and turns them into a 2D grid of coordinates 
fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, slice_fit_ind, cmap='plasma', alpha=0.9)
ax.set_xlabel('Induction')
ax.set_ylabel('Fitness')
ax.set_zlabel('Log-likelihood')
ax.set_title(f'Likelihood surface: fitness x induction\n(decay fixed at {best_decay})')
fig.colorbar(surf, shrink=0.5)
plt.savefig(os.path.join(docs_dir, 'exp_tp53_3d_surface_fit_ind.png'), dpi=150, bbox_inches='tight')
plt.show()


# Fitness X Decay

slice_fit_dec = result_3d[:, best_j, :]
X, Z = np.meshgrid(decay_vals, fitness_vals)
fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Z, slice_fit_dec, cmap='plasma', alpha=0.9)
ax.set_xlabel('Decay')
ax.set_ylabel('Fitness')
ax.set_zlabel('Log-likelihood')
ax.set_title(f'Likelihood surface: fitness x decay\n(induction fixed at {best_induction})')
fig.colorbar(surf, shrink=0.5)
plt.savefig(os.path.join(docs_dir, 'exp_tp53_3d_surface_fit_dec.png'), dpi=150, bbox_inches='tight')
plt.show()

# Induction X Decay

slice_ind_dec = result_3d[best_i, :, :]
Y, Z = np.meshgrid(decay_vals, induction_vals) 
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(Y, Z, slice_ind_dec, cmap='plasma', alpha=0.9)
ax.set_xlabel('Induction')
ax.set_ylabel('Decay')
ax.set_zlabel('Log-likelihood')
ax.set_title(f'Likelihood surface: induction x decay\n(fitness fixed at {best_fitness})')
fig.colorbar(surf, shrink=0.5)
plt.savefig(os.path.join(docs_dir, 'exp_tp53_3d_surface_ind_dec.png'), dpi=150, bbox_inches='tight')
plt.show()

TP53 = load_data_tp53(DATA_FILE)
PARAMS_3D = {'fitness': best_fitness, 'induction': best_induction, 'decay': best_decay}
tp53_means_3d, tp53_high_3d, tp53_low_3d = get_sim_means(PARAMS_3D, 10, TP53)

plt.figure(figsize=(6.5, 5))
times = [0] + [t for t in TP53]
plt.plot(times, tp53_means_3d, label='Time-varying fit', c='blue')
plt.fill_between(times, tp53_high_3d, tp53_low_3d, alpha=0.3, color='blue')
plot_data(TP53, label='TP53 data', colour='black')
plt.xlim(left=0) # sets starting limit of x-axis
plt.legend()
plt.ylabel('Mutant Takeover Fraction')
plt.xlabel('Time (days)')
plt.ylim(bottom=0, top=0.3)
plt.title(f'Best fit: fitness={best_fitness}, decay={best_decay}')
plt.savefig(os.path.join(docs_dir, 'exp_tp53_best_fit_3d.png'), dpi=150, bbox_inches='tight')


# Plotting both best fit lines on the same graph 
plt.figure(figsize=(8, 5))
plt.plot(times, tp53_means, label='Without fitness feedbacks', c='blue')
plt.plot(times, tp53_means_3d, label='With fitness feedbacks', c='red') 
plot_data(TP53, label='TP53 data', colour='black')
plt.xlim(left=0)
plt.ylabel('Mutant Takeover Fraction')
plt.xlabel('Time (days)')
plt.legend()
plt.title('Best fit lines plotted with and without decreasing fitness')
plt.savefig(os.path.join(docs_dir, 'exp_best_fits_with_and_without_feedbacks.png'), dpi=150, bbox_inches='tight')
plt.show()


# Want to plot marginals as well:
prob_3d =  np.exp(result_3d - np.nanmax(result_3d)) # finds largest least -ve value and subtracts so we get a value of 0 for the 
                                                    # best probability (ensures probability of 1 for the best value since exp(0)=1)
prob_3d /= prob_3d.sum()                      # normalises so that all values sum to 1
marginal_fitness = prob_3d.sum(axis=(1,2)) # sums over induction and decay
marginal_induction = prob_3d.sum(axis=(0,2))  # sums over fitness and decay
marginal_decay = prob_3d.sum(axis=(0,1))      # sums over fitness and induction

fitness_vals = np.array([fg3*i + fi3 for i in range(steps_3d)])
induction_vals = np.array([ig3*j + ii3 for j in range(steps_3d)])
decay_vals = np.array([dg3*k + di3 for k in range(steps_3d)]) # converts grid indices back into parameter values

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,4)) # plot all 3 marginal plots side by side 
for ax, vals, marginal, label in zip(
    [ax1, ax2, ax3],
    [fitness_vals, induction_vals, decay_vals],
    [marginal_fitness, marginal_induction, marginal_decay],
    ['Fitness', 'Induction', 'Decay']):
    ax.plot(vals, marginal)
    ax.set_xlabel(label)
    ax.set_ylabel("Marginal Probability")
    ax.set_title(f'Marginal: {label}')

plt.tight_layout()
plt.savefig(os.path.join(docs_dir, 'exp_tp53_3d_marginals.png'), dpi=150, bbox_inches='tight')


plt.show()

print(f"Best: fitness={best_fitness:.4f}, induction={best_induction:.6f}, decay={best_decay:.6f}")