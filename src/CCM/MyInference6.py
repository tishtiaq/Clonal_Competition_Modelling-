 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 11:33:19 2026

@author: syedtariqishtiaq
"""

# My code (wrote new) for running inference data. This is for All 6 timepoints. 
# Based off abc_fitting.py (Abby et al)


from pyabc import (ABCSMC, RV, Distribution, PNormDistance)
# Here we're importing classes to use later. pyabc is used for 
# approximate bayesian computation. 
from pyabc.sampler import MulticoreEvalParallelSampler, SingleCoreSampler
from functools import partial, update_wrapper
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
# This is one of the tests used to compare whether two samples come 
# from the same distribution 
from collections import OrderedDict
from clone_competition_simulation.parameters import Parameters
import pyabc.visualization.credible as credible
import sys

from clone_competition_simulation.parameters import (Parameters, 
    PopulationParameters, TimeParameters, FitnessParameters, LabelParameters)                                                 
# Above, we also have to import the parameters as separate classes
# due to the update in the code


# DATA_FILE = "09-03-21 Final clonal counting dataset.xlsx"
# This tells us which data file we're going to be using in the code,
# to plot our inference data over. This is for the NOTCH1 code inference.

SHAPE_OF_GRID = (500,500)
DIVISION_RATE = 0.27
CELLS = SHAPE_OF_GRID[0]*SHAPE_OF_GRID[1]
# These are the fixed parameters for the simulation. Shape tells us
# that we're working in a grid that is 500 by 500; Division Rate tells
# us the rate of division per time step; Cells looks for the 0th
# element of the list SHAPE_OF_GRID and the 1st element of it. These
# are 500 and 500 respectively. It multiples them together, and says 
# the total number of cells is 250,000, their product

ERROR_OBJECT = {'distance': 100000}
# Or any large number will work: the pyabc algorithm will use this as
# an upper bound

LOOP_LIMITS = 50
# 4 grids is a full mouse oesophagus: 50 grids is then 12.5 mice worth
# Therefore if fitness or induction is very low early on, clones will
# die out. So, 50 is a safe cap of, try 50 times, and at least some 
# will survive

def get_grid(fitness, induction, grid_shape, cells):
    initial_grid = np.zeros(grid_shape, dtype = int)
# Sets up initial grid with zeros, in the shape specified in the 
# grid_shape function (defines array dimensions) and data type (dtype) 
# tells the elements to be integers
    total_mutants = int(induction*cells)
# Gives how many total mutants there are to start with. Induction is 
# the percentage of initially infected cells
    mutant_locations = np.random.choice(grid_shape[0]*grid_shape[1], total_mutants, replace=False)
# Randomly scatters a number (depending on induction) of mutant cells 
# across the grid by picking random numbers from the grid. Replace 
# means whether the sample is with or without replacment.
    mutant_locations = [(m // grid_shape[1], m % grid_shape[1]) for m in mutant_locations]
# m // grid_shape[1] tells you which row you're in,
# m % grid_shape[1] tells you how far along that row you are. So we
# convert the exact position of the element in an array rather than 
# just a number between 0 and 249,000
    count = 0
    for i in range(total_mutants):
        initial_grid[mutant_locations[count]] = i + 1
        count += 1 
# This labels each mutant cell with a unique clone ID. It goes through
# the initial grid which is full of zeroes and replaces the coordinates
# given by the mutant_location function with the clone ID
    fitness_array = [1] + [fitness]*total_mutants
# Fitness given in the function, then this creates a list with first
# element 1 and the other elements are the fitness that has been 
# listed, repeated the number of times of how many mutants there are
    label_array = [0] + [1]*total_mutants
# Creates another list with first element 0 and then 1 repeated for
# as many times as there are total mutants
    return initial_grid, fitness_array, label_array
# This is to close the function: this function returns the matrix
# with clone IDs, and a list containing every mutant with the same
# fitness, and each labelled 1

def get_mutant_takeover(sim):
    mutants = sim.population_array.toarray()[1:]
# sim.population_array means to open the sim object, and hand over
# whatever is inside the population_array folder. Then we slice the 
# list to ignore the wild-type cell and only get the mutants.
    mutant_pop = mutants.sum(axis=0)
# This sums DOWN each column. So we get total number of mutant cells
# combined across all mutants clones, at each time point 
    if mutant_pop.max() == sim.total_pop:
            mutant_pop[np.argmax(mutant_pop):] = sim.total_pop
            # returns index of max value in an array
    return mutant_pop/CELLS
# Hence this function gives the percentage of mutated cells in 
# the array. Note that mutant_pop is an array, and we are 
# dividing each element in the array by the float that is CELLS

def distance_ks(target, sim_results):
    return ks_2samp(target, sim_results).statistic
    # This creates a wrapper function to rename the ks_2samp function

def distance_ir2(target, sim_results):
    residuals = target - sim_results
    sumSquareResiduals = np.sum(residuals ** 2)
    diffMean = target - np.mean(target)
    sumSquareDiff = np.sum(diffMean ** 2)
    return((sumSquareResiduals/sumSquareDiff))
# Another tracker of how far simulation is from real results using
# squares instead

def transformBySampling(takeover, samples):
    return([t for (t, individualSample) in zip(takeover, samples) for i in range(individualSample) ])
# The simulation returns 1 mouse per timepoint: yet the real 
# experiment had several mice per timepoint. To be able to compare
# them, both lists need to be the same length so that is what this
# function does

def run_sim(parameters, times, samplesPerTimepoint, target_data, return_takeover=False):
# This is what pyabc actually calls when it wants to test parameters
    fitness, induction = parameters['fitness'], parameters['induction']
    # parameters is a dictionary so this calls certain keys in it
    try: 
        initial_grid, fitness_array, label_array = get_grid(fitness, induction, SHAPE_OF_GRID, CELLS)
        if len(fitness_array) == 1:  # i.e no mutants, induction rate is too low
            return ERROR_OBJECT
        
        p = Parameters(
            algorithm='WF2D',
            population=PopulationParameters(initial_grid=initial_grid, cell_in_own_neighbourhood=True),
            times=TimeParameters(times=times, division_rate=DIVISION_RATE),
            fitness=FitnessParameters(
                initial_fitness_array=fitness_array,
                ), labels=LabelParameters(initial_label_array=label_array))
        # We had to change to using TimeParameters, FitnessParameters,
        # by line 102 of TakeoverFitting code

        full_results = []
# We removed the if here
        for loop in range(LOOP_LIMITS):
            s = p.get_simulator()
            s.run_sim()
            full_results.append(get_mutant_takeover(s))
        takeover = np.mean(full_results, axis=0)

# Runs simulation up to 50 times, takes average of each timepoint.
# This bit of code was changed from the initial to average the 
# timepoints instead of taking just the final one of 50
    
        if return_takeover:
            return transformBySampling(takeover, samplesPerTimepoint)
        
        result = transformBySampling(takeover, samplesPerTimepoint)
        total_distance = distance_ir2(target_data, result)
# A value for target_data was taken in the initial function. This
# just calculates how far apart the result and the target_data are,
# using one of the specified metrics
        return {'distance': total_distance} 
        # returns a dictionary with one key
        
    except (Exception, SystemExit) as e:
        # tells the code which errors to catch and give the caught
        # error a name: "e"
        print('Error!')
        print(e)
        return ERROR_OBJECT # returns the defined upper bound
        
if __name__ == "__main__":
    priors = Distribution(
        fitness = RV("uniform", 0, 50), # can tighten upper bound for greater efficiency
        induction = RV("uniform", 0, 0.1)
        )
# Sets up initial beliefs about the data to apply ABC on (gives
# very wide initial berth. RV is imported from pyabc.)

    DATA_FILE = "41467_2022_33945_MOESM5_ESM.xlsx"
    # load the data file
    times = np.array([7*i for i in [1.5, 3, 6, 12, 24, 52]])          
    # Multiplying by 7 gives us the timepoints in days instead
    # weeks, and transforming it a numpy array makes it easier
    # to work with 
    dataset = pd.read_excel(DATA_FILE, sheet_name="Supplementary Data 5", skiprows=5, skipfooter=1,
                           usecols="E", header=None, engine='openpyxl').to_numpy()[:,0]
    # Open an excel file and read it into in a table pandas 
    # can work with. sheet_name picks the exact sheet to work 
    # with. skiprows says skip the first 5 rows. skipfooter
    # tells us to ignore the last row. usecols says to only use
    # column E which shows % GFP area: the  fraction of the 
    # mouse's tissue taken over by labelled clone. engine says
    # which pandas library to open. to_numpy says convert to
    # a numpy array. [:,0] converts from a 2D object to a 1D array

    samplesPerTimepoint = [4,4,3,4,4,3]
    # number of mice sampled at each timepoint
    distance = PNormDistance()
    # tells pyabc to interpret what its been given as a distance
    sampler = MulticoreEvalParallelSampler()
    # tells pyabc to split up parameter guesses across machine
    # and run tests simulataneously
    
    f = partial(run_sim, times=times, samplesPerTimepoint=samplesPerTimepoint, target_data=dataset)
    # creates a partial function of the run_sim function so that
    # only parameters is still open
    update_wrapper(f, run_sim)
    # a partial object has no __name__: this changes that so 
    # f's __name__ is now run_sim
    abc = ABCSMC(f, priors, distance, population_size=100, sampler=sampler)
    db_path_all = ("sqlite:///" + "TP53All"+'_pyabc.db')
    # constructs address of where the database should live
    
    r = abc.new(db_path_all, {'distance': 0})
    # creates new database file. {'distance': 0} tells pyabc 
    # what the data looks like, in the same shape that run_sim
    # returns its results in
    print("RunID:", r.id)
    hist_all = abc.run(minimum_epsilon=0.1, max_nr_populations=15)
    # gives the 2 stopping conditions on the simulation: either
    # the simulation continues until distance has gone below 0.1
    # or we hit 15 generations
    
    def get_estimate_and_ci_for_param(param, df, w, confidence=0.95):
        vals = np.array(df[param])
        # vals becomes a flat array of 100 fitness values.
        # df[param] is a column in the table with heading 'param' 
        lb, ub = credible.compute_credible_interval(vals, w, confidence)
        median = credible.compute_quantile(vals, w, 0.5)
        return {'median': median, 'CI_lower_bound': lb, 'CI_upper_bound':ub}
    
    df, w = hist_all.get_distribution() 
    # df is the table and w are the weights of each value
    for p in ['fitness', 'induction']:
        print(p, get_estimate_and_ci_for_param(p, df, w))



















