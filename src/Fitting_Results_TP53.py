#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 16:18:31 2026

@author: syedtariqishtiaq
"""

# FOR PLOTTING GRAPHS
# My code (wrote new) for PLOTTING inference data
# Code based off of Fitting_results (Abby et al)

import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
import os
from pyabc.visualization import plot_kde_2d, plot_kde_matrix, plot_kde_1d, plot_credible_intervals
from pyabc import History

from pyabc.visualization.credible import compute_credible_interval, compute_kde_max, compute_quantile
from collections import OrderedDict


def get_inferred_fit(df,w):
    # datafile, weights
    median_fitness = compute_quantile(np.array(df['fitness']), w, alpha=0.5)
    median_induction = compute_quantile(np.array(df['induction']), w, alpha=0.5)
    fitness_confidence_interval = compute_credible_interval(np.array(df['fitness']), w)
    induction_confidence_interval = compute_credible_interval(np.array(df['induction']), w)
    print('median fitness:', median_fitness)
    print('median induction:', median_induction)
    print('95CI fitness:', fitness_confidence_interval)
    print('95CI induction:', induction_confidence_interval)
# Works out medians and confidence intervals for the 
# fitness and induction. Creates numpy arrays for both

    return {'fitness': median_fitness, 'induction': median_induction,
            'fitness_confidence_interval': fitness_confidence_interval, 
            'induction_confidence_interval': induction_confidence_interval}

def plot_results(all_results):
    
    ticks = []
    # the labels that end up printed under each bar on the x-axis
    bottoms = []
    # the lower edge of each bar
    tops = []
    # the upper edge of each bar
    medians = []
    # the median of each bar
    for k, (df, w) in all_results.items():
        ci_l, ci_u = compute_credible_interval(np.array(df['fitness']), w)
        # create tuple of upper and lower bounds of confidence interval 
        median = compute_quantile(np.array(df['fitness']), w, alpha=0.5)
        ticks.append(k)
        # k is the keys in all_results. This appends the empty list ticks with 'TP53First2', etc...
        bottoms.append(ci_l)
        # append the bottoms list with the lower bound for confidence intervals
        tops.append(ci_u)
        # ""
        medians.append(median)
        # each of these lists wil end up with 3 entries by the end. 
    
    plt.figure(figsize=(4.2,4))
    # also may need changing to better fit data
    plt.bar(range(3), np.array(tops)-np.array(bottoms), bottom=bottoms, facecolor='w', edgecolor='k', linestyle='-')
    for i, m in enumerate(medians):
        # enumerate allows you to loop through an iterable with access to the index and element itself
        plt.plot([i-0.4, i+0.4], [m, m], c='k', linestyle='--', linewidth=1)
        
    plt.ylim(bottom=0)
    # sets max and min of y-axis
    plt.yticks(range(10))
    # forces y-axis to show a labelled tick mark at every whole number in the range 
    plt.hlines(1,-1, 10, linestyles='--', linewidth=0.5)
    # plots horizontal line at y=1 from x=-1 to x=10. y=1 is showing wild-type.
    plt.xlim([-0.5, 2.5])
    # sets max and min of x-axis
    
    # TODO: Change these axis ranges to be suitable for TP53 
    
    plt.ylabel('Inferred Fitness relative to wild-type')
    plt.tight_layout()
    # measures what is currently in the figure and adjusts internal spacing
    
    plt.xticks(range(3), ['First 2 timepoints (Weeks 1.5, 3)', 'First 3 timepoints (Weeks 1.5, 3, 6)',
                          'All 6 timepoints (Weeks 1.5, 3, 6, 12, 24, 52)'])
    
    # TODO: Add something which saves graph to files 
    
def main():
    # Setting up databases
    db_path_first2 = "sqlite:///" + "TP53First2_pyabc.db"
    # tells the next line where to open
    hist_first2 = History(db_path_first2)
    # creates an object in the class History (part of pyabc)
    # it is an object storing everything about the data file excel_sheet
    df_first2, w_first2 = hist_first2.get_distribution(m=0)
    # hist_first is an object, get_distribution is a function you can run on it.
    # m is model index, so this asks for model 0 (WF2D)
    FIRST2_PARAMS = get_inferred_fit(df_first2, w_first2)
    # This is a dictionary containing 4 items, which are specifically the 4 items in 
    # 'return' on line 38. Just adapted for each dataset. 
    
    db_path_first3 = "sqlite:///" + "TP53First3_pyabc.db"
    hist_first3 = History(db_path_first3)
    df_first3, w_first3 = hist_first3.get_distribution(m=0)
    FIRST3_PARAMS = get_inferred_fit(df_first3, w_first3)
    # creates file for the first 3 timepoints
    
    db_path_all = "sqlite:///" + "TP53All_pyabc.db"
    hist_all = History(db_path_all)
    df_all, w_all = hist_all.get_distribution(m=0)
    ALL_PARAMS = get_inferred_fit(df_all, w_all)
    # creates file for all timepoints
    
    
    all_results = OrderedDict([('TP53First2', (df_first2, w_first2)),
                                ('TP53First3', (df_first3, w_first3)),
                                ('TP53All', (df_all, w_all))])
    # Creates list of tuples containing the data sets. OrderedDict preserves the order that key-value pairs (items)
    # are stored in the dictionary
    plot_results(all_results)



if __name__ == "__main__":
    main()
    
    



























