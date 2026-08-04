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
from pyabc import History
from pyabc.visualization import plot_kde_2d, plot_kde_matrix, plot_kde_1d, plot_credible_intervals

from pyabc.visualization.credible import compute_credible_interval, compute_kde_max, compute_quantile
from collections import OrderedDict


def get_inferred_fit(df,w):
    # datafile, weights
   
    median_a_intercept = compute_quantile(np.array(df['a_intercept']), w, alpha=0.5)
    median_b_slope = compute_quantile(np.array(df['b_slope']), w, alpha=0.5)
    a_intercept_confidence_interval = compute_credible_interval(np.array(df['a_intercept']), w)
    b_slope_confidence_interval = compute_credible_interval(np.array(df['b_slope']), w)
    print('median a_intercept', median_a_intercept)
    print('median b_slope', median_b_slope)
    print('95CI a_intercept', a_intercept_confidence_interval)
    print('95CI b_slope', b_slope_confidence_interval)
# Works out medians and confidence intervals for the 
# fitness and induction. Creates numpy arrays for both

    return {'a_intercept': median_a_intercept, 'b_slope': median_b_slope,
            'a_intercept_confidence_interval': a_intercept_confidence_interval, 
            'b_slope_confidence_interval': b_slope_confidence_interval}

def plot_results(all_results):
    
    ticks = []
    # the labels that end up printed under each bar on the x-axis
    a_intercept_bottoms, a_intercept_tops, a_intercept_medians = [], [], []
    b_slope_bottoms, b_slope_tops, b_slope_medians = [], [], []

    # the median of each bar
    for k, (df, w) in all_results.items():
        ticks.append(k)

        # a_intercept
        ci_l, ci_u = compute_credible_interval(np.array(df['a_intercept']), w)    
        median = compute_quantile(np.array(df['a_intercept']), w, alpha=0.5)
        a_intercept_bottoms.append(ci_l)
        a_intercept_tops.append(ci_u)
        a_intercept_medians.append(median)
        
        # b_slope
        ci_l, ci_u = compute_credible_interval(np.array(df['b_slope']), w)    
        median = compute_quantile(np.array(df['b_slope']), w, alpha=0.5)
        b_slope_bottoms.append(ci_l)
        b_slope_tops.append(ci_u)
        b_slope_medians.append(median)
 
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))
    
    # Left plot: a_intercept
    ax1.bar(range(3), np.array(a_intercept_tops)-np.array(a_intercept_bottoms), bottom=a_intercept_bottoms, 
                                facecolor='w', edgecolor='k', linestyle='-')
    
    for i, m in enumerate(a_intercept_medians):
        ax1.plot([i-0.4, i+0.4], [m, m], c='k', linestyle='--', linewidth=1)
    ax1.set_ylim(bottom=0.9, top=1.9)
    ax1.set_yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8])
    ax1.axhline(1, color='grey', linestyle='--', linewidth=0.5) # Adds horizontal line spanning axis
    ax1.set_xlim([-0.5, 2.5])
    ax1.set_ylabel('Inferred starting fitness (a_intercept)')
    ax1.set_xticks(range(3))

    ax1.set_xticklabels(['First 2\n(Wks 1.5, 3)', 'First 3\n(Wks 1.5, 3, 6)', 
                         'All 6\n(Wks 1.5-52)'], fontsize=8)
    ax1.set_title('Starting Fitness (a_intercept)')


    # Right plot — b_slope
    ax2.bar(range(3), np.array(b_slope_tops)-np.array(b_slope_bottoms),
            bottom=b_slope_bottoms, facecolor='w', edgecolor='k', linestyle='-')
    for i, m in enumerate(b_slope_medians):
        ax2.plot([i-0.4, i+0.4], [m, m], c='k', linestyle='--', linewidth=1)
    ax2.set_xlim([-0.5, 2.5])
    ax2.set_ylabel('Inferred decay rate (b_slope)')
    ax2.set_xticks(range(3))
    ax2.set_xticklabels(['First 2\n(Wks 1.5, 3)', 'First 3\n(Wks 1.5, 3, 6)', 
                         'All 6\n(Wks 1.5-52)'], fontsize=8)
    ax2.set_title('Decay Rate (b_slope)')

    fig.suptitle('Time-Varying Fitness Inference: Constant vs Linearly Decreasing',
                 fontsize=11)
    plt.tight_layout()
    plt.savefig('docs/combined_fitting_results_timevarying.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    
    
def main():
    # Setting up databases
    
    db_path_tv_first2 = "sqlite:///" + "TP53First2_TimeVarying_pyabc.db"
    # tells the next line where to open
    hist_tv_first2 = History(db_path_tv_first2)
    # creates an object in the class History (part of pyabc)
    # it is an object storing everything about the data file excel_sheet
    df_tv_first2, w_tv_first2 = hist_tv_first2.get_distribution(m=0)
    # hist_first is an object, get_distribution is a function you can run on it.
    # m is model index, so this asks for model 0 (WF2D)
    FIRST2_PARAMS = get_inferred_fit(df_tv_first2, w_tv_first2)
    # This is a dictionary containing 4 items, which are specifically the 4 items in 
    #  return' on line 38. Just adapted for each dataset. 
    
    db_path_tv_first3 = "sqlite:///" + "TP53First3_TimeVarying_pyabc.db"
    hist_tv_first3 = History(db_path_tv_first3)
    df_tv_first3, w_tv_first3 = hist_tv_first3.get_distribution(m=0)
    FIRST3_PARAMS = get_inferred_fit(df_tv_first3, w_tv_first3)
    # creates file for the first 3 timepoints
   
    db_path_tv_all = "sqlite:///" + "TP53All_TimeVarying_pyabc.db"
    hist_tv_all = History(db_path_tv_all)
    df_tv_all, w_tv_all = hist_tv_all.get_distribution(m=0)
    ALL_PARAMS = get_inferred_fit(df_tv_all, w_tv_all)
    # creates file for all timepoints
    

    all_results = OrderedDict([('TP53First2_TimeVarying', (df_tv_first2, w_tv_first2)),('TP53First3_TimeVarying',(df_tv_first3, w_tv_first3)),
                               ('TP53All_TimeVarying', (df_tv_all, w_tv_all))])
    # Creates list of tuples containing the data sets. OrderedDict preserves the order that key-value pairs (items)
    # are stored in the dictionary
    plot_results(all_results)


if __name__ == "__main__":
    main()
    
    