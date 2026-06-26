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
import matplotlib

import os
from pyabc.visualization import plot_kde_2d, plot_kde_matrix, plot_kde_1d, plot_credible_intervals
from pyabc import History

from pyabc.visualization.credible import compute_credible_interval, compute_kde_max, compute_quantile
from collections import OrderedDict


def get_inferred_fit(df,w):
    median_fitness = compute_quantile(np.array(df['fitness']), w, alpha=0.5)
    median_induction = compute_quantile(np.array(df['induction']), w, alpha=0.5)
    print('median fitness:', median_fitness)
    print('median induction:', median_induction)
    print('95CI fitness:', compute_credible_interval(np.array(df['fitness']), w))
    print('95CI induction:', compute_credible_interval(np.array(df['induction']), w))
# Works out medians and confidence intervals for the 
# fitness and induction. Creates numpy arrays for the both

    return {'fitness': median_fitness, 'induction': median_induction}

db_path = "sqlite:///" + "TP53First_pyabc.db"
# tells the next line where to open
hist_first = History(db_path)
# 
df_first, w_first = hist_first.get_distribution(m=0)

# hist_first is an object, get_distribution is a function you can run on it. m is model index, so 
# this asks for model 0 (our setup only ever uses 1 model anyway)
FIRST_PARAMS = get_inferred_fit(df_first, w_first)

















