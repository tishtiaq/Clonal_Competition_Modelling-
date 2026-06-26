import sys
db = sys.argv[1]

import numpy as np

import os
from pyabc.visualization import plot_kde_2d, plot_kde_matrix, plot_kde_1d, plot_credible_intervals
from pyabc import History

from pyabc.visualization.credible import compute_credible_interval, compute_kde_max, compute_quantile

from collections import OrderedDict

def get_inferred_fit(df, w):
    median_fitness = compute_quantile(np.array(df['fitness']), w, alpha=0.5)
    median_induction = compute_quantile(np.array(df['induction']), w, alpha=0.5)
    print('median fitness', median_fitness)
    print('median induction', median_induction)
    print('95CI fitness', compute_credible_interval(np.array(df['fitness']), w))
    print('95CI induction', compute_credible_interval(np.array(df['induction']), w))
    return {'fitness': median_fitness, 'induction': median_induction}

limits = {"fitness": [0, 50], "induction": [0, 0.1]}  # These were the ranges of the prior parameter distributions

db_path = "sqlite:///" + db
hist_tp53 = History(db_path)
df_tp53, w_tp53 = hist_tp53.get_distribution(m=0)

TP53_PARAMS = get_inferred_fit(df_tp53, w_tp53)
