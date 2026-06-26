import numpy as np
import matplotlib.pyplot as plt
import matplotlib
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

db_path = "sqlite:///" + "het_pyabc.db"
hist_het = History(db_path)
df_het, w_het = hist_het.get_distribution(m=0)
HET_PARAMS = get_inferred_fit(df_het, w_het);
median fitness 2.2798065488953827
median induction 0.02274515298934408
95CI fitness (2.020980478606173, 2.613100253787201)
95CI induction (0.017229010816853157, 0.028727685574963524)


db_path = "sqlite:///" + "hom_pyabc.db"
hist_hom = History(db_path)
df_hom, w_hom = hist_hom.get_distribution(m=0)
HOM_PARAMS = get_inferred_fit(df_hom, w_hom)
median fitness 7.047314812900756
median induction 0.004549076782422122
95CI fitness (6.188934222147853, 8.616281278944859)
95CI induction (0.000132129100465068, 0.00843781923578664)


db_path = "sqlite:///" + "hom_ctl_pyabc.db"
hist_hom_ctl = History(db_path)
df_hom_ctl, w_hom_ctl = hist_hom_ctl.get_distribution(m=0)
get_inferred_fit(df_hom_ctl, w_hom_ctl);
median fitness 1.0302071592122177
median induction 0.0519616963914189
95CI fitness (0.9588132360712948, 1.1242175608983702)
95CI induction (0.002760051329693719, 0.09703246284647551)


db_path = "sqlite:///" + "het_ctl_pyabc.db"
hist_het_ctl = History(db_path)
df_het_ctl, w_het_ctl = hist_het_ctl.get_distribution(m=0)
get_inferred_fit(df_het_ctl, w_het_ctl);
median fitness 0.9557602570802904
median induction 0.05072329352166699
95CI fitness (0.9252462787104764, 0.9945615082393384)
95CI induction (0.003530570512338681, 0.09802850248573641)


db_path = "sqlite:///" + "wt_pyabc.db"
hist_wt = History(db_path)
df_wt, w_wt = hist_wt.get_distribution(m=0)
WT_PARAMS = get_inferred_fit(df_wt, w_wt)
median fitness 0.9898450258100082
median induction 0.053996727633463
95CI fitness (0.9576270223344407, 1.0252912041178648)
95CI induction (0.001658416740796023, 0.09896343959995683)


all_results = OrderedDict([
    ('WT', (df_wt, w_wt)),
    ('iHET_ctl', (df_het_ctl, w_het_ctl)),
    ('iHOM_ctl', (df_hom_ctl, w_hom_ctl)),
    ('iHET', (df_het, w_het)),
    ('iHOM', (df_hom, w_hom))]
)


ticks = []
bottoms = []
tops = []
medians = []
for k, (df, w) in all_results.items():
    ci_l, ci_u = compute_credible_interval(np.array(df['fitness']), w)
    median = compute_quantile(np.array(df['fitness']), w, alpha=0.5)
    ticks.append(k)
    bottoms.append(ci_l)
    tops.append(ci_u)
    medians.append(median)
plt.figure(figsize=(4.2, 4))
plt.bar(range(5), np.array(tops)-np.array(bottoms), bottom=bottoms, facecolor='w', edgecolor='k', linestyle='-')
for i, m in enumerate(medians):
    plt.plot([i-0.4, i+0.4], [m, m], c='k', linestyle='--', linewidth=1)
plt.ylim(bottom=0)
plt.yticks(range(10))
plt.hlines(1, -1, 10, linestyles='--', linewidth=0.5)
plt.xlim([-0.5, 4.5])
plt.ylim(bottom=0)
plt.ylabel('Fitness')
plt.tight_layout()
plt.xticks(range(5), ['+/+', '+/+\n(iHET)', '+/+\n(iHOM)', '+/-', '-/-']);

# Version with just WT, iHET and iHOM
plt.figure(figsize=(3, 4))
tops2 = np.array(tops)[[0, 3, 4]]
bottoms2 = np.array(bottoms)[[0, 3, 4]]
medians2 = np.array(medians)[[0, 3, 4]]

plt.bar(range(3), np.array(tops2)-np.array(bottoms2), bottom=bottoms2, facecolor='w', edgecolor='k', linestyle='-')
for i, m in enumerate(medians2):
    plt.plot([i-0.4, i+0.4], [m, m], c='k', linestyle='--', linewidth=1)
plt.ylim(bottom=0)
plt.yticks(range(10))
plt.hlines(1, -1, 10, linestyles='--', linewidth=0.5)
plt.xlim([-0.5, 2.5])
plt.ylim(bottom=0)
plt.ylabel('Fitness')
plt.xticks(range(3), ['+/+', '+/-', '-/-'])
plt.tight_layout()
plt.savefig("ABC_inferred_fitness.pdf");

# All results on the fitness-induction plane
plt.figure(figsize=(3, 7))
plt.scatter(df_het['induction'], df_het['fitness'], s=w_het*1000, alpha=0.3, c='m')
plt.scatter(df_hom['induction'], df_hom['fitness'], s=w_hom*1000, alpha=0.3, c='r')
plt.scatter(df_hom_ctl['induction'], df_hom_ctl['fitness'], s=w_hom_ctl*1000, alpha=0.1, c='b')
plt.scatter(df_het_ctl['induction'], df_het_ctl['fitness'], s=w_het_ctl*1000, alpha=0.1, c='c')
plt.scatter(df_wt['induction'], df_wt['fitness'], s=w_wt*1000, alpha=0.1, c='g')
plt.ylim([0, 50])
plt.xlim([0, 0.1])
plt.xlabel('induction') 
plt.ylabel('fitness');

# Zoom in and only show the +/+, +/- and -/- cases
plt.figure(figsize=(3, 5))
plt.scatter(df_het['induction'], df_het['fitness'], s=w_het*1000, alpha=0.3, c='m')
plt.scatter(df_hom['induction'], df_hom['fitness'], s=w_hom*1000, alpha=0.3, c='r')
plt.scatter(df_wt['induction'], df_wt['fitness'], s=w_wt*1000, alpha=0.1, c='k')
plt.ylim([0, 10])
plt.xlim([0, 0.1])
plt.yticks(range(11))
plt.xlabel('induction') 
plt.ylabel('fitness')
plt.tight_layout()
plt.savefig("ABC_inferred_fitness_induction.pdf");

# Best fitting simulations
# Run 100 simulations with the median of the inferred fitness and induction values for the WT, het and hom data.
# Use 70 clones from each simulation to roughly match the experimental number of clones.

from abc_fitting import *
NUM_CLONES = 70
def get_sim_means(parameters, num, data):
    # Take many samples of N clones. 
    # Get mean clone size in each simulation.
    # Get mean of the means and the 95% interval. 
    
    res = [[] for t in data]
    for i in range(num):
        print(i, end=' ')
        np.random.seed(i)
        clone_sizes = run_sim(parameters=parameters, target_data=data, return_clone_sizes=True)
        for j, (t, clones) in enumerate(clone_sizes.items()):
            res[j].append(clones.mean())
        
    
    # Start with single cell clones at time zero
    intervals_high = [1]
    intervals_low = [1]
    means = [1]
    for r in res:
        intervals_high.append(np.quantile(r, 0.975))
        intervals_low.append(np.quantile(r, 0.025))
        means.append(np.mean(r))
    
    return means, intervals_high, intervals_low
from scipy.stats import sem

def plot_data(results, colour=None, label=None, err_stat=sem, elinewidth=2, capsize=4, markersize=3):
    x = results.keys()
    y = [results[k].mean() for k in results]
    yerr = [err_stat(results[k]) for k in results]
    plt.errorbar(x, y, yerr=yerr, label=label, c=colour, 
                 elinewidth=elinewidth, 
                 capsize=capsize, capthick=elinewidth,
                 markersize=markersize,
                 fmt='o')
# Load the data
HET = load_data('het', DATA_FILE)
HOM = load_data('hom', DATA_FILE)
WT = load_data('wt', DATA_FILE)
# Run the simulations
het_means, het_interval_high, het_interval_low = get_sim_means(HET_PARAMS, 100, HET)
0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 
wt_means, wt_interval_high, wt_interval_low = get_sim_means(WT_PARAMS, 100, WT)
0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 
hom_means, hom_interval_high, hom_interval_low = get_sim_means(HOM_PARAMS, 100, HOM)
0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 
plt.figure(figsize=(6.5, 5))

data = HET
colour_sim = 'm'
colour_data = 'm'
times = [0] + [t for t in data]
plt.plot(times, het_means, label='+/- fit', c=colour_sim)
plt.fill_between(times, het_interval_high, het_interval_low, alpha=0.3, color=colour_sim)
plot_data(data, label='+/- data', colour=colour_data, markersize=3)


data = HOM
colour_sim = 'r'
colour_data = 'r'
times = [0] + [t for t in data]
plt.plot(times, hom_means, label='-/- fit', c=colour_sim)
plt.fill_between(times, hom_interval_high, hom_interval_low, alpha=0.3, color=colour_sim)
plot_data(data, label='-/- data', colour=colour_data)

data = WT
colour_sim = 'k'
colour_data = 'k'
times = [0] + [t for t in data]
plt.plot(times, wt_means, label='+/+ fit', c=colour_sim)
plt.fill_between(times, wt_interval_high, wt_interval_low, alpha=0.3, color=colour_sim)
plot_data(data, label='+/+ data', colour=colour_data)

plt.xlim(left=0)
plt.legend(bbox_to_anchor=(1.03, 1))
plt.yscale('log')
plt.ylabel('Mean clone size')
plt.xlabel('Time (days)')
plt.ylim(bottom=1)
plt.yticks([1, 10, 100], [1, 10, 100]);
plt.tight_layout()
plt.savefig('Simulations_mean_clone_size.pdf');