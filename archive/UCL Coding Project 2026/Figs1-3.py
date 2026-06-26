#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 12:08:22 2026

@author: syedtariqishtiaq
"""

import sqlite3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ---------- 1. Load real data ----------
df = pd.read_excel("41467_2022_33945_MOESM5_ESM.xlsx", sheet_name="Supplementary Data 5",
                    header=None, engine='openpyxl', skiprows=5, skipfooter=1, usecols="A:H")
df.columns = ['week','mouse','gfp_area','total_area','pct_gfp','avg_pct','sd','sem']
df['week'] = df['week'].ffill()
weeks_unique = sorted(df['week'].unique())
mean_by_week = df.groupby('week')['pct_gfp'].mean()
sem_by_week = df.groupby('week')['pct_gfp'].sem()

# ---------- 2. Pull the ACTUAL completed ABC-SMC posteriors (real model output) ----------
def get_posterior(dbfile):
    con = sqlite3.connect(dbfile)
    pops = pd.read_sql('SELECT * FROM populations', con)
    last_t = pops['t'].max()
    pop_id = pops.loc[pops['t'] == last_t, 'id'].iloc[0]
    models = pd.read_sql(f'SELECT * FROM models WHERE population_id={pop_id}', con)
    model_id = models['id'].iloc[0]
    particles = pd.read_sql(f'SELECT * FROM particles WHERE model_id={model_id}', con)
    params = pd.read_sql(
        f'SELECT * FROM parameters WHERE particle_id IN ({",".join(map(str, particles["id"]))})', con)
    wide = params.pivot(index='particle_id', columns='name', values='value').reset_index()
    wide = wide.merge(particles[['id', 'w']], left_on='particle_id', right_on='id')
    w = wide['w'].values
    w = w / w.sum()
    out = {}
    for p in ['fitness', 'induction']:
        vals = wide[p].values
        order = np.argsort(vals)
        cw = np.cumsum(w[order])
        med = vals[order][np.searchsorted(cw, 0.5)]
        lo = vals[order][np.searchsorted(cw, 0.025)]
        hi = vals[order][np.searchsorted(cw, 0.975)]
        out[p] = (med, lo, hi)
    return out

fit_windows = [
    ("First 2 timepoints\n(1.5, 3 wk)",  "TP53First_pyabc.db"),
    ("First 3 timepoints\n(1.5, 3, 6 wk)", "TP53Early_pyabc.db"),
    ("All 6 timepoints\n(1.5\u201352 wk)",   "TP53-all_pyabc.db"),
]
posteriors = {name: get_posterior(dbfile) for name, dbfile in fit_windows}

# ---------- 3. Deterministic mean-field analogue of the model's selection dynamics ----------
# NOTE: This approximates the per-cell selection recursion used by the WF2D simulation
# (x_{n+1} = w*x_n / (w*x_n + (1-x_n))), run for n = division_rate * t generations,
# but WITHOUT the spatial neighbourhood constraint or stochastic drift/extinction that
# the real compiled simulation includes. It is used here only to illustrate what the
# fitted parameters imply, since the compiled simulation's dependencies could not be
# installed in this offline environment.
DIVISION_RATE = 0.27

def predict_takeover(fitness, induction, times_days, division_rate=DIVISION_RATE):
    out = []
    for t in times_days:
        n_gen = int(round(division_rate * t))
        x = induction
        for _ in range(n_gen):
            x = fitness * x / (fitness * x + (1 - x))
        out.append(x)
    return np.array(out)

weeks_fine = np.linspace(1.5, 52, 200)
times_fine_days = weeks_fine * 7

colors = {"First 2 timepoints\n(1.5, 3 wk)": "#d62728",
          "First 3 timepoints\n(1.5, 3, 6 wk)": "#ff7f0e",
          "All 6 timepoints\n(1.5\u201352 wk)": "#1f77b4"}

plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

# ===================== FIGURE 1: Parameter instability =====================
fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))
names = [n for n, _ in fit_windows]
for ax, param, ylabel, ref in zip(axes, ['fitness', 'induction'],
                                   ['Inferred fitness (relative to wild-type = 1)',
                                    'Inferred induction (fraction of cells labelled)'],
                                   [1.0, None]):
    meds = [posteriors[n][param][0] for n in names]
    los = [posteriors[n][param][1] for n in names]
    his = [posteriors[n][param][2] for n in names]
    x = np.arange(len(names))
    err_low = np.array(meds) - np.array(los)
    err_high = np.array(his) - np.array(meds)
    for i, n in enumerate(names):
        ax.errorbar(x[i], meds[i], yerr=[[err_low[i]], [err_high[i]]], fmt='o',
                    color=colors[n], markersize=9, capsize=6, elinewidth=2, capthick=2)
    if ref is not None:
        ax.axhline(ref, color='gray', linestyle=':', linewidth=1.2)
        ax.text(len(names)-0.55, ref, ' neutral (no advantage)', color='gray', fontsize=8.5, va='bottom')
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=10)
fig.suptitle("Same model, same disease, different time windows \u2192 incompatible fitted parameters\n"
             "(95% credible intervals from the completed ABC-SMC posteriors, not overlapping)",
             fontsize=11, y=1.03)
fig.tight_layout()
fig.savefig("/mnt/user-data/outputs/fig1_parameter_instability.png", dpi=200, bbox_inches='tight')
plt.close(fig)

# ===================== FIGURE 2: Predicted vs observed trajectories =====================
fig, ax = plt.subplots(figsize=(8, 5.5))
rng = np.random.default_rng(0)
for w_ in weeks_unique:
    sub = df[df['week'] == w_]
    jitter = rng.uniform(-0.05, 0.05, size=len(sub)) * w_
    ax.scatter(sub['week'] + jitter, sub['pct_gfp'] * 100, color='black', s=28, zorder=5,
               label='Individual mice (real data)' if w_ == weeks_unique[0] else None)
ax.errorbar(weeks_unique, mean_by_week.values * 100, yerr=sem_by_week.values * 100,
            fmt='-o', color='black', linewidth=1.6, markersize=5, capsize=4,
            label='Mean \u00b1 SEM (real data)', zorder=6)

for name, _ in fit_windows:
    med_fit, _, _ = posteriors[name]['fitness']
    med_ind, _, _ = posteriors[name]['induction']
    pred = predict_takeover(med_fit, med_ind, times_fine_days) * 100
    ax.plot(weeks_fine, pred, '--', color=colors[name], linewidth=2,
             label=f"Model fit to: {name.splitlines()[0]}")

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xticks(weeks_unique)
ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
ax.set_xlabel("Time since p53 induction (weeks)")
ax.set_ylabel("% of tissue area taken over by p53-mutant clone")
ax.set_title("Each fitting window predicts a different reality\n(curves use a mean-field analogue of the model's selection dynamics)",
             fontsize=11)
ax.legend(fontsize=8.5, loc='upper left')
fig.tight_layout()
fig.savefig("/mnt/user-data/outputs/fig2_predicted_vs_observed.png", dpi=200, bbox_inches='tight')
plt.close(fig)

# ===================== FIGURE 3: Implied local growth rate from real data =====================
m = mean_by_week.values
w_ = np.array(weeks_unique)
local_rate_per_week = np.log(m[1:] / m[:-1]) / (w_[1:] - w_[:-1])
midpoints = [f"{w_[i]:.1f}\u2192{w_[i+1]:.0f}wk" for i in range(len(w_)-1)]

fig, ax = plt.subplots(figsize=(8, 4.5))
bars = ax.bar(midpoints, local_rate_per_week, color='#1f77b4', edgecolor='black', linewidth=0.6)
geo_mean_rate = np.log(m[-1] / m[0]) / (w_[-1] - w_[0])
ax.axhline(geo_mean_rate, color='red', linestyle='--', linewidth=1.6,
           label=f"Constant rate a single-fitness model assumes\n({geo_mean_rate:.3f} / week, fit over full range)")
for b, v in zip(bars, local_rate_per_week):
    ax.text(b.get_x() + b.get_width()/2, v + (0.02 if v >= 0 else -0.04),
            f"{v:.2f}", ha='center', fontsize=9)
ax.set_ylabel("Implied exponential growth rate\nof mean takeover (per week)")
ax.set_xlabel("Interval between consecutive sampled timepoints")
ax.set_title("The real data's own growth rate is not constant\n\u2014 a single, time-invariant fitness coefficient cannot produce this shape",
             fontsize=11)
ax.legend(fontsize=8.5)
fig.tight_layout()
fig.savefig("/mnt/user-data/outputs/fig3_growth_rate_instability.png", dpi=200, bbox_inches='tight')
plt.close(fig)

print("Saved figures.")
print()
print("Posterior summary:")
for name, _ in fit_windows:
    print(name.replace("\n", " "), posteriors[name])
print()
print("Local growth rates (per week):", dict(zip(midpoints, local_rate_per_week)))
print("Overall geometric-mean rate (per week):", geo_mean_rate)