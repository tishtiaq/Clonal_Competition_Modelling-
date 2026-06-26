# Bug Fix Log — clone-competition-simulation

## Bug: `get_vafs_for_all_biopsies` crashes with NaN gene-mutated-id

**Date found:** During tutorial 10 (SimulatingBiopsies), while running the
larger WF2D example with three genes.

**File affected:**
`src/clone_competition_simulation/tissue_sampling/sim_sampling.py`, line 73
(inside the `get_vafs_for_all_biopsies` function).

**Error produced:**
```
ValueError: cannot convert float NaN to integer
```

**Root cause:**

Every simulation starts with one or more "initial clones" — cells present at
time 0 that have not yet acquired any mutation. These clones have `NaN`
stored in their `gene_mutated_id` field (in `clones_array`), since no gene
was responsible for creating them — they were never created by a mutation
at all.

The original line of code:
```python
mutant_gene_map = {i: int(clone[sim.gene_mutated_idx]) for i, clone in enumerate(sim.clones_array)}
```
attempts to convert *every* clone's gene-mutated-id to an integer, including
these initial clones. Converting `NaN` to `int` is undefined and raises a
`ValueError`, regardless of the `remove_initial_clones` setting, because
this filtering only happens later in the function — after this line has
already failed.

This affects ANY simulation that:
1. Has at least one initial clone (i.e. almost every simulation), AND
2. Has a non-zero `mutation_rates` with a `FitnessCalculator` defined.

Switching from `initial_cells` to `initial_grid` (giving every starting cell
a unique clone ID) does NOT avoid the bug — those clones are still "initial"
clones with no mutation, so they still carry `NaN`.

## Fix applied

Changed line 73 to skip the integer conversion for `NaN` values, storing
`None` instead:

```python
# Before:
mutant_gene_map = {i: int(clone[sim.gene_mutated_idx]) for i, clone in enumerate(sim.clones_array)}

# After:
mutant_gene_map = {i: (int(clone[sim.gene_mutated_idx]) if not np.isnan(clone[sim.gene_mutated_idx]) else None) for i, clone in enumerate(sim.clones_array)}
```

`None` is a safe placeholder here because initial clones are filtered out
by `remove_initial_clones` before `mutant_gene_map` is ever looked up for
them later in the function (around line 88, where it fetches the gene name
for display).

## How the fix was applied

Edited directly via Python in Spyder (not by hand in a text editor), using:

```python
filepath = ".../sim_sampling.py"
with open(filepath, "r") as f:
    content = f.read()
old_line = "    mutant_gene_map = {i: int(clone[sim.gene_mutated_idx]) for i, clone in enumerate(sim.clones_array)}"
new_line = "    mutant_gene_map = {i: (int(clone[sim.gene_mutated_idx]) if not np.isnan(clone[sim.gene_mutated_idx]) else None) for i, clone in enumerate(sim.clones_array)}"
content = content.replace(old_line, new_line)
with open(filepath, "w") as f:
    f.write(content)
```

After applying the fix, the Spyder kernel was restarted to clear any cached
imports before re-running the simulation.

## Verification

Re-ran the tutorial 10 WF2D example (10,000 cells, 3 genes, mutation_rates=0.01,
max_time=200) successfully after the fix:
- 1025 clones in the exact ground truth (`biopsy_sample`)
- 29 clones detected after simulating realistic sequencing
  (`detection_limit=5, coverage=100`)

## Notes for the future

- If the package is ever re-downloaded fresh from GitHub (e.g. a new clone
  of the repo, or `pip install` from PyPI), this fix will need to be
  re-applied, since it only exists in this local copy of the source code.
- Worth considering reporting this as a bug to the package's author
  (Michael Hall, Wellcome Sanger Institute) via the GitHub repository's
  issues page: https://github.com/michaelhall28/clone-competition-simulation/issues
- This is exactly the kind of fix that should be mentioned in any
  methods/appendix section of the eventual project writeup, since it
  affects reproducibility for anyone else trying to run the same analysis.
