

# READ ME markdown file for my Clonal_Competition_Modelling- github repo


This repository improves on the code from clone_competition_simulation, a repository by Michael Hall, and code from EvolutionaryModellingPractical-Teacher, a repository by Ben Hall, by introducing time-based feedbacks to the fitness of TP53 cells.

Previously, in clonal simulations, the fitness of TP53 was modelled as a constant float, and this fitted the mouse data from (41467_2022_33945_MOESM5_ESM.xlsx) very well at the early timepoints, but the simulations hugely overpredicted total mutant takeover to what was actually found at the end of the experiment. This led us to believe that there was something slowing down the competitiveness of TP53 cells over time: whether that be a simple fitness function which linearly decreased with time, or something as complex as a response to the fitness of neighbouring cells (i.e if neighbouring cells have higher fitness, your fitness will increase, and vice versa).


Guide to getting around:
1. *src*
The folder CCM contains all the code I have written across the project. ABC_Inference contains all the python code which uses ABC-SMC as a way of inferring fitness. Custom_Classes contains py files of all the custom classes written to override how next generations are generated.
   
2. *docs*
a) ProbFittingResults.md
This contains all the final results of my project. The heatmaps from the 2D grid searches and from the 3D cube searches are contained within here. Also the marginal plots showing the likelihood functions for fitness, induction and decay. Moreover this contains the residuals plots and mean clone size plots as a way to compare the models. The greatly improved best fit of the exponential and linear functions can be seen, and the poor fit of the step function can also be seen. 
b) PlotsOfNonNeutral.md
Contains Muller plots, mean clone size plots and survival rates, using the WF algorithm. The second table uses the WF2D algorithm and also shows the spatial grid plots, which is a plot displaying the wild-types vs mutants at the end of the simulation. The pixels individually show the cells, and so you can see the patterns with which they takeover and spread. 
c) ProbabilisticPlotsTable.md
This is the first introduction to the effect of using the grid search algorithm. It shows the log-likelihoods of fitness and induction as well as a best fit plot, before introducing a fitness feedback and showing how this improves the best fit, whilst also calculating the optimal values for fitness, induction and decay. (By optimal values, I mean the values for which the simulation uses that best fit the actual mouse data). In both decay 3D heatmaps we can see some redundancy where the yellow line (showing highest probability) works for many values of induction at a fixed decay and fitness value. This reflects what we see also in the marginal plots for decay where there are many spikes rather than one clean bell curve: the grid search may struggle to decipher between the parameters.  
However, the final plot shows cleanly what is it we aimed to do with the project: improve the fit of the simulation to the data, especially the later timepoints.
d) ComparisonWithDecreasingFitness.md
Shows simple muller plots, mean clone sizes, survival rates of WF and Moran, varying number of clones, fitness of non-wild-type, number of wild-types, models. 
e) combined_fitting_results.png shows the inferred fitness from the ABC-SMC simulations using 3 sets of timepoints (first2, first3, all6). The box plots are non-overlapping which shows that there must be something decreasing fitness over time, which is not yet in place in the code. This is an excellent plot that shows motivation for my project. 

All other subfolders in the docs section show a range of plots and heatmaps, and are organised based off which model they use. The names for all plots are intuitive.

3. *data*
Contains the excel spreadsheet of the mouse data from the experiment. 















