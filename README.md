

# READ ME markdown file for my Clonal_Competition_Modelling- github repo


This repository improves on the code from clone_competition_simulation, a repository by Michael Hall, and code from EvolutionaryModellingPractical-Teacher, a repository by Ben Hall, by introducing time-based feedbacks to the fitness of TP53 cells.

Previously, in clonal simulations, the fitness of TP53 was modelled as a constant float, and this fitted the mouse data from *GIVE FILE NAME OF MOUSE DATA* very well at the early timepoints, but the simulations hugely overpredicted total mutant takeover to what was actually found at the end of the experiment. This led us to believe that there was something slowing down the competitiveness of TP53 cells over time, whether that be a simple fitness function which linearly decreased with time, or 


Guide to getting around:
1. src:
   
2. docs
a) ProbFittingResults.md
This contains all the final results of my project. The heatmaps from the 2D grid searches and from the 3D cube searches are contained within here. Also the marginal plots showing the likelihood functions for fitness, induction and decay. Moreover this contains the residuals plots and mean clone size plots as a way to compare the models. The greatly improved best fit of the exponential and linear functions can be seen, and the poor fit of the step function can also be seen. 
















