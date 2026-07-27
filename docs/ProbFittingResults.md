
# Storing the results and plots from the exponential, linear and step function runs

Using the grid search sweep algorithm to find fitness and induction ONLY, we have shown that this results in a bad fit for the TP53 mouse data. In the following simulations, we have added in the search for an extra variable besides just fitness and induction. In the step function, we search for a time at which the fitness drops immediately from its starting value to neutral. In exponential, we search for a value for the exponential coefficient of fitness and in the linear function, we search for a value for the linear coefficient of fitness.

In each of the sections below you can see how these additional searches have sometimes improved the fit of the simulation to the real data points, and how sometimes it has not. You can see the probability distributions we have found for the different parameters, and both 2D and 3D plots showing heatmaps for maximum likelihood for the parameters. 


## Step Function Simulation
| best fit w/o feedbacks | best fit with feedbacks | best fit with and w/o feedbacks | marginals coarse 2d | marginals tight 2d | marginals 3d | surface plot fit-ind | surface plot fit-time | surface plot ind-time |
|---|---|---|---|---|---|---|---|---|
| ![](step_sim_tp53_best_fit.png)| ![](step_sim_tp53_best_fit_3d.png) | ![](step_sim_best_fits_with_and_without_feedbacks.png) | ![](step_sim_tp53_marginals_coarse.png) | ![](step_sim_tp53_marginals_tight.png) | ![](step_sim_tp53_3d_marginals.png) | ![](step_sim_tp53_3d_surface_fit_ind.png) | ![](step_sim_tp53_3d_surface_fit_dec.png) | ![](step_sim_tp53_3d_surface_ind_dec.png) |

## Exponential Function Simulation 
| best fit w/o feedbacks | best fit with feedbacks | best fit with and w/o feedbacks | marginals coarse 2d | marginals tight 2d | marginals 3d | surface plot fit-ind | surface plot fit-dec | surface plot ind-dec |
|---|---|---|---|---|---|---|---|---|
| ![](step_sim_tp53_best_fit.png)| ![](exp_tp53_best_fit_3d.png) | ![](exp_best_fits_with_and_without_feedbacks.png) | ![](tp53_marginals_coarse.png) | ![](tp53_marginals_tight.png) | ![](exp_tp53_3d_marginals.png) | ![](exp_tp53_3d_surface_fit_ind.png) | ![](exp_tp53_3d_surface_fit_dec.png) | ![](exp_tp53_3d_surface_ind_dec.png) | 

## Linear Function Simulation 
| best fit w/o feedbacks | best fit with feedbacks | best fit with and w/o feedbacks | surface plot fit-ind | surface plot fit-dec | surface plot ind-dec | marginals 3d |
|---|---|---|---|---|---|---|
| ![](tp53_best_fit.png)| ![](lin_tp53_best_fit_3d.png) | ![](lin_best_fits_with_and_without_feedbacks.png) | ![](lin_tp53_3d_surface_fit_ind.png) | ![](lin_tp53_3d_surface_fit_dec.png) | ![](lin_tp53_3d_surface_ind_dec.png) | ![](lin_tp53_3d_marginals.png) |


## Linear Heatmaps and marginals without using feedbacks:
| heatmap loglikelihood tight | heatmap probability tight | marginals coarse 2d | marginals tight 2d |
|---|---|---|---|
|  ![](tp53_heatmap_tight_loglik.png)  |  ![](tp53_heatmap_tight_prob.png)  | ![](tp53_marginals_coarse.png) | ![](tp53_marginals_tight.png) |


## Showing the difference between the models

| % Wild-type | Fitness of Mutants | Muller Plots | Mean Clone Size | Survival Rate |
|---|---|---|---|---|
| 90 | 1.05 |  ![](models_comparison_muller_plot.png)  | ![](models_comparison_mean_clone_size_combined.png) |  ![](models_comparison_survival.png) |
| 90 | 1.15 |  ![](1.15_models_comparison_muller_plot.png)  | ![](1.15_models_comparison_mean_clone_size_combined.png) |  ![](1.15_models_comparison_survival.png) |


# Residuals
| Residuals Plot |
|---|
| ![](residuals_all_models.png) |













