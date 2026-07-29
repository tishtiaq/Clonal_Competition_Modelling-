
# Storing the results and plots from the exponential, linear and step function runs

Using the grid search sweep algorithm to find fitness and induction ONLY, we have shown that this results in a bad fit for the TP53 mouse data. In the following simulations, we have added in the search for an extra variable besides just fitness and induction. In the step function, we search for a time at which the fitness drops immediately from its starting value to neutral. In exponential, we search for a value for the exponential coefficient of fitness and in the linear function, we search for a value for the linear coefficient of fitness.

In each of the sections below you can see how these additional searches have sometimes improved the fit of the simulation to the real data points, and how sometimes it has not. You can see the probability distributions we have found for the different parameters, and both 2D and 3D plots showing heatmaps for maximum likelihood for the parameters. 


## Step Function Simulation
### Searched for t_time in range 50-350
| best fit w/o feedbacks | best fit with feedbacks | best fit with and w/o feedbacks | marginals coarse 2d | marginals tight 2d | marginals 3d | surface plot fit-ind | surface plot fit-time | surface plot ind-time |
|---|---|---|---|---|---|---|---|---|
| ![](StepPlots/step_sim_tp53_best_fit.png)| ![](StepPlots/step_sim_tp53_best_fit_3d.png) | ![](StepPlots/step_sim_best_fits_with_and_without_feedbacks.png) | ![](StepPlots/step_sim_tp53_marginals_coarse.png) | ![](StepPlots/step_sim_tp53_marginals_tight.png) | ![](StepPlots/step_sim_tp53_3d_marginals.png) | ![](StepPlots/step_sim_tp53_3d_surface_fit_ind.png) | ![](StepPlots/step_sim_tp53_3d_surface_fit_dec.png) | ![](StepPlots/step_sim_tp53_3d_surface_ind_dec.png) |

### Searched for t_time in range 0-50
| best fit w/o feedbacks | best fit with feedbacks | best fit with and w/o feedbacks | marginals 3d | surface plot fit-ind | surface plot fit-time | surface plot ind-time |
|---|---|---|---|---|---|---|
| ![](StepPlots/0-50_step_sim_tp53_best_fit.png)| ![](StepPlots/0-50_step_sim_tp53_best_fit_3d.png) | ![](StepPlots/0-50_step_sim_best_fits_with_and_without_feedbacks.png) | ![](StepPlots/0-50_step_sim_tp53_3d_marginals.png) | ![](StepPlots/0-50_step_sim_tp53_3d_surface_fit_ind.png) | ![](StepPlots/0-50_step_sim_tp53_3d_surface_fit_dec.png) | ![](StepPlots/0-50_step_sim_tp53_3d_surface_ind_dec.png) | 

## Exponential Function Simulation 
| best fit w/o feedbacks | best fit with feedbacks | best fit with and w/o feedbacks | marginals coarse 2d | marginals tight 2d | marginals 3d | surface plot fit-ind | surface plot fit-dec | surface plot ind-dec |
|---|---|---|---|---|---|---|---|---|
| ![](ExponentialPlots/step_sim_tp53_best_fit.png)| ![](ExponentialPlots/exp_tp53_best_fit_3d.png) | ![](ExponentialPlots/exp_best_fits_with_and_without_feedbacks.png) | ![](Old+Repeats/tp53_marginals_coarse.png) | ![](Old+Repeats/tp53_marginals_tight.png) | ![](ExponentialPlots/exp_tp53_3d_marginals.png) | ![](ExponentialPlots/exp_tp53_3d_surface_fit_ind.png) | ![](ExponentialPlots/exp_tp53_3d_surface_fit_dec.png) | ![](ExponentialPlots/exp_tp53_3d_surface_ind_dec.png) | 

## Linear Function Simulation 
| best fit w/o feedbacks | best fit with feedbacks | best fit with and w/o feedbacks | surface plot fit-ind | surface plot fit-dec | surface plot ind-dec | marginals 3d |
|---|---|---|---|---|---|---|
| ![](Old+Repeats/tp53_best_fit.png)| ![](LinearPlots/lin_tp53_best_fit_3d.png) | ![](LinearPlots/lin_best_fits_with_and_without_feedbacks.png) | ![](LinearPlots/lin_tp53_3d_surface_fit_ind.png) | ![](LinearPlots/lin_tp53_3d_surface_fit_dec.png) | ![](LinearPlots/lin_tp53_3d_surface_ind_dec.png) | ![](LinearPlots/lin_tp53_3d_marginals.png) |


## Linear Heatmaps and marginals without using feedbacks:
| heatmap loglikelihood tight | heatmap probability tight | marginals coarse 2d | marginals tight 2d |
|---|---|---|---|
|  ![](Old+Repeats/tp53_heatmap_tight_loglik.png)  |  ![](Old+Repeats/tp53_heatmap_tight_prob.png)  | ![](Old+Repeats/tp53_marginals_coarse.png) | ![](Old+Repeats/tp53_marginals_tight.png) |


## Showing the difference between the models

| % Wild-type | Fitness of Mutants | Muller Plots | Mean Clone Size | Survival Rate |
|---|---|---|---|---|
| 90 | 1.05 |  ![](ModelComparisons/models_comparison_muller_plot.png)  | ![](ModelComparisons/models_comparison_mean_clone_size_combined.png) |  ![](ModelComparisons/models_comparison_survival.png) |
| 90 | 1.15 |  ![](ModelComparisons/1.15_models_comparison_muller_plot.png)  | ![](ModelComparisons/1.15_models_comparison_mean_clone_size_combined.png) |  ![](ModelComparisons/1.15_models_comparison_survival.png) |


# Residuals
| Residuals Plot |
|---|
| ![](ModelComparisons/residuals_all_models.png) |













