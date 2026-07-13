# Non-Neutral Simulations: Modelling the difference between a large clone of wild-types and a variety of non-neutral single-cell clones, all of the same fitness. In the two simulations, we show the difference between non-neutral clones which hold a constant fitness and those which have linearly decreasing fitness (ending the simulation at a wild-type  fitness of 1)

Using WF algorithm:

| Total cells | No. of cells in wild-type clone | No. of non-neutral cells | Fitness of non-neutral cells | Muller Plot | Mean Clone Size | Clone Survival |
|---|---|---|---|---|---|---| 
| 10000 |     9950      |    50    |     1.3     |  ![](/docs/9950v50_1.3_muller_comparison_wf.png)     |  ![](/docs/9950v50_1.3_size_comparison_wf.png)   |   ![](/docs/9950v50_1.3_survival_comparison_wf.png)   | 
| 10000 |     9900      |    100   |     1.3     |  ![](/docs/9900v100_1.3_muller_comparison_wf.png)    |  ![](/docs/9900v100_1.3_size_comparison_wf.png)  |   ![](/docs/9900v100_1.3_survival_comparison_wf.png)  | 
| 10000 |     9900      |    100   |     1.15     |  ![](/docs/9900v100_1.15_muller_comparison_wf.png)    |  ![](/docs/9900v100_1.15_size_comparison_wf.png)  |   ![](/docs/9900v100_1.15_survival_comparison_wf.png)  | 
| 10000 |     9900      |    100   |     1.05     |  ![](/docs/9900v100_1.05_muller_comparison_wf.png)    |  ![](/docs/9900v100_1.05_size_comparison_wf.png)  |   ![](/docs/9900v100_1.05_survival_comparison_wf.png)  | 


| Total cells | No. of cells in wild-type clone | No. of non-neutral cells | Fitness of non-neutral cells | Muller Plot | Spatial Grid | 
|---|---|---|---|---|---|
| 100 |   50   |   50  |   1.3   |    ![](WF2D_1.3_muller_comparison.png)    |   ![](WF2D_1.3_spatial_comparison.png)    |
| 100 |   50   |   50  |   1.15  |    ![](WF2D_1.15_muller_comparison.png)   |   ![](WF2D_1.15_spatial_comparison.png)   |
| 100 |   50   |   50  |   1.05  |    ![](WF2D_1.05_muller_comparison.png)   |   ![](WF2D_1.05_spatial_comparison.png)   |

Using Moran algorithm: 

| Total cells | No. of cells in wild-type clone | No. of non-neutral cells | Fitness of non-neutral cells | Muller Plot | Spatial Grid | 
|---|---|---|---|---|---|
| 100 |   50   |   50  |   1.3   |    ![]()    |   ![]()    |
| 100 |   50   |   50  |   1.15  |    ![]()    |   ![]()    |

