# PH-FCox_model
Tumor shape plays a critical role in influencing both growth and metastasis. We introduce a novel topological radiomic feature derived from persistent homology to characterize tumor shape, focusing on its association with time-to-event outcomes in gliomas. These features effectively capture diverse tumor shape patterns that are not represented by conventional radiomic measures. To incorporate these features into survival analysis, we employ a functional Cox regression model in which the topological features are represented in a functional space. We further include interaction terms between shape features and tumor location to capture lobe-specific effects. This approach enables interpretable assessment of how tumor morphology relates to survival risk. We evaluate the proposed method in two case studies using radiomic images of high-grade and low-grade gliomas. The findings suggest that the topological features serve as strong predictors of survival prognosis, remaining significant after adjusting for clinical variables, and provide additional clinically meaningful insights into tumor behavior.

# Simulation
## "./Simulation" folder

<code>simul1_data.R</code> : Generates 300 datasets of the sample size 140 that mimic the real data for simulation.

<code>simul1_data.RData</code> : Contains 300 simulated datasets of the sample size 140.

<code>simul1_simulation.R</code> : Fitting the PH-FCox model and the existing alternatives to the simulated datasets.

<code>simul1_simulation.RData</code> : Contains the model fitting results for the 300 simulated data sets.

