# PH-FCox_model
Tumor shape plays a critical role in influencing both growth and metastasis. We introduce a novel topological radiomic feature derived from persistent homology to characterize tumor shape, focusing on its association with time-to-event outcomes in gliomas. These features effectively capture diverse tumor shape patterns that are not represented by conventional radiomic measures. To incorporate these features into survival analysis, we employ a functional Cox regression model in which the topological features are represented in a functional space. We further include interaction terms between shape features and tumor location to capture lobe-specific effects. This approach enables interpretable assessment of how tumor morphology relates to survival risk. We evaluate the proposed method in two case studies using radiomic images of high-grade and low-grade gliomas. The findings suggest that the topological features serve as strong predictors of survival prognosis, remaining significant after adjusting for clinical variables, and provide additional clinically meaningful insights into tumor behavior.

# Application to GBM
## "./GBM" folder
<code>GBM_functions.R</code> : Contains R functions used throughout this research.

<code>GBM_load_PD.R</code> : Loads text files of persistent diagrams (PD) extracted by GUDHI library.

<code>GBM_PH-FCox_loocv.R</code> : Obtain the optimal hyper-parameters for Gaussian smoothing and regularization.

<code>GBM_PH_FCox_loocv.RData</code> : Contains the predictive risks calculated under each hyper-parameter vector.

<code>GBM_PH-FCox_fitting.R</code> : Fits the proposed PH-FCox model to the GBM data set.

<code>GBM_PH_FCox_fitting.RData</code> : Contains final result of fitting the selected PH-FCox model to GBM dataset.

<code>GBM_radiomic-Cox.R</code> : Fits the penalized Cox regression model to the GBM radiomic data set.

<code>GBM_radiomic-Cox.RData</code> : Contains the result of fitting penalized Cox regression model to the GBM radiomic data set.

<code>GBM_SECT-FCox.R</code> : Fits the SECT-FCox model to the GBM radiomic data set.

<code>GBM_SECT-FCox.RData</code> : Contains the result of fitting SECT-FCox model to the GBM radiomic data set.

<code>GBM_dataset.RData</code> : Contains the clinical variables.

<code>gbm_tumor_locs.txt</code> : Contains the xyz-coordinates of tumor centers and corresponding cortical lobes.

# Application to LGG
## "./LGG" folder

<code>LGG_load_PD.R</code> : Loads text files of persistent diagrams (PD) extracted by GUDHI library.

<code>LGG_PH-FCox_loocv.R</code> : Obtain the optimal hyper-parameters for Gaussian smoothing and regularization.

<code>LGG_PH-FCox_fitting.R</code> : Fits the proposed PH-FCox model to the GBM data set.

<code>LGG_radiomic-Cox.R</code> : Fits the penalized Cox regression model to the GBM radiomic data set.

<code>LGG_dataset.RData</code> : Contains the clinical variables.

<code>LGG_PH_FCox_loocv.RData</code> : Contains the predictive risks calculated under each hyper-parameter vector.

<code>LGG_PH_FCox_fitting.RData</code> : Contains final result of fitting the selected PH-FCox model to GBM dataset.

<code>LGG_SECT-FCox.R</code> : Fits the SECT-FCox model to the GBM radiomic data set.

<code>LGG_SECT-FCox.RData</code> : Contains the result of fitting SECT-FCox model to the GBM radiomic data set.

<code>lgg_tumor_locs.txt</code> : contains the xyz-coordinates of tumor centers and corresponding cortical lobes.

# Simulation
## "./Simulation" folder

<code>simul1_data.R</code> : Generates 300 datasets of the sample size 140 that mimic the real data for simulation.

<code>simul1_data.RData</code> : Contains 300 simulated datasets of the sample size 140.

<code>simul1_simulation.R</code> : Fitting the PH-FCox model and the existing alternatives to the simulated datasets.

<code>simul1_simulation.RData</code> : Contains the model fitting results for the 300 simulated data sets.

