# PH-FCox_model
Tumor shape plays a critical role in influencing both growth and metastasis. We introduce a novel topological feature obtained using persistent homology, to characterize tumor shape with a focus on its relationship to time-to-event data. These features, which are invariant to scale-preserving transformations, effectively capture diverse tumor shape patterns. To utilize these features in survival analysis, we adopt a functional Cox regression model that takes the shape features represented in a functional space as predictors. Furthermore, interactions between shape features and tumor location are incorporated to reflect location-specific effects. This approach enables interpretable analysis of the association between tumor shape characteristics and survival risks. Two case studies were conducted using radiomic images of high-grade and low-grade gliomas. The findings suggest that the topological features serve as strong predictors of survival prognosis, remaining significant after adjusting for clinical variables, and contribute further clinically relevant understanding.

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

<code>gbm_tumor_locs.txt</code> : contains the xyz-coordinates of tumor centers and corresponding cortical lobes.

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

