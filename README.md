# FIRO_synthetic-ensemble-forecasts
Repository to support WRR manuscript 'Synthetic forecast ensembles for evaluating Forecast Informed Reservoir Operations (FIRO)'
Submitted 15 March 2023
## Description
The code below supports the data processing, model fitting, generation, and plotting routines to support the aforementioned manuscript.
## Getting started
### Dependencies
Raw data to support this code in in the following Zenodo repository: https://doi.org/10.5281/zenodo.7688974
### Installing
Requires following R packages:
* BigVAR
### Executing program
The workflow below is configured to run from the file configuration when the Zenodo repository is unzipped and stored in a repository named 'data'
#### syn-HEFS synthetic forecast generation - hindcast period
Numbering indicates order in which scripts must be run  
Runtimes (in parentheses at end) are estimated with parallelization where applicable on an HPC resource 

1) data_process.R: Processes raw forecast data from individual .csv files in data repository for 3 sites associated with Lake Mendocino (10 min)
2) lamc_init-fit-model.R: Fits initial global parameters for synthetic forecast model, parameter arrays saved in 'fit' repository (10 min)
3) lamc_fit-model.R: Fits model parameters to each ensemble member to enable synthetic forecast generation, parameter arrays saved in 'fit' repository (12 hours)
4) lamc_synthetic-gen_hc.R: Generates specified number of synthetic ensemble samples saved in 'out' repository (6 hrs per 100 samples)
5) error_check-remove_hc.R: Synthetic ensemble post-processing script to remove/replace errant members (10 min per 100 samples)
6) syn-hefs_py-transfer_feather.R: Transfers R array ensemble output to .feather files for compatibility with Python (30 min per 100 samples)
7) ens_sample_npz-transfer.py: Arranges .feather files to Numpy array and zips individual synthetic ensembles for compatibility with EFO model (30 min per 100 samples)

#### syn-HEFS synthetic forecast generation - pre-hindcast period
As above but for out-of-sample pre-hindcast generation

1) lamc_synthetic-gen_pre-hc.R: Generates specified number of pre-hindcast synthetic ensemble samples saved in 'out' repository (6 hrs per 100 samples)
2) error_check-remove_pre-hc.R: Synthetic ensemble post-processing script to remove/replace errant members, pre-hindcast period (10 min)
3) syn-hefs_py-transfer_feather_pre-hc.R: Transfers R array ensemble output to .feather files for compatibility with Python
4) ens_sample_npz-transfer.py: Arranges .feather files to Numpy array and zips individual synthetic ensembles for compatibility with EFO model

#### syn-GEFS synthetic forecast generation
Scripts to generate synthetic meteorological forecasts to support the syn-GEFS generation methodology in the manuscript

1) temp_precip_data-process.R: Processes raw GEFS forecast and observed MAT/MAP data for model fitting (30 min)
2) temp_precip_model-fit.R: Fits synthetic forecast model to meteorological data from step 1 (1 hour)
3) temp_precip_synthetic-gen.R: Generates specified number of syn-GEFS samples to be input into HEFS (FEWS) (1 hour per 100 samples)
4) temp_precip_fews_process.R: Arranges syn-GEFS samples to be processed by external CNRFC HEFS architecture

#### Plotting routines

- main_plot.R: Main plotting script for manuscript figures
- calc_ensemble_stats_top-10.R: Calculates cumulative ensemble statistics for top-10 inflow events
- calc_ensemble_stats_top-100.R: Calculates cumulative ensemble statistics for top-10 inflow events

#### Miscellaneous

#### Contact
Zach Brodeur, zpb4@cornell.edu
