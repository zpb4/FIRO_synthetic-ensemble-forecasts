# FIRO_synthetic-ensemble-forecasts
Repository to support WRR manuscript 'Synthetic forecast ensembles for evaluating Forecast Informed Reservoir Operations (FIRO)'
Submitted 15 March 2023
## Description
The code below supports the data processing, model fitting, generation, and plotting routines to support the aforementioned manuscript.
## Getting started
### Dependencies
Raw data to support this code are stored here: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7688974.svg)](https://doi.org/10.5281/zenodo.7688974)  
Releases of this software are stored permanently here:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10019063.svg)](https://doi.org/10.5281/zenodo.10019063)

Developmental versions of synthetic forecast algorithm are located in the following public GitHub repos:   
https://github.com/zpb4/Synthetic-Forecast-v1-FIRO-DISES   
https://github.com/zpb4/Synthetic-Forecast-v2-FIRO-DISES   
We recommend that interested parties consult these repos for more generalizable implementations of this modeling approach   
### Installing
Requires following R packages:
* BigVAR
* fGarch
### Executing program
The workflow below is configured to run from the file configuration when the Zenodo repository is unzipped and stored in a repository named 'data'
#### syn-HEFS synthetic forecast generation - hindcast period
Numbering indicates order in which scripts must be run  
Runtimes (in parentheses at end) are estimated with parallelization where applicable on an HPC resource 

1) data_process_syn-hefs.R: Processes raw forecast data from individual .csv files in data repository for 3 sites associated with Lake Mendocino (20 min)
2) lamc_init-fit-model.R: Fits initial global parameters for synthetic forecast model, parameter arrays saved in 'fit' repository (5 min)
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
5) ens-process_syn-gefs.R: Processes syn-GEFS .csv files from CNRFC to be R compatible for analysis

#### Plotting routines

- main_plot.R: Main plotting script for manuscript figures
- plot_supp-inf.R: Plotting script for supporting information figures
- calc_ensemble_stats_top-10.R: Calculates cumulative ensemble statistics for top-10 inflow events
- calc_ensemble_stats_top-100.R: Calculates cumulative ensemble statistics for top-10 inflow events
- forecast_verification_functions.R: Helper functions for ensemble forecast verification and plotting
- EFO-results_process.R: Process and plot EFO results for SI
- EFO-results_process_pre-hc.R: Process and plot EFO pre-hindcast results
- plot_cumul_ensembles_hopper.R: Processing and plotting functions for cumulative ensemble forecast plots
- plot_ensembles_hopper.R: Processing and plotting function for ensemble forecast plots
- plot_rank_hist.R: Processing and plotting function for rank histogram verification
- plot_cumul_rank_hist.R: Processing and plotting function for cumulative rank histogram plot
- plot_spread-skill_hopper.R: Processing and plotting function for binned spread error diagrams
- plot_ecrps_hopper.R: Processing and plotting function for ensemble CRPS (eCRPS) verification plots
- plot_top-10-events_hopper.R: Processing and plotting script for 'Top-10' cumulative ensemble statistics
- plot-top-100-events_hopper.R: Processing and plotting script for 'Top-100' cumulative ensemble statistics

#### Miscellaneous

#### Contact
Zach Brodeur, zpb4@cornell.edu
