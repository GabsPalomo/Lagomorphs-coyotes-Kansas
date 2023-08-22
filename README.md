# Mesopredators have differing influences on prey habitat use and diet activity in a multipredator landscape.

Gabriela Palomo-Munoz, Ty J. Werdel, Colleen W. Piper, Mason Fidino, Travis Gallo, Matthew S. Peek, Andrew M. Ricketts, Adam A. Ahlers. 

## Description

This file includes a description of all the directories and files in this project.

## Directories

1.  /data

This folder includes the following documents: - landcover_scales.csv: occupancy covariates for different scales. - bob_dh.csv: detection history for 28 days across all sites and years for bobcats. - sfox_dh.csv: detection history for 28 days across all sites and years for swift foxes. - badger_dh.csv: detection history for 28 days across all sites and years for badgers. - coyote_dh.csv: detection history for 28 days across all sites and years for coyotes. - BTJR_dh.csv: detection history for 28 days across all sites and years for Black-tailed jackrabbit. - ECTR_dh.csv: detection history for 28 days across all sites and years for Eastern cottontail rabbit. - days_active.csv: - scent.csv - doy.csv - det_altered.csv - fine_scale_habitat.csv:

2.  /figs

This folder contains all figures generated in this project. Stored as .png or .jpeg

3.  /rsd

4.  /functions

5. /tables 
This folder contains all tables generated in this project. Stored as .png or .jpeg

## Files

1.  big_pred_model.R

This script contains the multi-species model made by Mason Fidino.

2.  01_data_org.R

This document organizes the data to have it ready to create the objects that will go in the model.

3.  02_data_indices.R

This is where we format the data (data and indices) to fit into the model.

4.  03_run_model.R

-   Here is the code to run the Bayesian model in JAGS.
-   It also has the marginal occupancy and detection estimates of each species(prey and predator) organized in tables.
-   It sources /functions/functions.R and the 02_data_indices.R file.

5. 04_final_plots.R

This is where I coded the predictions and final graphs. It sources the results of the Bayesian analysis from the document 2022_09_15_fit.rds 