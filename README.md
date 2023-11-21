# Mesopredators have differing influences on prey habitat use and diel activity in a multipredator landscape.

Gabriela Palomo-Munoz, Mason Fidino, Ty J. Werdel, Colleen W. Piper, Travis Gallo, Matthew S. Peek, Andrew M. Ricketts, Adam A. Ahlers.

## Description

This `README` includes a description on the various scripts used for this analysis. However, please note that not all data files are saved in this repository because they may be too large to store here. I have the scripts numbered so it's easier to see the flow of the analysis and where to start.

## RDS file with results of the Bayesian model.

The Bayesian model took several weeks to run, but I saved the output in a .rds file. However, because it is too large it cannot be stored here but it can be pulled from [here](https://drive.google.com/file/d/1vbbBUzLehAv0XIIYJC7YBIW2b2YgAsCz/view?usp=sharing) as a zip file and read in the session and name it `fit` so that all analysis can be done. That way, anyone can corroborate my analysis.

## Directories

1.  `/data` This folder includes the following documents:

    -   landcover_scales.csv: occupancy covariates for different scales.
    -   bob_dh.csv: detection history for 28 days across all sites and years for bobcats.
    -   sfox_dh.csv: detection history for 28 days across all sites and years for swift foxes.
    -   badger_dh.csv: detection history for 28 days across all sites and years for badgers.
    -   coyote_dh.csv: detection history for 28 days across all sites and years for coyotes.
    -   BTJR_dh.csv: detection history for 28 days across all sites and years for Black-tailed jackrabbit.
    -   ECTR_dh.csv: detection history for 28 days across all sites and years for Eastern cottontail rabbit.
    -   days_active.csv:
    -   scent.csv
    -   doy.csv
    -   det_altered.csv
    -   fine_scale_habitat.csv

2.  `/figs` This folder contains all figures generated in this project. Stored as .png or .jpeg. NOT available in the repository but the figures are in the manuscript.

3.  `/functions` This folder contains some long functions necessary for some data cleaning and analysis, as well as reading the silhouettes.

4.  `/tables` This folder contains all tables generated in this project. Stored as .png or .jpeg. NOT available in the repository but the tables are in the manuscript.

## Files

1.  `big_pred_model.R`

This script contains the multi-species model made by Mason Fidino and Gabriela Palomo.

2.  `01_data_org.R`

This document organizes the data to have it ready to create the objects that will go in the model.

3.  `02_data_indices.R`

This is where we format the data (data and indices) to fit into the model.

4.  `03_run_model.R`

    -   Here is the code to run the Bayesian model in JAGS.

    -   It also has the marginal occupancy and detection estimates of each species(prey and predator) organized in tables.

    -   It sources /functions/functions.R and the 02_data_indices.R file.

5.  `04_final_plots.R`

This is where I coded the predictions and final graphs. It sources the results of the Bayesian analysis from the document 2022_09_15_fit.rds

All plots coded here are in the manuscript as figures.

6.  `05_1_summary_table.Rmd`

This script contains the code necessary to create a table with the mean and 95% credible intervals of each of the model parameters. Estimates are in the logit scale and appear in **Appendix S2** of the manuscript.
