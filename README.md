# Overview
This repository includes all relevant custom R scripts and associated configuration files for the imaging mass cytometry analysis of the J19113 clinical trial in "updated manuscript name".

# Input data
All raw MCD files and fully annotated data frames ("backup_output.rds" and "backup_output_P010.rds") are available on "zenodo link". 
The fully annotated data frames can be loaded onto the R script to generate the published figures in the manuscript. 
In `Config`, there is a metadata, panel, merge (annotation), and area file necessary to generate the heatmaps from Figure 2C.

# R scripts
Scripts stored in `Rscripts` need to be run in numerical order. 
The R script "IMC_J19113_manuscript.R" pertains to the entire clinical trial dataset while the R script "IMC_J19113_P010_manuscript.R" pertains to only P010 - one of the patients with stable disease.

# Output data
`Output_files` contain all figures generated by the R scripts in this repository.
Figures larger than 25MBs were omitted.
