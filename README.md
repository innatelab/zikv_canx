# Data Analysis of proteomics data for "Zika virus-induced remodeling of ER membranes" project


| DOI |
|:----|
| [![DOI](https://zenodo.org/badge/696955720.svg)](https://zenodo.org/badge/latestdoi/696955720) |

A collection of R and Julia scripts for the analysis of data for "SARS-CoV-2/SARS-CoV" project:
  * `prepare_data.R`: MaxQuant data import, normalization, MSGLM model setup
  * `protregroup.jl`: inference of protein groups
  * `effect_scales_normalize.stan`: inference of the scales of AP- and ZikV-associates MSGLM model effects for individual replicates
  * `msglm_fit_chunk.R`: per protein group MSGLM model fitting (run on the [LRZ](https://lrz.de) SLURM compute cluster)
  * `assemble_fits.R`: assembling the results of MSGLM model fitting and report generation
  * `msglm_plots.R`: volcano and scatter plots for MSGLM model fitting results
