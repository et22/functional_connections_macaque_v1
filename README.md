# Functional Connections Among Neurons within Single Columns of Macaque V1
## Overview
This repository contains code for analyzing cross-correlations between the spikes trains of pairs of neurons in macaque V1 and generating the figures for the paper, [Functional Connections Among Neurons within Single Columns of Macaque V1](https://www.biorxiv.org/content/10.1101/2022.02.18.481095v1).

## Data format 
Raw and intermediate data is not currently included in the repository but will be added and described soon. Once data is included, the following directions can be used to replicate the results in the paper. 


## Code organization 
The main directory contains all scripts necessary to compute, cluster, postprocess, and plot CCGs. The `run_all` function executes each of these scrips in the order in which they should be executed, and executing `run_all` will generate all the results of the paper from the raw data. However, the computation of CCGs can take a substantial amount of time, so CCG analyses were run in parallel on multiple CPUs prior to running the other scripts. Here, we describe each of the scripts in the main directory: 
* `config` stores constants for the analyses
* `compute_ccgs` computes all cross-correlograms for a given session. 
* `postprocess` extracts relevant CCG attributes such as the peak, lag, and neuronal layer pairing for each CCG
* `compute_clusters` clusters the CCGs using parameters and methods defined in config
* `plot_figureX` plots the indicated main figure in the paper
* `plot_distribution` plots the supplemental figure showing the distribution of clusters across layers
* `plot_figures123.ipynb` is a jupyter notebook that contains code necessary for plotting the first three main figures and various supplementary figures

The remaining folders contain various functions organized as follows: 

├───analysis - *helper functions for computing CCGs* \
├───plotting - *helper functions for plotting* \
├───helpers - *generic helper functions* \
├───data - *spike trains from recordings (not currently included)* \
├───output - *intermediate output from analyses (not currently included)*
