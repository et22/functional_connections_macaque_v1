# Functional Connections Among Neurons within Single Columns of Macaque V1
## Overview
This repository contains code for analyzing cross-correlations between the spikes trains of pairs of neurons in macaque V1 and generating the figures for the paper, [Functional Connections Among Neurons within Single Columns of Macaque V1](https://www.biorxiv.org/content/10.1101/2022.02.18.481095v1).

## Data format 
Here, we include intermediate data files in the `output` directory that can be used to replicate all figures. The intermediate files contain the cross-correlograms for all pairs of neurons and the output of clustering these cross-correlograms. 

## Reproducing figures
To reproduce all figures, clone the repo and run `plot_all.m` in MATLAB and `plot_figures123.ipynb` in Jupyter Notebook or VS code.

## Reproducing all analyses from raw data
Raw data has been deposited as described in the publication, and can be used to reproduce the results.

The main directory contains all scripts necessary to compute, cluster, postprocess, and plot CCGs. The `run_all` function executes each of these scrips in the order in which they should be executed. To reproduce all figures from the raw data, place raw data in a directory named `data` and then run `run_all.m`  Note that the computation of jitter-corrected CCGs can take a substantial amount of time, so CCG analyses were run in parallel on multiple CPUs prior to running the other scripts. 

## Code organization 
Here, we describe each of the scripts in the main directory: 
* `config` stores constants for the analyses
* `compute_ccgs` computes all cross-correlograms for a given session. 
* `postprocess` extracts relevant CCG attributes such as the peak, lag, and neuronal layer pairing for each CCG
* `compute_clusters` clusters the CCGs using parameters and methods defined in config
* `plot_figureX` plots the indicated main figure in the paper
* `plot_distribution` plots the supplemental figure showing the distribution of clusters across layers
* `plot_figures123.ipynb` is a jupyter notebook that contains code necessary for plotting the first three main figures and various supplementary figures

The remaining folders are organized as follows: 
├───analysis - *helper functions for computing CCGs* \
├───helpers - *generic helper functions* \
├───data - *spike trains from recordings (not currently included)* \
├───output - *intermediate output from analyses*
