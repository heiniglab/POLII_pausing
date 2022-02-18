# POLII_pausing
This git repository contains the code accompanying the manuscript

**Predictive model of transcriptional elongation control identifies trans regulatory factors from chromatin signatures**

*Toray S. Akcan*, *Matthias Heinig*

## Table of Contents
  * [Script Overview](#script-overview)
  * [Dependencies](#dependencies)
  * [Data Availability](#data-availability)
  * [Computational Requirements](#computational-requirements)
  * [Script Execution](#script-execution)
  * [Key Data Structures](#key-data-structures)
  * [R Session Info](#r-session-info)
  
# Script Overview

[**Transcriptional_Pausing.R:**](src/Transcriptional_Pausing.R) Main script to perform analyses; Uses all other R-script files.

[**data_preprocessing.R:**](src/data_preprocessing.R) Functions for preprocessing (parsing) of all data sets for usage.

[**data_processing.R:**](src/data_processing.R) Functions for processing of all data sets for downstream analyses.

[**data_analyses.R:**](src/data_analyses.R) Actual analyses of processed datasets.

[**helper_functions.R:**](src/helper_functions.R) Generic helper functions.

# Dependencies
**R-version:** 4.0.3

**Bioconductor-Version:** 3.12 (BiocManager 1.30.15), R 4.0.3 (2020-10-10)

**CRAN-Packages:** 
readr, fastcluster, pvclust, scales, reticulate", parallel, 
foreach, doParallel,  LSD, plyr, feather, msigdbr, 
log4r,  plyr, ggplot2, optparse, tools, DBI, 
VennDiagram, bedr, Rcpp, rlang, tidyr,  stringi,
rlang, magrittr, tidyverse, viridis, dplyr, magrittr,
circlize, dynamicTreeCut, h2o, reshape2,
protr, gridExtra, caret, log4r, optparse, 
fitdistrplus, ROCR, reticulate, xgboost, data.table, here,
jsonlite, mclust, igraph, quantreg, stringr, SHAPforxgboost, 
anchors, seqinr, ape, ggpubr, ggforce, reconPlots, cowplot,
ggsci, gridGraphics, WriteXLS
                 
**Bioconductor-Packages:** 
msa, GenomicFeatures,CAGEr, GenomicRanges, biomaRt,  Biostrings, topGO,  goSTAG,
tracktables, GenomicAlignments, Rsamtools, BSgenome.Hsapiens.UCSC.hg19,
groHMM, HiTC, rtracklayer, RCAS, TFutils, universalmotif, org.Hs.eg.db, GOSim
                 
**External-Packages:** reconPlots

All neccessary packages will be installed automatically upon sourcing [init.R](src/init.R), however system-level dependencies required by specific packages might exist.

# Data Availability
The data folder structure has been replicated in the [data](data) folder and each subfolder contains a *README.txt* with detailed steps to obtain relevant data sets pertaining to each subfolder. This repository is accompanied by a *Zenodo repository* that hosts relevant data sets available at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5236311.svg)](https://doi.org/10.5281/zenodo.5236311)
Please refer to [**resources/folder-structure.txt**](resources/folder-structure.txt) for an overview of the data folder structure.

All necessary data sets with associated accession numbers are also listed in individual xlsx-sheets in file [**data-accessions.xlsx**](resources/data-accessions.xlsx) located in the resources folder.

# Computational Requirements
The following computational resources are recommended to execute the whole script

* 360 GB RAM
* 16 Cores

# Script Execution
1) Clone this repository with **git clone https://github.com/heiniglab/POLII_pausing.git**
2) Specify working directory in [Transcriptional_Pausing.R](src/Transcriptional_Pausing.R) (line 2); Specify number of available cores for low, average and high load computations (line 9); Defaults to 6, 12, 18 cores respectively
3) Obtain all relevant data, see section [**Data Availability**](#data-availability)
4) Run/Source [Transcriptional_Pausing.R](src/Transcriptional_Pausing.R)
5) All plots will be available under [**src/plots**](src/plots) and resulting R-data structures under the [**results**](results) folder

# Key Data Structures
* File **/results/Predictions/model_data/feature.vectors.RDS** contains feature matrices for each cell line
* File **/results/Predictions/model_data/model.matrices.RDS** contains matrices with features (and feature sub-spaces) and targets for each cell line to train predictive models
* File **/results/Predictions/model_evaluation/model.training.results.RDS** contains all model results (incl. model performance table) obtained by training XGB tree models
* File **/results/chipseq.peaks.RDS** contains a genomic ranges object of protein coding transcripts bound by specific proteins on the DNA
* File **/results/eclipseq.peaks.RDS** contains a genomic ranges object of protein coding transcripts bound by specific proteins on the RNA
* File **/results/rn7sk.binders.RDS** contains a genomic ranges object of 7SK transcript variants bound by proteins for each cell line identified from eCLIP-seq data
* File **/results/traveling_ratio.RDS** contains a genomic ranges object with pausing indices/traveling ratios for protein coding transcripts of each cell line

# R Session Info
![R Session info](https://user-images.githubusercontent.com/15715191/137335979-478dca86-5ec5-475c-bd58-79644c9213b7.png)
