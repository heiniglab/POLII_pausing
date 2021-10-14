# POLII_pausing
This git repository contains the code accompanying the manuscript

**Predictive model of transcriptional elongation control identifies trans regulatory factors from chromatin signatures**

Toray S. Akcan, Matthias Heinig

## Table of Contents
  * [Script Outline](#script-outline)
  * [Dependencies](#dependencies)
  * [Data Availability](#data-availability)
  * [Script Execution](#script-execution)
  * [Key Data Structures](#key-data-structures)

# Script Outline 

**Transcriptional_Pausing_MAIN.R:** Main script to perform analyses; Uses all other R-script files.

**data_preprocessing.R:** Functions for preprocessing (parsing) of all data sets for usage.

**data_processing.R:** Functions for processing of all data sets for downstream analyses.

**data_analyses.R:** Actual analyses of processed datasets.

**helper_functions.R:** Generic helper functions.

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

All neccessary packages will be installed automatically upon sourcing "init.R", however system-level dependencies might exist required by specific packages.

# Data Availability
All necessary data sets with associated accession numbers are listed in individual xlsx-sheets in file **"Data Acessions.xlsx"** located in the resources folder. To increase replicability the data folder structure has been replicated in the "data" folder and each subfolder contains a "README.txt" with detailed steps to obtain relevant data sets pertaining each subfolder. 

This git repository is accompanied by a Zenodo repository that hosts relevant data sets available at: **10.5281/zenodo.5236311**

# Script Execution
1) Clone this repository with **git clone https://github.com/heiniglab/POLII_pausing.git**
2) Specify working directory in Transcriptional_Pausing.R (line 2)
3) Specify number of available cores for low, average and high load computations (line 9); Defaults are 6, 12, 16 cores respectively
4) Obtain all relevant data, see previous section **Data Availability**
5) Run/Source Transcriptional_Pausing.R
6) All plots will be avaialbe under **src/plots** and resulting R-data structures under **results**

# Key Data Structures
