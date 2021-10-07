## Define function to initialize session
init <- function(r_version, install=T){
  
  ## Source functions
  # Generic functions
  #tryCatch(source("./gfun.R"), error = function(e) stop("Cannot source generic_functions.R"))
  # Plotting functions
  #tryCatch(source("./pfun.R"), error = function(e) stop("Cannot source plotting_functions.R"))
  # Project specifgic functions
  #tryCatch(source("./rn7sk_fun.R"), error = function(e) stop("Cannot sourcern7sk_fun.R"))
  
  ## Set library locations
  lib_path <- paste("../libs/", r_version, "/", sep ="")
  non_shared <<- lib_path
  if(dir.exists(lib_path)) 
  {non_shared <<- normalizePath(lib_path) } else{ 
    cat("Creating library directory \"../libs\" ") 
    if(!dir.create(lib_path)){ stop("Could not create library directory \"../libs/\" ") }
  }
  .libPaths(c(non_shared))
  #x shared_epigenreg = normalizePath("/storage/groups/groups_epigenereg/packages/2017/R/3.4")
  #x shared_server = normalizePath("/home/software/icb/fedora25/r_library")
  
  ## Build folder structure
  # Check if data directory is provided
  if(dir.exists("../data")) { INPUT <<- "../data/" } else{ stop("Input data directory \"data\" does not exist. ") }
  # Set output directory for script outputs
  if(!dir.exists("../results")){ 
    cat("Creating result output directory \"../results\" ") 
    if(!dir.create("../results/")){ stop("Could not create output directory \"../results/\" ") }
  } 
  OUTPUT <<- "../results/"
  # Set logging directory for script logs
  if(!dir.exists("../logs")){ 
    cat("Creating result output directory \"../logs\" ") 
    if(!dir.create("../logs")){ stop("Could not create log directory \"../logs/\" ") }
  }
  LOGS <<- "../logs/"
  # Set logging directory for script logs
  if(!dir.exists("../temp/")){ 
    cat("Creating temp directory \"./temp\" ") 
    if(!dir.create("../temp/")){ stop("Could not create log directory \"../logs/\" ") }
  }
  TEMP <<- "../temp/"
  
  ## Define data sub-directories of input data sets
  eCLIP <<- "eClip/"
  CLIP <<- "CLIP/"
  CHIP <<- "CHIPseq/"
  GENEANNOT <<- "GeneAnnotations/"
  CAGE <<- "CTSS/"
  RN7SK <<- "7SK/"
  RNA <<- "RNAseq/"
  GRO <<- "GROseq/"
  GENOME <<- "Genome/"
  STRING <<- "STRING/"
  CHROMHMM <<- "chromHMM"
  CDATA <<- "CCC/"
  MODEL <<- "Predictions/"
  DNASE <<- "DNAse/"
  SHRNA <<- "shRNAseq/"
  ## Define output directories for results
  PLOTS <<- "Plots/"
  HK <<- "HousekeepingGenes/"
  
  
  #!!! This script uses shared AND non-shared libraries.
  # Shared libraries are those as defined by the EpiGenReg group.
  # Non-shared libraries are the following.
  #installed.packages(lib.loc = non_shared)[, 1] #x
  
  ## Load/install required libraries
  #"FitHiC"
  #"acepack", "Hmisc", 
  # "doParallel", 
  # "parallel",
  # "foreach",
  # GDCRNATools", 
  # LncPath
  # "mvtnorm", "iml",
  # "FactoMineR", 
  #"biomartr",
  
  remotes::install_github("andrewheiss/reconPlots", upgrade = "never")
  
  cran.libs = c( "readr", "fastcluster", "pvclust",  "scales", "reticulate", "parallel", 
                 "foreach", "doParallel",  "LSD", "plyr", "feather", "msigdbr", 
                 "log4r",  "plyr", "ggplot2", "optparse", "tools", "DBI", 
                 "VennDiagram", "bedr", "Rcpp", "rlang", "tidyr",  "stringi",
                 "rlang", "magrittr", "tidyverse", "viridis", "dplyr", "magrittr",
                 "circlize", "dynamicTreeCut", "h2o", "reshape2",
                 "protr", "gridExtra", "caret", "log4r", "optparse", 
                 "fitdistrplus", "ROCR", "reticulate", "xgboost", "data.table", "here",
                 "jsonlite", "mclust", "igraph", "quantreg", "stringr", "SHAPforxgboost", 
                 "anchors", "seqinr", "ape", "ggpubr", "ggforce", "reconPlots", "cowplot",
                 "ggsci", "gridGraphics", "WriteXLS")
  
  
  # "GDCRNATools",
  bcon.libs = c( "msa", "GenomicFeatures","CAGEr", "GenomicRanges", "biomaRt",  "Biostrings", "topGO",  "goSTAG",
                 "tracktables", "GenomicAlignments", "Rsamtools", "BSgenome.Hsapiens.UCSC.hg19",
                 "groHMM", "HiTC", "rtracklayer", "RCAS", "TFutils", "universalmotif", "org.Hs.eg.db", "GOSim")#, , "SGSeq", "SPLINTER")
  
  #"bamsignals"
  
  for( library in c(cran.libs, bcon.libs)){
    # Load package if possible
    if( !require(library, character.only = TRUE)){
      # Install required package
      if(library %in% cran.libs & install){
        install.packages(library, lib = non_shared, repos = "http://cran.us.r-project.org")
      }else{
        if(install){
          #tryCatch(source("https://bioconductor.org/biocLite.R"), error = function(e) stop("Could not connect to Bioconductor.") )
          BiocManager::install(library)
          #biocLite(library, lib = non_shared, lib.loc = non_shared)
        }
      }
      #  Load package after installing
      require(library, character.only = TRUE)
    }
  }
  
  # Negation
  "%ni%" <- Negate("%in%")
  
  ## Setup python interpreter
  #reticulate::use_python("/naslx/projects/t1172/gu27nif2/home/conda3/bin/python", required = T)
  #reticulate::py_config()
  
  ## Initialize script level logging of project specific content
  LOG <<- create.logger()
  logfile(LOG) <- paste(LOGS, Sys.info()["nodename"], ".slog", sep = "")
  level(LOG) <- "DEBUG"
  # log4r::info(LOG,paste("<> Initialized Session <>\n\n", 
  #                       "Library Paths:\n", non_shared, "\n", shared_epigenreg,"\n", shared_server, "\n\n",
  #                       "Session Info:\n", sessionInfo(), "\n\n",
  #                       "Data Input Path:\n", INPUT, "\n\n",
  #                       "Data Output Path:\n", OUTPUT,"\n\n",
  #                       "Notebook Files:\n", NBCACHE,"\n\n",
  #                       "Temporary Files Path:\n", TEMP, "\n\n\n\n",
  #                       sep = ""))
  log4r::info(LOG,"<> Start Session <>\n")
  return(lib_path)
}
## Define function to parse command line arguments
parse.cargs <- function(){
  # Define list of parameter options
  option_list <- list(
    make_option(c("-m", "--host"), type="character", default="remote", 
                help="First argument to the script. Indicates whether the scripts runs locally or on a remote host [default= %default]", metavar="character"),
    make_option(c("-n", "--new"), type="logical", default=F, 
                help="Perform a new run, i.e. discard all previously generated data [default= %default]", metavar="logical"),
    make_option(c("-c", "--cell"), type="character", default=c("K562", "HepG2"), 
                help="Target cell line [default= %default]", metavar="vector"),
    make_option(c("-w", "--workers"), type="numeric", default=c(6,12,18), 
                help="Number of cpu workers per cluster with a low, average and high expected load [default= %default]", metavar="vector")
  )
  ## Actually parse the command line input
  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
  
  # Append names attribute to workers vector to enable indexing by name
  workers <- opt$workers
  names(workers) <- c("low", "avg", "high")
  return(list(cell.line = opt$cell, 
              new       = opt$new,
              workers   = workers))
}
###
# Initialize parallel computation capabilities with machine level logging
###
init.grid <- function(cores){
  # makeForkCluster(nnodes = cores, 
  #                 outfile = paste("./", Sys.info()["nodename"], "-", cores, "x.mlog", sep = ""), 
  #                 type = "PSOCK")
  
  makePSOCKcluster( names = cores,
                    outfile = paste("./", Sys.info()["nodename"], "-", cores, "x.mlog", sep = ""))
}
###############################################################################

#! Recommended/Tested r_version = 4.0.1
# Retrieve current R version 
r_version <- paste(R.version$major, R.version$minor, sep=".")
# Initialize session - Install/Load packages; set global variables; source relevant files
init(r_version, install = T)
CARGS <<- parse.cargs()

