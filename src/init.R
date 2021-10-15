## Define function to initialize session
init <- function(r_version, install=T){

  ## Set library locations
  lib_path <- paste("../libs/", r_version, "/", sep ="")
  non_shared <<- lib_path
  if(dir.exists(lib_path)){non_shared <<- normalizePath(lib_path) } else{ 
    cat(paste0("Creating library directory ", normalizePath(lib_path))) 
    if(!dir.create(lib_path)){ stop("Could not create library directory \"../libs/\" ") }else{
      cat(paste0("Created library path ", normalizePath(lib_path))) 
    }
  }
  .libPaths(c(non_shared))
  
  dir.create(paste("./plots/", sep = ""))
  dir.create(paste("../data", sep = ""))
  dir.create(paste("../results", sep = ""))
  dir.create(paste("../logs", sep = ""))
  dir.create(paste("../temp", sep = ""))
  
  ## Define global variables of data directories
  INPUT <<- "../data/"
  OUTPUT <<- "../results/"
  LOGS <<- "../logs/"
  TEMP <<- "../temp/"
  eCLIP <<- "eClip/"
  CHIP <<- "CHIPseq/"
  GENEANNOT <<- "GeneAnnotations/"
  CAGE <<- "Cage/"
  RNA <<- "RNAseq/"
  GRO <<- "GROseq/"
  MODEL <<- "Predictions/"
  PLOTS <<- "Plots/"
  HK <<- "HousekeepingGenes/"
  CPG <<- "CpG/"
  
  ## Package lists
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
                 "ggsci", "gridGraphics", "WriteXLS", "png")
  
  bcon.libs = c( "msa", "GenomicFeatures","CAGEr", "GenomicRanges", "biomaRt",  "Biostrings", "topGO",  "goSTAG",
                 "tracktables", "GenomicAlignments", "Rsamtools", "BSgenome.Hsapiens.UCSC.hg19",
                 "groHMM", "HiTC", "rtracklayer", "RCAS", "TFutils", "universalmotif", "org.Hs.eg.db", "GOSim")#, , "SGSeq", "SPLINTER")
  
  ## Load or install packages if neccessary
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
        }
      }
      #  Load package after installing
      require(library, character.only = TRUE)
    }
  }
  
  # Define new "Negation" operator 
  "%ni%" <- Negate("%in%")
  
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

  makePSOCKcluster( names = cores,
                    outfile = paste("../logs/", Sys.info()["nodename"], "-", cores, "x.mlog", sep = ""))
}
###############################################################################

#! Recommended/Tested r_version = 4.0.1
# Retrieve current R version 
r_version <- paste(R.version$major, R.version$minor, sep=".")
# Initialize session - Install/Load packages; set global variables; source relevant files
init(r_version, install = T)
CARGS <<- parse.cargs()

