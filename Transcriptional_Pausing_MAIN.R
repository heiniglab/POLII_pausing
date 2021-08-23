## Set session's working directory to project directory "src" containing source files
#! Modify according to own environment
src.dir <- "/dss/dssfs02/lwp-dss-0001/pn73zi/pn73zi-dss-0000/gu27nif3/Workspace/gene-reg/src/"
if(is.null(setwd(src.dir))){ stop( "Cannot enter src directory.") }
## Initialize session
source("init.R")
source("data_preprocessing.R")
source("data_processing.R")
source("data_analysis.R")
source("helper_functions.R")

# Parrallel execution environment
CARGS$workers <- c(low = 6, avg = 12, high = 16)
WORKERS <- list(low = NULL, avg = NULL, high = NULL)
WORKERS$low <- init.grid(CARGS$workers["low"])
WORKERS$avg <- init.grid(CARGS$workers["avg"])
WORKERS$high <- init.grid(CARGS$workers["high"])

# TODO download metadata files # xargs -L 1 curl -O -J -L < encode_chip_urls_HepG2.txt
# TODO download filtered experiments
#!!! Since the ENCODE narrowPeak format is 0-based per chromsome (https://genome.ucsc.edu/FAQ/FAQformat.html#format12),
# we have to shorten the Gencode gene annotation ranges by 1, since these are 1-based https://www.gencodegenes.org/data_format.html
# DONE To be able to better compare models of pausing vs expression, also restrict expression models to only genes with only one transcript
# already done by subsetting transcripts for transcript expression target by cell line specific traveling ratio target
# TODO mapping of gene domain indices to their associated functions
# TODO Kncokdown data accessions

                                                      #### Data Pre-Processing
# Process gene annotations
GENE.ANNOT <- parse.gene.annotations()
# Process gene quantifications
GENE.QUANT <- load.RNAseq(filter = F, min.reads = 0)
# Process CAGE transcription start sites [unfinished]
CTSS <- load.CAGE.tss(ncores = CARGS$workers["avg"])
# Process CHIP-seq data sets
chipseq.files <- retrieve.CHIPseq.experiments()
CHIPseq       <- parse.CHIPseq.experiments(chipseq.files)
# Process eCLIP-seq data sets
eclipseq.files <- retrieve.eCLIPseq.experiments()
eCLIPseq       <- parse.eCLIPseq.experiments(eclipseq.files)
# Process shared protein domain annotations from HGNC
SDOMAINS <- retrieve.shared.protein.domains()
# Process housekeeping gene annotations
HK.GENES <- load.house.keeping.genes()
# Retrieve sequence specific factors
SEQ.SPEC <- retrieve.sequence.specific.factors()
# Retrieve knockdown gene expressions data
KO.GENE.QUANT <- load.shRNAseq.data()

                                                      #### Data Processing
# Retrieve valid coding transcripts
valid.coding.transcripts <- retrieve.valid.coding.transcripts()
# Retrieve valid non-coding transcripts
valid.non.coding.transcripts <- retrieve.valid.non.coding.transcripts()
# Calculate traveling ratios for protein coding transcripts
traveling.ratios <- calculate.traveling.ratio()
# Retrieve 7SK-ncRNA transcripts
RN7SK.ANNOT <- load.7SKncRNA.data()
# Identify novel 7SK-ncRNA binding proteins
novel.rn7sk.binders <- identify.7SK.binders(RN7SK.ANNOT, eCLIPseq)
# Group transcripts by bound RNA-binding proteins (RBPs)
tx.by.rbp <- get.target.transcripts(eCLIPseq, name = "tx.by.rbp.eclip")
# Group transcripts by bound DNA-binding proteins (DBPs)
tx.by.dbp <- get.target.transcripts(CHIPseq, name = "tx.by.dbp.chip")

                                                    #### Data Analysis
# Build feature vectors for predictive models
feature.vectors <- build.feature.vector()

# Build feature matrices for predictive models
model.matrices <- build.model.matrices(feature.vectors, exclude = c(".*Quant$"))
# Submit models to slurm grid
#source("submit_models.R")

# Retrain models on local machine with optimized parameter setup
model.results <- train.xgb.models(model.matrices, append = "temp")
random.model.results <- train.randomized.xgb.models()

                                                    #### Model interpretations
# Retrieve feature importances of predictive models
evaluate.feature.effects(model.results = model.results,
                         feature.subspace = "All",
                         model.target = "target_traveling.ratio",
                         plot.base.size =  16)

prepare.paper.figures()

plot.base.size = 10                                            #### Descriptive analysis
#visualize.transcript.quantifications()
# visualize.traveling.ratio(plot.base.size)
# visualize.traveling.ratio.transcript.expression.distribution(plot.base.size)
# visualize.7SK.binders()
# visualize.features()



# model.results$performances[model.results$performances$cline=="K562" & 
#                              model.results$performances$modeltype == "individual.model.matrices" &
#                              model.results$performances$target == "target_traveling.ratio", ]
# 
# model.results$performances[model.results$performances$cline=="K562" & 
#                              model.results$performances$modeltype == "synchronised.model.matrices" &
#                              model.results$performances$target == "target_traveling.ratio", ]
# 
# model.results$performances[model.results$performances$cline=="HepG2" & 
#                              model.results$performances$modeltype == "individual.model.matrices" &
#                              model.results$performances$target == "target_traveling.ratio", ]
# 
# model.results$performances[model.results$performances$cline=="HepG2" & 
#                              model.results$performances$modeltype == "synchronised.model.matrices" &
#                              model.results$performances$target == "target_traveling.ratio", ]


## CpG island annotations
# The remaining columns are island length, number of CpGs in the island, the number 
# of C and G in the island, the percentage of island that is CpG, the percentage of 
# island that is C or G, and the ratio of observed(cpgNum) to expected(numC*numG/length) CpG in island.
cnames <- c("chr", "start", "end", "name", "length","cpg.counts","cg.num", "percent.cpg", "percent.cg", "obs.v.exp")
cpg.island.annotations <- read.table(paste0(INPUT, "CpG/cpgIslandExt.txt"), fill = T, sep = "\t", header = F)
cpg.island.annotations <- cpg.island.annotations[, -1]
colnames(cpg.island.annotations) <- cnames
fixed.entries <- which(grepl("fix", cpg.island.annotations$chr))
cpg.island.annotations$chr[fixed.entries] <- gsub("_.*", "", cpg.island.annotations$chr[fixed.entries])
cpg.island.annotations <- 
makeGRangesFromDataFrame(cpg.island.annotations,
                         keep.extra.columns=T,
                         ignore.strand=T,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         starts.in.df.are.0based=FALSE)

## Timecourse data
#timecourse.cata <- import(paste0(INPUT, "TimeCoursePausing/GSE123980_transcript.annotation.timecourse.refined.TSS.corrected.gtf"))



##############################################################################

#X
# model.matrices <- model.matrices$individual.model.matrices$All
# saveRDS(model.matrices, "./feature.matrix.RDS")
# static.features <- model.matrices[, !grepl("^chip|clip", colnames(model.matrices))]
# static.features <- static.features[, !grepl("target_traveling.ratio", colnames(static.features))]
# pol.chip <- model.matrices[, grepl("POL", colnames(model.matrices))]

