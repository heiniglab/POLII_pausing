# Set current working directory
setwd("YOUR/PATH/TO/SRC/WORKING/DIRECTORY")
# Initialize all accompanying scripts
source("init.R")
source("data_preprocessing.R")
source("data_processing.R")
source("data_analysis.R")
source("helper_functions.R")

# Specify number of available cores for low, average and high load computations
CARGS$workers <- c(low = 6, avg = 12, high = 18)
WORKERS <- list(low = NULL, avg = NULL, high = NULL)
# Initialize parallel workers
WORKERS$low <- init.grid(CARGS$workers["low"])
WORKERS$avg <- init.grid(CARGS$workers["avg"])
WORKERS$high <- init.grid(CARGS$workers["high"])

                                                      #### Data Pre-Processing
# Process gene annotations
GENE.ANNOT <- parse.gene.annotations()
# Process gene quantifications
GENE.QUANT <- load.RNAseq()
# Process CAGE transcription start sites [unfinished]
CTSS <- load.CAGE.tss(ncores = CARGS$workers["avg"])
# Process CHIP-seq data sets
chipseq.files <- retrieve.CHIPseq.experiments()
CHIPseq       <- parse.CHIPseq.experiments(chipseq.files)
# Process eCLIP-seq data sets
eclipseq.files <- retrieve.eCLIPseq.experiments()
eCLIPseq       <- parse.eCLIPseq.experiments(eclipseq.files)
# Process housekeeping gene annotations
HK.GENES <- load.house.keeping.genes()
# Retrieve sequence specific factors
SEQ.SPEC <- retrieve.sequence.specific.factors()
# Load CpG Island Annotations
CPG.ISLANDS <- load.cpg.island()

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
# Retrain models on local machine with optimized parameter setup
model.results <- train.xgb.models(model.matrices)
# Run random models
random.model.results <- train.randomized.xgb.models()

                                                    #### Model interpretations
# Retrieve feature importances of predictive models
evaluate.feature.effects(model.results = model.results,
                         feature.subspace = "All",
                         model.target = "target_traveling.ratio",
                         plot.base.size =  16)




#prepare.paper.figures()

plot.base.size = 10                                
#visualize.transcript.quantifications()
# visualize.traveling.ratio(plot.base.size)
# visualize.traveling.ratio.transcript.expression.distribution(plot.base.size)
# visualize.7SK.binders()
# visualize.features()
