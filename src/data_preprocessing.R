###
# This function loads the Gencode comprehensive gene annotation file for hg19 
# genome build (GrCH37) and subdivides the transcripts into intervals of genomic
# features of the transcripts (i.e. exons, introns, utr etc.) as well as provide 
# a mapping of different ids (ensemble, refseq). 
###
parse.gene.annotations <- function(){
  # Check if a precalculated version exists
  target.file <- paste(INPUT, GENEANNOT, "Gencode/comprehensive_gene_annotation_hg19.RDS", sep="")
  if(!(CARGS$new) & file.exists(target.file)){
    gene.annot <- readRDS(target.file)
    return(gene.annot)
  }
  ## IMPORTANT
  # The following code is deprecated due to a package update, however we provide the 
  # function specific result as a seperate R-Data structure (see Data Availability)
  ##
  if(F){
    # Build list of annotations by annotation type
    annotations <- vector("list", 5)
    names(annotations) <- c("five_prime", 
                            "exons", 
                            "introns", 
                            "coding.exons", 
                            "three_prime")
    
    # Create transcript database object from genecode gff file
    TxDb <- makeTxDbFromGFF(paste(INPUT, "GeneAnnotations/Gencode/", "gencode.v19.annotation.gff3", sep =""))
    genecode.annot <- import(paste(INPUT, "GeneAnnotations/Gencode/", "gencode.v19.annotation.gff3", sep =""))
    
    # Retrieve genes
    genes.data <- transcripts(TxDb, columns='TXNAME')
    # Change naming of col "tx_name" to "transcript id
    names(mcols(genes.data))[1] <- "transcript_id"
    # Append gene_id column
    mcols(genes.data)["gene_id"] <- mcols(transcripts(TxDb, columns='GENEID'))["GENEID"] 
    # Remove version numbers from ids
    mcols(genes.data)["transcript_id"] <- strip.version(mcols(genes.data)[["transcript_id"]], names = T)
    names(genes.data) <- unlist(mcols(genes.data)["gene_id"])
    # Append gene type
    mcols(genes.data)[, "gene_type"] <- mcols(genecode.annot)[ match(names(genes.data), mcols(genecode.annot)[, "gene_id"]),"gene_type"]
    mcols(genes.data)["gene_id"] <-  strip.version(mcols(genes.data)[["gene_id"]])
    # Reformat missing gene type annotations (NA to "none")
    genes.data$gene_type[which(is.na(genes.data$gene_type))] <- "none"
    # Append gene and transcript names
    mcols(genes.data)[, "gene_name"] <-  mcols(genecode.annot)[ match(names(genes.data), mcols(genecode.annot)[, "gene_id"]),"gene_name"]
    mcols(genes.data)[, "transcript_name"] <-  mcols(genecode.annot)[ match(names(genes.data), mcols(genecode.annot)[, "gene_id"]),"transcript_name"]
    
    # Retrieve transcripts
    transcripts.data <- transcripts(TxDb, columns=c("tx_name", "gene_id"))
    # Change naming of col "tx_name" to "transcript id
    names(mcols(transcripts.data))[1] <- "transcript_id"
    # Append transcript type
    mcols(transcripts.data)[, "transcript_type"] <- mcols(genecode.annot)[ match(transcripts.data$transcript_id, mcols(genecode.annot)[, "transcript_id"]),"transcript_type"]
    # Append gene and transcript names
    mcols(transcripts.data)[, "gene_name"] <-  mcols(genecode.annot)[ match(transcripts.data$transcript_id, mcols(genecode.annot)[, "transcript_id"]),"gene_name"]
    mcols(transcripts.data)[, "transcript_name"] <-  mcols(genecode.annot)[ match(transcripts.data$transcript_id, mcols(genecode.annot)[, "transcript_id"]),"transcript_name"]
    # Append gene ids to names attribute
    #names(transcripts.data) <- mcols(transcripts.data)[, "gene_id"]
    # Remove version numbers from ids
    mcols(transcripts.data)["transcript_id"] <- strip.version(mcols(transcripts.data)[["transcript_id"]])
    mcols(transcripts.data)["gene_id"] <- strip.version(mcols(transcripts.data)[["gene_id"]], names = T)
    # Mapping of transcripts to genes
    tx2gene <- unlist(elementMetadata(transcripts.data)$gene_id)
    names(tx2gene) <- elementMetadata(transcripts.data)$transcript_id
    
    ## Extract relevant annotation types
    # Retrieve 3' UTR
    threeUTRs.data <- threeUTRsByTranscript(TxDb, use.names=T)
    names(threeUTRs.data) <- strip.version(names(threeUTRs.data))
    #names(threeUTRs.data) <- tx2gene[names(threeUTRs.data)]
    annotations$three_prime <- threeUTRs.data
    # Retrieve 5' UTR
    fiveUTRs.data <- fiveUTRsByTranscript(TxDb, use.names=T)
    names(fiveUTRs.data) <- strip.version(names(fiveUTRs.data))
    #names(fiveUTRs.data) <- tx2gene[names(fiveUTRs.data)]
    annotations$five_prime <- fiveUTRs.data
    ## Retrieve exons
    exons.data <- exonsBy(TxDb, by='tx', use.names=T)
    names(exons.data) <- strip.version(names(exons.data))
    #names(exons.data) <- tx2gene[names(exons.data)]
    annotations$exons <- exons.data
    ## Retrieve introns
    introns.data <- intronsByTranscript(TxDb, use.names=T)
    ## Append intron identifiers and ranks (analogous to exon_name)
    num.introns.per.tx <- lengths(introns.data)
    intron.ranks <- unlist(lapply(num.introns.per.tx, function(num.introns){
      if(num.introns == 0){return(0)}
      1:num.introns
    }))
    # Exclude transcripts with no intron
    intron.ranks <- intron.ranks[intron.ranks !=0]
    # Append intron id (Analogous to exon id)
    temp <- unlist(introns.data)
    mcols(temp)["intron_id"] <- 1:length(temp)
    mcols(temp)["intron_name"] <- paste(names(temp), intron.ranks, sep = ":" )
    mcols(temp)["intron_rank"] <- intron.ranks
    # Split intron data by transcripts
    introns.data <- split(temp, names(temp))
    # Strip version number from transcripts ids
    names(introns.data) <- strip.version(names(introns.data))
    #names(introns.data) <- tx2gene[names(introns.data)]
    annotations$introns <- introns.data
    ## Retrieve coding sequences
    cds.data <- cdsBy(TxDb, by='tx', use.names=T)
    names(cds.data) <- strip.version(names(cds.data))
    #names(cds.data) <- tx2gene[names(cds.data)]
    annotations$coding.exons <- cds.data
    # Reformat annotations to GRanges
    annotations <- lapply(annotations, function(annot){
      annot.reformatted <- unlist(annot, use.names=FALSE)
      names(annot.reformatted) <- rep(names(annot), elementNROWS(annot))
      annot.reformatted
    })
    # Strip off version numbers from sub region annotations
    annotations <- lapply(annotations, function(annot){
      mcols(annot)[, 2] <- strip.version(mcols(annot)[, 2], pattern = "\\.[0-9]+:", replacement = ":", names = F)
      annot
    })
    # Retrieve genes
    annotations$genes <- genes.data
    # Retrieve transcripts
    annotations$transcripts <- transcripts.data
    # Remove version numbers from all gene ids in each GRanges object
    annotations <- lapply(annotations, function(annot){
      names(annot) <- strip.version(names(annot))
      annot
    })
    names(annotations$transcripts) <- NULL
    names(annotations$genes) <- NULL
    
    
    ## Create mapping of transcript ids to gene ids (along with ref seq ids)
    # Mapping of transcripts to genes
    id.map <- unlist(elementMetadata(annotations$transcripts)$gene_id)
    names(id.map) <- elementMetadata(annotations$transcripts)$transcript_id
    id.map <- as.data.frame(id.map)
    # Append column for transcript id
    id.map$ensembl_transcript_id <- rownames(id.map)
    # Append column for protein id
    id.map$ensembl_protein_id <- genecode.annot$protein_id[match(id.map$ensembl_transcript_id, strip.version(genecode.annot$transcript_id))]
    id.map$ensembl_protein_id <- strip.version(id.map$ensembl_protein_id)
    #  Load refseq ids
    ensemble.refseq.mapping <- read.table(paste(INPUT, "GeneAnnotations/Gencode/", "gencode.v19.metadata.RefSeq",
                                                sep =""), sep = "\t")
    # Delete version numbers from ids
    ensemble.refseq.mapping <- as.data.frame(apply(ensemble.refseq.mapping, 2, function(ids){gsub("\\..*", "", ids)}))
    colnames(ensemble.refseq.mapping) <- c("ensembl_transcript_id", "refseq_transcript_id", "refseq_protein_id")
    ## Extend target df by refseq annotaitons
    id.map$refseq_transcript_id <- ensemble.refseq.mapping$refseq_transcript_id[match(id.map$ensembl_transcript_id, ensemble.refseq.mapping$ensembl_transcript_id)]
    id.map$refseq_protein_id <-ensemble.refseq.mapping$refseq_protein_id[match(id.map$ensembl_transcript_id, ensemble.refseq.mapping$ensembl_transcript_id)]
    colnames(id.map) <- c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_protein_id", "refseq_transcript_id", "refseq_protein_id")
    # Retrieve and append hgnc symbols
    ensembl2hgnc <- read.table(paste(INPUT, "GeneAnnotations/Gencode/", "gencode.v19.metadata.HGNC", sep =""), sep = "\t")
    colnames(ensembl2hgnc) <- c("ensemble_transcript_id", "hgnc")
    ensembl2hgnc$ensemble_transcript_id <- strip.version(ensembl2hgnc$ensemble_transcript_id)
    id.map$hgnc_symbol <- ensembl2hgnc$hgnc[match(id.map$ensembl_transcript_id, ensembl2hgnc$ensemble_transcript_id)]
    # Save id mapping to annotations
    annotations$id_map <- id.map
    
    # Append list of refseq transcript identifier for quick access
    refseq.transcripts <- annotations$id_map$ensembl_transcript_id[!is.na(GENE.ANNOT$id_map$refseq_transcript_id)]
    refseq.supported.ensbl.transcripts <- as.character(refseq.transcripts)
    annotations$refseq_transcript_support <- refseq.supported.ensbl.transcripts
    
    ## Save data sets
    # Define output directory for sub data sets
    out.dir <- paste(INPUT, GENEANNOT, "Gencode", sep="")
    lapply(names(annotations), function(annot){
      saveRDS(annotations[[annot]], paste(out.dir, "/", annot, ".RDS", sep = ""))
      return()
    })
    saveRDS(annotations, target.file)
  }
  return(annotations)
}
###
# This function loads the transcript/gene quantifications from ENCODE for the
# K562 and HepG2 cell lines from different cell fractions and averages the total
# quantifications per transcript and identifies nuclear expressed transcripts 
# based on the cell fraction transcript expressions
###
load.RNAseq <- function(){
  
  # Check if a precalculated version exists
  target.file <- paste(OUTPUT,  "GENE.QUANT.RDS", sep = "")
  if(!(CARGS$new) & file.exists(target.file)){
    GENE.QUANT <- readRDS(target.file)
    return(GENE.QUANT)
  }
  
  # Helper function to filter by valid transcripts
  filter.quantifications <- function(sub.quants){
    
    # Subset quantifications by only those which have a valid ensembl id
    sub.quants$transcript_id <- as.character(sub.quants$transcript_id)
    sub.quants <- sub.quants[startsWith(sub.quants$transcript_id, "ENST"), ]
    # Only retain transcripts which are expressed
    sub.quants <- sub.quants[sub.quants$FPKM > 0, ]
    ## Remove versions from transcript/gene identifiers 
    sub.quants$transcript_id <- strip.version(sub.quants$transcript_id)
    sub.quants$gene_id <- strip.version(sub.quants$gene_id)
    ## Remove transcripts without annotations
    sub.quants <- sub.quants[sub.quants$transcript_id %in% GENE.ANNOT$transcripts$transcript_id, ]
    return(sub.quants)
  }
  
  ## Prepare vector to save transcript quantifications
  cline.quants <- vector("list", 2)
  names(cline.quants) <- CARGS$cell.line
  
  ###Load RNAseq data for the K562 cell line
  ## Load gene quantifications for each cell fraction present in both cell lines
  total.quants.rep1 <- read_tsv(paste(INPUT, RNA, CARGS$cell.line[1], "/ENCFF424CXV.tsv", sep =""))
  total.quants.rep2 <- read_tsv(paste(INPUT, RNA, CARGS$cell.line[1], "/ENCFF073NHK.tsv", sep =""))
  cline.quants$K562 <- list(total.quants.rep1 = total.quants.rep1, 
                            total.quants.rep2 = total.quants.rep2)
  # Filter all quantifications from all cell fractions by valid ones
  cline.quants$K562 <- lapply(cline.quants$K562, filter.quantifications)
  
  ###Load RNAseq data for the HepG2 cell line
  ## Load gene quantifications for each cell fraction present in both cell lines
  total.quants.rep1 <- read_tsv(paste(INPUT, RNA, CARGS$cell.line[2], "/ENCFF205WUQ.tsv", sep =""))
  total.quants.rep2 <- read_tsv(paste(INPUT, RNA, CARGS$cell.line[2], "/ENCFF915JUZ.tsv", sep =""))
  cline.quants$HepG2 <- list(total.quants.rep1 = total.quants.rep1, 
                             total.quants.rep2 = total.quants.rep2)
  # Filter all quantifications from all cell fractions by valid ones
  cline.quants$HepG2 <- lapply(cline.quants$HepG2, filter.quantifications)
  
  # Process quantifications of both cell lines
  quantifications <- 
    lapply(cline.quants, function(quants){
      ## Filter by transcripts which are supported in both replicates and average quantifications
      quantifications <- quants$total.quants.rep1[which(quants$total.quants.rep1$transcript_id %in% quants$total.quants.rep2$transcript_id), ]
      quantifications$FPKM <- colMeans(rbind(quantifications$FPKM, 
                                             quants$total.quants.rep2$FPKM[match(quantifications$transcript_id, quants$total.quants.rep2$transcript_id)]))
      ## Collapse transcripts to genes and take sum of transcript quantificaations as gene quantification
      # Get unique gene ids
      genes <- unique(quantifications$gene_id)
      # Sum up the transcript quantifications for each gene
      clusterExport(WORKERS$avg, c("quantifications"), envir = environment())
      gene.quants <- parSapply(WORKERS$avg, genes, function(gene.id){
        sum(quantifications$FPKM[which(quantifications$gene_id == gene.id)])
      })
      quantifications$gene_FPKM <- 0
      quantifications$gene_FPKM <- gene.quants[match(quantifications$gene_id, names(gene.quants))]
      quantifications$gene_FPKM <- log10(quantifications$gene_FPKM)
      quantifications$orig_FPKM <- quantifications$FPKM
      quantifications$FPKM <- log10(quantifications$FPKM)
      
      # Append transcript type
      quantifications$transcript_type <- GENE.ANNOT$transcripts$transcript_type[match(quantifications$transcript_id, GENE.ANNOT$transcripts$transcript_id)]
      
      # Append variable to label non-coding transcripts
      quantifications$non_coding <- 0
      quantifications[quantifications$transcript_type %in% c("misc_RNA", "lincRNA", 
                                                             "miRNA", "snoRNA", "snRNA"), "non_coding"] <- 1 
      quantifications$non_coding <- factor(quantifications$non_coding)
      return(quantifications)
    })
  names(quantifications) <- CARGS$cell.line
  
  saveRDS(quantifications, target.file)
  return(quantifications)
}
###
# This function loads the cell type specific cage tags, normalizes and clusters these tags
# to yield transcription start sites
###
load.CAGE.tss <- function(ncores = CARGS$workers["avg"]){
  
  # Check if a precalculated version exists
  target.file <- paste(INPUT, CAGE,  "CTSS.RDS", sep = "")
  if(!(CARGS$new) & file.exists(target.file)){
    ctss <- readRDS(target.file)
    return(ctss)
  }
  dir.create(paste(OUTPUT, "CTSS/", sep = ""))
  
  ## IMPORTANT
  # The following code is deprecated due to a package update, however we provide the 
  # function specific result as a seperate R-Data structure (see Data Availability)
  ##
  if(F){
    ctss <- 
      lapply(CARGS$cell.line, function(cline){
        
        # Retrieve ctss annotations
        cage.data <- load(paste(INPUT, "Cage/ENCODEprojectCAGE/data/", cline, ".RData", sep = ""))
        data <- get(cage.data)
        cell.fractions <- names(data)
        
        ## Prepare input for "Cageset" coercion
        data <- lapply(data, function(cline.tss.fraction){
          # Only index columns with tag data
          cols <- colnames(cline.tss.fraction)[-c(1:3)]
          # Shift tag data indices appropriately, since first 3 cols are skipped
          indices <- 1:length(cols) + 3
          # Transform columns with tss tag counts to integers
          cline.tss.fraction[ , indices] <- lapply(indices, function(col.index){
            as.integer(cline.tss.fraction[, col.index])
          })
          # Return data frame
          cline.tss.fraction
        })
        
        ## Coerce data frames to "CAGEsets"
        data <- lapply(data, function(cline.tss.fraction){
          # Coerce data frames to CAGEsets
          cline.tss.fraction <- as(cline.tss.fraction, "CAGEset")
        })
        
        # Normalize data sets
        data <- lapply(data, function(cell.fraction.ctss){
          # Normalize samples
          normalizeTagCount(cell.fraction.ctss,
                            method = "simpleTpm")
          cell.fraction.ctss
        })
        
        # Identify which cell fractions have multiple replicates
        multi.samples <- which(sapply(lapply(data, sampleLabels), length) > 1)
        corr <- lapply(multi.samples, function(cell.fraction.index){
          # Retrieve sample names of replicates
          sample.names <- sampleLabels(data[[cell.fraction.index]])
          # Make a pairwise correlatin plot between all samples per fraction
          sample.corr <- plotCorrelation(data[[cell.fraction.index]], 
                                         samples = sample.names[grepl(".*rep.*",sample.names)], 
                                         method = "pearson")
          # Return the sample correlations per fraction
          sample.corr
        })
        
        ## Merge samples with the highest correlation
        lapply(names(corr), function(cell.fraction){
          print(cell.fraction)
          # Retrieve sample correlations for a particular cell fraction
          cell.fraction.sample.corrs <- corr[[cell.fraction]]
          # Set the diagonal to 0 in order to be able to apply the "max" function
          diag(cell.fraction.sample.corrs) <- 0
          # Identify the indices of samples which have the maximum correlations
          max.corr.index <- which(cell.fraction.sample.corrs == max(cell.fraction.sample.corrs), arr.ind = TRUE)[1, ]
          # Retrieve the sample names based on retrieved indices
          s1 <- row.names(cell.fraction.sample.corrs)[max.corr.index["row"]]
          s2 <- row.names(cell.fraction.sample.corrs)[max.corr.index["col"]]
          # Merge identified samples with maximum correlation and drop all other samples
          sample.labels <- sampleLabels(data[[cell.fraction]])
          merge.index <- 1:length(sample.labels)
          merge.index[which(sample.labels %in% c(s1, s2))] <- max(merge.index)-2+1
          new.labels <- c(sample.labels[-c(which(sample.labels %in% c(s1, s2)))], "merged")
          # Extract CAGEset from the list, so subsequent methods can apply changes
          temp.CAGEset <- data[[cell.fraction]]
          mergeSamples(temp.CAGEset, 
                       mergeIndex = merge.index,
                       mergedSampleLabels = new.labels)
          data[[cell.fraction]] <<- temp.CAGEset
          return()
        })
        
        # Normalize data sets
        data <- 
          lapply(data, function(cell.fraction){
            # Normalize samples
            normalizeTagCount(cell.fraction,
                              method = "simpleTpm")
            cell.fraction
          })
        
        # Create new CAGEset with only merged samples
        data <- lapply(names(data), function(cell.fraction){
          fractional.ctss <- CTSStagCount(data[[cell.fraction]])
          if(!dim(fractional.ctss)[2] > 4){
            colnames(fractional.ctss)[4] <- "merged"
            fractional.ctss$merged <- as.integer(fractional.ctss$merged)
            return(as(fractional.ctss, "CAGEset"))
          }else{
            fractional.ctss <- fractional.ctss[ , c("chr", "pos", "strand", "merged")]
            fractional.ctss$merged <- as.integer(fractional.ctss$merged)
            return(as(fractional.ctss, "CAGEset"))
          }
        })
        
        # Normalize data sets
        data <- 
          lapply(data, function(cell.fraction){
            # Normalize samples
            normalizeTagCount(cell.fraction,
                              method = "simpleTpm")
            cell.fraction
          })
        names(data) <- cell.fractions
        
        log4r::info(LOG, paste("\t Generate clustered ctss data..", cline, sep = ""))
        # Generate clustered ctss data
        ctss <- lapply(names(data), function(cell.fraction){
          
          # Retrieve individual Cageset
          fractional.ctss <- data[[cell.fraction]]
          # Cluster ctss
          clusterCTSS(fractional.ctss, 
                      method = "paraclu", 
                      threshold = 0.01,
                      thresholdIsTpm = T,
                      removeSingletons = TRUE, 
                      keepSingletonsAbove = 0.1,
                      reduceToNonoverlapping = T, 
                      useMulticore = T, 
                      nrCores = ncores)
          # ## Divide signals into two regions, by quantiles
          # At the 5??? end the position of the ???lower??? quantile qLow is determined,
          # which is defined as the point that divides the cluster into two parts,
          # such that the 5??? part contains < qLow * 100% of the CAGE signal of that cluster.
          # Accordingly, position of the ???upper??? quantile qUp is determined near the 3??? end,
          # which is defined as the point that divides the cluster into two parts such
          # that the 5??? part contains >= qUp * 100% of the CAGE signal of that cluster.
          # cumulativeCTSSdistribution(fractional.ctss,
          #                            clusters = "tagClusters",
          #                            useMulticore = T,
          #                            nrCores = ncores)
          # # Determine ctss cluster positions such that 95% of ctss are at the 5'end
          # quantilePositions(fractional.ctss,
          #                   qUp = 0.9,
          #                   clusters = "tagClusters",
          #                   useMulticore = T,
          #                   nrCores = ncores)
          
          # Update CAGEsets
          data[[cell.fraction]] <<- fractional.ctss
          ## Reformat ctss clusters with divisions as genomic ranges object
          tag.clusters <- tagClusters(fractional.ctss, 
                                      "merged")
          
          # Reformat data as GRanges object for efficient range overlap operations
          GRanges(seqnames = tag.clusters$chr,
                  ranges = IRanges(start = tag.clusters$start, tag.clusters$end),
                  strand =  tag.clusters$strand,
                  score = tag.clusters$tpm,
                  tss_count=tag.clusters$nr_ctss,
                  dominant_ctss=tag.clusters$dominant_ctss,
                  tpm_dominant_ctss=tag.clusters$tpm.dominant_ctss,
                  min_density = tag.clusters$min_density,
                  max_density = tag.clusters$max_density
                  #qUpper = tag.clusters$`q_0.9`,
                  #interquantile_width = tag.clusters$interquantile_width
          )
          
        })
        names(ctss) <- cell.fractions
        # Retrieve and reformat individual ctss
        # individual.ctss <- lapply(data, function(cell.fraction.ctss){
        #   cell.fraction.ctss <- CTSSnormalizedTpm(cell.fraction.ctss)
        #   GRanges(seqnames = cell.fraction.ctss$chr,
        #           ranges = IRanges(start = cell.fraction.ctss$pos, end = cell.fraction.ctss$pos),
        #           strand = cell.fraction.ctss$strand,
        #           score = cell.fraction.ctss$merged)
        # })
        
        #return(list(cage.sets = data, clustered.ctss = ctss, individual.ctss <- individual.ctss))
        return(ctss)
      })
    names(ctss) <- CARGS$cell.line
    ## Save data sets
    # Save data set in original CAGEset format
    saveRDS(ctss, target.file)
  }
  return(ctss)

}
###
# This function extract chip-seq experiments from the ENCODE metadata file of 
# all available chip-seq experiments for the K562 and HepG2 cell line.
###
retrieve.CHIPseq.experiments <- function(){
  # For each cell line identify target experiments
  targeted.chip.experiments <- 
    lapply(CARGS$cell.line, function(cline){
      # Load meta data file with file accessions
      metadata.orig <- read.table(paste(INPUT, CHIP, cline, "/bed/metadata_CHIPseq_", cline, ".tsv", sep = ""), header = T, sep = "\t")
      metadata <- metadata.orig
      # Filter for hg19 assambly based files
      metadata <- metadata[which(metadata$Assembly == "hg19"), ]
      # Filter for optimal idr thresholded peaks
      metadata <- metadata[which(metadata$Output.type == "optimal idr thresholded peaks"), ]
      # Filter out perturbation experiments
      metadata <- metadata[which(metadata$Biosample.treatments.duration == ""), ]
      # Filter for released files
      metadata <- metadata[which(metadata$File.Status == "released"), ]
      # Filter for files with replicates
      metadata <- metadata[which(metadata$Biological.replicate.s. == "1, 2"), ]
      # Filter for newer releases
      metadata <- 
        lapply(unique(metadata$Experiment.target), function(target){
          # Retrieve target protein's meta data
          file.meta <- metadata[metadata$Experiment.target==target, ]
          # Take newer release if multiple releases are available
          if(dim(file.meta)[1] != 1){
            release <- as.Date(file.meta$Experiment.date.released)
            file.meta <- file.meta[which.max(release), ]
            return(file.meta)
          }
          return(file.meta)
        })
      metadata <- do.call(rbind, metadata)
      ## Prioritize untagged experiments of target proteins
      # Identify untagged experiments
      targets <- metadata$Experiment.target
      untagged.experiments <- metadata[which(!startsWith(targets,"eGFP" ) & !startsWith(targets,"3xFLAG" )), ]
      # Identify tagged experiments
      tagged.experiments <- metadata[which(startsWith(targets,"eGFP" ) | startsWith(targets,"3xFLAG" )), ]
      tagged.experiments$Experiment.target <- gsub("eGFP-", "", tagged.experiments$Experiment.target)
      tagged.experiments$Experiment.target <- gsub("3xFLAG-", "", tagged.experiments$Experiment.target)
      # Identify tagged experiments for which no untagged experiments exists
      tagged.missing.experiments <- tagged.experiments[(!(tagged.experiments$Experiment.target %in% untagged.experiments$Experiment.target)), ]
      # Filter tagged experiments for newer releases
      tagged.missing.experiments <- 
        lapply(unique(tagged.missing.experiments$Experiment.target), function(target){
          # Retrieve target protein's meta data
          file.meta <- tagged.missing.experiments[tagged.missing.experiments$Experiment.target==target, ]
          # Take newer release if multiple releases are available
          if(dim(file.meta)[1] != 1){
            release <- as.Date(file.meta$Experiment.date.released)
            file.meta <- file.meta[which.max(release), ]
            return(file.meta)
          }
          return(file.meta)
        })
      tagged.missing.experiments <- do.call(rbind, tagged.missing.experiments)
      
      # Aggregate tagged and untagged experiments
      metadata <- rbind(untagged.experiments, tagged.missing.experiments)
      # Rename experiment targets
      metadata$Experiment.target <- gsub("-human", "", metadata$Experiment.target)
      
      # Save metadata file with file accessions
      write.table(metadata, paste0(OUTPUT, "filtered.metadata.chipseq.", cline, ".txt"),  row.names = F, quote = F,col.names = F)
      # Generate a table with file accessions for pubplication purposes
      encode.table <- metadata[, c("Experiment.accession", "File.accession", "Experiment.target" )]
      # Reconstruct original experiment targets names
      encode.table$Experiment.target <- metadata.orig$Experiment.target[match(encode.table$File.accession, metadata.orig$File.accession)]
      # Save encode table
      write.table(encode.table, paste0(OUTPUT, "filtered.encode.chipseq.accessions.table.", cline, ".txt"), row.names = F, quote = F,col.names = F)
      
      #ä Save table with download URLS
      file.access  <- metadata.orig$File.download.URL[match(encode.table$File.accession, metadata.orig$File.accession)]
      write.table(file.access, paste0(OUTPUT, "filtered.encode.chipseq.download.urls.", cline, ".txt"), row.names = F, quote = F, col.names = F)
      
      # Return full paths of files to be included in further analyses
      target.files <- paste(INPUT, CHIP, cline, "/bed/",  metadata$File.accession,".bed.gz", sep = "")
      target.files <- normalizePath(target.files)
      names(target.files) <- metadata$Experiment.target
      return(target.files)
    })
  names(targeted.chip.experiments) <- CARGS$cell.line
  return(targeted.chip.experiments)
}
###
# This function parses the CHIP-seq experiments bed peak files for each 
# cell line (K562, HepG2)
###
parse.CHIPseq.experiments <- function(chipseq.files){
  
  # Check if a precalculated version exists
  target.file <- paste(OUTPUT, "chipseq.peaks.RDS", sep="")
  if(!(CARGS$new) & file.exists(target.file)){
    chip.peaks <- readRDS(target.file)
    return(chip.peaks)
  }
  
  # Specify type of extraCols
  extraCols <- c(V7 = "numeric",
                 V8 = "numeric",
                 V9 = "numeric",
                 V10 = "numeric",
                 V11 = "numeric",
                 V12 = "numeric",
                 V13 = "numeric")
  
  # Parse chip-seq peaks from each cell line as GRanges objects
  chip.peaks <- lapply(CARGS$cell.line, function(cline){
    # Retrieve list of chip-seq experiments in the targeted cell line
    chip.bed.file.paths <- chipseq.files[[cline]]
    # Load peaks for each experiment target 
    chip.peaks <- 
      lapply(names(chip.bed.file.paths), function(dbp){
        # Retrieve bed peak file path for experiment target
        path <- chip.bed.file.paths[dbp]
        ## Load the peak calls
        # Define a helper function to load peak calls and handle exceptions where files cannot be loaded
        peak.parser <- function(path, extraCols) {return(tryCatch(rtracklayer::import(gzfile(path), extraCols = extraCols), error=function(e) NULL))}
        peaks <- peak.parser(path, extraCols)
        if(is.null(peaks)){return(NULL)}
        # Discard peak metadata
        mcols(peaks) <- NULL
        # Annotate the experiment target
        peaks$hgnc_symbol <- dbp
        # Annotate experiment cell line
        peaks$cell_line <- cline
        # Annotate ensemble id
        peaks$ensembl_id <- GENE.ANNOT$id_map$ensembl_gene_id[match(peaks$hgnc_symbol, GENE.ANNOT$id_map$hgnc_symbol)]
        return(peaks)
      })
    chip.peaks <- do.call(c, chip.peaks)
    return(chip.peaks)
  })
  names(chip.peaks) <- CARGS$cell.line
  
  # Save CHIP-seq peak calls
  saveRDS(chip.peaks, target.file)
  return(chip.peaks)
}
###
# This function extract eclip-seq experiments from the encode metadata file of 
# all available chip-seq experiments for the K562 and HepG2 cell line.
###
retrieve.eCLIPseq.experiments <- function(){
  
  # For each cell line identify target experiments
  targeted.clip.experiments <- 
    lapply(CARGS$cell.line, function(cline){
      # Load meta data file with file accessions
      metadata.orig <- read.table(paste(INPUT, eCLIP, cline, "/bed/metadata_eCLIP_", cline, ".tsv", sep = ""), header = T, sep = "\t")
      metadata <- metadata.orig
      # Filter for hg19 assambly based files
      metadata <- metadata[which(metadata$File.assembly == "hg19"), ]
      # Filter out perturbation experiments
      metadata <- metadata[which(is.na(metadata$Biosample.treatments.duration)), ]
      # Filter for released files
      metadata <- metadata[which(metadata$File.Status == "released"), ]
      # Filter for files with replicates
      metadata <- metadata[which(metadata$Biological.replicate.s. == "1, 2"), ]
      # Filter for newer releases
      metadata <- 
        lapply(unique(metadata$Experiment.target), function(target){
          # Retrieve target protein's meta data
          file.meta <- metadata[metadata$Experiment.target==target, ]
          # Take newer release if multiple releases are available
          if(dim(file.meta)[1] != 1){
            release <- as.Date(file.meta$Experiment.date.released)
            file.meta <- file.meta[which.max(release), ]
            return(file.meta)
          }
          return(file.meta)
        })
      metadata <- do.call(rbind, metadata)
      
      # Rename experiment targets
      metadata$Experiment.target <- gsub("-human", "", metadata$Experiment.target)
      # Save metadata file with file accessions
      write.table(metadata, paste0(OUTPUT, "filtered.metadata.eclipseq.", cline, ".txt"),  row.names = F, quote = F,col.names = F)
      # Generate a table with file accessions for publication purposes
      encode.table <- metadata[, c("Experiment.accession", "File.accession", "Experiment.target" )]
      # Save encode table
      write.table(encode.table, paste0(OUTPUT, "filtered.encode.eclipseq.accessions.table.", cline, ".txt"), row.names = F, quote = F,col.names = F)
      
      #ä Save table with download URLS
      file.access  <- metadata.orig$File.download.URL[match(encode.table$File.accession, metadata.orig$File.accession)]
      write.table(file.access, paste0(OUTPUT, "filtered.encode.eclipseq.download.urls.", cline, ".txt"), row.names = F, quote = F, col.names = F)
      
      # Return full paths of files to be included in further analyses
      target.files <- paste(INPUT, eCLIP, cline,"/bed/", metadata$File.accession,".bed.gz", sep = "")
      target.files <- normalizePath(target.files)
      names(target.files) <- metadata$Experiment.target
      return(target.files)
    })
  names(targeted.clip.experiments) <- CARGS$cell.line
  return(targeted.clip.experiments)
}
###
# This function parses the eCLIP-seq experiments bed peak files for each 
# cell line (K562, HepG2)
###
parse.eCLIPseq.experiments <- function(eclipseq.files){
  
  # Check if a precalculated version exists
  target.file <- paste(OUTPUT, "eclipseq.peaks.RDS", sep="")
  if(!(CARGS$new) & file.exists(target.file)){
    chip.peaks <- readRDS(target.file)
    return(chip.peaks)
  }
  
  # Col definitions of eCLIP bed files
  cnames <- c("chrom", "start", "end", "target", 
              "score", "strand", "l2_eClip_VS_SMinput_enrichment", 
              "l10_eClip_VS_SMinput_fisher", "V9", "V10")
  
  # Parse chip-seq peaks from each cell line as GRanges objects
  eclip.peaks <- lapply(CARGS$cell.line, function(cline){
    # Retrieve list of chip-seq experiments in the targeted cell line
    eclip.bed.file.paths <- eclipseq.files[[cline]]
    # Load peaks for each experiment target 
    eclip.peaks <- 
      lapply(names(eclip.bed.file.paths), function(dbp){
        # Retrieve bed peak file path for experiment target
        path <- eclip.bed.file.paths[dbp]
        ## Load the peak calls
        # Define a helper function to load peak calls and handle exceptions where files cannot be loaded
        peak.parser <- function(path, cnames) {return(tryCatch( read.table(gzfile(path), col.names = cnames, header = F), error=function(e) NULL))}
        peaks <- peak.parser(path, cnames)
        if(is.null(peaks)){return(NULL)}
        
        # Rename experiment target
        peaks$target <- gsub("_.*", "", peaks$target)
        # Transform peak calls into a GRanges object
        peaks <- 
          GRanges(seqnames = peaks$chrom, 
                  strand = peaks$strand, 
                  ranges = IRanges(start = peaks$start, end = peaks$end), 
                  hgnc_symbol = peaks$target,
                  cell_line = cline,
                  enrichment = peaks$l2_eClip_VS_SMinput_enrichment,
                  significance = peaks$l10_eClip_VS_SMinput_fisher)
        
        # Discard peak metadata
        mcols(peaks) <- NULL
        # Annotate the experiment target
        peaks$hgnc_symbol <- dbp
        # Annotate experiment cell line
        peaks$cell_line <- cline
        # Annotate ensemble id
        peaks$ensembl_id <- GENE.ANNOT$id_map$ensembl_gene_id[match(peaks$hgnc_symbol, GENE.ANNOT$id_map$hgnc_symbol)]
        return(peaks)
      })
    eclip.peaks <- do.call(c, eclip.peaks)
    return(eclip.peaks)
  })
  names(eclip.peaks) <- CARGS$cell.line
  
  # Save CHIP-seq peak calls
  saveRDS(eclip.peaks, target.file)
  return(eclip.peaks)
}
###
# This function parses a housekeeping gene annotation data set
###
load.house.keeping.genes <- function(){
  
  # Check if a precalculated version exists
  target.file <- paste(OUTPUT, "housekeeping.genes.RDS", sep="")
  if(!(CARGS$new) & file.exists(target.file)){
    housekeeping.genes <- readRDS(target.file)
    return(housekeeping.genes)
  }
  
  # Parse data
  hk.genes <- read.delim(paste(INPUT, HK, "house_keeping.txt", sep = ""), header = F)
  colnames(hk.genes) <- c("hgnc_symbol", "refseq_transcript_id")
  ## Append ensemble transcript ids
  # Match by refseq id
  hk.genes$ensembl_transcript_id <- GENE.ANNOT$id_map$ensembl_transcript_id[match(hk.genes$refseq_transcript_id, GENE.ANNOT$id_map$refseq_transcript_id)]
  # Match by hgnc symbol
  hk.genes$ensembl_transcript_id[which(is.na(hk.genes$ensembl_transcript_id))] <- GENE.ANNOT$id_map$ensembl_transcript_id[match(hk.genes$hgnc_symbol[which(is.na(hk.genes$ensembl_transcript_id))], GENE.ANNOT$id_map$hgnc_symbol)]
  # Save annotation file
  saveRDS(hk.genes, target.file)
  return(hk.genes)
}
###
# This function parses CpG island annotations from UCSC
###
load.cpg.island <- function(){
  
  # Check if a precalculated version exists
  target.file <- paste0(OUTPUT,"cpg.islands.RDS", sep="")
  if(!(CARGS$new) & file.exists(target.file)){
    cpg.island.annotations <- readRDS(target.file)
    return(cpg.island.annotations)
  }
  
  ## CpG island annotations
  # The remaining columns are island length, number of CpGs in the island, the number 
  # of C and G in the island, the percentage of island that is CpG, the percentage of 
  # island that is C or G, and the ratio of observed(cpgNum) to expected(numC*numG/length) CpG in island.
  cnames <- c("chr", "start", "end", "name", "length","cpg.counts","cg.num", "percent.cpg", "percent.cg", "obs.v.exp")
  cpg.island.annotations <- read.table(paste0(INPUT, CPG, "cpgIslandExt.txt"), fill = T, sep = "\t", header = F)
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
  saveRDS(cpg.island.annotations, target.file)
  return(cpg.island.annotations)
}





