##
# This function filters GENCODE annotated transcripts for protein coding, 
# expressed, mappable and CTSS and RefSeq supported transcripts
##
retrieve.valid.coding.transcripts <- function(){
  
  # Check if a precalculated version exists
  target.file <- paste(OUTPUT,  "valid.coding.transcripts.RDS", sep = "")
  if(!(CARGS$new) & file.exists(target.file)){
    valid.transcripts <- readRDS(target.file)
    return(valid.transcripts)
  }
  
  valid.transcripts <- 
    lapply(CARGS$cell.line, function(cline){
      # Prepare target vector of valid transcripts ids that are CTSS supported
      valid.transcripts <- c()
      # Merge ctss data from different cell fractions
      merged.ctss <- CTSS[[cline]]
      names(merged.ctss) <- NULL
      merged.ctss <- do.call(c, merged.ctss)
      
      # Retrieve all protein coding transcripts which are expressed in the specific cell line
      all.tx <- GENE.ANNOT$transcripts
      valid.transcripts <- all.tx
      # Subset by only protein coding transcripts
      valid.transcripts <- valid.transcripts[which(valid.transcripts$transcript_type == "protein_coding")]
      # Subset by those which are RefSeq supported
      valid.transcripts <- valid.transcripts[valid.transcripts$transcript_id %in% GENE.ANNOT$refseq_transcript_support]
      # Subset by those which are expressed
      valid.transcripts <- valid.transcripts[valid.transcripts$transcript_id %in% GENE.QUANT[[cline]]$transcript_id]
      
      ## Subset by CAGE dominant CTSS supported transcripts
      start.range.valid.transcripts <- valid.transcripts
      end(start.range.valid.transcripts) <- start(start.range.valid.transcripts)
      dominant.ctss <- merged.ctss
      start(dominant.ctss) <- merged.ctss$dominant_ctss
      end(dominant.ctss) <- merged.ctss$dominant_ctss
      nearest.tss <- distanceToNearest(valid.transcripts, dominant.ctss)
      valid.clusters <- which(mcols(nearest.tss)[ , "distance"] == 0)
      valid.transcripts <- valid.transcripts[queryHits(nearest.tss[valid.clusters])]
      ctss.clusters <- merged.ctss[subjectHits(nearest.tss[valid.clusters])]
      # Assign cluster ids to each transcript
      valid.transcripts$ctss_cluster_id <- subjectHits(nearest.tss[valid.clusters])
      valid.transcripts$ctss_cluster_dominant_tss <- ctss.clusters$dominant_ctss
      valid.transcripts$ctss_cluster_width <- width(ctss.clusters)
      
      # Average exon width per transcript
      annotated.exons <- GENE.ANNOT$exons
      exon.widths <- width(annotated.exons)
      names(exon.widths) <- names(annotated.exons)
      clusterExport(WORKERS$avg, "exon.widths", envir = environment())
      mean.exon.width.per.transcript <- 
        parLapply(WORKERS$avg, unique(names(exon.widths)), 
                  function(transcript.id){mean(exon.widths[which(names(exon.widths) == transcript.id)])})
      mean.exon.width.per.transcript <- unlist(mean.exon.width.per.transcript)
      names(mean.exon.width.per.transcript) <- unique(names(exon.widths))
      mean.exon.width.per.transcript <- mean.exon.width.per.transcript[names(mean.exon.width.per.transcript) %in% valid.transcripts$transcript_id]
      mcols(valid.transcripts)[match(names(mean.exon.width.per.transcript), valid.transcripts$transcript_id), "avg_exon_width"] <- 
        mean.exon.width.per.transcript
      
      # Retrieve sequences of transcripts
      transcript.sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg19, valid.transcripts)
      # Check if only base letters are present in the transcript
      registerDoParallel(CARGS$workers["low"])
      transcript.letter.ambiguity <-
        foreach(index = 1:length(transcript.sequences)) %dopar% {
          hasOnlyBaseLetters(transcript.sequences[index])
        }
      transcript.letter.ambiguity <- unlist(transcript.letter.ambiguity)
      valid.transcripts <- valid.transcripts[transcript.letter.ambiguity, ]
      registerDoSEQ()
      
      return(valid.transcripts)
    })
  names(valid.transcripts) <- CARGS$cell.line
  ## Save data sets
  # Save data set in original CAGEset format
  saveRDS(valid.transcripts, target.file)
  return(valid.transcripts)
}
##
# This function identifies expressed non-coding transcripts
##
retrieve.valid.non.coding.transcripts <- function(){
  
  # Check if a precalculated version exists
  target.file <- paste(OUTPUT,  "valid.non.coding.transcripts.RDS", sep = "")
  if(!(CARGS$new) & file.exists(target.file)){
    valid.transcripts <- readRDS(target.file)
    return(valid.transcripts)
  }
  
  valid.non.coding.transcripts <- 
    lapply(CARGS$cell.line, function(cline){
      # Retrieve expressed transcripts by cell type specific expression
      expressed.transcripts <- GENE.QUANT[[cline]]
      # Subset by non-coding transcripts
      non.coding.transcripts.idx <- expressed.transcripts$transcript_id[which(expressed.transcripts$non_coding == 1)]
      # Subset by those which are also annotated in GENCODE
      target.transcripts <- GENE.ANNOT$transcripts[match(non.coding.transcripts.idx, GENE.ANNOT$transcripts$transcript_id)]
      # Retrieve cellular fraction in which non-coding transcripts are expressed
      target.transcripts$cellular_fraction <- GENE.QUANT[[cline]]$cellular_fraction[which(expressed.transcripts$non_coding == 1)]
      
      # Average exon width per transcript
      annotated.exons <- GENE.ANNOT$exons
      exon.widths <- width(annotated.exons)
      names(exon.widths) <- names(annotated.exons)
      clusterExport(WORKERS$avg, "exon.widths", envir = environment())
      mean.exon.width.per.transcript <- 
        parLapply(WORKERS$avg, unique(names(exon.widths)), 
                  function(transcript.id){mean(exon.widths[which(names(exon.widths) == transcript.id)])})
      mean.exon.width.per.transcript <- unlist(mean.exon.width.per.transcript)
      names(mean.exon.width.per.transcript) <- unique(names(exon.widths))
      mean.exon.width.per.transcript <- mean.exon.width.per.transcript[names(mean.exon.width.per.transcript) %in% target.transcripts$transcript_id]
      mcols(target.transcripts)[match(names(mean.exon.width.per.transcript), target.transcripts$transcript_id), "avg_exon_width"] <- mean.exon.width.per.transcript
      
      target.transcripts
    })
  names(valid.non.coding.transcripts) <- CARGS$cell.line
  ## Save data sets
  # Save data set in original CAGEset format
  saveRDS(valid.non.coding.transcripts, target.file)
  return(valid.non.coding.transcripts)
}
###
# This function calculates the traveling ratio per valid coding transcript
###
calculate.traveling.ratio <- function(){
  
  # Check if a precalculated version exists
  target.file <- paste(OUTPUT, "traveling_ratio.RDS", sep = "")
  if(!(CARGS$new) & file.exists(target.file)){
    processed.transcripts <- readRDS(target.file)
    return(processed.transcripts)
  }
  
  # Set the window boundaries to define TSS and genic regions
  tss.window.boundaries <- setNames(vector("list", 2), CARGS$cell.line)
  tss.window.boundaries$K562$upstream <-  1
  tss.window.boundaries$K562$downstream <-  1
  tss.window.boundaries$K562$offset <- round((tss.window.boundaries$K562$upstream + tss.window.boundaries$K562$downstream)/2, digits = -2)
  tss.window.boundaries$HepG2$upstream <-  1
  tss.window.boundaries$HepG2$downstream <-  1
  tss.window.boundaries$HepG2$offset <- round((tss.window.boundaries$HepG2$upstream + tss.window.boundaries$HepG2$downstream)/2, digits = -2)
  
  processed.transcripts <- 
    lapply(CARGS$cell.line, function(cline){
      print(cline)
      ## Load Gro-seq data, correct strand and score information and combine data sets
      gro.plus <- import.bw(paste(INPUT, GRO, cline, "/", cline, "_GROseq_plus.bigWig", sep = ""))
      strand(gro.plus) <- "+"
      gro.minus <- import.bw(paste(INPUT, GRO, cline, "/", cline, "_GROseq_minus.bigWig", sep = ""))
      strand(gro.minus) <- "-"
      score(gro.minus) <- abs(score(gro.minus))
      gro.data <- list(gro.plus, gro.minus)
      names(gro.data) <- c("plus", "minus")
      
      # Only calculate traveling ratio for coding transcripts
      coding.transcripts <- valid.coding.transcripts[[cline]]
      # Remove transcripts shorter than 1kb
      #coding.transcripts <- coding.transcripts[width(coding.transcripts) > 1e3]
      # Split valid protein coding transcripts by strand
      coding.transcripts <- split(coding.transcripts, strand(coding.transcripts))
      names(coding.transcripts) <- c("plus", "minus", "undef")
      processed.transcripts <- 
        lapply(c("plus", "minus"), function(strand){
          print(strand)
          # Subset by transcript on a specific strand
          strand.coding.transcripts <- coding.transcripts[[strand]]
          # Remove genes with multiple transcripts
          invalid.genes <- names(which(table(unlist(strand.coding.transcripts$gene_id)) >= 2))
          strand.coding.transcripts <- strand.coding.transcripts[!(unlist(strand.coding.transcripts$gene_id) %in% invalid.genes), ]
          # Make list of transcripts per gene
          transcripts.by.gene <- split(strand.coding.transcripts, unlist(strand.coding.transcripts$gene_id))
          clusterExport(WORKERS$avg, c("gro.data", "cline", "strand", "tss.window.boundaries"), envir = environment())
          processed.transcripts <- 
            parLapply(WORKERS$avg, transcripts.by.gene, function(gene.transcripts){
            #lapply(transcripts.by.gene, function(gene.transcripts){
              ## Set tss region boundaries as identified by heatmap plots
              transcript.tss.region <- gene.transcripts
              end(transcript.tss.region) <- start(transcript.tss.region) + tss.window.boundaries[[cline]]$downstream
              start(transcript.tss.region) <- start(transcript.tss.region) - tss.window.boundaries[[cline]]$upstream + 1

              transcript.body.region <- gene.transcripts
              start(transcript.body.region) <- end(transcript.tss.region) + 1
              #end(transcript.body.region) <-  end(transcript.body.region)+3e3
              ### Get the number of GROseq reads mapping into the tss and body region regions
              ## TSS
              tss.ovs <- findOverlaps(transcript.tss.region, gro.data[[strand]])
              # Includes pseudo count
              transcript.tss.region.reads <- 1 + sum(score(gro.data[[strand]][subjectHits(tss.ovs)]))
              gene.transcripts$tss_region_width <- width(transcript.tss.region)
              gene.transcripts$avg_tss_gro_count <- transcript.tss.region.reads/width(transcript.tss.region)
              # Body
              body.ovs <- findOverlaps(transcript.body.region, gro.data[[strand]])
              transcript.body.region.reads <- 1 + sum(score(gro.data[[strand]][subjectHits(body.ovs)]))
              gene.transcripts$body_region_width <- width(transcript.body.region)
              gene.transcripts$avg_body_gro_count <- transcript.body.region.reads/width(transcript.body.region)
              
              ## Calculate the "traveling ratio", i.e. log2((tss_gro_reads/tss_width)/(body_gro_reads/body_width)) 
              ## quantifying the polymerase pausing
              gene.transcripts$traveling.ratio <- log2((transcript.tss.region.reads/width(transcript.tss.region))/
                                                         (transcript.body.region.reads/width(transcript.body.region)))
              
              # Note if tss OR body region had a zero count, i.e. only the pseudo count
              gene.transcripts$unbalanced <- F
              gene.transcripts$zero_count <- F
              if(transcript.tss.region.reads == 1 | transcript.body.region.reads == 1){
                gene.transcripts$unbalanced <- T
                # Note if tss AND body region had a zero count, i.e. only the pseudo count
                if(transcript.tss.region.reads == 1 & transcript.body.region.reads == 1){
                  gene.transcripts$zero_count <- T
                }
              }
              # Annotate that this gene has only one transcript
              gene.transcripts$singleton <- T
              gene.transcripts$pseudo_singleton <- F
              return(gene.transcripts)
            })
          # Aggrgate results
          names(processed.transcripts) <- NULL
          processed.transcripts <- do.call(base::c, processed.transcripts)
          return(processed.transcripts)
        })
      # Aggregate results
      names(processed.transcripts) <- NULL
      processed.transcripts <- do.call(base::c, processed.transcripts)
      return(processed.transcripts)
    })
  names(processed.transcripts) <- CARGS$cell.line
  
  saveRDS(processed.transcripts, target.file)
  return(processed.transcripts)
}
###
# This function identifies 7SK related non-coding trasncripts
###
load.7SKncRNA.data <- function(){
  
  # Check if a precalculated version exists
  target.file <- paste(OUTPUT,"RN7SK.RDS", sep = "")
  if(!(CARGS$new) & file.exists(target.file)){
    RN7SK <- readRDS(target.file)
    return(RN7SK)
  }
  
  # Get all unique 7SK entries including pseudo transcripts
  # Non-pseudo entry has name "RN7SK" and is on chromosome 6 with range [52860418, 52860748]
  RN7SK.related <- 
    lapply(CARGS$cell.line, function(cline){
      # Subset gene annotations by transcript and gene name matching "7SK"
      by.transcipt.name <- grepl("7SK", GENE.ANNOT$transcripts$transcript_name)
      by.gene.name <- grepl("7SK", GENE.ANNOT$transcripts$gene_name)
      RN7SK.related <- GENE.ANNOT$transcripts[by.transcipt.name & by.gene.name]
      # Subset by those which are expressed 
      RN7SK.related <- RN7SK.related[which(unlist(RN7SK.related$gene_id) %in% GENE.QUANT[[cline]]$gene_id)]
      # Subset by those which are expressed above mean ncRNA expression levels
      RN7SK.related <- RN7SK.related[which(GENE.QUANT[[cline]]$FPKM[match(RN7SK.related$transcript_id, GENE.QUANT[[cline]]$transcript_id)] >= 
                                             median(GENE.QUANT[[cline]]$FPKM[match(valid.non.coding.transcripts[[cline]]$transcript_id, GENE.QUANT[[cline]]$transcript_id)]))]
      # Append column to differentiate between true and pseudo 7SK variants
      mcols(RN7SK.related)["is_RN7SK"] <- F 
      # Indicate "true" 7SK annotation
      mcols(RN7SK.related[which(RN7SK.related$gene_name == "RN7SK")])["is_RN7SK"] <- T
      ## Remove version numbers from ids 
      # Transcript ids
      ids <- RN7SK.related$transcript_id
      names(ids) <- gsub("\\..*" ,"", names(ids))
      RN7SK.related$transcript_id <- gsub("\\..*", "", ids)
      # Gene ids
      ids <- RN7SK.related$gene_id
      names(ids) <- gsub("\\..*" ,"", names(ids))
      RN7SK.related$gene_id <- gsub("\\..*", "", ids)
      return(RN7SK.related)
    })
  # Set cell line identifier
  names(RN7SK.related) <- CARGS$cell.line
  # Save data structure
  saveRDS(RN7SK.related, target.file)
  return(RN7SK.related)
}
###
# This function groups transcripts by bound factors for convenient down-stream analyes
###
get.target.transcripts <- function(peak.calls,
                                   targeted.annot = c("five_prime",
                                                      "introns",
                                                      "exons", 
                                                      "coding.exons",
                                                      "three_prime"),
                                   filter = NULL, 
                                   exon.rank = NULL, 
                                   name = "tx.by.factors"){
  
  # Check if a precalculated version exists
  target.file <- paste(OUTPUT, name, ".RDS", sep = "")
  if(!(CARGS$new) & file.exists(target.file)){
    bindings <- readRDS(target.file)
    return(bindings)
  }
  
  bindings <- 
    lapply(CARGS$cell.line, function(cline){
      # Subset by targeted annotations
      annotations <- GENE.ANNOT[targeted.annot]
      if(!is.null(exon.rank) & "exons" %in% targeted.annot){
        annotations[targeted.annot] <- annotations[[targeted.annot]][annotations[[targeted.annot]]$exon_rank == exon.rank]
      }
      if(!is.null(filter)){
        annotations <- lapply(annotations, function(annot){
          annot[names(annot) %in% filter[[cline]] ]
        })
      }
      # Retrieve cell line specific peak calls
      cline.peak.calls <- peak.calls[[cline]]
      # Split peak calls by targets
      cline.peak.calls <- split(cline.peak.calls, mcols(cline.peak.calls)[ ,"hgnc_symbol"])
      # Overlap annotations with target peak calls
      clusterExport(WORKERS$high, list("annotations"), environment())
      bindings <- 
        parLapply(WORKERS$high, cline.peak.calls, function(rbp.calls){
          bindings <- 
            lapply(annotations, function(annot){
              ovs <- findOverlaps(rbp.calls, annot)
              annot[subjectHits(ovs)]
            })
          names(bindings) <- names(annotations)
          bindings
        })
      names(bindings) <- names(bindings)
      bindings
    })
  names(bindings) <- names(peak.calls)
  # Save intermediate input data
  saveRDS(bindings, target.file)
  return(bindings)
}
###
# This function retrieves the GC content of a transcript
###
get.GC.content <- function(ranges){
  
  # Retrieve sequences of transcripts
  transcript.sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg19, ranges)
  # Check if only base letters are present in the transcript
  registerDoParallel(CARGS$workers["avg"])
  transcript.gc.content <- 
    foreach(index = 1:length(transcript.sequences)) %dopar% {
      gc.freq <- letterFrequency(transcript.sequences[index], letters="CG", as.prob = T)[1,1]
      names(gc.freq) <- NULL
      return(gc.freq)
    }
  transcript.gc.content <- unlist(transcript.gc.content)
  names(transcript.gc.content) <- ranges$transcript_id
  return(transcript.gc.content)
}
###
# This function builds features for predictive models to predict the traveling 
# ratio/expression profiles of transcripts
###
build.feature.vector <- function(targeted.annot = c("five_prime",
                                                    "introns",
                                                    "coding.exons",
                                                    "three_prime"),
                                 feature.set = "feature.vectors"){
  
  # Create target outpute directory which stores model spcific files
  dir.create(paste(OUTPUT, MODEL, "model_data", sep = ""))
  
  # Check if a precalculated version exists
  target.file <- paste(OUTPUT, MODEL,"model_data/", feature.set, ".RDS", sep = "")
  if(!(CARGS$new) & file.exists(target.file)){
    feature.vectors <- readRDS(target.file)
    return(feature.vectors)
  }
  
  feature.vectors <- 
    lapply(CARGS$cell.line, function(cline){
      
      # Feature:Transcrip strand
      transcript.strand <- as.numeric(strand(valid.coding.transcripts[[cline]]))
      # Feature:Transcrip GC Content
      transcript.gc.content <- get.GC.content(valid.coding.transcripts[[cline]])
      # Feature:Transcrip GC Content
      # transcript.terminator.gc.content <- get.GC.content(valid.coding.transcripts[[cline]])
      # transcript.terminator.regions <- valid.coding.transcripts[[cline]]
      # ranges(transcript.terminator.regions) <- ranges(GENE.ANNOT$three_prime[match(valid.coding.transcripts[[cline]]$transcript_id, names(GENE.ANNOT$three_prime)),])

      # Feature:Transcrip Length
      transcript.length <- rescale(width(valid.coding.transcripts[[cline]]))
      ## Feature:Number of exons per transcript
      targeted.transcripts <- GENE.ANNOT$transcripts$transcript_id[match(valid.coding.transcripts[[cline]]$transcript_id, 
                                                                         mcols(GENE.ANNOT$transcripts)[ ,"transcript_id"])]
      exon.counts <- table(names(GENE.ANNOT$exons))
      transcript.exon.count <- rescale(as.vector(exon.counts[targeted.transcripts]))
      # Feature:Ratio of exon count to transcript length
      exon.count.ratio <- rescale(transcript.length/transcript.exon.count)
      ## Feature:Fraction of exonic sequence
      clusterExport(WORKERS$high, c("GENE.ANNOT"), envir = environment())
      build.feature.exonic.share <- function(){
        targeted.transcripts <- GENE.ANNOT$transcripts$transcript_id[match(valid.coding.transcripts[[cline]]$transcript_id, 
                                                                           mcols(GENE.ANNOT$transcripts)[ ,"transcript_id"])]
        exonic.share <- 
          parSapply(WORKERS$high, targeted.transcripts, function(transcript){
            transcript.exons <- GENE.ANNOT$exons[which(names(GENE.ANNOT$exons) == transcript), ]
            cum.exonic.sequence.length <- sum(width(transcript.exons))
            transcript.length <- width(GENE.ANNOT$transcripts[GENE.ANNOT$transcripts$transcript_id == transcript])
            exonic.share <- cum.exonic.sequence.length/transcript.length
          })
        exonic.share
      }
      exonic.share <- build.feature.exonic.share()
      # Feature:Average exon width
      avg.exon.width <- rescale(valid.coding.transcripts[[cline]]$avg_exon_width)
      # Feature:Genomic location of the transcript on the genome
      transcript.localization <- rescale(start(valid.coding.transcripts[[cline]]))
      # Feature:Chromosomal location of the transcript
      transcript.chromosome.localization <- as.numeric(seqnames(valid.coding.transcripts[[cline]]))
      ## Feature:RBP bindings on transcript sub regions (exon, intron etc.)
      build.feature.rbp.bindings <- function(){
        # Prepare indicator matrix
        rbp.bindings <- setNames(data.frame(matrix(0, ncol = length(names(tx.by.rbp[[cline]]))*length(targeted.annot), 
                                                   nrow = length(valid.coding.transcripts[[cline]]))), 
                                 paste("clip", rep(names(tx.by.rbp[[cline]]), each = length(targeted.annot)), targeted.annot, sep = "."))
        rownames(rbp.bindings) <- valid.coding.transcripts[[cline]]$transcript_id
        # Fill indicator matrix per sub region
        for(rbp in names(tx.by.rbp[[cline]])){
          for(sub.region in names(tx.by.rbp[[cline]][[rbp]])){
            if(sub.region %in% targeted.annot){
              bound.transcripts <- unique(names(tx.by.rbp[[cline]][[rbp]][[sub.region]]))
              bound.transcripts <- bound.transcripts[bound.transcripts %in% valid.coding.transcripts[[cline]]$transcript_id]
              rbp.bindings[bound.transcripts, paste("clip", rbp, sub.region, sep = ".")] <- 1
            }
          }
        }
        rbp.bindings <- apply(rbp.bindings, 2, as.factor)
      }
      rbp.bindings <- build.feature.rbp.bindings()
      ## Feature:DBP bindings on transcript sub regions (exon, intron etc.)
      build.feature.dna.bindings <- function(){
        # Prepare indicator matrix
        dna.bindings <- setNames(data.frame(matrix(0, ncol = length(names(tx.by.dbp[[cline]]))*length(targeted.annot), 
                                                   nrow = length(valid.coding.transcripts[[cline]]))), 
                                 paste("chip", rep(names(tx.by.dbp[[cline]]), 
                                                   each = length(targeted.annot)), targeted.annot , sep = "."))
        rownames(dna.bindings) <- valid.coding.transcripts[[cline]]$transcript_id
        for(dbp in names(tx.by.dbp[[cline]])){
          for(sub.region in names(tx.by.dbp[[cline]][[dbp]])){
            if(sub.region %in% targeted.annot){
              bound.transcripts <- unique(names(tx.by.dbp[[cline]][[dbp]][[sub.region]]))
              bound.transcripts <- bound.transcripts[bound.transcripts %in% rownames(dna.bindings)]
              dna.bindings[bound.transcripts, paste("chip", dbp, sub.region, sep = ".")] <- 1
            }
          }
        }
        dna.bindings <- apply(dna.bindings, 2, factor)
      }
      dna.bindings <- build.feature.dna.bindings()
      # Remove POLRII signals from dna bindings since it should be correlated with the targets quite naturally
      dna.bindings <- dna.bindings[, -which(grepl( "POLR2", colnames(dna.bindings)))]
      # ## Feature:Expression of nearest RN7SK
      # nearest.RN7SK.quant <- data.frame(matrix(0, ncol = 1, nrow = length(valid.coding.transcripts[[cline]])), 
      #                                   row.names = valid.coding.transcripts[[cline]]$transcript_id)
      # colnames(nearest.RN7SK.quant) <- "nearest.RN7SK.quant"
      # # Retreive transcript ranges
      # tx.ranges <- valid.coding.transcripts[[cline]]
      # # Retreive the nearest RN7SK annotations per transcript
      # nearest.RN7SK.quant.annotations <- nearest(tx.ranges, RN7SK.ANNOT[[cline]])
      # names(nearest.RN7SK.quant.annotations) <- tx.ranges$transcript_id
      # nearest.RN7SK.quant.annotations <- nearest.RN7SK.quant.annotations[!is.na(nearest.RN7SK.quant.annotations)]
      # nearest.RN7SK.quant[names(nearest.RN7SK.quant.annotations), "nearest.RN7SK.quant"] <- 
      #   GENE.QUANT[[cline]]$FPKM[match(RN7SK.ANNOT[[cline]][nearest.RN7SK.quant.annotations]$transcript_id, GENE.QUANT[[cline]]$transcript_id)]
      # ## Feature:Distance to nearest 7SK
      # distance.nearest.RN7SK <- data.frame(matrix(-1, ncol = 1, nrow = length(valid.coding.transcripts[[cline]])), 
      #                                      row.names = valid.coding.transcripts[[cline]]$transcript_id)
      # colnames(distance.nearest.RN7SK) <- "distance.nearest.RN7SK"
      # names(tx.ranges) <- tx.ranges$transcript_id
      # distance.nearest.RN7SK[names(nearest.RN7SK.quant.annotations), "distance.nearest.RN7SK"] <-
      #   distance(Pairs(tx.ranges[names(nearest.RN7SK.quant.annotations)], 
      #                  RN7SK.ANNOT[[cline]][nearest.RN7SK.quant.annotations]))
      
      ## Feature:Distance to and quantification of n nearest lincRNAs
      linc.RNA.transcripts <- valid.non.coding.transcripts[[cline]]
      names(linc.RNA.transcripts) <- mcols(linc.RNA.transcripts)[, "transcript_id"]
      # Exclude 7SK annotatons from lincRNAs, as these are seperate predictors
      linc.RNA.transcripts <- linc.RNA.transcripts[!names(linc.RNA.transcripts) %in% as.vector(RN7SK.ANNOT[[cline]]$transcript_id)]
      nearest.linc.RNA.quant <- data.frame(matrix(0, ncol = 1, nrow = length(valid.coding.transcripts[[cline]])), 
                                           row.names = valid.coding.transcripts[[cline]]$transcript_id)
      colnames(nearest.linc.RNA.quant) <- "nearest.linc.RNA.quant"
      tx.ranges <- valid.coding.transcripts[[cline]]
      num.of.nearest.lincRNAs <- 2
      clusterExport(WORKERS$high, list("GENE.QUANT","tx.ranges", "num.of.nearest.lincRNAs", 
                                       "linc.RNA.transcripts", "cline"),envir = environment())
      nearest.lincRNA.metrics <- 
        parLapply(WORKERS$high, seq_along(tx.ranges), function(index){
          transcript <- tx.ranges[index, ]
          distances.to.linc.RNAs <- distance(transcript, linc.RNA.transcripts)
          names(distances.to.linc.RNAs) <- linc.RNA.transcripts$transcript_id
          distances.to.linc.RNAs <- distances.to.linc.RNAs[!is.na(distances.to.linc.RNAs)]
          distances.to.linc.RNAs <- sort(distances.to.linc.RNAs, decreasing = F)
          nearest.linc.RNAs.quants <- GENE.QUANT[[cline]]$FPKM[match(names(distances.to.linc.RNAs), GENE.QUANT[[cline]]$transcript_id)]
          n.nearest.linc.RNA.quant <- vector("numeric", num.of.nearest.lincRNAs)
          n.nearest.linc.RNA.quant[1:num.of.nearest.lincRNAs] <- -1
          n.nearest.linc.RNA.quant[1:num.of.nearest.lincRNAs] <- nearest.linc.RNAs.quants[1:num.of.nearest.lincRNAs]
          names(n.nearest.linc.RNA.quant) <- names(distances.to.linc.RNAs)[1:num.of.nearest.lincRNAs]
          data.frame(id = names(n.nearest.linc.RNA.quant)[1:num.of.nearest.lincRNAs], 
                     quant = n.nearest.linc.RNA.quant[1:num.of.nearest.lincRNAs],
                     distance = distances.to.linc.RNAs[1:num.of.nearest.lincRNAs], stringsAsFactors = F, row.names = NULL)
        })
      
      # ncRNA IDs
      nearest.lincRNA.ids <- 
        lapply(nearest.lincRNA.metrics, function(nearest.lincRNAs){nearest.lincRNAs[, "id"]})
      nearest.lincRNA.ids <- do.call(rbind, nearest.lincRNA.ids)
      colnames(nearest.lincRNA.ids) <- paste("Proximal.ncRNA.", 1:num.of.nearest.lincRNAs, ".ID", sep = "")
      # ncRNA Quantifications
      nearest.lincRNA.quant <- 
        lapply(nearest.lincRNA.metrics, function(nearest.lincRNAs){nearest.lincRNAs[, "quant"]})
      nearest.lincRNA.quant <- do.call(rbind, nearest.lincRNA.quant)
      colnames(nearest.lincRNA.quant) <- paste("Proximal.ncRNA.", 1:num.of.nearest.lincRNAs, ".Quant", sep = "")
      # # ncRNA Distances
      # nearest.lincRNA.dist <- 
      #   lapply(nearest.lincRNA.metrics, function(nearest.lincRNAs){rescale(nearest.lincRNAs[, "distance"])})
      # nearest.lincRNA.dist <- do.call(rbind, nearest.lincRNA.dist)
      # colnames(nearest.lincRNA.dist) <- paste("Proximal.ncRNA.", 1:num.of.nearest.lincRNAs, ".Distance", sep = "")
      
      # Feature:Chip signals on ncRNAs
      build.feature.ncRNA.dna.bindings <- function(ncRNAs, index){
        # Prepare indicator matrix
        dna.bindings <- setNames(data.frame(matrix(0, ncol = length(names(tx.by.dbp[[cline]]))*length(targeted.annot),
                                                   nrow = length(valid.coding.transcripts[[cline]]))),
                                 paste("chip", rep(names(tx.by.dbp[[cline]]),
                                                   each = length(targeted.annot)), targeted.annot, "Proximal.ncRNA", index, sep = "."))
        for(dbp in names(tx.by.dbp[[cline]])){
          # Fill indicator matrix per sub region
          for(sub.region in names(tx.by.dbp[[cline]][[dbp]])){
            if(sub.region %in% targeted.annot){
              bound.transcripts <- unique(names(tx.by.dbp[[cline]][[dbp]][[sub.region]]))
              bound.transcripts <- bound.transcripts[bound.transcripts %in% ncRNAs]
              if(length(bound.transcripts)>0){
                dna.bindings[match(bound.transcripts, ncRNAs), paste("chip", dbp, sub.region,"Proximal.ncRNA", index,  sep = ".")] <- 1
              }
            }
          }
        }
        dna.bindings <- apply(dna.bindings, 2, factor)
        rownames(dna.bindings) <- valid.coding.transcripts[[cline]]$transcript_id
        dna.bindings
      }
      ncRNA.dbp.bindings <-
        lapply(1:dim(nearest.lincRNA.ids)[2], function(index){
          ncRNA.ids <- nearest.lincRNA.ids[, index]
          build.feature.ncRNA.dna.bindings(ncRNA.ids, index)
        })
      ncRNA.dbp.bindings <- do.call(cbind, ncRNA.dbp.bindings)
      exclude <- which(grepl( "POLR2", colnames(ncRNA.dbp.bindings)))
      if(length(exclude)>0){
        ncRNA.dbp.bindings <- ncRNA.dbp.bindings[, -exclude]
      }
      
      # Feature:Clip signals on ncRNAs
      build.feature.ncRNA.rbp.bindings <- function(ncRNAs, index){
        # Prepare indicator matrix
        rbp.bindings <- setNames(data.frame(matrix(0, ncol = length(names(tx.by.rbp[[cline]]))*length(targeted.annot),
                                                   nrow = length(valid.coding.transcripts[[cline]]))),
                                 paste("clip", rep(names(tx.by.rbp[[cline]]), each = length(targeted.annot)), targeted.annot,"Proximal.ncRNA", index, sep = "."))
        rownames(rbp.bindings) <- valid.coding.transcripts[[cline]]$transcript_id
        # Fill indicator matrix per sub region
        for(rbp in names(tx.by.rbp[[cline]])){
          for(sub.region in names(tx.by.rbp[[cline]][[rbp]])){
            if(sub.region %in% targeted.annot){
              bound.transcripts <- unique(names(tx.by.rbp[[cline]][[rbp]][[sub.region]]))
              bound.transcripts <- bound.transcripts[bound.transcripts %in% ncRNAs]
              if(length(bound.transcripts)>0){
                rbp.bindings[match(bound.transcripts, ncRNAs), paste("clip", rbp, sub.region,"Proximal.ncRNA" ,index, sep = ".")] <- 1
              }
            }
          }
        }
        rbp.bindings <- apply(rbp.bindings, 2, as.factor)
      }
      ncRNA.rbp.bindings <-
        lapply(1:dim(nearest.lincRNA.ids)[2], function(index){
          print(index)
          ncRNA.ids <- nearest.lincRNA.ids[, index]
          build.feature.ncRNA.rbp.bindings(ncRNA.ids, index)
        })
      ncRNA.rbp.bindings <- do.call(cbind, ncRNA.rbp.bindings)
      
      # Feature:TSS Cluster width
      tss.cluster.width <- rescale(valid.coding.transcripts[[cline]]$ctss_cluster_width)
      
      # Feature:TSS cluster AT content
      merged.ctss <- CTSS[[cline]]
      names(merged.ctss) <- NULL
      merged.ctss <- do.call(c, merged.ctss)
      ctss.cluster.ranges <- merged.ctss[valid.coding.transcripts[[cline]]$ctss_cluster_id]
      ctss.cluster.ranges$transcript_id <- valid.coding.transcripts[[cline]]$transcript_id
      ctss.cluster.at.content <- 1-get.GC.content(ctss.cluster.ranges)
      
      ## Housekeeping gene annotations
      hk.gene <- revalue(factor(valid.coding.transcripts[[cline]]$transcript_id %in% HK.GENES$ensembl_transcript_id), c("TRUE"=1, "FALSE"=0))
      
      ## Retrieve cpg island annotations
      tx.ranges <- valid.coding.transcripts[[cline]]
      dist.nearest.cpg.island <- distanceToNearest(tx.ranges, CPG.ISLANDS)
      nearest.cpg.islands <- CPG.ISLANDS[subjectHits(dist.nearest.cpg.island)]
      cpg.island.dist <- rescale(dist.nearest.cpg.island@elementMetadata$distance)
      cpg.island.length <- rescale(nearest.cpg.islands$length)
      cpg.island.count <- rescale(nearest.cpg.islands$cpg.counts)
      cpg.island.percent.cpg <- nearest.cpg.islands$percent.cpg
      cpg.island.percent.cg <- nearest.cpg.islands$percent.cg
      cpg.island.cg.num <-  nearest.cpg.islands$cg.num
      cpg.island.exp.obs <- nearest.cpg.islands$obs.v.exp
      
      # Aggregate features into a feature vector
      feature.vectors <- data.frame(housekeeping = hk.gene,
                                    cpg.island.dist = cpg.island.dist,
                                    cpg.island.length = cpg.island.length,
                                    cpg.island.count = cpg.island.count,
                                    cpg.island.percent.cpg = cpg.island.percent.cpg,
                                    cpg.island.percent.cg = cpg.island.percent.cg,
                                    cpg.island.exp.obs = cpg.island.exp.obs,
                                    tx.strand = transcript.strand,
                                    tx.gc.seq = transcript.gc.content,
                                    tx.len = transcript.length,
                                    tx.ex.count = transcript.exon.count,
                                    tx.ex.ratio = exon.count.ratio,
                                    tx.ex.width = avg.exon.width,
                                    tx.ex.seq = exonic.share,
                                    tx.loc = transcript.localization,
                                    tx.chr.loc = transcript.chromosome.localization,
                                    tx.tss.width = tss.cluster.width,
                                    tx.tss.at.cont = ctss.cluster.at.content,
                                    rbp.bindings,
                                    dna.bindings,
                                    nearest.lincRNA.ids,
                                    ncRNA.dbp.bindings,
                                    ncRNA.rbp.bindings)
      
      # Convert CHIP/eCLIP features to factors
      cnames <- colnames(feature.vectors)[grepl("^chip|clip", colnames(feature.vectors))]
      feature.vectors[cnames] <- lapply(feature.vectors[cnames], factor)
      
      return(feature.vectors)
    })
  
  # Set names of design matrices to cell line names
  names(feature.vectors) <- CARGS$cell.line
  # Save data sets
  saveRDS(feature.vectors, target.file)
  return(feature.vectors)
}
###
# This function removes features from a feature vector
###
exclude.features <- function(feature.vectors, patterns){
  
  ## Exclude features which correspond to given patterns
  feature.vectors <- 
    lapply(feature.vectors, function(feature.vector){
      exclude <- 
        lapply(patterns, function(p){
          exclude <- grep(p,colnames(feature.vector))
        })
      exclude <- unique(unlist(exclude))
      if(length(exclude)>0){
        feature.vector <- feature.vector[, -exclude]
      }
      feature.vector
    })
  return(feature.vectors)
}
### 
# This function builds sub feature sets based on prior knowledge about the 7SK ncRNA
###
build.pausing.associated.factor.sub.feature.spaces <- function(){
  
  # Convert ensemble ids to hgnc symbols
  convert.ids <- function(ids, reverse = F){
    if(reverse){
      return(as.character(GENE.ANNOT$id_map$ensembl_gene_id[match(ids, GENE.ANNOT$id_map$hgnc_symbol)]))
    }
    return(as.character(GENE.ANNOT$id_map$hgnc_symbol[match(ids, GENE.ANNOT$id_map$ensembl_gene_id)]))
  }
  
  # # Function to retrieve genes which share a domain with query genes
  # sdom.factors <- function(targets){
  #   shared.domain.factors <- 
  #     lapply(targets, function(factor){
  #       # Retrieve factor's domains
  #       domains <- SDOMAINS$members[[factor]]
  #       # For each factor's domain, retrieve associated factors
  #       related.factors <- 
  #         sapply(domains, function(domain){
  #           sdom.all.targets[sdom.all.targets %in% SDOMAINS$groups[[domain]]]
  #         })
  #       unique(unlist(related.factors, use.names = F))
  #     }) 
  #   shared.domain.factors <- c(shared.domain.factors, targets) %>% unlist() %>% unique()
  # }
  
  # Exclude factors from set B from set A
  set.diff.factors <- function(A, B){
    new.list <- 
      lapply(names(A), function(element){
        setdiff(A[[element]], B[[element]])
      })
    names(new.list) <- names(A)
    new.list
  }
  
  ## Retrieve all targets from all assays
  clip.targets <- lapply(eCLIPseq, function(peaks){
    as.character(unique(peaks$ensembl_id))
  })
  chip.targets <- lapply(CHIPseq, function(peaks){
    as.character(unique(peaks$ensembl_id))
  })
  all.targets  <- lapply(CARGS$cell.line, function(cline){
    all.targets <- unique(c(clip.targets[[cline]], chip.targets[[cline]]))
    all.targets[!is.na(all.targets)]
  })
  names(all.targets) <- CARGS$cell.line
  full <- all.targets
  all.targets <- unique(unlist(all.targets))
  
  # # Remove RNA bninding containg motif grouped shared domain factos group
  # SDOMAINS$groups$`725` <- c()
  # 
  # # Identify factors which have shared domains
  # sdom.all.targets <- all.targets[all.targets %in% names(SDOMAINS$members)]
  
  ##  Define (A)  the set of known pausing factors from literature
  known.pausing.factors <- c("NELFE", "SUPT5H", 
                             "MLLT1", "LARP7", 
                             "BRD4", 
                             "MYC", "TAF1", 
                             "TBP", "PAF1", 
                             "SUPT16H", "SUPT6H",
                             "SUPT4H1", "NELFB", 
                             "NELFA", "NELFCD", 
                             "CDK9", "HEXIM1", 
                             "HEXIM2", "MEPCE",
                             "CCNT1", "CCNT2", 
                             "ELL", "ELL2", 
                             "ELL3", "AFF1",
                             "AFF4", "MLLT1", 
                             "MLLT3")
  known.pausing.factors.ids <- c("ENSG00000204356", "ENSG00000196235", 
                                 "ENSG00000130382", "ENSG00000174720",
                                 "ENSG00000141867", 
                                 "ENSG00000136997", "ENSG00000147133", 
                                 "ENSG00000112592", "ENSG00000006712",
                                 "ENSG00000092201", "ENSG00000109111", 
                                 "ENSG00000213246", "ENSG00000188986",
                                 "ENSG00000185049", "ENSG00000101158", 
                                 "ENSG00000136807", "ENSG00000186834",
                                 "ENSG00000168517", "ENSG00000146834",
                                 "ENSG00000129315", "ENSG00000082258",
                                 "ENSG00000105656", "ENSG00000118985",
                                 "ENSG00000128886", "ENSG00000172493",
                                 "ENSG00000072364", "ENSG00000130382",
                                 "ENSG00000171843")
  A <- list(K562 = known.pausing.factors.ids,
            HepG2 = known.pausing.factors.ids)    
  
  # ## Define (B) the set of factors which share a binding domain with known pausing factors
  # B <- sdom.factors(known.pausing.factors.ids)
  # B <- list(K562 = B, HepG2 = B)
  
  ## Define (C) the set of true 7SK binders
  rn7sk.targets <- 
    lapply(CARGS$cell.line, function(cline){
      as.character(unique(novel.rn7sk.binders[[cline]]$ensembl_id[novel.rn7sk.binders[[cline]]$is_RN7SK]))
    })
  names(rn7sk.targets) <- CARGS$cell.line
  C <- rn7sk.targets
  
  ## Define (D) the set of all true and pseudo 7SK binding targets
  pseudo.rn7sk.targets <- 
    lapply(CARGS$cell.line, function(cline){
      as.character(unique(novel.rn7sk.binders[[cline]]$ensembl_id[!novel.rn7sk.binders[[cline]]$is_RN7SK]))
    })
  names(pseudo.rn7sk.targets) <- CARGS$cell.line
  D <- combine.lists(C, pseudo.rn7sk.targets)
  
  ## Define (E) the set of all true 7SK targets and shared domain associated factors
  #E <- lapply(C, sdom.factors)
  
  ## Define (F) the set of true and pseudo 7SK binding targets and shared domain factors
  #FF <- lapply(D, sdom.factors)
  
  ## Define (G) the set of pausing related factors, i.e. known, known sdom, 7SK binding, pseudo 7SK binding and all their sdom targets
  #G <- combine.lists(B, FF)
  
  ## Define (H) the set of all factors excluding pausing associated factors
  #H <- set.diff.factors(full, G)
  
  #  Define (I) the set of all known and novel pausing factors excluding sdom factors
  I <- combine.lists(A, D)
  
  # Remove 7SK binders from established pausing factor's shared domain factors
  #B <- set.diff.factors(B, FF)
  
  # Remove established pausing factors from 7SK binder's associated shared domain factors
  #E <- set.diff.factors(E, B)
  #FF <- set.diff.factors(FF, B)
  
  ## Convert transcript ids to hgnc symbols for subsetting feature matrices
  full <-  lapply(full, convert.ids)
  A <-  lapply(A, convert.ids)
  #B <-  lapply(B, convert.ids)
  C <-  lapply(C, convert.ids)
  D <-  lapply(D, convert.ids)
  #E <-  lapply(E, convert.ids)
  #FF <-  lapply(FF, convert.ids)
  #G <-  lapply(G, convert.ids)
  #H <-  lapply(H, convert.ids)
  I <-  lapply(I, convert.ids)
  
  feature.sets <- list(All = full,
                       Known = A,
                       #Known.sdom = B,
                       RN7SK = C,
                       all.RN7SK = D,
                       #RN7SK.sdom= E,
                       #all.RN7SK.sdom = FF,
                       #pausing.sdom = G,
                       #non.pausing= H,
                       pausing = I)
  
  names(feature.sets) <- c("All", 
                           "Known", 
                           #"Known.sdom", 
                           "7SK.Binding*", 
                           "7SK.Binding", 
                           #"7SK.sdom", 
                           #"All.7SK.sdom", 
                           #"Pausing.sdom", 
                           #"Non-Pausing", 
                           "Pausing" )
  
  ## (K) Random sets of factors
  #size.FF <- lengths(FF)
  #size.G <- lengths(G)
  # selected.feature.sets <- list(known.pausing.factors_setA = A,
  #                               known.pausing.w.sdom.factors_setB = B,
  #                               true.7sk.binders_setC = C,
  #                               true.pseudo.7sk.binders_setD = D,
  #                               true.7sk.binders.w.sdom_setE = E,
  #                               true.pseudo.7sk.binders.w.sdom_setFF = FF,
  #                               known.pausing.true.pseudo.7sk.w.sdom_setG = G)
  # random.samples<- 
  #   lapply(names(selected.feature.sets), function(fset){
  #     random.samples <- 
  #       lapply(1:100, function(iteration){
  #         cline.set <- 
  #           lapply(CARGS$cell.line, function(cline){
  #             available.factors <- setdiff(full[[cline]], selected.feature.sets[[fset]][[cline]])
  #             sample(available.factors, size = length(selected.feature.sets[[fset]][[cline]]), replace = F)
  #           })
  #         names(cline.set) <- CARGS$cell.line
  #         cline.set
  #       })
  #     names(random.samples) <- as.character(paste0(fset,":random",":",1:100))
  #     feature.sets <<- c(feature.sets, random.samples)
  #     #random.samples
  #     return()
  #   })
  return(feature.sets)
}
### 
# This function identifies sequence specific factors
###
retrieve.sequence.specific.factors <- function(){
  
  # Check if a precalculated version exists
  target.file <- paste0(OUTPUT,"factor.sequence.specificity.RDS", sep="")
  if(!(CARGS$new) & file.exists(target.file)){
    factor.sequence.specificity <- readRDS(target.file)
    return(factor.sequence.specificity)
  }
  
  ## Retrieve all available DNA binding and RNA binding factors from both cell lines
  # Retrieve all DNA binding factors
  dbp.factors <- 
    lapply(CHIPseq, function(peaks){
      unique(peaks$hgnc_symbol)
    })
  # Retrieve all RNA binding factors
  rbp.factors <- 
    lapply(eCLIPseq, function(peaks){
      unique(peaks$hgnc_symbol)
    })
  # Aggregate all factors
  factors <- unique((c(unlist(rbp.factors), unlist(dbp.factors))))
  
  ## Retrieve data on sequence specific binding of factors based on external databases
  # Msig database
  data(tftColl)
  data(tftCollMap)
  TFs_MSIG = TFCatalog(name="MsigDb.TFT", nativeIds=names(tftColl),
                       HGNCmap=data.frame(tftCollMap,stringsAsFactors=FALSE))
  msig.factors <- HGNCmap(TFs_MSIG)[, "tftname"]
  msig.factors <- unique(c(msig.factors, HGNCmap(TFs_MSIG)[, "hgnc.heur"]))
  # CisBP database
  data(cisbpTFcat)
  TFs_CISBP = TFCatalog(name="CISBP.info", nativeIds=cisbpTFcat[,1],
                        HGNCmap = cisbpTFcat)
  cisbp.factors <- HGNCmap(TFs_CISBP)[, "HGNC"]
  ## Hocomoco database
  data(hocomoco.mono.sep2018)
  TFs_HOCO = TFCatalog(name="hocomoco11", nativeIds=hocomoco.mono.sep2018[,1],
                       HGNCmap=hocomoco.mono.sep2018)
  hocomoco.factors <- HGNCmap(TFs_HOCO)[, "HGNC"]
  # Aggregate all sequence specific factors from various sources (msigdb, hocomoco, cisbp)
  all.sequence.specific.factors <- unique(c(msig.factors, cisbp.factors, hocomoco.factors))
  
  ## Build matrix of factor categorizations
  factor.classes <- data.frame(factor = factors)
  # Classifiy by target sequence binding
  factor.classes$binding_sequence[match(unique(unlist(dbp.factors)), factor.classes$factor)] <- "dbp"
  factor.classes$binding_sequence[match(unique(unlist(rbp.factors)), factor.classes$factor)] <- "rbp"
  # Annotate whether sequence specific binding
  factor.classes$sequence_specific = 0
  factor.classes$sequence_specific[which(factor.classes$factor %in% all.sequence.specific.factors)] = 1
  # Subdivide into two groups, sequence specific and non-sequence specific factors
  sequence.specific <- factor.classes$factor[factor.classes$sequence_specific==1]
  nonsequence.specific <- factor.classes$factor[factor.classes$sequence_specific==0]
  
  factor.sequence.specificity <- list(sequence.specific = sequence.specific, 
                                      nonsequence.specific = nonsequence.specific, 
                                      dbp.factors = dbp.factors, 
                                      rbp.factors = rbp.factors)
  
  saveRDS(factor.sequence.specificity, target.file)
  return(factor.sequence.specificity)
}
###
# This function build sub feature spaces based on prior knowledge of sequence 
# specificity and biological functions of factors (DBPs, RBPs)
###
build.sequence.specific.biologically.functional.feature.sub.spaces <- function(stratify = NULL){
  
  # Retrieve sequence specificity of factors
  sequence.specific <- SEQ.SPEC$sequence.specific
  nonsequence.specific <- SEQ.SPEC$nonsequence.specific
  dbp.factors <- SEQ.SPEC$dbp.factors
  rbp.factors <- SEQ.SPEC$rbp.factors
  factors <- unique(c(unlist(dbp.factors), unlist(rbp.factors)))
  
  ## Identify go terms the set of all DBP and RBPfactors are associated with 
  # Retrieve entrez ids for querying GO DB
  entrez.factor.ids <- mapIds(org.Hs.eg.db, factors, 'ENTREZID', 'SYMBOL')
  
  # Retrieve associated go terms
  setOntology(ont="BP")
  go.terms <- getGOInfo(entrez.factor.ids)
  colnames(go.terms) <- names(entrez.factor.ids)[match(colnames(go.terms), entrez.factor.ids)]
  
  ## Build sets of DBP/RBP factors that share a common functional biological process
  # Define functional groups based on broad go terms around the transcriptional process
  target.terms <- 
    list(Chromatin=c("chromosome organization",
                     "chromatin organization",
                     "chromatin remodeling"), 
         
         Initiation = c("RNA polymerase II preinitiation complex assembly",
                        "transcription initiation from RNA polymerase II promoter"),
         
         #pos_regulation=c("positive regulation of transcription by RNA polymerase II",
         #                 "positive regulation of transcription, DNA-templated"),
         
         #neg_regulation=c("negative regulation of transcription by RNA polymerase II",
         #                 "negative regulation of transcription, DNA-templated"),
         
         Elongation=c("transcription elongation from RNA polymerase II promoter"),
         
         Termination=c("termination of RNA polymerase II transcription"),
         
         #Translation=c("translational initiation" ,"translation"), 
         
         Splicing=c("mRNA splicing, via spliceosome",
                    "regulation of alternative mRNA splicing, via spliceosome"),
         
         Processing=c("mRNA export from nucleus", "mRNA 3'-end processing"))
  
  ## For each broad GO term identify which DBP/RBP are associated with it
  factor.groups <- 
    lapply(target.terms, function(terms){
      term.associated.factors <- 
        lapply(terms, function(term){
          term.group <- 
            lapply(colnames(go.terms), function(factor){
              desc <- go.terms[2, factor]
              if(grepl(term, desc, ignore.case = T)){return(factor)}
            })
          term.group <- Filter(Negate(is.null), term.group)
        })
      term.associated.factors <- unlist(term.associated.factors)
    })
  names(factor.groups) <- names(target.terms)
  
  if(!is.null(stratify)){
    # Combine known factors with elongation factors
    factor.groups$Elongation <- unique(c(factor.groups$Elongation, unlist(stratify$Known)))
    stratify <- lapply(lapply(stratify, unlist, use.names = F), unique)
    factor.groups <- append(stratify, factor.groups)
    # Build factor set with union of elongation and 7SK binding
    factor.groups$`Elongation+7SK` <- unique(c(factor.groups$Elongation, factor.groups$`7SK.Binding`))
 }
  
  ## Sub-divide factor lists into sequence specific and non-specific factor sets
  ss.factor.groups <- 
    lapply(names(factor.groups), function(group.name){
      group.factors <- factor.groups[[group.name]]
      sub.groups <- paste(group.name, c("ss", "nss"), sep = ":")
      ll <- list(group.factors[group.factors %in% sequence.specific], 
                 group.factors[group.factors %in% nonsequence.specific])
      names(ll) <- sub.groups
      ll
    })
  ss.factor.groups <- unlist(ss.factor.groups, recursive = F)
  
  # Extend factor groups by sequence specificity defined groups
  factor.groups <- append(ss.factor.groups, factor.groups)
  
  ## Group functional factor groups by individual cell lines
  # Aggregate DPBs/RBPs into cell line specific factor sets
  cell.line.factor <- setNames(vector("list", 2), CARGS$cell.line)
  cell.line.factor$K562 <- unique(c(dbp.factors$K562, rbp.factors$K562))
  cell.line.factor$HepG2 <- unique(c(dbp.factors$HepG2, rbp.factors$HepG2))
  # Stratify functional factor groups by cell line
  lapply(names(factor.groups), function(functional.group.name){
    functional.group <- factor.groups[[functional.group.name]]
    cline.functional.groups <- 
      lapply(CARGS$cell.line, function(cline){
        functional.group[functional.group %in% cell.line.factor[[cline]]]
      })
    names(cline.functional.groups) <- CARGS$cell.line
    factor.groups[[functional.group.name]] <<- cline.functional.groups
    return()
  })
  
  return(factor.groups)
}
###
# This function builds feature matrices to train predictive models
###
build.model.matrices <- function(feature.vectors, exclude = NULL){
  
  ## Create output directories
  dir.create(paste(OUTPUT, MODEL, "model_data/", sep = ""))
  dir.create(paste(OUTPUT, MODEL, "model_evaluation/", sep = ""))
  
  ## Check if a precalculated version exists
  target.file <- paste(OUTPUT, MODEL, "model_data/", "model.matrices.RDS", sep = "")
  if(!(CARGS$new) & file.exists(target.file)){
    model.matrices <- readRDS(target.file)
    return(model.matrices)
  }
  

  ## Get common transcripts across cell lines
  # cline.transcript.overlap <- intersect(intersect(traveling.ratios$K562$transcript_id,
  #                                                 traveling.ratios$HepG2$transcript_id), 
  #                                       intersect(GENE.QUANT$K562$transcript_id, 
  #                                                 GENE.QUANT$HepG2$transcript_id))

  
  valid.traveling.ratios <- 
    lapply(traveling.ratios, function(cline.traveling.ratios){
      names(cline.traveling.ratios) <- cline.traveling.ratios$transcript_id
      cline.traveling.ratios <- cline.traveling.ratios[!cline.traveling.ratios$zero_count, ]
      #cline.traveling.ratios[match(cline.transcript.overlap, cline.traveling.ratios$transcript_id)]
    })
  
  cline.transcript.overlap <- intersect(valid.traveling.ratios$K562$transcript_id,
                                        valid.traveling.ratios$HepG2$transcript_id)

  
  # Build synchronised traveling ratio and expression target
  synchronised.targets <- 
    lapply(CARGS$cell.line, function(cline){
      
      ## Process raw traveling ratios as a regression target
      # Prepare vector with traveling ratio target targets
      targets <- data.frame(row.names = valid.traveling.ratios[[cline]][cline.transcript.overlap]$transcript_id, 
                            target_traveling.ratio = valid.traveling.ratios[[cline]][match(cline.transcript.overlap, valid.traveling.ratios[[cline]]$transcript_id)]$traveling.ratio)
      
      ## Make target with raw gene expressions
      #targets$target_transcript.expression <- GENE.QUANT[[cline]]$FPKM[match(rownames(targets), GENE.QUANT[[cline]]$transcript_id)]
      
      return(targets)
    })
  names(synchronised.targets) <- CARGS$cell.line
  
  # Build individual traveling ratio and expression target
  individual.targets <-
    lapply(CARGS$cell.line, function(cline){
      
      # # Retrieve samples for which there is a GRO-seq and RNA-seq signal
      # common.samples <- intersect(traveling.ratios[[cline]]$transcript_id,
      #                             GENE.QUANT[[cline]]$transcript_id)
      # 
      # 
      # ## Process raw traveling ratios as a regression target
      # # Prepare vector with traveling ratio target targets
      # targets <- data.frame(row.names = common.samples,
      #                       target_traveling.ratio = traveling.ratios[[cline]]$traveling.ratio[match(common.samples, traveling.ratios[[cline]]$transcript_id)])

      ## Make target with raw gene expressions
      #targets$target_transcript.expression <- GENE.QUANT[[cline]]$FPKM[match(rownames(targets), GENE.QUANT[[cline]]$transcript_id)]
      
      targets <- data.frame(row.names = valid.traveling.ratios[[cline]]$transcript_id,
                            target_traveling.ratio = valid.traveling.ratios[[cline]]$traveling.ratio)
      return(targets)
    })
  names(individual.targets) <- CARGS$cell.line
  
  ## Exclude ID features which identify proximal ncRNAs, which are only used 
  ## for back-tracing during model evaluation
  reduced.feature.vectors <- exclude.features(feature.vectors, patterns=c(
    "^chip.*Proximal\\.ncRNA\\.[0-9]+\\.ID$", 
    "^clip.*Proximal\\.ncRNA\\.[0-9]+\\.ID$",
    "^Proximal\\.ncRNA\\.[0-9]+\\.ID"))
  
  cline.feature.overlap <- intersect(colnames(reduced.feature.vectors$K562), 
                                     colnames(reduced.feature.vectors$HepG2))
  
  ## Exclude specific features if given
  if(!is.null(exclude)){ reduced.feature.vectors <- exclude.features(reduced.feature.vectors, patterns = exclude)}
  
  # ## Build synchronised target model matrices
  # synchronised.model.matrices <- 
  #   lapply(CARGS$cell.line, function(cline){
  #     # Retrieve synchronised targets
  #     targets <- synchronised.targets[[cline]]
  #     # Retrive feature.vectors
  #     feat.mat <- reduced.feature.vectors[[cline]]
  #     # Subset by synchonised target's samples
  #     feat.mat <- feat.mat[match(rownames(targets), rownames(feat.mat)), ]
  #     # Map targets to samples
  #     cbind(feat.mat, targets)
  #   })
  # names(synchronised.model.matrices) <- CARGS$cell.line
  
  ## Build synchronised target model matrices
  synchronised.model.matrices <- 
    lapply(CARGS$cell.line, function(cline){
      # Retrieve synchronised targets
      targets <- individual.targets[[cline]]
      # Retrive feature.vectors
      feat.mat <- reduced.feature.vectors[[cline]]
      # Subset by synchonised target's samples
      feat.mat <- feat.mat[match(rownames(targets), rownames(feat.mat)), ]
      feat.mat <- feat.mat[, cline.feature.overlap]
      # Map targets to samples
      cbind(feat.mat, targets)
    })
  names(synchronised.model.matrices) <- CARGS$cell.line
  
  # Build individual model matrices
  individual.model.matrices <-
    lapply(CARGS$cell.line, function(cline){
      # Retrieve synchronised targets
      targets <- individual.targets[[cline]]
      # Retrive feature.vectors
      feat.mat <- reduced.feature.vectors[[cline]]
      # Subset by synchonised target's samples
      feat.mat <- feat.mat[match(rownames(targets), rownames(feat.mat)), ]
      # Map targets to samples
      cbind(feat.mat, targets)
    })
  names(individual.model.matrices) <- CARGS$cell.line
  
  # Aggregate models matrices
  model.matrices <- list(synchronised.model.matrices = synchronised.model.matrices,
                         individual.model.matrices = individual.model.matrices)
  
 # remove features with indiscrimiatory power
  model.matrices <-
    lapply(model.matrices, function(model.matrix.type){
      lapply(model.matrix.type, function(cline.model.matrix){
        non.binding.features <- colnames(cline.model.matrix)[!grepl("^chip|clip", colnames(cline.model.matrix))]
        binding.features <- colnames(cline.model.matrix)[grepl("^chip|clip", colnames(cline.model.matrix))]
        feature.dist <-
          lapply(cline.model.matrix[, binding.features], table)
        feature.dist <- lengths(feature.dist)
        feature.dist <- feature.dist[feature.dist>1]
        cline.model.matrix <- cline.model.matrix[, c(non.binding.features, names(feature.dist) )]
        return(cline.model.matrix)
      })
    })
  
  ## Build matrices with sub feature spaces
  # Retrieve sub features spaces based on factors to reduce feature matrices for model comparison
  sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
  sub.factor.feature.spaces <- 
    build.sequence.specific.biologically.functional.feature.sub.spaces(stratify = sub.factor.feature.spaces)
  #sub.factor.feature.spaces <- append(sub.factor.feature.spaces, sequence.specific.biologically.functional.feature.sub.spaces)
  
  # Subset matrices by sub feature spaces
  model.matrices <- 
    lapply(CARGS$cell.line, function(cline){
      matrix.types <- 
        lapply(model.matrices, function(matrix.subtype){
          matrix.subtypes <- 
            lapply(sub.factor.feature.spaces, function(feature.space.type){
              non.binding.signals <- colnames(matrix.subtype[[cline]])[!grepl("^chip|clip", colnames(matrix.subtype[[cline]]))]
              # Retrieve DNA or RNA binding factors from model matrix feature space
              binding.signals <- colnames(matrix.subtype[[cline]])[grepl("^chip|clip", colnames(matrix.subtype[[cline]]))]
              if(!length(feature.space.type[[cline]]) == 0){
                binding.signals <- binding.signals[which(grepl(paste0(feature.space.type[[cline]], collapse = "\\.|"), binding.signals))]
                # factors <- 
                #   lapply(str_split(binding.signals, pattern = "\\."), function(x){x[2]}) %>% 
                #   unlist()
                # factors <- toupper(factors)
                # # Extract non-binding signals
                # # Reduce features space by sub feature spaces
                # valid.factor.features <- factors %in% feature.space.type[[cline]]
                # 
                # valid.factor.features <- binding.signals[valid.factor.features]
                # matrix.subtype[[cline]] <- matrix.subtype[[cline]][ , c(non.binding.signals, valid.factor.features)]
                matrix.subtype[[cline]] <- matrix.subtype[[cline]][ , c(non.binding.signals, binding.signals)]
                return(matrix.subtype[[cline]])
              }
              matrix.subtype[[cline]] <- matrix.subtype[[cline]][ , non.binding.signals]
              return(matrix.subtype[[cline]])
            })
          names(matrix.subtypes) <- names(sub.factor.feature.spaces)
          matrix.subtypes
        })
      names(matrix.types) <- names(model.matrices)
      matrix.types
    })
  names(model.matrices) <- CARGS$cell.line
  
  # ## Synchronize feature matrices
  # clineA = CARGS$cell.line[1]
  # clineB = CARGS$cell.line[2]
  # mtype <- "synchronised.model.matrices"
  # feature.subspaces <- names(sub.factor.feature.spaces)
  # lapply(feature.subspaces, function(subspace){
  #   # Retrieve common samples
  #   samples.clineA <- rownames(model.matrices[[clineA]][[mtype]][[subspace]])
  #   samples.clineB <- rownames(model.matrices[[clineB]][[mtype]][[subspace]])
  #   common.samples <- intersect(samples.clineA, samples.clineB)
  #   
  #   # Retrieve common features
  #   features.clineA <- colnames(model.matrices[[clineA]][[mtype]][[subspace]])
  #   features.clineB <- colnames(model.matrices[[clineB]][[mtype]][[subspace]])
  #   common.features <- intersect(features.clineA, features.clineB)
  #   
  #   # Build col and rowwise synchronized versions of feature matrices
  #   model.matrices[[clineA]][[mtype]][[subspace]] <<- model.matrices[[clineA]][[mtype]][[subspace]][common.samples, common.features]
  #   model.matrices[[clineB]][[mtype]][[subspace]] <<- model.matrices[[clineB]][[mtype]][[subspace]][common.samples, common.features]
  #   return()
  # })
  
  ## Synchronize feature matrices
  clineA = CARGS$cell.line[1]
  clineB = CARGS$cell.line[2]
  mtype <- "synchronised.model.matrices"
  feature.subspaces <- names(sub.factor.feature.spaces)
  lapply(feature.subspaces, function(subspace){

    # Retrieve common features
    features.clineA <- colnames(model.matrices[[clineA]][[mtype]][[subspace]])
    features.clineB <- colnames(model.matrices[[clineB]][[mtype]][[subspace]])
    common.features <- intersect(features.clineA, features.clineB)
    
    # Build col and rowwise synchronized versions of feature matrices
    model.matrices[[clineA]][[mtype]][[subspace]] <<- model.matrices[[clineA]][[mtype]][[subspace]][, common.features]
    model.matrices[[clineB]][[mtype]][[subspace]] <<- model.matrices[[clineB]][[mtype]][[subspace]][, common.features]
    return()
  })
  
  # Create models with only static features (non-binding features)
  # Subset matrices by sub feature spaces
  lapply(CARGS$cell.line, function(cline){
    matrix.types <- 
      lapply(names(model.matrices[[cline]]), function(matrix.subtype){
        model.matrices[[cline]][[matrix.subtype]]$Annotation <<- model.matrices[[cline]][[matrix.subtype]]$All[, !grepl("^chip|clip", colnames(model.matrices[[cline]][[matrix.subtype]]$All))]
        return()
      })
  })
  
  ## Save individual model matrices
  synchronised.model.matrices.nfo <- paste(OUTPUT, MODEL, "model_data/", "synchronised.model.matrices", ".txt", sep = "")
  individual.model.matrices.nfo <- paste(OUTPUT, MODEL, "model_data/", "individual.model.matrices", ".txt", sep = "")
  lapply(CARGS$cell.line, function(cline){
    lapply(names(model.matrices[[cline]]$synchronised.model.matrices), function(mnames){
      target.matrix <- paste(OUTPUT, MODEL, "model_data/", "synchronised.model.matrices.", mnames, ".", cline, ".RDS", sep = "")
      saveRDS(model.matrices[[cline]]$synchronised.model.matrices[[mnames]], target.matrix)
      target.matrix <- paste(OUTPUT, MODEL, "model_data/", "synchronised.model.matrices.", mnames, ".", cline, ".feather", sep = "")
      write_feather(model.matrices[[cline]]$synchronised.model.matrices[[mnames]], target.matrix)
      # Write out file with names of matrices as input for runner 
      write.table(file_path_as_absolute(target.matrix), synchronised.model.matrices.nfo, append = T, 
                  row.names = F, col.names = F)
      return()
    })
    lapply(names(model.matrices[[cline]]$individual.model.matrices), function(mnames){
      target.matrix <- paste(OUTPUT, MODEL, "model_data/", "individual.model.matrices.", mnames,".", cline, ".RDS", sep = "")
      saveRDS(model.matrices[[cline]]$individual.model.matrices[[mnames]], target.matrix)
      target.matrix <- paste(OUTPUT, MODEL, "model_data/", "individual.model.matrices.", mnames,".", cline, ".feather", sep = "")
      write_feather(model.matrices[[cline]]$individual.model.matrices[[mnames]], target.matrix)
      # Write out file with names of matrices as input for runner 
      write.table(file_path_as_absolute(target.matrix), individual.model.matrices.nfo, append = T, 
                  row.names = F, col.names = F)
      return()
    })
  })
  
  ## Save matrices list
  saveRDS(model.matrices, target.file)
  return(model.matrices)
}