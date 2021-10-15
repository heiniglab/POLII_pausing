###
# This function identifies 7SK related transcript binding RNA-binding proteins
###
identify.7SK.binders <- function(RN7SK, clip.calls){
  
  # Check if a precalculated version exists
  target.file <- paste(OUTPUT, "RN7SK.binders.RDS", sep = "")
  if(!(CARGS$new) & file.exists(target.file)){
    binders <- readRDS( target.file)
    return(binders)
  }
  
  binders <-
    lapply(names(clip.calls), function(cline){
      # Overlap clip calls with 7SK annotations
      ovs <- findOverlaps(clip.calls[[cline]], RN7SK.ANNOT[[cline]])
      # Get binders per 7SK related transcript
      binders <-
        lapply(unique(subjectHits(ovs)), function(rn7sk.index){
          # Get indices of overlap between rbp and 7SK related transcript
          bound <- which(subjectHits(ovs) == rn7sk.index)
          # Retrieve RBPs which bind a 7SK transcript
          binders <- clip.calls[[cline]][queryHits(ovs[bound])]
          # Annotate which 7SK transcript they bind
          binders$RN7SK.tx.id <- RN7SK.ANNOT[[cline]][rn7sk.index]$transcript_id
          # Annotate the type of 7SK transcripts bound (pseudo or non pseudo)
          binders$is_RN7SK <- RN7SK.ANNOT[[cline]]$is_RN7SK[match(binders$RN7SK.tx.id, RN7SK.ANNOT[[cline]]$transcript_id)]
          return(binders)
        })
      binders <- do.call(c, binders)
      return(binders)
    })
  names(binders) <- names(clip.calls)
  
  # Save intermediate input data
  saveRDS(binders, target.file)
  return(binders)
}
###
# REVISE
###
retrieve.shared.domains2 <- function(){
  ###
  # Overlap of 7SK binding RBPs between cell lines
  ###
  CLIPseq <- c(filtered.eClip.calls, HEK293=filtered.parClip.calls)
  
  
  
  # Check for 7SK binding RBPs from all three cell lines
  rn7sk.targets <- 
    lapply(novel.rn7sk.binders, function(calls){
      sort(unique(calls$ensembl_id))
    }) 
  clip.targets <- unique(unlist(rn7sk.targets))
  # Load gene group definitions from HGNC
  gene.groups <- read.table(paste(INPUT, GENEANNOT, "HGNC/", "hgnc_gene_groups.txt", sep = ""), 
                            header = T, sep = "\t", stringsAsFactors = F,quote = "\"",
                            fill = T)
  # Subset by 7SK binding RBPs
  gene.groups <- gene.groups[match(clip.targets, gene.groups$Ensembl.gene.ID), ]
  # Mark entries with no functional group assigned
  gene.groups <- gene.groups[which(!is.na(gene.groups$Gene.group.ID)),]
  gene.groups$Gene.group.ID[which(gene.groups$Gene.group.ID == "")] <- "-1"
  # Subset initial targeted ncRNAs by those which have functional group annotations
  rn7sk.targets <- 
    lapply(rn7sk.targets, function(binders){
      binders[binders %in% gene.groups$Ensembl.gene.ID ]
    })
  # Create list of functional group ids per gene
  unique.groups <- gene.groups$Gene.group.ID
  group.memberships <- sapply(unique.groups, function(group){strsplit(group,"\\|")})
  names(group.memberships) <-  gene.groups$Ensembl.gene.ID
  # Unique present groups
  unique.groups <- unique(unlist(group.memberships))
  # List of genes per group
  groups <- setNames(vector("list", length(unique.groups)), unique.groups)
  lapply(names(group.memberships), function(gene){
    sub.group <- group.memberships[[gene]]
    sapply(sub.group, function(sgroup){
      groups[[sgroup]] <<-  c(groups[[sgroup]], gene)
      return()
    })
  })
  
  
  
  ## Distribution of group sizes (singletons excluded)
  group.card <- sort(lengths(groups), decreasing = T)
  group.card <- group.card[group.card>2]
  pdf(paste(PLOTS, "rn7sk_binding_functional_group_size_dist.pdf", sep = ""))
  png(paste(PLOTS, "rn7sk_binding_functional_group_size_dist.png", sep = ""))
  barplot(group.card)
  dev.off()
  # With hgnc symbols for better interpretation
  groups.red <- groups[which(lengths(groups)>2)]
  groups.pretty <- 
    lapply(groups, function(group){
      as.character(GENE.ANNOT$id_map$hgnc_symbol[match(groups.red, GENE.ANNOT$id_map$ensembl_gene_id)])
    })
  # The group name of those groups, i.e. biological interpretation
  sapply(names(groups.red), function(group.id){
    gene.groups$Gene.group.name[match(group.id, gene.groups$Gene.group.ID)]
  })
  
  #x singletons <- groups[lengths(groups)==1]
  #x groups <- groups[lengths(groups)>2]
  groups.df.shared.domains <- 
    data.frame(ID=paste(unlist(groups), rep.int(names(groups), lengths(groups)), sep = ":"),
               genes=unlist(groups), 
               class = rep.int(names(groups), lengths(groups)),
               K562 = 0,
               HepG2 = 0,
               HEK293 = 0)
  lapply(names(rn7sk.targets), function(cline){
    groups.df.shared.domains[groups.df.shared.domains$genes %in% rn7sk.targets[[cline]], cline] <<- 1
    return()
  })
  
  
  groups.df.shared.domains <- split(groups.df.shared.domains, groups.df.shared.domains$class)
  shared.sdom.target.transcripts <- setNames(vector("list", 7), c("K562", "HepG2", "HEK293", "K562_HepG2", "K562_HEK293", "HepG2_HEK293", "K562_HepG2_HEK293"))
  for(domain.name in names(groups.df.shared.domains)){
    print(domain.name)
    present.clines <- c(any(groups.df.shared.domains[[domain.name]]$K562==1),
                        any(groups.df.shared.domains[[domain.name]]$HepG2==1),
                        any(groups.df.shared.domains[[domain.name]]$HEK293==1))
    names(present.clines) <- c("K562", "HepG2", "HEK293")
    sdom.factors <- as.character(groups.df.shared.domains[[domain.name]]$genes)
    
    for(cline in names(present.clines)[present.clines]){
      shared.sdom.target.transcripts[[cline]] <- c(shared.sdom.target.transcripts[[cline]], 
                                                   sdom.factors)
    }
    
    if(sum(present.clines)==2){
      paired.cline <- paste(names(present.clines[present.clines]), collapse = "_")
      shared.sdom.target.transcripts[[cline]] <- c(shared.sdom.target.transcripts[[paired.cline]], 
                                                   sdom.factors)
    }
    
    if(sum(present.clines)==3){
      
      shared.sdom.target.transcripts[["K562_HepG2"]] <- c(shared.sdom.target.transcripts[["K562_HepG2"]], 
                                                          sdom.factors)
      shared.sdom.target.transcripts[["K562_HEK293"]] <- c(shared.sdom.target.transcripts[["K562_HEK293"]], 
                                                           sdom.factors)
      shared.sdom.target.transcripts[["HepG2_HEK293"]] <- c(shared.sdom.target.transcripts[["HepG2_HEK293"]], 
                                                            sdom.factors)
      shared.sdom.target.transcripts[["K562_HepG2_HEK293"]] <- c(shared.sdom.target.transcripts[["K562_HepG2_HEK293"]], 
                                                                 sdom.factors)
    }
    
  }
  shared.sdom.target.transcripts <- lapply(shared.sdom.target.transcripts, unique)
  
  pdf(paste(PLOTS,"overlap_rn7sk_binding_rbps_all_clines_domain.genes.pdf", sep = ""))
  png(paste(PLOTS,"overlap_rn7sk_binding_rbps_all_clines_sdom.pseudo.png", sep = ""))
  venn.plot <- draw.triple.venn(
    area1 = length(shared.sdom.target.transcripts$K562),
    area2 = length(shared.sdom.target.transcripts$HepG2),
    area3 = length(shared.sdom.target.transcripts$HEK293),
    n12 = length(shared.sdom.target.transcripts$K562_HepG2),
    n13 = length(shared.sdom.target.transcripts$K562_HEK293),
    n23 = length(shared.sdom.target.transcripts$HepG2_HEK293),
    n123 = length(shared.sdom.target.transcripts$K562_HepG2_HEK293),
    category = c("K562" , "HepG2", "HEK293"),
    fill = viridis(3, end = 0.75),
    #fill = rainbow(4),
    #lty = "dashed",
    cex = 2,
    cat.cex = 2,
    cat.col = viridis(3, end = 0.75)
  )
  dev.off()
  
  # Domains shared bertween all pseudo 7SK binders between all cell lines
  unique(gene.groups$Gene.group.name[match(shared.sdom.target.transcripts$K562_HepG2_HEK293,gene.groups$Ensembl.gene.ID)])
  functional.domains <- sort(table(gene.groups$Gene.group.name[match(shared.sdom.target.transcripts$K562_HepG2_HEK293,gene.groups$Ensembl.gene.ID)]), decreasing = T)
  head(functional.domains)
  png(paste(PLOTS,"overlap_rn7sk_binding_rbps_all_clines_protein.domain.distribution.png", sep = ""))
  names(functional.domains) <- ""
  barplot(functional.domains, main = "Functional proteins domain count of 7SK binding proteins")
  dev.off()
  
  group.sizes <- lengths(groups)
  group.sizes <- group.sizes[group.sizes>1]
  group.sizes <- sort(group.sizes, decreasing = T)
  png(paste(PLOTS,"overlap_rn7sk_binding_rbps_all_clines_prots.per.domain.png", sep = ""))
  barplot(group.sizes, main = paste("Number of proteins per functional domain - Avg = ",format(mean(group.sizes),digits = 3)))
  dev.off()
  
  
  
  domains <- 
    lapply(shared.sdom.target.transcripts, function(group){
      unlist(sapply(group, function(member){
        group.memberships[[member]]
      } ))
    })
  domains <- 
    lapply(domains, unique)
  pdf(paste(PLOTS,"overlap_rn7sk_binding_rbps_all_clines_domain.genes.pdf", sep = ""))
  png(paste(PLOTS,"overlap_rn7sk_binding_rbps_prot_domains.png", sep = ""))
  venn.plot <- draw.triple.venn(
    area1 = length(domains$K562),
    area2 = length(domains$HepG2),
    area3 = length(domains$HEK293),
    n12 = length(domains$K562_HepG2),
    n13 = length(domains$K562_HEK293),
    n23 = length(domains$HepG2_HEK293),
    n123 = length(domains$K562_HepG2_HEK293),
    category = c("K562" , "HepG2", "HEK293"),
    fill = viridis(3, end = 0.75),
    #fill = rainbow(4),
    #lty = "dashed",
    cex = 2,
    cat.cex = 2,
    cat.col = viridis(3, end = 0.75)
  )
  dev.off()
  
  
  # K562.gene.classes <- as.character(unique(groups$class[which(groups$K562==1)]))
  # HepG2.gene.classes <- as.character(unique(groups$class[which(groups$HepG2==1)]))
  # HEK293.gene.classes <- as.character(unique(groups$class[which(groups$HEK293==1)]))
  # 
  # #pdf("")
  # # Draw venn diagram
  # pdf(paste(PLOTS,"overlap_rn7sk_binding_rbps_all_clines_functional_groups.pdf", sep = ""))
  # png(paste(PLOTS,"overlap_rn7sk_binding_rbps_all_clines_functional_groups.png", sep = ""))
  # venn.plot <- draw.triple.venn(
  #   area1 = length(K562.gene.classes),
  #   area2 = length(HepG2.gene.classes),
  #   area3 = length(HEK293.gene.classes),
  #   n12 = length(intersect(K562.gene.classes, HepG2.gene.classes )),
  #   n13 = length(intersect(K562.gene.classes, HEK293.gene.classes )),
  #   n23 = length(intersect(HepG2.gene.classes, HEK293.gene.classes )),
  #   n123 = length(intersect(intersect(K562.gene.classes, HepG2.gene.classes), HEK293.gene.classes )),
  #   category = c("K562" , "HepG2", "HEK293"),
  #   fill = viridis(3, end = 0.75),
  #   #fill = rainbow(4),
  #   #lty = "dashed",
  #   cex = 2,
  #   cat.cex = 2,
  #   cat.col = viridis(3, end = 0.75)
  # )
  # dev.off()
  
  
  
  # # Same for corresponding genes
  # true.binders.binders <- 
  #   lapply(novel.rn7sk.binders, function(bindings){
  #     unique(bindings$ensembl_id)
  #   })
  # true.rn7sk.bindings.K562 <- intersect(true.binders.binders$K562, true.binders.binders$HepG2)
  # true.rn7sk.bindings.HepG2 <- intersect(true.binders.binders$K562, true.binders.binders$HEK293)
  # true.rn7sk.bindings.HEK293 <- intersect(true.binders.binders$HepG2, true.binders.binders$HEK293)
  # 
  # 
  # pdf(paste(PLOTS,"overlap_rn7sk_binding_rbps_all_clines_functional_groups_genes.pdf", sep = ""))
  # png(paste(PLOTS,"overlap_rn7sk_binding_rbps_all_clines_functional_groups_genes.png", sep = ""))
  # venn.plot <- draw.triple.venn(
  #   area1 = length(groups.df$genes[groups.df$K562==1]),
  #   area2 = length(groups.df$genes[groups.df$HepG2==1]),
  #   area3 = length(groups.df$genes[groups.df$HEK293==1]),
  #   n12 = length(groups.df$genes[groups.df$K562==1 & groups.df$HepG2==1]),
  #   n13 = length(groups.df$genes[groups.df$K562==1 & groups.df$HEK293==1]),
  #   n23 = length(groups.df$genes[groups.df$HepG2==1 & groups.df$HEK293==1]),
  #   n123 = length(groups.df$genes[groups.df$K562==1 & groups.df$HepG2==1 & groups.df$HEK293==1]),
  #   category = c("K562" , "HepG2", "HEK293"),
  #   fill = viridis(3, end = 0.75),
  #   #fill = rainbow(4),
  #   #lty = "dashed",
  #   cex = 2,
  #   cat.cex = 2,
  #   cat.col = viridis(3, end = 0.75)
  # )
  # dev.off()
  
  
  
  
  
}
###
# DONE
###
rn7sk.phylogeny <- function(){
  
  rn7sk.metrics <- 
  lapply(RN7SK.ANNOT, function(rn7sk.variants){
    # Function to print sequence
    printSplitString <- function(x, width=getOption("width") - 1){
      starts <- seq(from=1, to=nchar(x), by=width) 
      for (i in 1:length(starts))cat(substr(x, starts[i], starts[i] + width - 1), "\n")
    }
    # Retrieve 7SK transcript's DNA sequences
    rn7sk.variant.sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg19, rn7sk.variants)
    # Perform a multiple sequence alignment of 7SK transcriptss
    multiple.align <- msa(rn7sk.variant.sequences, method = "ClustalW")
    # Calculate conservation scores
    data(PAM250)
    cons <- msaConservationScore(multiple.align, PAM250, gapVsGap=0, type="upperlower", thresh=c(40, 20))
    hist(cons, 50)
    mean.cons <- mean(cons)
    # Print consensus sequence
    printSplitString(multiple.align)
    
    # Print MSA to pdf
    msaPrettyPrint(multiple.align, output="asis", y=c(164, 213),subset=c(1:6), 
                   showNames="none", showLogo="top",logoColors="rasmol", 
                   shadingMode="similar",showLegend=FALSE, askForOverwrite=FALSE)
    #Pairwise sequence distance
    multiple.align.seqinr <- msaConvert(multiple.align, type="seqinr::alignment")
    seq.dist <- dist.alignment(multiple.align.seqinr, "identity")
    mean.similarity <- 1-mean(seq.dist)**2
    # Phylogenetic Tree
    phylo.tree <- nj(seq.dist)
    plot(phylo.tree, main="Phylogenetic Tree of 7SK sequences")
    list(sim = mean.similarity, 
           mean.cons = mean.cons)
  })
  names(rn7sk.metrics) <- names(RN7SK.ANNOT)
}
## 7SK binding ncRNA binders baserd on ENCODEand PARCLIP
# Retrieve 7SK associated factors across cell lines
# REVISE
shared.ncRNAs.7SK.binders <- function(){
  
  # Get additionally bound ncRNAs
  peak.signals <- c(eCLIPseq, HEK293 = PARCLIP)
  sub.peaks <- lapply(names(peak.signals), function(cline){peak.signals[[cline]][peak.signals[[cline]]$hgnc_symbol %in% unique(unlist(novel.rn7sk.binders[[cline]]$hgnc_symbol))]})
  names(sub.peaks) <- names(peak.signals)
  # Get unique 7SK binders
  new.rn7sk.binders <- 
    lapply(novel.rn7sk.binders, function(binders){
      sort(unique(binders$hgnc_symbol))
    })
  indicator.cols <- 
    lapply(names(new.rn7sk.binders), function(cline){
      paste(new.rn7sk.binders[[cline]], "_", cline, sep = "")
    })
  names(indicator.cols) <-names(new.rn7sk.binders)
  
  
  # X
  # valid.ncRNAs <- unique(unlist(lapply(valid.non.coding.transcripts, function(expressed.ncRNAs){
  #   expressed.ncRNAs$transcript_id
  # })))
  
  additonal.ncrnas <- 
    lapply(names(sub.peaks), function(cline){
      ovs <- findOverlaps(sub.peaks[[cline]], valid.non.coding.transcripts[[cline]], select = "all") 
      unique(valid.non.coding.transcripts[[cline]]$transcript_id[unique(subjectHits(ovs))])
    })
  names(additonal.ncrnas) <- names(sub.peaks)
  all.target.ncrnas <-  unique(unlist(additonal.ncrnas))
  # Indicator matrix of additonal ncRNAs that are bound by 7SK binding RBPs
  ncrna.indicator <- 
    matrix(0, ncol = length(unique(unlist(indicator.cols))), 
           nrow = length(all.target.ncrnas),
           dimnames = list(all.target.ncrnas, unique(unlist(indicator.cols))))
  # Fill indicator matrix
  lapply(names(additonal.ncrnas), function(cline){
    ncrna.indicator[additonal.ncrnas[[cline]], indicator.cols[[cline]]] <<- 1
    return()
  })
  # Identify transcripts that are expressed in all cell lines and bound in all cell lines by any cell line specific 7SK binders
  binding.order <- rowSums(ncrna.indicator)
  max.bindings <- table(binding.order)
  mostly.targeted <- max(as.numeric(names(max.bindings)))
  target.transcripts <- names(binding.order[which(binding.order == mostly.targeted)])
  target.transcripts <- GENE.ANNOT$transcripts[match(target.transcripts, GENE.ANNOT$transcripts$transcript_id)]
  table(target.transcripts$transcript_type)
  target.names <- target.transcripts$gene_name
  target.names <- sort(target.names, decreasing = F)
  
  shared.transcripts <- 
    lapply(names(valid.non.coding.transcripts), function(cline){
      dat <- 
        data.frame(id = valid.non.coding.transcripts[[cline]]$transcript_id, 
                   length = width(valid.non.coding.transcripts[[cline]])/1e6,
                   FPKM = GENE.QUANT[[cline]]$FPKM[match(valid.non.coding.transcripts[[cline]]$transcript_id, GENE.QUANT[[cline]]$transcript_id)], 
                   shared = 0,
                   cline = cline, stringsAsFactors = F)
      dat[!duplicated(dat),]
    })
  shared.transcripts <- do.call(rbind, shared.transcripts)
  shared.transcripts$shared[match(target.transcripts$transcript_id, shared.transcripts$id)] <- 1
  shared.transcripts$shared <- factor(shared.transcripts$shared)
  
  # med <- median(shared.transcripts$FPKM[shared.transcripts$shared==1])
  # shared.transcripts$above_med <- 0
  # shared.transcripts$above_med[shared.transcripts$shared == 0 &
  #                                shared.transcripts$FPKM > med] <- 1
  # shared.transcripts$above_med <- factor(shared.transcripts$above_med)
  
  # TODO indicate whether a transcript is bound or not
  shared.transcripts$is_bound <- 0
  shared.transcripts$is_bound[match(all.target.ncrnas, shared.transcripts$id)] <- 1
  shared.transcripts$is_bound <- factor(shared.transcripts$is_bound)
  #shared.transcripts$shared_above_med_bound <- interaction(shared.transcripts$shared, shared.transcripts$above_med, shared.transcripts$is_bound)
  
  #X
  # # shared vs non shared
  # ggplot(shared.transcripts, aes(x = shared_above_med_bound, y=FPKM, fill = shared_above_med_bound)) + 
  #   geom_boxplot()
  
  # Expression
  ggplot(shared.transcripts, aes(x=FPKM, color=is_bound)) +
    geom_freqpoly(bins = 20)
  temp <- shared.transcripts
  
  p1<- 
    ggplot(temp, aes(x=FPKM, color=is_bound)) +
    geom_histogram(aes(y = ..count../sum(..count..),fill = is_bound), position="identity", alpha=0.3, bins = 200)
  #scale_y_continuous(limits = c(0, 160))
  p2 <-
    ggplot(shared.transcripts, aes(x=FPKM, color=shared)) +
    geom_histogram(aes(y = ..count../sum(..count..), fill = shared), position="identity", alpha=0.3, bins = 200)
  #scale_y_continuous(limits = c(0, 160))
  grid.arrange(p1, p2, nrow = 2) 
  
  
  # Length
  temp <- shared.transcripts[shared.transcripts$shared == 0,]
  p1<- 
    ggplot(temp, aes(x=length, color=is_bound)) +
    geom_histogram(aes(y = ..count../sum(..count..),fill = is_bound), position="identity", alpha=0.3, bins = 200)
  #scale_y_continuous(limits = c(0, 160))
  p2 <-
    ggplot(shared.transcripts, aes(x=length, color=shared)) +
    geom_histogram(aes(y = ..count../sum(..count..), fill = shared), position="identity", alpha=0.3, bins = 200)
  #scale_y_continuous(limits = c(0, 160))
  grid.arrange(p1, p2, nrow = 2) 
  
  
  # X for visualization
  sort(GENE.ANNOT$transcripts$transcript_id[match(shared.transcripts$id[shared.transcripts$shared==1], GENE.ANNOT$transcripts$transcript_id)])
  
  # Retrieve ncRNAs that are shared between all cell lines
  shared.ncrnas <- GENE.ANNOT$transcripts[match(shared.transcripts$id[shared.transcripts$shared==1], GENE.ANNOT$transcripts$transcript_id),]
  #shared.ncrnas <- sort(unlist(shared.ncrnas))
  #shared.ncrnas
  
  #X
  # # Write out ncRNA sequences for mulötiple sequence alignment and secondary structure predition
  # ncrna.sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg19, shared.ncrnas)
  # names(ncrna.sequences) <- shared.ncrnas$transcript_id
  # ncrna.sequence.widths <- width(ncrna.sequences)
  # ncrna.sequences <- ncrna.sequences[c(97,103)]
  # ncrna.sequences <- ncrna.sequences[ncrna.sequence.widths < 7500]
  # sapply(ncrna.sequences, function(s){  write.table(as.character(s), 
  #                                                      paste(OUTPUT, "ncRNA_sequences2.fa", sep = ""), 
  #                                                      append = T, quote = F,  row.names = F, col.names = F)})
  # 
  # "ENST00000573866" 
  # "ENST00000574939"
  
  
  # Pathways synergistically regulated by the ncRNAs
  #PPI <- read.table(paste(INPUT, STRING, "/human_gene_hgnc_symbol.links.detailed.v10.txt", sep = ""), header = T)
  
  
  a1.clipped.prots.K562 <- additonal.ncrnas$K562
  a2.clipped.prots.HepG2 <- additonal.ncrnas$HepG2
  a3.clipped.prots.HEK293 <-additonal.ncrnas$HEK293
  
  #pdf("")
  # Draw venn diagram
  pdf(paste(PLOTS,"rn7sk_binding_rbps_all_clines.pdf", sep = ""))
  png(paste(PLOTS,"rn7sk_binding_rbps_shared.ncrnas.png", sep = ""))
  venn.plot <- draw.triple.venn(
    area1 = length(a1.clipped.prots.K562),
    area2 = length(a2.clipped.prots.HepG2),
    area3 = length(a3.clipped.prots.HEK293),
    n12 = length(intersect(a1.clipped.prots.K562, a2.clipped.prots.HepG2 )),
    n13 = length(intersect(a1.clipped.prots.K562, a3.clipped.prots.HEK293 )),
    n23 = length(intersect(a2.clipped.prots.HepG2, a3.clipped.prots.HEK293 )),
    n123 = length(intersect(intersect(a1.clipped.prots.K562, a2.clipped.prots.HepG2), a3.clipped.prots.HEK293 )),
    category = c("K562" , "HepG2", "HEK293"),
    fill = viridis(3, end = 0.75),
    #fill = rainbow(4),
    #lty = "dashed",
    cex = 2,
    cat.cex = 2,
    cat.col = viridis(3, end = 0.75)
  )
  dev.off()
  
  
  
  
  library(LncPath)
  background.lnc.mrna.net <- getNet()
  kegg.results <- lncPath(unlist(shared.ncrnas$gene_id), background.lnc.mrna.net,  Weighted = TRUE, PathwayDataSet = "KEGG", nperm = 100,
                          minPathSize = 0, maxPathSize = 500)
  reactome.results <- lncPath(unlist(shared.ncrnas$gene_id), background.lnc.mrna.net,  Weighted = TRUE, PathwayDataSet = "Reactome", nperm = 100,
                              minPathSize = 0, maxPathSize = 500)
  ## Print to table
  Table <- lncPath2Table(kegg.results)
  head(Table)
  Table[as.numeric(as.character(Table$`P Value`))<=0,]
  
  plotRunningES(Result, Name = "KEGG_PROTEASOME")
  plotRunningES(Result, Name = "KEGG_VEGF_SIGNALING_PATHWAY")
  plotRunningES(Result, Name = "KEGG_INSULIN_SIGNALING_PATHWAY")
  
  
  geneSetDetail(Result, Name = "KEGG_RIBOSOME")
  geneSetDetail(Result, Name = "KEGG_PROTEASOME")
  geneSetDetail(Result, Name = "KEGG_MELANOMA")
  
  # Profile <- getExampleData("Profile")
  # Labels <- getExampleData("Labels")
  # drawAHeatMap(Result, Name = "KEGG_RIBOSOME", PCExpr = Profile, Labels = Labels)
  # 
  
  # Global clustering by target ncRNA
  ncrna.indicator.dist <- dist(ncrna.indicator, diag = F, method = "binary")
  rownames(ncrna.indicator[binding.order[binding.order == 101],])
  ## Cluster transcripts (proteins) by semantic similarity
  # Compute pairwise rowdistances based on binary measure
  dist.mat <- dist(mat, method = "binary")
  # Hierarchically cluster matrix row-wise
  hc <- hclust(tx.go.sim.bp, method="ward.D2")
  # Dynamically cut hirarchy to yield clusters
  clusters <- cutreeHybrid(hc, distM = as.matrix(dist.mat), minClusterSize = 2)
  # Get ids of clusters
  cluster.ids <- sort(unique(clusters$labels), decreasing = T)
  
  
  # Retrieve overlapping 7SK targets between cell lines
  K562.HepG2 <- intersect(common.7sk.binders$K562, common.7sk.binders$HepG2)
  K562.HEK293 <- intersect(common.7sk.binders$K562, common.7sk.binders$HEK293)
  HepG2.HEK293 <- intersect(common.7sk.binders$HepG2, common.7sk.binders$HEK293)
  common.binders <- intersect(common.7sk.binders$K562, intersect(common.7sk.binders$HepG2, common.7sk.binders$HEK293))
  shared.factors <- list(K562.HepG2 = K562.HepG2,
                         K562.HEK293 = K562.HEK293,
                         HepG2.HEK293 = HepG2.HEK293,
                         K562.HepG2.HEK293 = common.binders)
  
  
  
  ##
  
  # X
  # # Retrieve other ncRNAs that are bound by 7SK binding factors
  # peak.signals <- c(eCLIPseq, HEK293 = PARCLIP)
  # # Factors that bind 7SK across cell lines
  # shared.7SK.factors <- K562.HepG2
  # sub.peaks <- lapply(peak.signals, function(peaks){peaks[peaks$hgnc_symbol %in% shared.7SK.factors]})
  # sub.peaks <- sub.peaks[lengths(sub.peaks) > 0]
  # #clines <- unlist(strsplit(sfactors, "\\."))
  # further.ncrnas <- 
  #   lapply(names(sub.peaks), function(cline){
  #     #ovs <- findOverlaps(peaks, GENE.ANNOT$transcripts)
  #     #unique(GENE.ANNOT$transcripts$transcript_id[unique(subjectHits(ovs))])
  #     ovs <- findOverlaps(sub.peaks[[cline]], valid.non.coding.transcripts[[cline]], select = "all") 
  #     unique(valid.non.coding.transcripts[[cline]]$transcript_id[unique(subjectHits(ovs))])
  #   })
  # names(further.ncrnas) <- further.ncrnas
  # # Get overlapping cnRNAs sets for shared factors
  # shared.ncrnas <- unique(unlist(further.ncrnas))
  # # Subset by those which are expressed in both cell lines
  # shared.ncrnas <- shared.ncrnas[shared.ncrnas %in% valid.non.coding.transcripts$K562$transcript_id &
  #                                  shared.ncrnas %in% valid.non.coding.transcripts$HepG2$transcript_id]
  # # Remove 7SK variants
  # by.transcipt.name <- grepl("7SK", GENE.ANNOT$transcripts$transcript_name)
  # by.gene.name <- grepl("7SK", GENE.ANNOT$transcripts$gene_name)
  # RN7SK.related <- GENE.ANNOT$transcripts[by.transcipt.name & by.gene.name]
  # shared.ncrnas <- shared.ncrnas[!shared.ncrnas %in% RN7SK.related$transcript_id]
  # saveRDS(shared.ncrnas, paste(OUTPUT, "shared.7SK.binders.ncRNAs.RDS", sep = ""))
  # 
  
}
###
# REIVISE
###
rn7sk.variant.expression <- function(){
  
  lapply(names(RN7SK.ANNOT), function(cline){
    
    # Retrieve average gene and ncRNA expression levels
    avg.gene.expr <- mean(GENE.QUANT[[cline]]$FPKM[match(valid.coding.transcripts[[cline]]$transcript_id, GENE.QUANT[[cline]]$transcript_id)])
    avg.ncrna.expr <- mean(GENE.QUANT[[cline]]$FPKM[match(valid.non.coding.transcripts[[cline]]$transcript_id, GENE.QUANT[[cline]]$transcript_id)])
    
    rn7sk.expressions <- GENE.QUANT[[cline]]$FPKM[match(unique(RN7SK.ANNOT[[cline]]$transcript_id), GENE.QUANT[[cline]]$transcript_id)]
    names(rn7sk.expressions) <- ""
    names(rn7sk.expressions)[which(RN7SK.ANNOT[[cline]]$is_RN7SK == T)] <- "RN7SK"
    rn7sk.expressions <- sort (rn7sk.expressions, decreasing = T)
    colors <- rep("black", length(rn7sk.expressions))
    colors[1] <- "red"
    #pdf(paste(PLOTS, "RN7SK_Expr_", cline,".pdf", sep = ""))
    png(paste(PLOTS, "RN7SK_Expr_", cline,".png", sep = ""))
    barplot(rn7sk.expressions, 
            col = colors, 
            main = paste("Expression levels of 7SK variants in the ", cline, " cell line", sep = ""), 
            xlab = "7SK Variants", 
            ylab = "Expression")
    abline(avg.gene.expr, 0, col = "blue")
    abline(avg.ncrna.expr, 0, col = "purple")
    legend(length(rn7sk.expressions)-7, max(rn7sk.expressions), 
           legend=c("Avg. Gene Expression", "Avg. ncRNA Expression"),
           col=c("blue", "purple"), 
           lty=1, cex=0.8)
    
    dev.off()
    return()
  })
  
}
###
# This function trains xgb models
###
apply.xgboost <- function(train.model.data, test.model.data, model.target){
  
  # Define function to prepare model data
  prepare.model.matrices <- function(model.data, model.target){
    # Training data set
    X = model.data
    Y = X[, model.target] %>% as.matrix()
    X <- within(X, rm(list=colnames(X)[grepl("^target_", colnames(X))]))
    dtypes <- sapply(X, class)
    X <- data.matrix(X)
    X[,dtypes=="factor"][X[,dtypes=="factor"] == 1]  <- 0
    X[,dtypes=="factor"][X[,dtypes=="factor"] == 2]  <- 1
  
    return(list(X=X, Y=Y))
  }
  
  #train.model.data <- model.matrices$HepG2$individual.model.matrices$Termination
  # 0.374
  # 0.377
  if(is.null(test.model.data)){
    
    # n <- ceiling(0.5*dim(train.model.data)[1]/100)
    # target <- train.model.data[, model.target]
    # names(target) <- rownames(train.model.data)
    # 
    # bin.members <- split(target, cut(target, quantile(target, prob = seq(0, 1, 0.01), names = F), include = T))
    # names(bin.members) <- NULL
    # 
    # test.data.samples <-
    # lapply(bin.members, function(bin){
    #   set.seed(7)
    #   sample(bin, n)
    # })
    # test.data.samples <- unlist(test.data.samples)
    # test.data.samples <- names(test.data.samples)
    

    n <- ceiling(0.5*dim(train.model.data)[1])
    set.seed(7)
    test.data.samples <- sample(1:dim(train.model.data)[1], n)
    test.model.data <- train.model.data[test.data.samples, ]
    train.model.data <- train.model.data[-test.data.samples, ]

    # test.model.data <- train.model.data[test.data.samples, ]
    # train.model.data <- train.model.data[!rownames(train.model.data) %in% test.data.samples, ]

  }
  
  # Prepare training data for model training
  train.data <- prepare.model.matrices(train.model.data, model.target)
  X <- train.data$X
  Y <- train.data$Y
  
  #Retrieve samples
  samples <- rownames(train.model.data)
  
  ## Retrain models
  param_list <- list(eta = 0.08,
                     max_depth = 7,
                     gamma = 0.5,
                     lambda = 0.1,
                     alpha = 0.01,
                     colsample_bytree = 0.7,
                     subsample = 0.7, 
                     min_child_weight = 50,
                     booster="gbtree")
  
  ## Retrieve model performances on cross validation hold out data sets 
  set.seed(7)
  cv.model <- xgboost::xgb.cv(data = X,
                              label = Y,
                              params = param_list,
                              nrounds = 120,
                              nfold = 5,
                              objective = "reg:squarederror",
                              prediction = T,
                              callbacks = list(cb.cv.predict(save_models = TRUE)),
                              nthread = CARGS$workers["high"])
  
  # Retrieve CV predictions
  train.preds <- cv.model$pred
  # Calculate root mean squared error of training data
  train.residuals <- as.vector(Y) - train.preds
  train.rmse <- sqrt(mean(train.residuals^2))
  ## Calculate the r-squared of the test data 
  # Calculate the mean of the test data observations
  mean.y.train <- mean(Y)
  # Calculate the total sum of squares
  train.tss =  sum((Y - mean.y.train)^2)
  # Calculate residual sum of squares
  train.rss =  sum(train.residuals^2)
  # Calculate R-squared
  rsq.train  =  1 - (train.rss/train.tss)
  
  #rsq.train <- cor(as.vector(Y), train.preds)^2
  
  
  ## Retrain model to be able to apply to independent test data set
  # Prepare test data for model evaluation
  test.data <- prepare.model.matrices(test.model.data, model.target)
  X.test <- test.data$X
  Y.test <- test.data$Y
  test.data.samples <- rownames(X.test)
  
  set.seed(7)
  model <- xgboost::xgboost(data = X, 
                            label = Y, 
                            params = param_list, 
                            nrounds = 120,
                            verbose = T,
                            objective = "reg:squarederror",
                            nthread = CARGS$workers["high"])
  
  # Make predictions on test data set
  test.preds <-  predict(model, newdata = X.test)
  ## Evaluate model performance
  # Calculate root mean squared error of test data
  test.residuals <- as.vector(Y.test)-test.preds
  test.rmse <- sqrt(mean(test.residuals^2))
  ## Calculate the r-squared of the test data 
  # Calculate the mean of the test data observations
  mean.y.test <- mean(Y.test)
  # Calculate the total sum of squares
  tss =  sum((Y.test - mean.y.test)^2)
  # Calculate residual sum of squares
  rss =  sum(test.residuals^2)
  # Calculate R-squared
  rsq.test  =  1 - (rss/tss)
  
  #rsq.test <- cor(as.vector(Y.test), test.preds)^2
  
  # Retrieve shap values
  shap.values <- shap.values(xgb_model = model, 
                             X_train = X)
  
  ## Set rownames for sub data sets
  names(train.preds) <- samples
  names(test.preds) <- test.data.samples
  rownames(Y) <- samples
  rownames(Y.test) <- test.data.samples
  shap.values$shap_score <- as.data.frame(shap.values$shap_score)
  rownames(shap.values$shap_score) <- samples
  
  
  return(list(model = model, 
              X = X, 
              Y = Y, 
              X.test = X.test, 
              Y.test=Y.test, 
              train.preds = train.preds,
              test.preds = test.preds,
              shap.values = shap.values,
              samples = samples,
              test.data.samples = test.data.samples, 
              train.rsqrd = rsq.train, 
              test.rsqrd = rsq.test))
}
###
# This function evalutes the features of a trained xgb model
###
train.xgb.models <- function(model.matrices, append = ""){
  
  # Check if a precalculated version exists
  target.file <- paste0(OUTPUT, MODEL, "model_evaluation/model.training.results", append,  ".RDS")
  if(!(CARGS$new) & file.exists(target.file)){
    model.results <- readRDS(target.file)
    return(model.results)
  }
  
  # Create table of model performances
  performances <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                           c("cline", "modeltype", "subspace", "target", 
                             "train.rsqrd", "test.rsqrd", "diff", "nfeatures", "nfactors"))
  
  cline.model.results <- 
    lapply(names(model.matrices), function(cline){
      cline.model.matrices <- model.matrices[[cline]]
      cline.model.matrices.results <- 
        lapply(names(cline.model.matrices), function(cline.model.type){
          model.type <- cline.model.matrices[[cline.model.type]]
          sub.space.results <- 
            lapply(names(model.type), function(model.feature.sub.space){
              train.model.data <- model.type[[model.feature.sub.space]]
              targets <- colnames(train.model.data)[grepl("^target_", colnames(train.model.data))]
              target.results <- 
                lapply(targets, function(model.target){
                  print(paste(cline, cline.model.type, model.feature.sub.space, model.target))
                  test.model.data <-  NULL
                  if(grepl("^synchronised", cline.model.type)){
                    cross.cline <- if(cline == "K562") "HepG2" else "K562"
                    test.model.data <- model.matrices[[cross.cline]][[cline.model.type]][[model.feature.sub.space]]
                  }
                  model.results <- apply.xgboost(train.model.data = train.model.data,
                                                 test.model.data = test.model.data,
                                                 model.target = model.target)
                  
                  binding.signals <- colnames(train.model.data)[grepl("^chip|clip", colnames(train.model.data))]
                  nfactors <-  gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*",
                                    "", binding.signals) %>% unique() %>% length()
                  
                  nfeatures <- dim(train.model.data)[2]
                  
                  mean.shap <- mean(model.results$shap.values$mean_shap_score[model.results$shap.values$mean_shap_score>0])
                  mean.shap2 <- mean(model.results$shap.values$mean_shap_score)
                  performances <<- rbind(performances, data.frame(cline = cline,
                                                                  modeltype = cline.model.type,
                                                                  subspace = model.feature.sub.space,
                                                                  target = model.target,
                                                                  train.rsqrd = model.results$train.rsqrd,
                                                                  test.rsqrd = model.results$test.rsqrd,
                                                                  diff = model.results$train.rsqrd - model.results$test.rsqrd,
                                                                  nfeatures = nfeatures,
                                                                  nfactors = nfactors,
                                                                  train.power = 100*(model.results$train.rsqrd/nfeatures),
                                                                  test.power = 100*(model.results$test.rsqrd/nfeatures),
                                                                  mean.shap = mean.shap,
                                                                  mean.shap2 = mean.shap2))
                  return(model.results)
                })
              names(target.results) <- targets
              target.results
            })
          names(sub.space.results) <- names(model.type)
          sub.space.results
        })
      names(cline.model.matrices.results) <- names(cline.model.matrices)
      cline.model.matrices.results
    })
  names(cline.model.results) <- names(model.matrices)
  
  results <- list(models = cline.model.results, performances = performances)
  saveRDS(results, target.file)
  saveRDS(performances, paste0(OUTPUT, MODEL, "model_evaluation/model.training.performances", append, ".RDS"))
  write.csv(performances, paste0(OUTPUT, MODEL, "model_evaluation/model.training.results", append, ".csv"))
  return(results)
}
###
# This function evalutes the features of a trained xgb model
###
train.randomized.xgb.models <- function(){
  
  # Check if a precalculated version exists
  target.file <- paste(OUTPUT, MODEL, "model_evaluation/randomized.model.perf.RDS", sep = "")
  if(!(CARGS$new) & file.exists(target.file)){
    performances <- readRDS(target.file)
    return(performances)
  }
  
  cline <- "K562"
  model.type <- "randomized"
  subspace <-  "All"
  mat <- model.matrices$K562$individual.model.matrices$All
  non.binding.signals <- colnames(mat)[!grepl("^chip|clip|target", colnames(mat))]
  
  binding.features <- colnames(mat)[grepl("^chip|clip", colnames(mat))]
  all.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*",
                      "", binding.features) %>% unique()
  model.target <- "target_traveling.ratio"
  
  performances <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), 
                           c("cline", "modeltype", "subspace", "target", 
                             "train.rsqrd", "test.rsqrd", "nfeatures", "nfactors"))
  
  # Draw random number of random factors and train model
  seeds <- 1:100
  model.results <- 
    lapply(seeds, function(seed){
      print(seed)
      set.seed(seed)
      random.factors <- sample(all.factors, sample(1:length(all.factors), 1), replace = F)
      random.factors.pattern <- paste0(random.factors, collapse="|")
      random.binding.features <- binding.features[grepl(random.factors.pattern, binding.features)]
      
      randomized.vectors <- 
        lapply(random.binding.features, function(binding.feature){
          freq <- table(mat[, binding.feature ])
          if(length(freq) != 2){return(NULL)}
          n <- sum(freq)
          random.binding <- vector("numeric", n)
          set.seed(seed)
          random.binding[sample(1:n, freq["1"])] <- 1
          random.binding <- data.frame(random.binding)
          colnames(random.binding) <- binding.feature
          random.binding
        })
      randomized.vectors <- do.call(cbind, randomized.vectors)
      
      # Randomized non binding features
      random.sequence.features <-
        lapply(non.binding.signals, function(feature){
          set.seed(seed)
          random.feat <- sample(mat[ ,feature], size = dim(mat)[1], replace = F)
          random.feat
        })
      random.features <- as.data.frame(do.call(cbind, random.sequence.features))
      colnames(random.features) <- non.binding.signals
      rownames(random.features) <- rownames(mat)
      
      train.model.data <- cbind( target_traveling.ratio=mat[, model.target], random.features, randomized.vectors)
      test.model.data <-  NULL
      model.results <- apply.xgboost(train.model.data = train.model.data, 
                                     test.model.data = test.model.data,
                                     model.target=model.target)
      
      binding.signals <- colnames(train.model.data)[grepl("^chip|clip", colnames(train.model.data))]
      nfactors <- length(random.factors)
      
      nfeatures <- dim(train.model.data)[2]
      
      #mean.shap <- mean(model.results$shap.values$mean_shap_score[model.results$shap.values$mean_shap_score>0])
      mean.shap <- mean(model.results$shap.values$mean_shap_score)
      performances <<- rbind(performances, data.frame(cline = cline, 
                                                      modeltype = model.type, 
                                                      subspace = seed, 
                                                      target = model.target, 
                                                      train.rsqrd = model.results$train.rsqrd, 
                                                      test.rsqrd = model.results$test.rsqrd, 
                                                      nfeatures = nfeatures,
                                                      nfactors = nfactors,
                                                      mean.shap = mean.shap))
      return(model.results)
      
    })
  
  saveRDS(performances, target.file)
  return(performances)
}

###
# This function evalutes the features of a trained xgb model
###
evaluate.feature.effects2 <- function(model.results, feature.subspace, model.target, plot.base.size =  16){
  
  # #X TEMP vars to test function
  model.target="target_traveling.ratio"
  matrix.type <- "individual.model.matrices"
  feature.subspace <- "All"
  model.target.type <- gsub("target_", "", model.target)
  train.cline <- "K562"
  test.cline <- "HepG2"
  plot.base.size =  16
  
  ## Draw venn for gene overlap between cell lines
  png(paste0("plots/", "paper_figure_", train.cline, "_and_", test.cline, "_gene_counts_", ".png"), res = 300)
  
  gene.venn <- 
  draw.pairwise.venn(length(model.matrices[[train.cline]][[matrix.type]][[feature.subspace]][[model.target]]),
                     length(model.matrices[[test.cline]][[matrix.type]][[feature.subspace]][[model.target]]),
                     length(intersect(rownames(model.matrices[[train.cline]][[matrix.type]][[feature.subspace]]),
                                      rownames(model.matrices[[test.cline]][[matrix.type]][[feature.subspace]]))),
                     category = c(train.cline, test.cline),
                     scaled = F,
                     cat.fontface = c("bold", "bold"),
                     fill = c(alpha("#440154ff",1), alpha('#21908dff',1)))
  dev.off()
  
  evaluation.results <- 
  lapply(CARGS$cell.line, function(train.cline){
     lapply(c("individual.model.matrices", "synchronised.model.matrices"), function(matrix.type){
      
  # train.cline <- "K562"
  # matrix.type <- "individual.model.matrices"
  # 
  # evaluation.results <-
  # lapply(c("K562"), function(train.cline){
  #   lapply(c("individual.model.matrices"), function(matrix.type){
      # Load relevant model results
      select.model.results <- model.results$models[[train.cline]][[matrix.type]][[feature.subspace]][[model.target]]
      
      # Extract relevant variables from model results
      X <- select.model.results$X
      Y <- select.model.results$Y
      X.test = select.model.results$X.test
      Y.test <- select.model.results$Y.test
      train.preds <- select.model.results$train.preds
      test.preds <- select.model.results$test.preds
      test.cline <- if(train.cline=="K562") "HepG2" else "K562"
      samples <- rownames(X)
      
      ## Individual hap scores
      # The sum of each row’s SHAP values (plus the BIAS column, which is like an intercept) 
      # is the predicted model output,i.e., the explanation’s attribution values sum up to the model output
      # but not true for CV models
      shap.values <- select.model.results$shap.values
      shap.scores <- as.data.frame(shap.values$shap_score)
      rownames(shap.scores) <- rownames(X)
      shap.bias <- shap.values$BIAS0
      # Mean |SHAP| scores for each feature over all samples
      mean.shap.scores <- shap.values$mean_shap_score
      
      # Plot obs vs pred prediction performances
      print(paste0("Model Performances ", train.cline, " ", matrix.type))
      model.performance <- 
        visualize.model.performance(select.model.results, 
                                    feature.subspace, 
                                    model.target.type, 
                                    plot.base.size = plot.base.size,
                                    train.cline = train.cline)
      
      # Plot obs target differences between cell lines against differences of predictions between cell lines
      print(paste0("Cross Cell Type Model Performances ", train.cline, " ", matrix.type))
      sync.select.model.results <- model.results$models[[train.cline]][["synchronised.model.matrices"]][[feature.subspace]][[model.target]]
      sync.reverse.model.results <- model.results$models[[test.cline]][["synchronised.model.matrices"]][[feature.subspace]][[model.target]]
      cell.type.specific.sample.model.performances <- 
        retrieve.cell.type.specific.sample.prediction.performances(sync.select.model.results, 
                                                                  sync.reverse.model.results,
                                                                  train.cline, 
                                                                  test.cline,
                                                                  matrix.type, 
                                                                  feature.subspace, 
                                                                  model.target.type, 
                                                                  plot.base.size)
      ## Model performance on cross cell type data
      print(paste0("Cross application performance ", train.cline, " ", matrix.type))
      sync.model.performance <- 
        visualize.model.performance(sync.select.model.results, 
                                    feature.subspace,
                                    model.target, 
                                    plot.base.size=plot.base.size,train.cline)
      
      ## Individual Feature distribution of top 25 features
      print(paste0("SHAP summary plot ", train.cline, " ", matrix.type))
      shap.summary.plot <- 
        visualize.shap.summary(shap.score.matrix = shap.scores, 
                               X= X, 
                               train.cline = train.cline, 
                               matrix.type = matrix.type, 
                               feature.subspace = feature.subspace, 
                               model.target = model.target,
                               plot.base.size = 12)
      
      ## DNA sequence feature distributions
      print(paste0("DNA sequence feature plots ", train.cline, " ", matrix.type))
      static.features.plots <- 
        visualize.static.features(shap.scores, 
                                  X= X, 
                                  train.cline = train.cline, 
                                  matrix.type = matrix.type, 
                                  feature.subspace = feature.subspace, 
                                  model.target = model.target,
                                  plot.base.size = plot.base.size)
    
      print(paste0("Binding feature plots ", train.cline, " ", matrix.type))
      dynamic.features.plots <- 
        visualize.dynamic.features(shap.scores = shap.scores, 
                                   X = X, 
                                   train.cline = train.cline, 
                                   matrix.type = matrix.type, 
                                   feature.subspace = feature.subspace, 
                                   model.target = model.target,
                                   plot.base.size = plot.base.size, 
                                   mean.shap.scores,
                                   append = "")
      
      dynamic.features.plots$aggregate.non.binding.feature.contrib.plot <- dynamic.features.plots$aggregate.non.binding.feature.contrib.plot+
        ggtitle("Contributions of gene annotation and composition features")
      cowplot::save_plot(paste0("plots/supplementary_figure_static_feature_contribs_", train.cline,"_", matrix.type, ".png"),
                         plot = dynamic.features.plots$aggregate.non.binding.feature.contrib.plot,
                         ncol = 1,
                         nrow = 1, 
                         dpi = 600)
      
      
      # Check enrichment of DNA intron binders
      # Retrieve top itron DNA binding proteins
      dna.intron.binding.contribs <- dynamic.features.plots$class.contributions$DNA_introns
      dna.intron.binding.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*",
                                          "", names(dna.intron.binding.contribs))
      names(dna.intron.binding.contribs) <- dna.intron.binding.factors
      # Retrieve functional sets
      sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
      sub.factor.feature.spaces <-
        build.sequence.specific.biologically.functional.feature.sub.spaces(stratify = sub.factor.feature.spaces)
      sub.factor.feature.spaces <- sub.factor.feature.spaces[!grepl("ss|All|sdom", names(sub.factor.feature.spaces))]
      sub.factor.feature.spaces <- sub.factor.feature.spaces[c("Chromatin", "Initiation", "Known", "All.7SK", "Pausing", "Elongation", "Termination","Splicing", "Processing", "Translation")]

      dna.intron.binders.functional.contribs <-
      lapply(sub.factor.feature.spaces, function(sub.space){
        factors <- names(dna.intron.binding.contribs)[names(dna.intron.binding.contribs) %in% sub.space[[train.cline]]]
        #sum(dna.intron.binding.contribs[factors])
      })
      #dna.intron.binders.functional.contribs <- dna.intron.binders.functional.contribs[c("Chromatin", "Initiation", "Known", "All.7SK", "Pausing", "Elongation", "Termination","Splicing", "Processing", "Translation")]
      dna.intron.binders.functional.contribs <- unlist(dna.intron.binders.functional.contribs)
      dna.intron.binders.functional.contribs <- data.frame(Process = names(dna.intron.binders.functional.contribs), Contribution = dna.intron.binders.functional.contribs)
      dna.intron.binders.functional.contribs <- na.omit(dna.intron.binders.functional.contribs)
      dna.intron.binders.functional.contribs <- dna.intron.binders.functional.contribs[order(dna.intron.binders.functional.contribs$Contribution, decreasing = F),]
      dna.intron.contrib.plot <-
        ggbarplot(dna.intron.binders.functional.contribs,
                  x = "Process",
                  y = "Contribution",
                  fill = rgb(0.2,0.4,0.6,0.6),
                  #label = paste0(paste0("Rsqrd. ", sub.model.performances$test.rsqrd), "\n#Factors ", sub.model.performances$nfactors, "\nContrib. ", format(round(sub.model.performances$mean.shap, 3), nsmall=3) ),
                  #label = paste0(sub.model.performances$nfactors),
                  rotate = T,
                  #lab.pos = "in",
                  lab.col = "black",
                  lab.font = "bold",
                  lab.size = 4,
                  lab.hjust = -0.5,
                  lab.vjust = 0.5,
                  position = position_dodge(0.9))+
        #facet_grid(~sequence.specific)+
        xlab("Process")+
        ylab("Contribution")+
        theme_pubr(base_size = plot.base.size, legend = "right")+
        #scale_fill_material("grey")+
        scale_fill_material("blue-grey")
      cowplot::save_plot(paste0("plots/supplementary_figure_dna_intronbindings_", train.cline,"_", matrix.type, ".png"),
                         plot = dna.intron.contrib.plot,
                         ncol = 1,
                         nrow = 1, 
                         dpi = 600)
      
      
      
      rna.intron.binding.contribs <- dynamic.features.plots$class.contributions$RNA_introns
      rna.intron.binding.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*",
                                         "", names(rna.intron.binding.contribs))
      names(rna.intron.binding.contribs) <- rna.intron.binding.factors
      # Retrieve functional sets
      sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
      sub.factor.feature.spaces <- 
        build.sequence.specific.biologically.functional.feature.sub.spaces(stratify = sub.factor.feature.spaces)
      sub.factor.feature.spaces <- sub.factor.feature.spaces[!grepl(":ss|All|sdom", names(sub.factor.feature.spaces))]
      sub.factor.feature.spaces <- sub.factor.feature.spaces[c("Chromatin", "Initiation", "Elongation", "7SK.Binding", "Elongation+7SK", "Termination","Splicing", "Processing", "Translation")]
      
      rna.intron.binders.functional.contribs <- 
        lapply(sub.factor.feature.spaces, function(sub.space){
          factors <- names(rna.intron.binding.contribs)[names(rna.intron.binding.contribs) %in% sub.space[[train.cline]]]
          sum(rna.intron.binding.contribs[factors])
        })
      #dna.intron.binders.functional.contribs <- dna.intron.binders.functional.contribs[c("Chromatin", "Initiation", "Known", "All.7SK", "Pausing", "Elongation", "Termination","Splicing", "Processing", "Translation")]
      rna.intron.binders.functional.contribs <- unlist(rna.intron.binders.functional.contribs)
      rna.intron.binders.functional.contribs <- data.frame(Process = names(rna.intron.binders.functional.contribs), Contribution = rna.intron.binders.functional.contribs)
      rna.intron.binders.functional.contribs <- na.omit(rna.intron.binders.functional.contribs)
      rna.intron.binders.functional.contribs <- rna.intron.binders.functional.contribs[order(rna.intron.binders.functional.contribs$Contribution, decreasing = F),]
      rna.intron.contrib.plot <- 
        ggbarplot(rna.intron.binders.functional.contribs, 
                  x = "Process", 
                  y = "Contribution",
                  fill = rgb(0.2,0.4,0.6,0.6), 
                  #label = paste0(paste0("Rsqrd. ", sub.model.performances$test.rsqrd), "\n#Factors ", sub.model.performances$nfactors, "\nContrib. ", format(round(sub.model.performances$mean.shap, 3), nsmall=3) ),
                  #label = paste0(sub.model.performances$nfactors),
                  rotate = T,
                  #lab.pos = "in",
                  lab.col = "black",
                  lab.font = "bold",
                  lab.size = 4,
                  lab.hjust = -0.5,
                  lab.vjust = 0.5,
                  position = position_dodge(0.9))+
        #facet_grid(~sequence.specific)+
        xlab("Process")+
        ylab("Contribution")+
        theme_pubr(base_size = plot.base.size, legend = "right")+
        #scale_fill_material("grey")+
        scale_fill_material("blue-grey")#+
      #scale_fill_uchicago()+
      #theme(panel.spacing = unit(2, "lines"))+
      #scale_y_continuous(limits = c(0, 0.8), breaks = seq(0,0.8, 0.1))
      cowplot::save_plot(paste0("plots/supplementary_figure_rna_intronbindings_", train.cline,"_", matrix.type, ".png"),
                         plot = rna.intron.contrib.plot,
                         ncol = 1,
                         nrow = 1, 
                         dpi = 600)
      
      
      all.rna.binding.factors <- lapply(sub.factor.feature.spaces, function(process){process[train.cline]})
      all.rna.binding.factors <- unique(unlist(all.rna.binding.factors))
      splicing.factors <- unique(sub.factor.feature.spaces[["Splicing"]][[train.cline]])
      rna.intron.binding.factors <-  
        lapply(sub.factor.feature.spaces, function(sub.space){
          factors <- names(rna.intron.binding.contribs)[names(rna.intron.binding.contribs) %in% sub.space[[train.cline]]]
        })
      rna.intron.binding.splicing.factors <- rna.intron.binding.factors[["Splicing"]]
      rna.intron.binding.factors <- unlist(rna.intron.binding.factors)

        
      non.splicing <- all.rna.binding.factors[!(all.rna.binding.factors %in% splicing.factors)]
      
      categories <- list(c("splicing", "non-splicing"), c("binding", "non-binding"))
      counts <- c(length(rna.intron.binding.splicing.factors),
                  sum(!splicing.factors %in% rna.intron.binding.factors),
                  sum(non.splicing %in% rna.intron.binding.factors),
                  sum(!non.splicing %in% rna.intron.binding.factors))
      rna.inton.splicing.enrichment.results <- perform.fisher.test(categories = categories, 
                                                                    counts = counts)
      
      
      # tob.intron.binders <- dna.intron.binding.contribs[dna.intron.binding.contribs >= quantile(dna.intron.binding.contribs, seq(0,1,0.05))["95%"]]
      # 
      # binding.features <- colnames(X)
      # all.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*",
      #                     "", binding.features) %>% unique()
      # universe <- factor(as.integer(all.factors %in% names(tob.intron.binders)))
      # names(universe) <- all.factors
      # dna.intron.binding.factor.gsea <- perform.gene.set.enrichment.analysis(universe, node.size = 50)
      
      # top.directional.effectors <- dynamic.features.plots$top.directional.effectors
      # individual.directional.factor.contrib.plots <- 
      #   visualize.individual.directional.dynamic.factor.effects(top.directional.effectors = top.directional.effectors,
      #                                                           shap.scores, 
      #                                                           X = X, 
      #                                                           train.cline = train.cline, 
      #                                                           matrix.type = matrix.type, 
      #                                                           feature.subspace = feature.subspace, 
      #                                                           model.target = model.target,
      #                                                           plot.base.size = plot.base.size)
      
      print(paste0("Force plots ", train.cline, " ", matrix.type))
      shap.force.plots <-
        visualize.shap.forces(shap.scores = shap.scores, 
                              mean.shap.scores = mean.shap.scores,
                              X = X,
                              train.cline = train.cline,
                              matrix.type = matrix.type,
                              feature.subspace = feature.subspace,
                              model.target = model.target,
                              plot.base.size = plot.base.size)
      
      
      
      ## Featre class and prior knowledge contributions
      prior.knowledge.contributions <-
        calculate.prior.knowledge.contributions(shap.scores = shap.scores,
                                                mean.shap.scores = mean.shap.scores,
                                                X = X,
                                                train.cline = train.cline,
                                                matrix.type = matrix.type,
                                                feature.subspace = feature.subspace,
                                                model.target = model.target,
                                                plot.base.size = plot.base.size)
      
      print(paste0("Cluster effectors plot ", train.cline, " ", matrix.type))
      shap.values <- select.model.results$shap.values
      shap.scores <- as.data.frame(shap.values$shap_score)
      rownames(shap.scores) <- rownames(X)
      shap.bias <- shap.values$BIAS0
      # Mean |SHAP| scores for each feature over all samples
      mean.shap.scores <- shap.values$mean_shap_score
      cluster.effectors <- visualize.cluster.factors(force.plot_data = shap.force.plots$force.plot_data,
                                                     shap.scores = shap.scores, 
                                                     X = X, 
                                                     Y=Y)
      #X
      # cluster.target.distribution.plots <- 
      #   visualize.shap.cluster.target.distributions(force.plot_data = shap.force.plots$force.plot_data,
      #                                               shap.scores, 
      #                                               X=X,
      #                                               train.cline = train.cline,
      #                                               matrix.type = matrix.type,
      #                                               feature.subspace = feature.subspace,
      #                                               model.target = model.target,
      #                                               plot.base.size = plot.base.size,
      #                                               optimal.factors = NULL)
      
      
      # print(paste0("Target deconvolution plot ", train.cline, " ", matrix.type))
      # deconvoluted.target.plots <- 
      #   visualize.target.deconvolution(force.plot_data = shap.force.plots$force.plot_data, 
      #                                  Y = Y,
      #                                  train.cline = train.cline,
      #                                  matrix.type = matrix.type,
      #                                  feature.subspace = feature.subspace,
      #                                  model.target = model.target,
      #                                  plot.base.size = plot.base.size, 
      #                                  optimal.factors = NULL)
      # 
      # print(paste0("Cluster 1 optimal factors plot ", train.cline, " ", matrix.type))
      # cluster1.optimal.factors.plot <- 
      #   retrieve.optimal.number.of.factors(cum.factor.contrib = cluster.effectors$undirectional.cluster.factor.contrib[[1]],
      #                                      train.cline = train.cline,
      #                                      matrix.type = matrix.type,
      #                                      feature.subspace = feature.subspace,
      #                                      model.target = model.target,
      #                                      plot.base.size = plot.base.size)
      # 
      # print(paste0("Cluster 2 optimal factors plot ", train.cline, " ", matrix.type))
      # cluster2.optimal.factors.plot <- 
      #   retrieve.optimal.number.of.factors(cum.factor.contrib = cluster.effectors$undirectional.cluster.factor.contrib[[2]],
      #                                      train.cline = train.cline,
      #                                      matrix.type = matrix.type,
      #                                      feature.subspace = feature.subspace,
      #                                      model.target = model.target,
      #                                      plot.base.size = plot.base.size)
      # 
      # #optimal.factors <- unique(c(cluster1.optimal.factors.plot$optimal.factors, cluster2.optimal.factors.plot$optimal.factors))
      # optimal.cluster.factors <- list(cluster1.optimal.factors.plot$optimal.factors, cluster2.optimal.factors.plot$optimal.factors)
      # all.optimal.factors <- unique(unlist(optimal.cluster.factors))
      # 
      # print(paste0("Optimal factor model ", train.cline, " ", matrix.type))
      # minimal.model.results <- train.minimal.model(all.optimal.factors)
      

      
      optimal.factors.plot <-
        retrieve.optimal.number.of.factors(cum.factor.contrib = dynamic.features.plots$cum.factor.contrib,
                                           train.cline = train.cline,
                                           matrix.type = matrix.type,
                                           feature.subspace = feature.subspace,
                                           model.target = model.target,
                                           plot.base.size = plot.base.size)
      all.optimal.factors <- optimal.factors.plot$optimal.factors
      minimal.model.results <- train.minimal.model(all.optimal.factors)
      
      cowplot::save_plot(paste0("plots/supplementary_minimal_model_", train.cline,"_", matrix.type, ".png"),
                         plot = minimal.model.results$minimal.model.perf,
                         ncol = 1,
                         nrow = 1, 
                         dpi = 600)
      
      
      
      shap.force.plots.minimal.model <-
        visualize.shap.forces(shap.scores = minimal.model.results$select.feature.model$shap.values$shap_score, 
                              mean.shap.scores = minimal.model.results$select.feature.model$shap.values$mean_shap_score,
                              X = minimal.model.results$select.feature.model$X,
                              train.cline = train.cline,
                              matrix.type = matrix.type,
                              feature.subspace = feature.subspace,
                              model.target = model.target,
                              plot.base.size = plot.base.size)
      
      print(paste0("Optimal factor target deconvolution ", train.cline, " ", matrix.type))
      optimal.factor.target.deconvolution <- 
        visualize.shap.cluster.target.distributions(force.plot_data = shap.force.plots.minimal.model$force.plot_data,
                                                    as.data.frame(minimal.model.results$select.feature.model$shap.values$shap_score), 
                                                    X=minimal.model.results$select.feature.model$X,
                                                    Y = minimal.model.results$select.feature.model$Y, 
                                                    train.cline = train.cline,
                                                    matrix.type = matrix.type,
                                                    feature.subspace = feature.subspace,
                                                    model.target = model.target,
                                                    plot.base.size = plot.base.size,
                                                    optimal.factors = all.optimal.factors)
      
      deconvoluted.target.plots.minimal.model <- 
        visualize.target.deconvolution(force.plot_data = shap.force.plots.minimal.model$force.plot_data, 
                                       Y = minimal.model.results$select.feature.model$Y,
                                       train.cline = train.cline,
                                       matrix.type = matrix.type,
                                       feature.subspace = feature.subspace,
                                       model.target = model.target,
                                       plot.base.size = plot.base.size, 
                                       optimal.factors = NULL)
      
      #visualize.shap.interactions(minimal.model.results$select.feature.model$model, minimal.model.results$select.feature.model$X)
      
      # optimal.factor.target.deconvolution <- 
      #   visualize.shap.cluster.target.distributions(force.plot_data = shap.force.plots$force.plot_data,
      #                                               shap.scores, 
      #                                               X=X,
      #                                               Y = Y, 
      #                                               train.cline = train.cline,
      #                                               matrix.type = matrix.type,
      #                                               feature.subspace = feature.subspace,
      #                                               model.target = model.target,
      #                                               plot.base.size = plot.base.size,
      #                                               optimal.factors = all.optimal.factors)
      
      print(paste0("Optimal factor directional contribution ", train.cline, " ", matrix.type))
      optimal.cluster.factor.directional.contrib <-
        visualize.optimal.factor.contributions( force.plot_data = shap.force.plots$force.plot_data,
                                                shap.scores = shap.scores,
                                               optimal.factors = all.optimal.factors,
                                               cluster.membership = all.optimal.factors ,
                                               X = X)
      
      print(paste0("Optimal factor individual contribution ", train.cline, " ", matrix.type))
      select.individual.directional.factor.contrib.plots <- 
        visualize.individual.directional.dynamic.factor.effects(top.directional.effectors = optimal.cluster.factor.directional.contrib$top.directional.binding.factors,
                                                                shap.scores, 
                                                                X = X, 
                                                                train.cline = train.cline, 
                                                                matrix.type = matrix.type, 
                                                                feature.subspace = feature.subspace, 
                                                                model.target = model.target,
                                                                plot.base.size = plot.base.size)

      
      #sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
      #all.pausing.related.factors <- unique(c(all.optimal.factors, unlist(sub.factor.feature.spaces$Pausing)))
      #final.pausing.model <- train.minimal.model(all.pausing.related.factors)
      
      print(paste0("Optimal factor enrichment sets ", train.cline, " ", matrix.type))
      contingency.test.results <- 
        lapply(optimal.cluster.factors, function(optimal.factors){
          ## Perform an enrichment test of optimal factors as NSS/SS factors
          binding.features <- names(shap.values$mean_shap_score)[grepl("^chip|clip", names(shap.values$mean_shap_score))]
          all.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*",
                              "", binding.features) %>% unique()
          
          sequence.specific <- SEQ.SPEC$sequence.specific[SEQ.SPEC$sequence.specific %in% all.factors]
          nonsequence.specific <- SEQ.SPEC$nonsequence.specific[SEQ.SPEC$nonsequence.specific %in% all.factors]
          non.optimal <- all.factors[!(all.factors %in% optimal.factors)]
          categories <- list(c("optimal", "non-optimal"), c("nss", "ss"))
          counts <- c(sum(optimal.factors %in% nonsequence.specific),
                      sum(optimal.factors %in% sequence.specific),
                      sum(non.optimal %in% nonsequence.specific),
                      sum(non.optimal %in% sequence.specific))
          sequence.speificity.enrichment.results <- perform.fisher.test(categories = categories, 
                                                                        counts = counts)
          
          ## Perform an enrichment test for optimal factors as 7SK binders/associated proteins
          sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
          rn7sk.assoc.factors <- sub.factor.feature.spaces$All.7SK[[train.cline]]
          rn7sk.assoc.factors <- rn7sk.assoc.factors[rn7sk.assoc.factors %in% all.factors]
          non.7sk.assoc.factors <- all.factors[!(all.factors %in% rn7sk.assoc.factors)]
          non.optimal <- all.factors[!(all.factors %in% optimal.factors)]
          categories <- list(c("optimal", "non-optimal"), c("7SK-associated", "non-7SK-associated"))
          counts <- c(sum(optimal.factors %in% rn7sk.assoc.factors),
                      sum(optimal.factors %in% non.7sk.assoc.factors),
                      sum(non.optimal %in% rn7sk.assoc.factors),
                      sum(non.optimal %in% non.7sk.assoc.factors))
          rn7sk.association.enrichment.results <- perform.fisher.test(categories = categories, 
                                                                      counts = counts)
          
          ## Perform an enrichment test for optimal factors as 7SK binders/associated proteins
          sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
          pausing.factors <- sub.factor.feature.spaces$Pausing[[train.cline]]
          pausing.factors <- pausing.factors[pausing.factors %in% all.factors]
          non.pausing <- all.factors[!(all.factors %in% pausing.factors)]
          non.optimal <- all.factors[!(all.factors %in% optimal.factors)]
          categories <- list(c("optimal", "non-optimal"), c("Pausing", "Non-Pausing"))
          counts <- c(sum(optimal.factors %in% pausing.factors),
                      sum(optimal.factors %in% non.pausing),
                      sum(non.optimal %in% pausing.factors),
                      sum(non.optimal %in% non.pausing))
          pausing.enrichment.results <- perform.fisher.test(categories = categories, 
                                                                      counts = counts)
          
          universe <- factor(as.integer(all.factors %in% optimal.factors))
          names(universe) <- all.factors
          optimal.factor.gsea <- perform.gene.set.enrichment.analysis(universe, node.size = 50)
          
          
          return(list(sequence.speificity.enrichment.results = sequence.speificity.enrichment.results,
                      rn7sk.association.enrichment.results = rn7sk.association.enrichment.results,
                      optimal.factor.gsea = optimal.factor.gsea))
        })
      names(contingency.test.results) <- c("cluster1", "cluster2")
      

      ## Enrichment test between clusters
      #cluster1.effectors.list <- cluster.effectors$undirectional.cluster.factor.contrib[[1]][1:30, "factor"]
      #cluster2.effectors.list <- cluster.effectors$undirectional.cluster.factor.contrib[[2]][1:30, "factor"]
      # cluster1.effectors.list
      
      
      ## Perform an enrichment test of optimal factors as NSS/SS factors
      # binding.features <- names(shap.values$mean_shap_score)[grepl("^chip|clip", names(shap.values$mean_shap_score))]
      # all.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*",
      #                     "", binding.features) %>% unique()
      # 
      # cluster1.effectors.list <- optimal.cluster.factors[[1]]
      # cluster2.effectors.list <- optimal.cluster.factors[[2]]
      # 
      # sequence.specific <- SEQ.SPEC$sequence.specific[SEQ.SPEC$sequence.specific %in% all.factors]
      # nonsequence.specific <- SEQ.SPEC$nonsequence.specific[SEQ.SPEC$nonsequence.specific %in% all.factors]
      # non.optimal <- all.factors[!(all.factors %in% optimal.factors)]
      # categories <- list(c("optimal", "non-optimal"), c("nss", "ss"))
      # counts <- c(sum(cluster1.effectors.list %in% nonsequence.specific),
      #             sum(cluster1.effectors.list %in% sequence.specific),
      #             sum(cluster2.effectors.list %in% nonsequence.specific),
      #             sum(cluster2.effectors.list %in% sequence.specific))
      # sequence.speificity.enrichment.results <- perform.fisher.test(categories = categories, 
      #                                                               counts = counts)
      # 
      # sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
      # rn7sk.assoc.factors <- sub.factor.feature.spaces$All.7SK[[train.cline]]
      # rn7sk.assoc.factors <- rn7sk.assoc.factors[rn7sk.assoc.factors %in% all.factors]
      # non.7sk.assoc.factors <- all.factors[!(all.factors %in% rn7sk.assoc.factors)]
      # categories <- list(c("cluster1", "cluster2"), c("7SK-associated", "non-7SK-associated"))
      # counts <- c(sum(cluster1.effectors.list %in% rn7sk.assoc.factors),
      #             sum(cluster1.effectors.list %in% non.7sk.assoc.factors),
      #             sum(cluster2.effectors.list %in% rn7sk.assoc.factors),
      #             sum(cluster2.effectors.list %in% non.7sk.assoc.factors))
      # rn7sk.association.cluster.enrichment.results <- perform.fisher.test(categories = categories, 
      #                                                             counts = counts)
      # 
      # sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
      # pausing.factors <- sub.factor.feature.spaces$Pausing[[train.cline]]
      # pausing.factors <- pausing.factors[pausing.factors %in% all.factors]
      # non.pausing <- all.factors[!(all.factors %in% pausing.factors)]
      # categories <- list(c("cluster1", "cluster2"), c("Pausing", "Non-Pausing"))
      # counts <- c(sum(cluster1.effectors.list %in% pausing.factors),
      #             sum(cluster1.effectors.list %in% non.pausing),
      #             sum(cluster2.effectors.list %in% pausing.factors),
      #             sum(cluster2.effectors.list %in% non.pausing))
      # pausing.enrichment.results <- perform.fisher.test(categories = categories, 
      #                                                   counts = counts)
      # 
      # dbps <- setdiff(SEQ.SPEC$dbp.factors[[train.cline]], SEQ.SPEC$rbp.factors[[train.cline]])
      # rbps <- setdiff(SEQ.SPEC$rbp.factors[[train.cline]], SEQ.SPEC$dbp.factors[[train.cline]])
      # dbps <- dbps[dbps %in% all.factors]
      # rbps <- rbps[rbps %in% all.factors]
      # categories <- list(c("cluster1", "cluster2"), c("RNA-binding", "Non-RNA-binding"))
      # counts <- c(sum(cluster1.effectors.list %in% rbps),
      #             sum(cluster1.effectors.list %in% dbps),
      #             sum(cluster2.effectors.list %in% rbps),
      #             sum(cluster2.effectors.list %in% dbps))
      # binding.mode.enrichment.results <- perform.fisher.test(categories = categories, 
      #                                                   counts = counts)
      
      
      ## Compare traveling ratio and transcript expression
      # target.model.results <- model.results$models[[train.cline]][[matrix.type]][[feature.subspace]]
      # target.correlations <- 
      #   visualize.target.correlations(target.model.results,
      #                                 train.cline = train.cline, 
      #                                 matrix.type = matrix.type, 
      #                                 feature.subspace = feature.subspace, 
      #                                 model.target = model.target,
      #                                 plot.base.size = plot.base.size)
      
      ## Compare forward and reverse model feature contributions
      # reverse.model.feature.correlations <- 
      #   visualize.reverse.model.feature.contributions(optimal.factors.plot$optimal.factors, 
      #                                                 dynamic.features.plots$directional.factor.contrib, 
      #                                                 model.results, X, train.cline, matrix.type, feature.subspace, model.target, plot.base.size)
      # 
      # 
      # reverse.target.model.dynamic.features.plots <- 
      #   visualize.dynamic.features(reverse.model.feature.correlations$reverse.mdoel$shap.values$shap_score, 
      #                              X = reverse.model.feature.correlations$reverse.mdoel$X, 
      #                              train.cline = train.cline, 
      #                              matrix.type = matrix.type, 
      #                              feature.subspace = feature.subspace, 
      #                              model.target = "target_transcript.expression",
      #                              plot.base.size = plot.base.size)
      
      # knockdown.correlations.plot <- 
      #   visualize.knockdown.correlations(optimal.factors.plot$optimal.factors, model.results)
      
      ## Perform gene set enrichment analysis of top factors
      # binding.features <- names(shap.values$mean_shap_score)[grepl("^chip|clip", names(shap.values$mean_shap_score))]
      # all.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*",
      #                     "", binding.features) %>% unique()
      # clip.targets <- lapply(eCLIPseq, function(peaks){
      #   as.character(unique(peaks$hgnc_symbol))
      # })
      # chip.targets <- lapply(CHIPseq, function(peaks){
      #   as.character(unique(peaks$hgnc_symbol))
      # })
      # all.factors <- unique(c((clip.targets[[train.cline]]), (chip.targets[[train.cline]])))
      
      print(paste0("Prior knowledge model array ", train.cline, " ", matrix.type))
      go.term.model.results.plot <- 
        visualize.go.term.model.results(model.results = model.results, 
                                        train.cline = train.cline, 
                                        matrix.type = matrix.type, 
                                        model.target = model.target, 
                                        plot.base.size = plot.base.size, 
                                        random.model.results = random.model.results)
      
      # go.term.model.results.plot$model.performance.plot <- 
      # go.term.model.results.plot$model.performance.plot+
      #   theme( legend.position = "top", plot.title = element_text(size=10))
      # prior.knowledge.contributions <- 
      # prior.knowledge.contributions+
      #   theme( legend.position = "top", plot.title = element_text(size=10))
      # 
      # 
      # model.perf.grid <-  
      #   plot_grid(go.term.model.results.plot$model.performance.plot,
      #             prior.knowledge.contributions,
      #             align = "vh",
      #             labels="AUTO",
      #             ncol = 2,
      #             nrow = 1,
      #             axis = "b")
      # 
      # cowplot::save_plot(paste0("plots/paper_figure_model_performances_", train.cline,"_", matrix.type, ".png"),
      #                    plot = model.perf.grid,
      #                    ncol = 2,
      #                    nrow = 3, 
      #                    dpi = 300,
      #                    base_height = 4.71)
      
      
      # custom.model.results.plot <- 
      #   visualize.custom.model.results(model.results = model.results, 
      #                                  train.cline = train.cline, 
      #                                  matrix.type = matrix.type, 
      #                                  model.target = model.target, 
      #                                  plot.base.size = plot.base.size)
      
      # custom.model.subspace.gsea <- feature.sub.space.gsea(model.results = model.results,
      #                                                      train.cline = train.cline,
      #                                                      matrix.type = matrix.type,
      #                                                      feature.subspace = feature.subspace, 
      #                                                      model.target = model.target)
      
      
      
      #####
      #  Model performances
      #####
      indiv.full.model.perf <- model.performance$obs.v.pred.plot+
        #ggtitle("Observed vs. Predicted Traveling Ratio\nof independent test data set (K562)")+
        xlab("Observed Traveling Ratio")+
        ylab("Predicted Traveling Ratio")+
        theme( legend.position = "right", plot.title = element_text(size=14))+
        theme_pubr( base_size = 14)
      
      sync.full.model.perf <- sync.model.performance$obs.v.pred.plot+
        #ggtitle("Observed vs. Predicted Traveling Ratio\nof independent cross cell type data set (HepG2)")+
        xlab("Observed Traveling Ratio")+
        ylab("Predicted Traveling Ratio")+
        theme( legend.position = "right", plot.title = element_text(size=14))+
        theme_pubr( base_size = 14)
      
      cell.type.specific.sample <- cell.type.specific.sample.model.performances$cline.specific.traveling.ratios+
        #ggtitle("Cell type specific samples")+
        xlab("Traveling Ratio (K562 cell line)")+
        ylab("Traveling Ratio (HepG2 cell line)")+
        theme( plot.title = element_text(size=14))+
        theme_pubr( base_size = 14)
      
      cell.type.specific.traveling.ratio.prediction.performances <- cell.type.specific.sample.model.performances$cline.specific.traveling.ratio.prediction.performances+
        #ggtitle("Cell type specific sample prediction performances")+
        xlab("log2(Obs_K562 / Obs_HepG2)")+
        ylab("log2(Pred_K562 / Pred_HepG2)")+
        theme( plot.title = element_text(size=14))+
        theme_pubr( base_size = 14)
      
      cell.type.specific.sample.prediction.performances <- cell.type.specific.sample.model.performances$cell.line.specific.gene.prediction.performances+
        #ggtitle("Cell type specific sample prediction performances")+
        xlab("Observed traveling ratio")+
        ylab("Predicted traveling ratio")+
        theme( plot.title = element_text(size=14))+
        theme_pubr( base_size = 14)
      
      model.perf.grid <-  
        plot_grid(indiv.full.model.perf,
                  sync.full.model.perf,
                  gene.venn,
                  cell.type.specific.sample.prediction.performances,
                  cell.type.specific.sample,
                  cell.type.specific.traveling.ratio.prediction.performances,
                  align = "vh",
                  labels="AUTO",
                  ncol = 2,
                  nrow = 3,
                  axis = "b")
      
      cowplot::save_plot(paste0("plots/paper_figure_model_performances_", train.cline,"_", matrix.type, ".png"),
                         plot = model.perf.grid,
                         ncol = 2,
                         nrow = 3, 
                         dpi = 300,
                         base_height = 4.71)
      
      
      #####
      # Model array
      #####
      # go.term.models <- go.term.model.results.plot$model.performance.plot+
      #   theme_pubr(base_size = 12)+
      #   theme( legend.position = "right")+
      #   xlab("")
      # cowplot::save_plot(paste0("plots/paper_figure_model_array_", train.cline,"_", matrix.type,".png"),
      #                    plot = go.term.models,
      #                    ncol = 2,
      #                    nrow = 2, 
      #                    dpi = 300,
      #                    base_height = 4.71)

      #####
      # Final model
      #####
      #####
      # Traveling Ratio Deconvolution
      #####
      forces.plot <- shap.force.plots$force.plot+
        xlab("Transcripts")+
        ylab("Feature Contribution")+
        ggtitle("")+
        scale_fill_jco()+
        theme(plot.margin=margin(r = 1, unit="cm"))
      
      cowplot::save_plot(paste0("plots/conference_figure_shap_focces_", train.cline,"_", matrix.type,".png"),
                         plot = forces.plot,
                         ncol = 1,
                         nrow = 1, 
                         dpi = 600)
      
      weighted.aggregate.class.contributions <- dynamic.features.plots$aggregate.weighted.class.contributions.plot+ggtitle("")+
        ylab("Overall aggregate absolute class contributions")+
        theme_pubr( base_size = 14)
      relative.aggregate.class.contributions <- dynamic.features.plots$relative.aggregate.class.contributions.plot+ggtitle("")+
        ylab("Relative aggregate absolute class contributions")+
        theme_pubr( base_size = 14)
      
      prior.knowledge.contributions <- prior.knowledge.contributions+
        theme( legend.position = "right", plot.title = element_text(size=10))+
        ylab("Aggregate absolute contributions")+
        xlab("Process")
      
      go.term.models <- 
      go.term.model.results.plot$model.performance.plot+
        theme( legend.position = "right", plot.title = element_text(size=10))+
        xlab("Process")
      
      
      
      # grid1 <- 
      #   plot_grid(force.cluster.plot, 
      #             deconvoluted.target,
      #             labels = c("B", "C"),
      #             align = "h",
      #             ncol = 2,
      #             nrow = 1,
      #             axis = "b")
      
      grid2 <- 
        plot_grid(weighted.aggregate.class.contributions,
                  relative.aggregate.class.contributions,
                  labels = c("D", "E" ),
                  align = "h",
                  ncol = 2,
                  nrow = 1,
                  axis = "b")
      
      grid3 <- 
        plot_grid(forces.plot, 
                  #grid1,
                  prior.knowledge.contributions,
                  go.term.models,
                  grid2,
                  labels = c("A", "B", "C"),
                  align = "h",
                  ncol = 1,
                  nrow = 4,
                  axis = "b")
      
      cowplot::save_plot(paste0("plots/paper_figure_feature_contrib_", train.cline,"_", matrix.type,".png"),
                         plot = grid3,
                         ncol = 1,
                         nrow = 5, 
                         dpi = 600, 
                         base_width = 14, 
                         base_height = 3.71)
      
      ### Factor Interpretations
      # all.cluster.factor.contribs.plots <- cluster.effectors$all.cluster.factor.contribs.plots
      # all.cluster.factor.contribs.plots.grid <-
      #   plot_grid(all.cluster.factor.contribs.plots[[1]],
      #             all.cluster.factor.contribs.plots[[2]],
      #             labels = "A",
      #             ncol = 2, 
      #             nrow = 1)

      force.cluster.plot <- shap.force.plots$force.plot.cluster.plot+
        ggtitle("")+
        xlab("Transcripts")+
        ylab("Feature Contribution")+
        scale_fill_jco()
      
      deconvoluted.target <- deconvoluted.target.plots.minimal.model$deconv.target.plot+
        ggtitle("")+
        xlab("Traveling Ratio")+
        ylab("Count")+
        theme(plot.margin=margin(t = 0, unit="cm"))

      library(png)
      novel.modulators <- readPNG("./plots/ModulatorsOfPausing.png")
      novel.modulators <- 
      ggplot() +
        annotation_raster(novel.modulators, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

      # optimal.factors.grid <- 
      #   plot_grid(optimal.cluster.factor.directional.contrib$optimal.factor.contrib.plot,
      #             novel.modulators,
      #             #optimal.factor.target.deconvolution,
      #             labels = c("B", "C"),
      #             align = "h",
      #             ncol = 2,
      #             nrow = 1,
      #             axis = "b", 
      #             rel_widths = c(0.7, 1.3))
      # 
      # final.model.results.grid <-
      #   plot_grid(all.cluster.factor.contribs.plots.grid,
      #             optimal.factors.grid,
      #             ncol = 1,
      #             nrow = 2,
      #             align = "h",
      #             axis = "b")
      
      final.model.results.grid <- 
        plot_grid(optimal.cluster.factor.directional.contrib$optimal.factor.contrib.plot, 
                  force.cluster.plot, 
                  deconvoluted.target,
                  novel.modulators,
                  labels = c("A", "B", "C", "D"),
                  align = "h",
                  ncol = 2,
                  nrow = 2,
                  axis = "b")
      
      cowplot::save_plot(paste0("plots/paper_figure_factor_selection_", train.cline,"_", matrix.type,".png"),
                         plot = final.model.results.grid,
                         ncol = 2,
                         nrow = 2, 
                         dpi = 600, 
                         base_height = 4)

      return(list(cluster.effectors, 
             cluster1.optimal.factors.plot = cluster1.optimal.factors.plot, 
             cluster2.optimal.factors.plot = cluster2.optimal.factors.plot,
             contingency.test.results = contingency.test.results,
             minimal.model.results = minimal.model.results))
      
      # # Prepare paper figures
      # prepare.paper.figures(model.performance = model.performance, 
      #                       sync.model.performance = sync.model.performance,
      #                       cell.type.specific.sample.model.performances = cell.type.specific.sample.model.performances, 
      #                       shap.summary.plot = shap.summary.plot, static.features.plots, 
      #                       dynamic.features.plots = dynamic.features.plots, 
      #                       individual.directional.factor.contrib.plots = individual.directional.factor.contrib.plots, 
      #                       shap.force.plots = shap.force.plots, 
      #                       cluster.target.distribution.plots = cluster.target.distribution.plots, 
      #                       ddeconvoluted.target.plots = econvoluted.target.plots, 
      #                       optimal.factors.plot = optimal.factors.plot, 
      #                       optimal.factor.gsea = optimal.factor.gsea, 
      #                       target.correlations = target.correlations, 
      #                       reverse.model.feature.correlations = reverse.model.feature.correlations, 
      #                       knockdown.correlations.plot = knockdown.correlations.plot, 
      #                       go.term.model.results.plot = go.term.model.results.plot, 
      #                       custom.model.results.plot = custom.model.results.plot)
    })
  })
  return()
}
evaluate.feature.effects <- function(model.results, feature.subspace, model.target, plot.base.size =  16){
  
  # #X TEMP vars to test function
  model.target="target_traveling.ratio"
  matrix.type <- "individual.model.matrices"
  feature.subspace <- "All"
  model.target.type <- gsub("target_", "", model.target)
  # train.cline <- "K562"
  # test.cline <- "HepG2"
  plot.base.size =  16

  
  K562.gene.count <- length(model.matrices[["K562"]][[matrix.type]][[feature.subspace]][[model.target]])
  HepG2.gene.count <- length(model.matrices[["HepG2"]][[matrix.type]][[feature.subspace]][[model.target]])
  cross <- length(intersect(rownames(model.matrices[["K562"]][[matrix.type]][[feature.subspace]]),
                            rownames(model.matrices[["HepG2"]][[matrix.type]][[feature.subspace]])))
  grid.newpage()
  gene.venn <- 
    draw.pairwise.venn(area1 = K562.gene.count,
                       area2 = HepG2.gene.count,
                       cross.area = cross,
                       category = c("K562", "HepG2"),
                       scaled = F, 
                       cat.fontface = c("bold", "bold"),
                       fill = c(alpha("#440154ff",1), alpha('#21908dff',1)))

  
  evaluation.results <- 
    lapply(CARGS$cell.line, function(train.cline){
      lapply(c("individual.model.matrices", "synchronised.model.matrices"), function(matrix.type){
        
        test.cline <- if(train.cline=="K562") "HepG2" else "K562"

        # Load relevant model results
        select.model.results <- model.results$models[[train.cline]][[matrix.type]][[feature.subspace]][[model.target]]
        
        # Extract relevant variables from model results
        X <- select.model.results$X
        Y <- select.model.results$Y
        X.test = select.model.results$X.test
        Y.test <- select.model.results$Y.test
        train.preds <- select.model.results$train.preds
        test.preds <- select.model.results$test.preds
        samples <- rownames(X)
        
        ## Individual hap scores
        # The sum of each row’s SHAP values (plus the BIAS column, which is like an intercept) 
        # is the predicted model output,i.e., the explanation’s attribution values sum up to the model output
        # but not true for CV models
        shap.values <- select.model.results$shap.values
        shap.scores <- as.data.frame(shap.values$shap_score)
        rownames(shap.scores) <- rownames(X)
        shap.bias <- shap.values$BIAS0
        # Mean |SHAP| scores for each feature over all samples
        mean.shap.scores <- shap.values$mean_shap_score
        
        # Plot obs vs pred prediction performances
        print(paste0("Model Performances ", train.cline, " ", matrix.type))
        model.performance <- 
          visualize.model.performance(select.model.results, 
                                      feature.subspace, 
                                      model.target.type, 
                                      plot.base.size = plot.base.size,
                                      train.cline = train.cline)
        
        # Plot obs target differences between cell lines against differences of predictions between cell lines
        print(paste0("Cross Cell Type Model Performances ", train.cline, " ", matrix.type))
        sync.select.model.results <- model.results$models[[train.cline]][["synchronised.model.matrices"]][[feature.subspace]][[model.target]]
        sync.reverse.model.results <- model.results$models[[test.cline]][["synchronised.model.matrices"]][[feature.subspace]][[model.target]]
        cell.type.specific.sample.model.performances <- 
          retrieve.cell.type.specific.sample.prediction.performances(sync.select.model.results, 
                                                                     sync.reverse.model.results,
                                                                     train.cline, 
                                                                     test.cline,
                                                                     matrix.type, 
                                                                     feature.subspace, 
                                                                     model.target.type, 
                                                                     plot.base.size)
        ## Model performance on cross cell type data
        print(paste0("Cross application performance ", train.cline, " ", matrix.type))
        sync.model.performance <- 
          visualize.model.performance(sync.select.model.results, 
                                      feature.subspace,
                                      model.target, 
                                      plot.base.size=plot.base.size,train.cline)
        
        ## Individual Feature distribution of top 25 features
        print(paste0("SHAP summary plot ", train.cline, " ", matrix.type))
        shap.summary.plot <- 
          visualize.shap.summary(shap.score.matrix = shap.scores, 
                                 X= X, 
                                 train.cline = train.cline, 
                                 matrix.type = matrix.type, 
                                 feature.subspace = feature.subspace, 
                                 model.target = model.target,
                                 plot.base.size = 12)
        
        ## DNA sequence feature distributions
        print(paste0("DNA sequence feature plots ", train.cline, " ", matrix.type))
        static.features.plots <- 
          visualize.static.features(shap.scores, 
                                    X= X, 
                                    train.cline = train.cline, 
                                    matrix.type = matrix.type, 
                                    feature.subspace = feature.subspace, 
                                    model.target = model.target,
                                    plot.base.size = plot.base.size)
        
        print(paste0("Binding feature plots ", train.cline, " ", matrix.type))
        dynamic.features.plots <- 
          visualize.dynamic.features(shap.scores = shap.scores, 
                                     X = X, 
                                     train.cline = train.cline, 
                                     matrix.type = matrix.type, 
                                     feature.subspace = feature.subspace, 
                                     model.target = model.target,
                                     plot.base.size = plot.base.size, 
                                     mean.shap.scores,
                                     append = "")
        
        dynamic.features.plots$aggregate.non.binding.feature.contrib.plot <- dynamic.features.plots$aggregate.non.binding.feature.contrib.plot+
          ggtitle("Contributions of gene annotation and composition features")
        cowplot::save_plot(paste0("plots/supplementary_figure_static_feature_contribs_", train.cline,"_", matrix.type, ".png"),
                           plot = dynamic.features.plots$aggregate.non.binding.feature.contrib.plot,
                           ncol = 1,
                           nrow = 1, 
                           dpi = 600)
        
        
        # Check enrichment of DNA intron binders
        # Retrieve top itron DNA binding proteins
        dna.intron.binding.contribs <- dynamic.features.plots$class.contributions$DNA_introns
        dna.intron.binding.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*",
                                           "", names(dna.intron.binding.contribs))
        names(dna.intron.binding.contribs) <- dna.intron.binding.factors
        # Retrieve functional sets
        sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
        sub.factor.feature.spaces <-
          build.sequence.specific.biologically.functional.feature.sub.spaces(stratify = sub.factor.feature.spaces)
        sub.factor.feature.spaces <- sub.factor.feature.spaces[!grepl(":ss|All|sdom", names(sub.factor.feature.spaces))]
        sub.factor.feature.spaces <- sub.factor.feature.spaces[c("Chromatin", "Initiation", "Elongation", "7SK.Binding", "Elongation+7SK", "Termination","Splicing", "Processing")]
        
        dna.intron.binders.functional.contribs <-
          lapply(sub.factor.feature.spaces, function(sub.space){
            factors <- names(dna.intron.binding.contribs)[names(dna.intron.binding.contribs) %in% sub.space[[train.cline]]]
            sum((dna.intron.binding.contribs[factors]))
          })
        dna.intron.binders.functional.contribs <- unlist(dna.intron.binders.functional.contribs)
        dna.intron.binders.functional.contribs <- data.frame(Process = names(dna.intron.binders.functional.contribs), Contribution = dna.intron.binders.functional.contribs)
        dna.intron.binders.functional.contribs <- na.omit(dna.intron.binders.functional.contribs)
        dna.intron.binders.functional.contribs <- dna.intron.binders.functional.contribs[order(dna.intron.binders.functional.contribs$Contribution, decreasing = F),]
        dna.intron.contrib.plot <-
          ggbarplot(dna.intron.binders.functional.contribs,
                    x = "Process",
                    y = "Contribution",
                    fill = rgb(0.2,0.4,0.6,0.6),
                    rotate = T,
                    #lab.pos = "in",
                    lab.col = "black",
                    lab.font = "bold",
                    lab.size = 4,
                    lab.hjust = -0.5,
                    lab.vjust = 0.5,
                    position = position_dodge(0.9))+
          xlab("Process")+
          ylab("Aggregate contribution")+
          theme_pubr(base_size = 12, legend = "right")+
          ggtitle("Process contributions on genomic intron regions")
        cowplot::save_plot(paste0("plots/supplementary_figure_dna_intronbindings_", train.cline,"_", matrix.type, ".png"),
                           plot = dna.intron.contrib.plot,
                           ncol = 1,
                           nrow = 1, 
                           dpi = 600)
        
        # DNA intron elongation factor enrichment results
        dbp.factors <- 
          lapply(CHIPseq, function(peaks){
            unique(peaks$hgnc_symbol)
          }) %>% unlist %>% unique 
        all.dna.binding.factors <- lapply(sub.factor.feature.spaces, function(process){process[train.cline]})
        all.dna.binding.factors <- unique(unlist(all.dna.binding.factors))
        elongation.7sk.factors <- unique(sub.factor.feature.spaces[["Elongation+7SK"]][[train.cline]])
        dna.intron.binding.factors <-  
          lapply(sub.factor.feature.spaces, function(sub.space){
            factors <- names(dna.intron.binding.contribs)[names(dna.intron.binding.contribs) %in% sub.space[[train.cline]]]
          })
        dna.intron.binding.elongation.7sk.factors <- dna.intron.binding.factors[["Elongation+7SK"]]
        dna.intron.binding.factors <- unlist(dna.intron.binding.factors)
        non.elongation.7sk.factors <- all.dna.binding.factors[!(all.dna.binding.factors %in% elongation.7sk.factors)]
        
        categories <- list(c("Elongation+7SK", "other"), c("binding", "non-binding"))
        counts <- c(length(dna.intron.binding.elongation.7sk.factors),
                    sum(non.elongation.7sk.factors %in% dna.intron.binding.factors),
                    sum(!elongation.7sk.factors %in% dna.intron.binding.factors),
                    sum(!non.elongation.7sk.factors %in% dna.intron.binding.factors))
        dna.intron.elongation.enrichment.results <- perform.fisher.test(categories = categories, 
                                                                       counts = counts)
        

        ## RNA intron bindings
        rna.intron.binding.contribs <- dynamic.features.plots$class.contributions$RNA_introns
        rna.intron.binding.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*",
                                           "", names(rna.intron.binding.contribs))
        names(rna.intron.binding.contribs) <- rna.intron.binding.factors
        # Retrieve functional sets
        sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
        sub.factor.feature.spaces <- 
          build.sequence.specific.biologically.functional.feature.sub.spaces(stratify = sub.factor.feature.spaces)
        sub.factor.feature.spaces <- sub.factor.feature.spaces[!grepl(":ss|All|sdom", names(sub.factor.feature.spaces))]
        sub.factor.feature.spaces <- sub.factor.feature.spaces[c("Chromatin", "Initiation", "Elongation", "7SK.Binding", "Elongation+7SK", "Termination","Splicing", "Processing")]
        
        rna.intron.binders.functional.contribs <- 
          lapply(sub.factor.feature.spaces, function(sub.space){
            factors <- names(rna.intron.binding.contribs)[names(rna.intron.binding.contribs) %in% sub.space[[train.cline]]]
            sum((rna.intron.binding.contribs[factors]))
          })
        rna.intron.binders.functional.contribs <- unlist(rna.intron.binders.functional.contribs)
        rna.intron.binders.functional.contribs <- data.frame(Process = names(rna.intron.binders.functional.contribs), Contribution = rna.intron.binders.functional.contribs)
        rna.intron.binders.functional.contribs <- na.omit(rna.intron.binders.functional.contribs)
        rna.intron.binders.functional.contribs <- rna.intron.binders.functional.contribs[order(rna.intron.binders.functional.contribs$Contribution, decreasing = F),]
        rna.intron.contrib.plot <- 
          ggbarplot(rna.intron.binders.functional.contribs, 
                    x = "Process", 
                    y = "Contribution",
                    fill = rgb(0.2,0.4,0.6,0.6), 
                    rotate = T,
                    #lab.pos = "in",
                    lab.col = "black",
                    lab.font = "bold",
                    lab.size = 4,
                    lab.hjust = -0.5,
                    lab.vjust = 0.5,
                    position = position_dodge(0.9))+
          xlab("Process")+
          ylab("Contribution")+
          theme_pubr(base_size = 11, legend = "right")+
          scale_fill_material("blue-grey")+
          ggtitle("Process contributions on transcriptomic intron regions")
        cowplot::save_plot(paste0("plots/supplementary_figure_rna_intronbindings_", train.cline,"_", matrix.type, ".png"),
                           plot = rna.intron.contrib.plot,
                           ncol = 1,
                           nrow = 1, 
                           dpi = 600)
        
        # RNA intron splicing factor enrichment results 
        # Retrieve all RNA binding factors
        rbp.factors <- 
          lapply(eCLIPseq, function(peaks){
            unique(peaks$hgnc_symbol)
          }) %>% unlist %>% unique 
        all.rna.binding.factors <- lapply(sub.factor.feature.spaces, function(process){process[train.cline]})
        all.rna.binding.factors <- unique(unlist(all.rna.binding.factors))
        all.rna.binding.factors <- all.rna.binding.factors[all.rna.binding.factors %in% unique(unlist(rbp.factors))]
        splicing.factors <- unique(sub.factor.feature.spaces[["Splicing"]][[train.cline]])
        rna.intron.binding.factors <-  
          lapply(sub.factor.feature.spaces, function(sub.space){
            factors <- names(rna.intron.binding.contribs)[names(rna.intron.binding.contribs) %in% sub.space[[train.cline]]]
          })
        rna.intron.binding.splicing.factors <- rna.intron.binding.factors[["Splicing"]]
        rna.intron.binding.factors <- unlist(rna.intron.binding.factors)
        non.splicing <- all.rna.binding.factors[!(all.rna.binding.factors %in% splicing.factors)]
        
        categories <- list(c("splicing", "non-splicing"), c("binding", "non-binding"))
        counts <- c(length(rna.intron.binding.splicing.factors),
                    sum(!splicing.factors %in% rna.intron.binding.factors),
                    sum(non.splicing %in% rna.intron.binding.factors),
                    sum(!non.splicing %in% rna.intron.binding.factors))
        rna.inton.splicing.enrichment.results <- perform.fisher.test(categories = categories, 
                                                                     counts = counts)

        # Shap force plots
        print(paste0("Force plots ", train.cline, " ", matrix.type))
        shap.force.plots <-
          visualize.shap.forces(shap.scores = shap.scores, 
                                mean.shap.scores = mean.shap.scores,
                                X = X,
                                train.cline = train.cline,
                                matrix.type = matrix.type,
                                feature.subspace = feature.subspace,
                                model.target = model.target,
                                plot.base.size = plot.base.size)
        
        ## Featre class and prior knowledge contributions
        prior.knowledge.contributions <-
          calculate.prior.knowledge.contributions(shap.scores = shap.scores,
                                                  mean.shap.scores = mean.shap.scores,
                                                  X = X,
                                                  train.cline = train.cline,
                                                  matrix.type = matrix.type,
                                                  feature.subspace = feature.subspace,
                                                  model.target = model.target,
                                                  plot.base.size = plot.base.size)

        # Reieve most influential factors
        optimal.factors.plot <-
          retrieve.optimal.number.of.factors(cum.factor.contrib = dynamic.features.plots$cum.factor.contrib,
                                             train.cline = train.cline,
                                             matrix.type = matrix.type,
                                             feature.subspace = feature.subspace,
                                             model.target = model.target,
                                             plot.base.size = plot.base.size)
        all.optimal.factors <- optimal.factors.plot$optimal.factors

        # Build minimal model with top influential factors
        minimal.model.results <- train.minimal.model(all.optimal.factors, train.cline, matrix.type, feature.subspace)
        minimal.model.results$minimal.model.perf <- 
          minimal.model.results$minimal.model.perf+
          xlab("Observed Pausing Index")+
          ylab("Predicted Pausing Index")
        cowplot::save_plot(paste0("plots/supplementary_minimal_model_", train.cline,"_", matrix.type, ".png"),
                           plot = minimal.model.results$minimal.model.perf,
                           ncol = 1,
                           nrow = 1, 
                           dpi = 600)
        
        print(paste0("Optimal factor directional contribution ", train.cline, " ", matrix.type))
        optimal.cluster.factor.directional.contrib <-
          visualize.optimal.factor.contributions( force.plot_data = shap.force.plots$force.plot_data,
                                                  shap.scores = shap.scores,
                                                  optimal.factors = all.optimal.factors,
                                                  cluster.membership = all.optimal.factors ,
                                                  X = X,
                                                  train.cline, 
                                                  matrix.type, 
                                                  feature.subspace,
                                                  model.target)

        
        print(paste0("Prior knowledge model array ", train.cline, " ", matrix.type))
        go.term.model.results.plot <- 
          visualize.go.term.model.results(model.results = model.results, 
                                          train.cline = train.cline, 
                                          matrix.type = matrix.type, 
                                          model.target = model.target, 
                                          plot.base.size = plot.base.size, 
                                          random.model.results = random.model.results)
        
        # Figure 2
        indiv.full.model.perf <- model.performance$obs.v.pred.plot+
          #ggtitle("Observed vs. Predicted Traveling Ratio\nof independent test data set (K562)")+
          xlab("Observed Pausing Index")+
          ylab("Predicted Pausing Index")+
          theme( legend.position = "right", plot.title = element_text(size=14))+
          theme_pubr( base_size = 14)
        
        sync.full.model.perf <- sync.model.performance$obs.v.pred.plot+
          #ggtitle("Observed vs. Predicted Traveling Ratio\nof independent cross cell type data set (HepG2)")+
          xlab("Observed Pausing Index")+
          ylab("Predicted Pausing Index")+
          theme( legend.position = "right", plot.title = element_text(size=14))+
          theme_pubr( base_size = 14)
        
        cell.type.specific.sample <- cell.type.specific.sample.model.performances$cline.specific.traveling.ratios+
          #ggtitle("Cell type specific samples")+
          xlab(paste0("Pausing Index (", train.cline, " cell line)"))+
          ylab(paste0("Pausing Index (", test.cline, " cell line)"))+
          theme( plot.title = element_text(size=14))+
          theme_pubr( base_size = 14)
        
        cell.type.specific.traveling.ratio.prediction.performances <- cell.type.specific.sample.model.performances$cline.specific.traveling.ratio.prediction.performances+
          #ggtitle("Cell type specific sample prediction performances")+
          xlab("Obs_K562_PI - Obs_HepG2_PI")+
          ylab("Pred_K562_PI - Pred_HepG2_PI")+
          theme( plot.title = element_text(size=14))+
          theme_pubr( base_size = 14)
        
        cell.type.specific.sample.prediction.performances <- cell.type.specific.sample.model.performances$cell.line.specific.gene.prediction.performances+
          #ggtitle("Cell type specific sample prediction performances")+
          xlab("Observed Pausing Index")+
          ylab("Predicted Pausing Index")+
          theme( plot.title = element_text(size=14))+
          theme_pubr( base_size = 14)
        
        model.perf.grid <-  
          plot_grid(indiv.full.model.perf,
                    sync.full.model.perf,
                    gene.venn,
                    cell.type.specific.sample.prediction.performances,
                    cell.type.specific.sample,
                    cell.type.specific.traveling.ratio.prediction.performances,
                    align = "vh",
                    labels="AUTO",
                    ncol = 2,
                    nrow = 3,
                    axis = "b")
        
        cowplot::save_plot(paste0("plots/paper_figure_model_performances_", train.cline,"_", matrix.type, ".png"),
                           plot = model.perf.grid,
                           ncol = 2,
                           nrow = 3, 
                           dpi = 300,
                           base_height = 4.71)
        
        #Figure 3
        forces.plot <- shap.force.plots$force.plot+
          xlab("Transcripts")+
          ylab("Feature Contribution")+
          ggtitle("")+
          scale_fill_jco()+
          theme(plot.margin=margin(r = 1, unit="cm"))
        
        cowplot::save_plot(paste0("plots/conference_figure_shap_focces_", train.cline,"_", matrix.type,".png"),
                           plot = forces.plot,
                           ncol = 1,
                           nrow = 1, 
                           dpi = 600)
        
        weighted.aggregate.class.contributions <- dynamic.features.plots$aggregate.weighted.class.contributions.plot+ggtitle("")+
          ylab("Overall aggregate absolute class contributions")+
          theme_pubr( base_size = 14)
        relative.aggregate.class.contributions <- dynamic.features.plots$relative.aggregate.class.contributions.plot+ggtitle("")+
          ylab("Relative aggregate absolute class contributions")+
          theme_pubr( base_size = 14)
        
        prior.knowledge.contributions <- prior.knowledge.contributions+
          theme( legend.position = "right", plot.title = element_text(size=10))+
          ylab("Aggregate absolute contributions")+
          xlab("Process")
        
        go.term.models <- 
          go.term.model.results.plot$model.performance.plot+
          theme( legend.position = "right", plot.title = element_text(size=10))+
          xlab("Process")
        
        # Figure 4
        grid1 <- 
          plot_grid(#weighted.aggregate.class.contributions,
                    relative.aggregate.class.contributions,
                    labels = c("D", "E" ),
                    align = "h",
                    ncol = 1,
                    nrow = 1,
                    axis = "b")
        
        grid2 <- 
          plot_grid(forces.plot, 
                    #grid1,
                    prior.knowledge.contributions,
                    go.term.models,
                    grid1,
                    labels = c("A", "B", "C"),
                    align = "h",
                    ncol = 1,
                    nrow = 4,
                    axis = "b")
        
        cowplot::save_plot(paste0("plots/paper_figure_feature_contrib_", train.cline,"_", matrix.type,".png"),
                           plot = grid2,
                           ncol = 1,
                           nrow = 5, 
                           dpi = 600, 
                           base_width = 14, 
                           base_height = 3.71)
        
        ### Figure 5
        novel.modulators <- readPNG("./plots/ModulatorsOfPausing.png")
        novel.modulators <- 
          ggplot() +
          annotation_raster(novel.modulators, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
        
        final.model.results.grid <- 
          plot_grid(optimal.cluster.factor.directional.contrib$optimal.factor.contrib.plot, 
                    novel.modulators,
                    labels = c("A", "B", "C", "D"),
                    align = "h",
                    ncol = 2,
                    nrow = 1,
                    axis = "b")
        
        cowplot::save_plot(paste0("plots/paper_figure_factor_selection_", train.cline,"_", matrix.type,".png"),
                           plot = final.model.results.grid,
                           ncol = 2,
                           nrow = 1, 
                           dpi = 600, 
                           base_height = 4)
        
        
        ### Supplementary figures
        sequence.feature.distributions <- visualize.features()
        cowplot::save_plot(paste0("plots/suppl_figure_sequence_feature_distribution_", train.cline, "_", matrix.type, ".png"),
                           plot = sequence.feature.distributions[[train.cline]],
                           ncol = 2,
                           nrow = 8, 
                           dpi = 600,  
                           base_width = 14, 
                           base_height = 5)
        
        ## SHAP contributions
        shap.summary.plot <- shap.summary.plot+
          theme( axis.text.y = element_text(size=12))
        
        cowplot::save_plot(paste0("plots/suppl_figure_shap_contributions_", train.cline, "_", matrix.type, ".png"),
                           plot = shap.summary.plot,
                           ncol = 2,
                           nrow = 2, 
                           dpi = 600)
        
        
        cowplot::save_plot(paste0("plots/suppl_figure_shap_static_feature_distributions_", train.cline, "_", matrix.type, ".png"),
                           plot = static.features.plots,
                           ncol = 2,
                           nrow = 8, 
                           dpi = 600,
                           base_width = 14, 
                           base_height = 5)
        
        
        cowplot::save_plot(paste0("plots/suppl_figure_shap_non_binding_feature_contributions_", train.cline, "_", matrix.type, ".png"),
                           plot = dynamic.features.plots$aggregate.non.binding.feature.contrib.plot,
                           ncol = 1,
                           nrow = 1, 
                           dpi = 600)
        
        
        #visualize.7SK.binders()
        
        
        return(list(dna.intron.elongation.enrichment.results = dna.intron.elongation.enrichment.results,
                    rna.inton.splicing.enrichment.results = rna.inton.splicing.enrichment.results,
                    all.optimal.factors = all.optimal.factors))
      })
    })
  names(evaluation.results) <- CARGS$cell.line
  
  ## Inverse correlation of transcript traveling ratios with transcript expressions 
  inv.corr.tr.exp <- visualize.traveling.ratio.transcript.expression.distribution(plot.base.size)
  inv.corr.tr.exp.grid <- 
    plot_grid(inv.corr.tr.exp[[1]], 
              inv.corr.tr.exp[[2]],
              ncol = 1, 
              nrow = 2,
              align = "b",
              labels = "AUTO")
  cowplot::save_plot(paste0("plots/suppl_figure_inverse_corr_tr_exp.png"),
                     plot = inv.corr.tr.exp.grid,
                     ncol = 1,
                     nrow = 2, 
                     dpi = 600)
  
  ## Traveling ratio distribution
  traveling.ratio.distributions <- visualize.traveling.ratio(plot.base.size)
  traveling.ratio.distributions <- 
    plot_grid(traveling.ratio.distributions[[1]], 
              traveling.ratio.distributions[[2]],
              ncol = 1, 
              nrow = 2,
              align = "b",
              labels = "AUTO")
  cowplot::save_plot(paste0("plots/suppl_figure_traveling_ratios.png"),
                     plot = traveling.ratio.distributions,
                     ncol = 1,
                     nrow = 2, 
                     dpi = 600)
  return()
}
###
# This function visualizes the performance of model predictions
###
visualize.model.performance <- function(select.model.results, 
                                        feature.subspace, 
                                        model.target, 
                                        plot.base.size, 
                                        train.cline){
  
  # Extract relevant variables from model results
  X <- select.model.results$X
  Y <- select.model.results$Y
  X.test = select.model.results$X.test
  Y.test <- select.model.results$Y.test
  train.preds <- select.model.results$train.preds
  test.preds <- select.model.results$test.preds
  shap.values <- select.model.results$shap.values
  samples <- select.model.results$samples
  
  # Identify test cell line
  test.cline <- if(train.cline=="K562") "HepG2" else "K562"
  
  ## Evaluate model performance
  # Calculate root mean squared error of test data
  test.residuals <- Y.test - test.preds
  test.rmse <- sqrt(mean(test.residuals^2))
  ## Calculate the r-squared of the test data 
  # Calculate the mean of the test data observations
  mean.y.test <- mean(Y.test)
  # Calculate the total sum of squares
  tss =  sum((Y.test - mean.y.test)^2)
  # Calculate residual sum of squares
  rss =  sum(test.residuals^2)
  # Calculate R-squared
  rsq.test  =  1 - (rss/tss)
  
  # Observed vs predicted values with residuals
  obs.v.pred = data.frame(obs = Y.test, pred = test.preds, resid = abs(test.residuals))
  # Plot observed vs predicted targets without residual errors
  obs.v.pred.plot <- 
    ggscatter(obs.v.pred, 
              x = "obs", 
              y = "pred",
              color = "resid",
              add ="reg.line",
              conf.int = TRUE,
              add.params = list(color = "red",
                                fill = "lightgray",
                                size = 0.1),
              #palette = pal_uchicago("dark")(9)[9],
              alpha = 1/10
    ) +
    theme_pubr(base_size = plot.base.size)+
    #ggtitle(paste0("Observed vs. Predicted target \n", "(", train.cline, " / ", model.target.type,"/", feature.subspace, ")"))+
    xlab("Observed test data target")+
    ylab("Predicted test data target")+
    ggpubr::stat_cor(cor.coef.name = c("rho"))+
    #stat_binhex(color = "white", binwidth = 0.25)+
    #scale_fill_gradient(low = "grey", high = "darkblue")+
    gradient_color(c("black", pal_uchicago("light")(9)[1]))
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS,"obs_vs_pred_target_",gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target),"_", train.cline, ".pdf"), 
                     plot = obs.v.pred.plot)
  
  # Plot distribution of residuals
  residuals <- data.frame(residuals = test.residuals)
  residual.plot <- 
    ggdensity(residuals, 
              x = "residuals", 
              alpha = 0.7, 
              fill = pal_uchicago("light")(9)[2],
              add = "mean") +
    theme_pubr(base_size = plot.base.size)+
    #ggtitle(paste0("Residual distribution \n", "(", train.cline, " / ", model.target.type,"/", feature.subspace, ")"))+
    xlab("Residual")
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS,"prediction_residuals_",gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target),"_", train.cline, ".pdf"), 
                     plot = residual.plot)
  
  return(list(obs.v.pred.plot = obs.v.pred.plot, 
              residual.dist.plot = residual.plot))
}
### 
# This function visualizes the model performances on cross cell type specific samples
###
retrieve.cell.type.specific.sample.prediction.performances <- function(sync.select.model.results, 
                                                                      sync.reverse.model.results, 
                                                                      train.cline, 
                                                                      test.cline,
                                                                      matrix.type, 
                                                                      feature.subspace, 
                                                                      model.target, 
                                                                      plot.base.size){
  
  # sync.select.model.results, 
  # sync.reverse.model.results,
  # samples,
  # train.cline, matrix.type, 
  # feature.subspace, model.target, 
  # plot.base.size = plot.base.size
  
  ### Evaluate cell type specific samples prediction performances
  ## Retrieve transcript quantifications 
  train.cline.target <- sync.select.model.results$Y
  test.cline.target <- sync.reverse.model.results$Y
  common.samples <- intersect(rownames(train.cline.target), rownames(test.cline.target))
  train.cline.target <- train.cline.target[match(common.samples, rownames(train.cline.target)),]
  test.cline.target <- test.cline.target[match(common.samples, rownames(test.cline.target)),]
  observed.target.differences <-  train.cline.target-test.cline.target
  
  #observed.target.differences.test <- (test.cline.target-(train.cline.target))/abs(train.cline.target)*100
  #observed.target.differences.train<- (train.cline.target-(test.cline.target))/abs(test.cline.target)*100

  
  # Retrieve cell type specific expressions via cutoff on the traveling ratio
  if(model.target.type == "traveling.ratio"){
    train.cline.specific <- which(observed.target.differences >= 1)
    test.cline.specific <- which(observed.target.differences <= -1)
  }else{
    train.cline.specific <- which(observed.target.differences >= log10(2))
    test.cline.specific <- which(observed.target.differences <= -log10(2))
  }

  ## Genes with cell type specific traveling ratio distributions
  target.profiles <- data.frame( sample = common.samples,
                                 train.cline.target = train.cline.target,
                                 test.cline.target = test.cline.target,
                                 specific = "none")
  target.profiles$specific[train.cline.specific] <- train.cline
  target.profiles$specific[test.cline.specific] <- test.cline  
  
  # target.profiles$specific[observed.target.differences.test >= 100] <- test.cline 
  # target.profiles$specific[observed.target.differences.test <= -100] <- train.cline 
  # 
  # target.profiles$specific[observed.target.differences.train >= 100] <- train.cline
  # target.profiles$specific[observed.target.differences.train <= -100] <- test.cline
  target.profiles$specific <- as.factor(target.profiles$specific)
  
  cline.specific.traveling.ratios <- 
    ggpubr::ggscatter(target.profiles, 
                      x = "train.cline.target", 
                      y = "test.cline.target",
                      color = "specific", 
                      palette = pal_uchicago("dark")(9)[c(4:6)],           # Color by groups "cyl"
                      #shape = "specific",                        # Change point shape by groups "cyl"
                      fullrange = TRUE,                         # Extending the regression line
                      rug = F,                                # Add marginal rug
                      alpha = 1/5,
                      show.legend.text = FALSE,
                      xlim = c(-5, 20),
                      ylim = c(-5, 20),
                      #xticks.by = 1,
                      #yticks.by = 1
    )+
    xlab(paste0(train.cline," ", "traveling ratio"))+
    ylab(paste0(test.cline," ",  "traveling ratio"))+
    scale_color_manual(name = "At least twice as high pausing index in ",
                       breaks = c("HepG2", "K562", "none"),
                       values = c("HepG2" = pal_uchicago("dark")(9)[c(4:6)][1],
                                  "K562" = pal_uchicago("dark")(9)[c(4:6)][2],
                                  "none" = pal_uchicago("dark")(9)[c(4:6)][3]),
                       guide = guide_legend(override.aes = list(alpha = 1, size = 6)))+
    ggpubr::stat_cor(aes(color = specific), cor.coef.name = c("rho"))+
    theme_pubr( base_size = plot.base.size)
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, train.cline, "_and_", test.cline, "_specific_traveling_ratio_distributions_", gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target),"_", train.cline, ".pdf"), 
                     plot = cline.specific.traveling.ratios)
  
  
  ## Prediction performances of a model of genes with cross cell type specific traveling ratio distributions
  train.pred.mat <- data.frame(sample = target.profiles$sample[match(common.samples,target.profiles$sample)], 
                               predictions = sync.select.model.results$test.preds[match(common.samples, names(sync.select.model.results$test.preds))], 
                               observed.cross.cline = sync.select.model.results$Y.test[match(common.samples, rownames(sync.select.model.results$Y.test))], 
                               cell.type = train.cline, 
                               specific = "both", 
                               row.names = NULL)
  test.pred.mat <- data.frame(sample = target.profiles$sample[match(common.samples,target.profiles$sample)], 
                              predictions = sync.reverse.model.results$test.preds[match(common.samples, names(sync.reverse.model.results$test.preds))], 
                              observed.cross.cline = sync.reverse.model.results$Y.test[match(common.samples, rownames(sync.reverse.model.results$Y.test))], 
                              cell.type = test.cline, 
                              specific = "both", 
                              row.names = NULL)
  
  all.sample.predictions <- data.frame(samples = common.samples,
                                       pred.diff = train.pred.mat$predictions-test.pred.mat$predictions,
                                       #target.diff =  observed.target.differences,
                                       target.diff = train.pred.mat$observed.cross.cline - test.pred.mat$observed.cross.cline,
                                       expressed = "both")
  
  all.sample.predictions$expressed[all.sample.predictions$sample %in% target.profiles$sample[target.profiles$specific==train.cline]] <- train.cline
  all.sample.predictions$expressed[all.sample.predictions$sample %in% target.profiles$sample[target.profiles$specific==test.cline]] <- test.cline
  all.sample.predictions <- all.sample.predictions[!all.sample.predictions$expressed=="both",]
  #all.sample.predictions$direction <- revalue(as.factor(all.sample.predictions$target.diff>0), c("TRUE" = "K562", "FALSE" = "HepG2"))
  
  cline.specific.traveling.ratio.prediction.performances <- 
    ggpubr::ggscatter(all.sample.predictions, 
                      x = "target.diff", 
                      y = "pred.diff",
                      color ="expressed",
                      #palette = pal_uchicago("dark")(9)[c(4:6)],           
                      fullrange = TRUE,                         
                      rug = F,
                      alpha = 1/5
    )+
    ggpubr::stat_cor(aes(color = expressed), cor.coef.name = c("rho"))+
    #ggpubr::stat_cor(cor.coef.name = c("rho"))+
    xlab("Observed traveling ratio difference")+
    ylab("Predicted traveling ratio difference")+
    scale_color_manual(name = "Genes with higher pausing indices in ",
                       breaks = c("HepG2", "K562"),
                       values = c("HepG2" = pal_uchicago("dark")(9)[c(4:6)][1],
                                  "K562" = pal_uchicago("dark")(9)[c(4:6)][2]),
                       guide = guide_legend(override.aes = list(alpha = 1, size = 6)))+
    theme_pubr( base_size = plot.base.size)
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, train.cline, "_and_", test.cline, "_specific_traveling_ratio_sample_predictions_", gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target),"_", train.cline, ".pdf"), 
                     plot = cline.specific.traveling.ratio.prediction.performances)
  
  # # Add genes that are exclusively present in one cell line but not the other
  # Train cline exclusively expressed genes
  additional.train.cline.genes <- sync.select.model.results$Y[!rownames(sync.select.model.results$Y) %in% rownames(sync.reverse.model.results$Y),]
  # Test cline exclusively expressed genes
  additional.test.cline.genes <-  sync.reverse.model.results$Y[!rownames(sync.reverse.model.results$Y) %in% rownames(sync.select.model.results$Y),]
  additional.train.cline.genes <- data.frame(sample = names(additional.train.cline.genes),
                                             observed.target = additional.train.cline.genes,
                                             cross.cline.preds = sync.reverse.model.results$test.preds[match(names(additional.train.cline.genes), names(sync.reverse.model.results$test.preds))],
                                             specific = train.cline)
  rownames(additional.train.cline.genes) <- NULL
  additional.test.cline.genes <- data.frame(sample = names(additional.test.cline.genes),
                                             observed.target = additional.test.cline.genes,
                                             cross.cline.preds = sync.select.model.results$test.preds[match(names(additional.test.cline.genes), names(sync.select.model.results$test.preds))],
                                             specific = test.cline)
  rownames(additional.test.cline.genes) <- NULL
  additional.target.profiles <- rbind(additional.train.cline.genes,
                                      additional.test.cline.genes)
  
  cell.line.specific.gene.prediction.performances <- 
    ggpubr::ggscatter(additional.target.profiles, 
                      x = "observed.target", 
                      y = "cross.cline.preds",
                      color = "specific", 
                      #palette = pal_uchicago("dark")(9)[c(4:6)][c(1,3)],           # Color by groups "cyl"
                      #shape = "specific",                        # Change point shape by groups "cyl"
                      fullrange = TRUE,                         # Extending the regression line
                      rug = F,                                # Add marginal rug
                      alpha = 1/5,
                      add = "reg.line",
                      show.legend.text = FALSE
    )+
    xlab("Observed traveling ratio of cell line specific genes")+
    ylab("Cross cell line model prediction")+
    scale_color_manual(name = "Expressed only in",
                       breaks = c("HepG2", "K562"),
                       # values = c("HepG2" = pal_uchicago("dark")(9)[c(4:6)][1],
                       #            "K562" = pal_uchicago("dark")(9)[c(4:6)][2]),
                       values = c("HepG2" = alpha('#21908dff',1),
                                  "K562" = alpha("#440154ff",1)),
                       guide = guide_legend(override.aes = list(alpha = 1, size = 6)))+
    ggpubr::stat_cor(aes(color = specific), cor.coef.name = c("rho"), digits = 3)+
    theme_pubr( base_size = plot.base.size)


    cowplot::save_plot(paste0(OUTPUT, PLOTS, train.cline, "_and_", test.cline, "_specific_sample_predictions_", gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target),"_", train.cline, ".pdf"), 
                     plot = cell.line.specific.gene.prediction.performances)
  
  return(list(cline.specific.traveling.ratios = cline.specific.traveling.ratios,
              cline.specific.traveling.ratio.prediction.performances = cline.specific.traveling.ratio.prediction.performances,
              cell.line.specific.gene.prediction.performances = cell.line.specific.gene.prediction.performances))
}
###
# This function generates a shap summary plot
###
visualize.shap.summary <- function(shap.score.matrix, X, train.cline, matrix.type, feature.subspace, model.target, plot.base.size){
  
  # Shap plot of top50 feautures
  shap.summary.plot <- 
    shap.plot.summary.wrap2(shap_score = shap.score.matrix, X = X, top_n=25)+
    #ggtitle(paste0("Shap Summary plot of top 50 features (", train.cline, ")"))+
    theme_pubr(base_size = plot.base.size)+
    scale_color_gradient(low = pal_uchicago()(6)[4], high = pal_uchicago()(6)[5], 
                         breaks = c(0, 1), labels = c("Low", "High "), guide = guide_colorbar(barwidth = 12, 
                                                                                               barheight = 0.3))+
    theme(axis.text.x = element_text(size=9))
  
  cowplot::save_plot(paste0("./plots/shap_summary_plot_",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), 
                     shap.summary.plot)
  return(shap.summary.plot)
} 
###
# This function visualizes the importances scores of static features 
###
visualize.static.features <- function(shap.scores, X, train.cline, matrix.type, feature.subspace, model.target, plot.base.size){
  
  # Shap dependence plots of continous features
  plot.base.size <- 16
  x.axis.text.size <- 12
  shap_long <- shap.prep(shap_contrib = shap.scores, X_train = X)
  
  p1 <- shap.plot.dependence(data_long = shap_long, x = 'housekeeping', y = 'housekeeping', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of houskeeeping width vs. \nhouskeeeping")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Housekeeping gene annotation")
  
  p2 <- shap.plot.dependence(data_long = shap_long, x = 'tx.chr.loc', y = 'tx.chr.loc', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of chrom. location vs. \nchrom. location")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    scale_x_continuous(breaks=1:23,
                     labels=c(paste0("chr", 1:22), "X"))+
    ggtitle("Chromosome location")
    
  p3 <- shap.plot.dependence(data_long = shap_long, x = 'tx.strand', y = 'tx.strand', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of chrom. location vs. \nchrom. location")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Strand specification")
  
  p4 <- shap.plot.dependence(data_long = shap_long, x = 'tx.loc', y = 'tx.loc', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of tx location vs. \ntx location")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Position on the linear genome")
  
  p5 <- shap.plot.dependence(data_long = shap_long, x = 'tx.len', y = 'tx.len', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of transcript length vs. \ntranscript length")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Transcript length")
  
  p6 <- shap.plot.dependence(data_long = shap_long, x = 'tx.tss.width', y = 'tx.tss.width', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of TSS width vs. \nTSS width")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("CAGE tss cluster width")
  
  p7 <- shap.plot.dependence(data_long = shap_long, x = 'tx.tss.at.cont', y = 'tx.tss.at.cont', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of TSS width vs. \nTSS width")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("CAGE tss cluster AT frequency")
  
  p8 <- shap.plot.dependence(data_long = shap_long, x = 'tx.gc.seq', y = 'tx.gc.seq', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of GC content vs. \nGC content")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Transcript GC frequency")
  
  p9 <- shap.plot.dependence(data_long = shap_long, x = 'tx.ex.count', y = 'tx.ex.count', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of exon count vs. \nexon count")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Transcript exon count")
  
  p10 <- shap.plot.dependence(data_long = shap_long, x = 'tx.ex.ratio', y = 'tx.ex.ratio', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of exon ratio vs. \nexon ratio")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Transcript exon density")
  
  p11 <- shap.plot.dependence(data_long = shap_long, x = 'tx.ex.width', y = 'tx.ex.width', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of exon width vs. \nexon width")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Average exon width in transcript")
  
  p12 <- shap.plot.dependence(data_long = shap_long, x = 'tx.ex.seq', y = 'tx.ex.seq', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of exonic sequence vs. \nexonic sequence")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Proportion of exonic sequence per transcript")
  
  p13 <- shap.plot.dependence(data_long = shap_long, x = 'cpg.island.dist', y = 'cpg.island.dist', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of exonic sequence vs. \nexonic sequence")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Distance to nearest CpG island")
  
  p14 <- shap.plot.dependence(data_long = shap_long, x = 'cpg.island.length', y = 'cpg.island.length', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of exonic sequence vs. \nexonic sequence")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Length of nearest CpG island")

  p15 <- shap.plot.dependence(data_long = shap_long, x = 'cpg.island.count', y = 'cpg.island.count', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of exonic sequence vs. \nexonic sequence")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Number of CpG islands")
  
  p16 <- shap.plot.dependence(data_long = shap_long, x = 'cpg.island.percent.cpg', y = 'cpg.island.percent.cpg', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of exonic sequence vs. \nexonic sequence")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Percentage of CpG in the nearest island")
  
  p17 <- shap.plot.dependence(data_long = shap_long, x = 'cpg.island.percent.cg', y = 'cpg.island.percent.cg', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of exonic sequence vs. \nexonic sequence")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Percentage of CG in the nearest island")
  
  p18 <- shap.plot.dependence(data_long = shap_long, x = 'cpg.island.exp.obs', y = 'cpg.island.exp.obs', color_feature = 'Column_WV') +  
    #ggtitle("SHAP values of exonic sequence vs. \nexonic sequence")+
    ylab("Shap value")+
    theme_pubr(base_size = plot.base.size)+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme( axis.text.x = element_text(size=x.axis.text.size))+
    ggtitle("Ratio of observed (CpG_num) to expected ((numC * numG) / IslandLength) CpG in island")
  
  #x gridExtra::grid.arrange(g1, g2, ncol = 2)
  static.feature.plots <- 
    plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,p11, p12, p13, p14, p15, p16, p17, p18,
              labels = "AUTO", 
              #label_size = 10, 
              align = "v", 
              ncol = 2,
              nrow = 9,
              axis = "b", 
              base_width = 14, 
              base_height = 5)
  cowplot::save_plot(paste0("./plots/shap_feature_contributions_", train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".png"), 
                     plot = static.feature.plots, 
                     ncol = 2, 
                     nrow = 9,
                     dpi = 600,
                     base_width = 14, 
                     base_height = 5)
  return(static.feature.plots)
}
###
# This function visualizes the dynamic features
###
visualize.dynamic.features <- function(shap.scores, X, train.cline, matrix.type, feature.subspace, model.target, mean.shap.scores, plot.base.size, append = ""){
  
  # Temp vars to test function
  # shap.scores = shap.scores
  # X = X
  # train.cline = train.cline
  # matrix.type = matrix.type
  # feature.subspace = feature.subspace
  # model.target = model.target
  # plot.base.size = plot.base.size
  
  
  ### Feature class (DBP, RBP etc.) ranking with absolute shap values
  # Retrieve indices of binary binding features of factor bindings
  binding.signals <- colnames(shap.scores)[grepl("^chip|clip", colnames(shap.scores))]
  factor.binding.contrib <- as.data.frame(shap.scores)[ , binding.signals]
  # Make shap values absolute to make contributions undirectional
  absolute.feature.contrib <- as.data.frame(apply(factor.binding.contrib, 2, abs))
  # Filter for non-zero contributions
  aggregate.feature.contrib <- colSums(absolute.feature.contrib)
  aggregate.feature.contrib <- aggregate.feature.contrib[aggregate.feature.contrib>0]
  # Calculate contributions for each class of features (DNA-binding, RNA-binding, ncRNA-binding)
  valid.features <- names(aggregate.feature.contrib)
  # General DNA/RNA bindings
  dna.binding.contrib <- aggregate.feature.contrib[grepl("^chip", valid.features) & !grepl("Proximal", valid.features)]
  rna.binding.contrib <- aggregate.feature.contrib[grepl("^clip", valid.features) & !grepl("Proximal", valid.features)]
  # General ncRNA bindings
  ncrna.dna.binding.contrib <- aggregate.feature.contrib[grepl("^chip", valid.features) & grepl("Proximal", valid.features)]
  ncrna.rna.binding.contrib <- aggregate.feature.contrib[grepl("^clip", valid.features) & grepl("Proximal", valid.features)]
  # Specific DNA bindings on genomic features
  dna.exon.binding.contrib <- aggregate.feature.contrib[grepl("^chip.*exons", valid.features) & !grepl("Proximal", valid.features)]
  dna.intron.binding.contrib <- aggregate.feature.contrib[grepl("^chip.*introns", valid.features) & !grepl("Proximal", valid.features)]
  dna.five.prime.binding.contrib <- aggregate.feature.contrib[grepl("^chip.*five_prime", valid.features) & !grepl("Proximal", valid.features)]
  dna.three.prime.binding.contrib <- aggregate.feature.contrib[grepl("^chip.*three_prime", valid.features) & !grepl("Proximal", valid.features)]
  # Specific RNA bindings on genomic features
  rna.exon.binding.contrib <- aggregate.feature.contrib[grepl("^clip.*exons", valid.features) & !grepl("Proximal", valid.features)]
  rna.intron.binding.contrib <- aggregate.feature.contrib[grepl("^clip.*introns", valid.features) & !grepl("Proximal", valid.features)]
  rna.five.prime.binding.contrib <- aggregate.feature.contrib[grepl("^clip.*five_prime", valid.features) & !grepl("Proximal", valid.features)]
  rna.three.prime.binding.contrib <- aggregate.feature.contrib[grepl("^clip.*three_prime", valid.features) & !grepl("Proximal", valid.features)]
  
  non.binding.signals <- colnames(shap.scores)[!grepl("^chip|clip", colnames(shap.scores))]
  non.binding.signals.contrib <- as.data.frame(shap.scores)[ , non.binding.signals]
  absolute.non.binding.feature.contrib <- as.data.frame(apply(non.binding.signals.contrib, 2, abs))
  aggregate.non.binding.feature.contrib <- colSums(absolute.non.binding.feature.contrib)
  
  # Aggregate feature contributions per class
  class.contributions <- 
    list(DNA = dna.binding.contrib,
         RNA = rna.binding.contrib,
         DNA_ncRNA = ncrna.dna.binding.contrib,
         RNA_ncRNA = ncrna.rna.binding.contrib,
         DNA_exons = dna.exon.binding.contrib,
         DNA_introns = dna.intron.binding.contrib,
         DNA_five_prime = dna.five.prime.binding.contrib,
         DNA_three_prime = dna.three.prime.binding.contrib,
         RNA_exons = rna.exon.binding.contrib,
         RNA_introns = rna.intron.binding.contrib,
         RNA_five_prime = rna.five.prime.binding.contrib,
         RNA_three_prime = rna.three.prime.binding.contrib)
  
  aggregate.class.contributions <- lapply(class.contributions,sum)
  aggregate.class.contributions <- unlist(aggregate.class.contributions)
  aggregate.class.contributions <- sort(aggregate.class.contributions, decreasing = T)
  
  percent.dna.sequence.contrib <- sum(aggregate.non.binding.feature.contrib[-1]) / sum(aggregate.class.contributions["DNA"], aggregate.class.contributions["RNA"])*100
  percent.dna.sequence.contrib.housekeeping.incl <- sum(aggregate.non.binding.feature.contrib) / sum(aggregate.class.contributions["DNA"], aggregate.class.contributions["RNA"])*100
  
  ## Normalize by number of binding sites in each class
  # General DNA/RNA bindings
  dna.bindings <- X[grepl("^chip", valid.features) & !grepl("Proximal", valid.features)]
  dna.bindings <- sum(dna.bindings == 1, na.rm = T)
  rna.bindings <- X[grepl("^clip", valid.features) & !grepl("Proximal", valid.features)]
  rna.bindings <- sum(rna.bindings == 1, na.rm = T)
  
  # General ncRNA bindings
  ncrna.dna.bindings <- X[grepl("^chip", valid.features) & grepl("Proximal", valid.features)]
  ncrna.dna.bindings <- sum(ncrna.dna.bindings == 1, na.rm = T)
  ncrna.rna.bindings <- X[grepl("^clip", valid.features) & grepl("Proximal", valid.features)]
  ncrna.rna.bindings <- sum(ncrna.rna.bindings == 1, na.rm = T)
  
  # Specific DNA bindings on genomic features
  dna.exon.bindings <- X[grepl("^chip.*exons", valid.features) & !grepl("Proximal", valid.features)]
  dna.exon.bindings <- sum(dna.exon.bindings == 1, na.rm = T)
  dna.intron.bindings <- X[grepl("^chip.*introns", valid.features) & !grepl("Proximal", valid.features)]
  dna.intron.bindings <- sum(dna.intron.bindings == 1, na.rm = T)
  dna.five.prime.bindings<- X[grepl("^chip.*five_prime", valid.features) & !grepl("Proximal", valid.features)]
  dna.five.prime.bindings <- sum(dna.five.prime.bindings == 1, na.rm = T)
  dna.three.prime.bindings <- X[grepl("^chip.*three_prime", valid.features) & !grepl("Proximal", valid.features)]
  dna.three.prime.bindings <- sum(dna.three.prime.bindings == 1, na.rm = T)
  # Specific RNA bindings on genomic features
  rna.exon.bindings <- X[grepl("^clip.*exons", valid.features) & !grepl("Proximal", valid.features)]
  rna.exon.bindings <- sum(rna.exon.bindings == 1, na.rm = T)
  rna.intron.bindings <- X[grepl("^clip.*introns", valid.features) & !grepl("Proximal", valid.features)]
  rna.intron.bindings <- sum(rna.intron.bindings == 1, na.rm = T)
  rna.five.prime.bindings <- X[grepl("^clip.*five_prime", valid.features) & !grepl("Proximal", valid.features)]
  rna.five.prime.bindings <- sum(rna.five.prime.bindings == 1, na.rm = T)
  rna.three.prime.bindings <- X[grepl("^clip.*three_prime", valid.features) & !grepl("Proximal", valid.features)]
  rna.three.prime.bindings <- sum(rna.three.prime.bindings == 1, na.rm = T)
  
  total.binding.events <- dna.bindings+rna.bindings
  weighted.aggregate.class.contributions <- vector()
  weighted.aggregate.class.contributions["DNA"] <- (aggregate.class.contributions["DNA"]*dna.bindings)/total.binding.events
  weighted.aggregate.class.contributions["RNA"] <- (aggregate.class.contributions["RNA"]*rna.bindings)/total.binding.events
  weighted.aggregate.class.contributions["DNA_five_prime"] <- (aggregate.class.contributions["DNA_five_prime"]*dna.five.prime.bindings)/total.binding.events
  weighted.aggregate.class.contributions["DNA_introns"] <- (aggregate.class.contributions["DNA_introns"]*dna.intron.bindings)/total.binding.events
  weighted.aggregate.class.contributions["DNA_exons"] <- (aggregate.class.contributions["DNA_exons"]*dna.exon.bindings)/total.binding.events
  weighted.aggregate.class.contributions["DNA_three_prime"] <- (aggregate.class.contributions["DNA_three_prime"]*dna.three.prime.bindings)/total.binding.events
  weighted.aggregate.class.contributions["RNA_introns"] <- (aggregate.class.contributions["RNA_introns"]*rna.intron.bindings)/total.binding.events
  weighted.aggregate.class.contributions["RNA_exons"] <- (aggregate.class.contributions["RNA_exons"]*rna.exon.bindings)/total.binding.events
  weighted.aggregate.class.contributions["RNA_three_prime"] <- (aggregate.class.contributions["RNA_three_prime"]*rna.three.prime.bindings)/total.binding.events
  weighted.aggregate.class.contributions["RNA_five_prime"] <- (aggregate.class.contributions["RNA_five_prime"]*rna.five.prime.bindings)/total.binding.events
  weighted.aggregate.class.contributions["DNA_ncRNA"] <- (aggregate.class.contributions["DNA_ncRNA"]*ncrna.dna.bindings)/total.binding.events
  weighted.aggregate.class.contributions["RNA_ncRNA"] <- (aggregate.class.contributions["RNA_ncRNA"]*ncrna.rna.bindings)/total.binding.events
  weighted.aggregate.class.contributions[is.nan(weighted.aggregate.class.contributions)] <- 0
  weighted.aggregate.class.contributions <- weighted.aggregate.class.contributions[!names(weighted.aggregate.class.contributions) %in% c("DNA", "RNA")]
  
  weighted.aggregate.class.contributions <- as.data.frame(weighted.aggregate.class.contributions)
  weighted.aggregate.class.contributions <- cbind(weighted.aggregate.class.contributions, class = rownames(weighted.aggregate.class.contributions))
  weighted.aggregate.class.contributions$aggregate.class.contributions <- rescale(weighted.aggregate.class.contributions$aggregate.class.contributions)
  aggregate.weighted.class.contributions.plot <- 
    ggpubr::ggbarplot(weighted.aggregate.class.contributions, 
                      x="class",
                      y = "weighted.aggregate.class.contributions",
                      color = "black",
                      fill = rgb(0.2,0.4,0.6,1),
                      palette = "uchicago",           
                      sort.val = "asc",          
                      sort.by.groups = FALSE,     
                      x.text.angle = 90,         
                      xlab = "Feature Class",
                      ylab = "Aggregate SHAP Feature Class Contributions (rescaled)",
                      rotate = TRUE,
                      main = "DNA & RNA Binding Feature Class Contributions"
    )+
    theme_pubr(base_size = 10)+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  
  relative.aggregate.class.contributions <- vector()
  relative.aggregate.class.contributions["DNA"] <- aggregate.class.contributions["DNA"]/dna.bindings
  relative.aggregate.class.contributions["RNA"] <- aggregate.class.contributions["RNA"]/rna.bindings
  relative.aggregate.class.contributions["DNA_five_prime"] <- aggregate.class.contributions["DNA_five_prime"]/dna.five.prime.bindings
  relative.aggregate.class.contributions["DNA_introns"] <- aggregate.class.contributions["DNA_introns"]/dna.intron.bindings
  relative.aggregate.class.contributions["DNA_exons"] <- aggregate.class.contributions["DNA_exons"]/dna.exon.bindings
  relative.aggregate.class.contributions["DNA_three_prime"] <- aggregate.class.contributions["DNA_three_prime"]/dna.three.prime.bindings
  relative.aggregate.class.contributions["RNA_introns"] <- aggregate.class.contributions["RNA_introns"]/rna.intron.bindings
  relative.aggregate.class.contributions["RNA_exons"] <- aggregate.class.contributions["RNA_exons"]/rna.exon.bindings
  relative.aggregate.class.contributions["RNA_three_prime"] <- aggregate.class.contributions["RNA_three_prime"]/rna.three.prime.bindings
  relative.aggregate.class.contributions["RNA_five_prime"] <- aggregate.class.contributions["RNA_five_prime"]/rna.five.prime.bindings
  relative.aggregate.class.contributions["DNA_ncRNA"] <- aggregate.class.contributions["DNA_ncRNA"]/ncrna.dna.bindings
  relative.aggregate.class.contributions["RNA_ncRNA"] <- aggregate.class.contributions["RNA_ncRNA"]/ncrna.rna.bindings
  relative.aggregate.class.contributions[is.nan(relative.aggregate.class.contributions)] <- 0
  relative.aggregate.class.contributions <- relative.aggregate.class.contributions[!names(relative.aggregate.class.contributions) %in% c("DNA", "RNA")]
  
  relative.aggregate.class.contributions <- as.data.frame(relative.aggregate.class.contributions)
  relative.aggregate.class.contributions <- cbind(relative.aggregate.class.contributions, class = rownames(relative.aggregate.class.contributions))
  relative.aggregate.class.contributions$aggregate.class.contributions <- rescale(relative.aggregate.class.contributions$aggregate.class.contributions)
  relative.aggregate.class.contributions.plot <- 
    ggpubr::ggbarplot(relative.aggregate.class.contributions, 
                      x="class",
                      y = "relative.aggregate.class.contributions",
                      color = "black",
                      fill = rgb(0.2,0.4,0.6,1),
                      palette = "uchicago",           
                      sort.val = "asc",          
                      sort.by.groups = FALSE,     
                      x.text.angle = 90,         
                      xlab = "Feature Class",
                      ylab = "Aggregate SHAP Feature Class Contributions (rescaled)",
                      rotate = TRUE,
                      main = "DNA & RNA Binding Feature Class Contributions"
    )+
    theme_pubr(base_size = 10)+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  
  ### Plot class contributions
  aggregate.non.binding.feature.contrib <- as.data.frame(aggregate.non.binding.feature.contrib)
  aggregate.non.binding.feature.contrib <- cbind(aggregate.non.binding.feature.contrib, class = rownames(aggregate.non.binding.feature.contrib))
  aggregate.non.binding.feature.contrib$aggregate.non.binding.feature.contrib <- rescale(aggregate.non.binding.feature.contrib$aggregate.non.binding.feature.contrib)
  aggregate.non.binding.feature.contrib.plot <- 
    ggpubr::ggbarplot(aggregate.non.binding.feature.contrib, 
                      x="class",
                      y = "aggregate.non.binding.feature.contrib",
                      color = "black", 
                      fill = rgb(0.2,0.4,0.6,0.6),
                      palette = "uchicago",           
                      sort.val = "asc",          
                      sort.by.groups = FALSE,     
                      x.text.angle = 90,         
                      xlab = "Feature",
                      ylab = "Aggregate SHAP Feature Contributions (rescaled)",
                      rotate = TRUE,
                      main = "DNA sequence feature contributions (incl. housekeeping feature)"
    )+
    theme_pubr(base_size = 10)+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "shap_DNA_sequence_contributions_",append, "_", train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), aggregate.non.binding.feature.contrib.plot)
  
  
  #### Binding features contributions in either direction with absolute shap values, to show binding features with biggest effect
  factor.contrib <- sort(aggregate.feature.contrib, decreasing = T)
  factor.contrib <- data.frame(contrib = factor.contrib, factor.binding = names(factor.contrib))
  factor.contrib$contrib.rescaled <- rescale(factor.contrib$factor.contrib)
  factor.contrib.plot <- 
    ggpubr::ggbarplot(head(factor.contrib, 25), 
                      x="factor.binding",
                      y = "contrib",
                      color = "black",            # Set bar border colors to white
                      palette = "uchicago",            # jco journal color palett. see ?ggpar
                      sort.val = "asc",          # Sort the value in descending order
                      sort.by.groups = FALSE,     # Don't sort inside each group
                      x.text.angle = 90,          # Rotate vertically x axis texts
                      ylab = "Model contributions (rescaled)",
                      rotate = TRUE,
                      ggtheme = theme_pubr(base_size = plot.base.size),
                      main = "Factor's binding contributions"
    )
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "shap_total_factor_feature_contributions_",append, "_", train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), 
                     factor.contrib.plot)
  
  ### Directional Factor contributions direction with mean shap values, to show factor's effects
  avg.factor.contrib <- mean.shap.scores
  avg.factor.contrib <- avg.factor.contrib[grepl("^chip|clip", names(avg.factor.contrib))]
  avg.factor.contrib <- avg.factor.contrib[avg.factor.contrib>0]
  avg.factor.contrib <- data.frame(contrib = avg.factor.contrib, factor.binding = names(avg.factor.contrib))
  avg.factor.contrib$factor.contrib <- rescale(avg.factor.contrib$factor.contrib)
  avg.factor.contrib.plot <- 
    ggpubr::ggbarplot(head(avg.factor.contrib, 25),
                      x="factor.binding",
                      y = "contrib",
                      color = "black",            # Set bar border colors to white
                      palette = "jco",            # jco journal color palett. see ?ggpar
                      sort.val = "asc",          # Sort the value in descending order
                      sort.by.groups = FALSE,     # Don't sort inside each group
                      x.text.angle = 90,          # Rotate vertically x axis texts
                      ylab = "Average model contributions (rescaled)",
                      rotate = TRUE,
                      ggtheme = theme_minimal(),
                      main = "Average factor's binding contributions"
    )
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "shap_avg_factor_feature_contributions_",append, "_", train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), 
                     avg.factor.contrib.plot)
  
  
  ### Binding features contributions in either direction with absolute shap values, to show binding features with biggest effect
  contributing.factors <- factor.contrib$factor.binding
  contributing.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*", 
                               "", contributing.factors)
  contributing.factors <- unique(contributing.factors)
  cum.factor.contrib <- 
    lapply(contributing.factors, function(factor){
      cum.factor.contrib <- sum(factor.contrib$contrib[grepl(factor, factor.contrib$factor.binding)])
      data.frame(factor = factor, contrib = cum.factor.contrib)
    })
  cum.factor.contrib <- do.call(rbind, cum.factor.contrib)
  cum.factor.contrib <- cum.factor.contrib[order(cum.factor.contrib$contrib, decreasing = T), ]
  cum.factor.contrib$contrib.rescaled <- rescale(cum.factor.contrib$contrib)
  cum.factor.contrib.plot <- 
    ggpubr::ggbarplot(head(cum.factor.contrib, 50), 
                      x="factor",
                      y = "contrib.rescaled",
                      color = "black",            # Set bar border colors to white
                      palette = "uchicago",            # jco journal color palett. see ?ggpar
                      sort.val = "asc",          # Sort the value in descending order
                      sort.by.groups = FALSE,     # Don't sort inside each group
                      x.text.angle = 90,          # Rotate vertically x axis texts
                      xlab = "DNA & RNA binding factors",
                      ylab = "Cumulative factor contribution (rescaled)",
                      rotate = TRUE,
                      ggtheme = theme_pubr(base_size = plot.base.size),
                      main = "Factor model contributions"
    )
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "shap_total_factor_contributions_",append, "_", train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), 
                     cum.factor.contrib.plot)
  
  
  ### Directional factor contributions
  directional.factor.contrib <- 
    lapply(contributing.factors, function(factor){
      factor.binding.patterns <- X[, grepl(factor, colnames(X))]
      factor.features <- colnames(X)[grepl(factor, colnames(X))]
      factor.feature.contrib <- 
        lapply(factor.features, function(factor.binding.feature){
          bound.samples <- which(factor.binding.patterns[, factor.binding.feature]==1)
          sum(factor.binding.contrib[bound.samples, factor.binding.feature])
        })
      factor.feature.contrib <- sum(unlist(factor.feature.contrib))
      data.frame(factor = factor, factor.contrib = factor.feature.contrib, direction = if (factor.feature.contrib<0) "decrease" else "increase" )
    })
  directional.factor.contrib <- do.call(rbind, directional.factor.contrib)
  directional.factor.contrib <- directional.factor.contrib[order(abs(directional.factor.contrib$factor.contrib), decreasing = F),]
  directional.factor.contrib$direction <- as.factor(directional.factor.contrib$direction)
  directional.factor.contrib.plot <- 
    ggpubr::ggbarplot(tail(directional.factor.contrib, 50), 
                      x="factor",
                      y = "factor.contrib",
                      fill = "direction",
                      color = "black",            # Set bar border colors to white
                      palette = "uchicago",            # jco journal color palett. see ?ggpar
                      #sort.val = "asc",          # Sort the value in descending order
                      legend.title = "Impact on target",
                      sort.by.groups = FALSE,     # Don't sort inside each group
                      x.text.angle = 90,          # Rotate vertically x axis texts
                      xlab = "DNA & RNA binding factors",
                      ylab = "Cumulative directional factor contribution",
                      rotate = TRUE,
                      ggtheme = theme_pubr(base_size = plot.base.size),
                      main = "Directional factor model contributions"
    )
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "shap_total_diretional_factor_contributions_", append, "_", train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), 
                     directional.factor.contrib.plot)
  
  
  # Seperate by positive and negative regulators
  directional.factor.contrib <- split(directional.factor.contrib, directional.factor.contrib$direction)
  directional.factor.contrib <- 
    lapply(directional.factor.contrib, function(directional.effects){directional.effects[order(abs(directional.effects$factor.contrib), decreasing = F),]})
  #directional.factor.contrib$decrease$factor.contrib <- rescale(directional.factor.contrib$decrease$factor.contrib, c(-1,0))
  #directional.factor.contrib$increase$factor.contrib <- rescale(directional.factor.contrib$increase$factor.contrib, c(0,1))
  
  dbp.factors <- 
    lapply(CHIPseq, function(peaks){
      unique(peaks$hgnc_symbol)
    }) %>% unlist %>% unique 
  
  # Retrieve all RNA binding factors
  rbp.factors <- 
    lapply(eCLIPseq, function(peaks){
      unique(peaks$hgnc_symbol)
    }) %>% unlist %>% unique 
  
  top.directional.effectors <- 
    lapply(directional.factor.contrib, function(directional.effects){
      tail(directional.effects, 10)
    })
  top.directional.effectors <- do.call(rbind, top.directional.effectors)
  
  top.directional.effectors$factor.type <- "factor"
  top.directional.effectors$factor.type[top.directional.effectors$factor %in% dbp.factors] <- "DBP"
  top.directional.effectors$factor.type[top.directional.effectors$factor %in% rbp.factors] <- "RBP"
  top.directional.effectors$factor.type[top.directional.effectors$factor %in% rbp.factors &
                                          top.directional.effectors$factor %in% dbp.factors ] <- "DBP/RBP"
  
  top.directional.effectors.plot <- 
    ggpubr::ggbarplot(top.directional.effectors, 
                      x="factor",
                      y = "factor.contrib",
                      fill = "factor.type",
                      color = "black",            # Set bar border colors to white
                      palette = "uchicago",            # jco journal color palett. see ?ggpar
                      sort.val = "desc",          # Sort the value in descending order
                      legend.title = "Factor Type",
                      sort.by.groups = FALSE,     # Don't sort inside each group
                      x.text.angle = 90,          # Rotate vertically x axis texts
                      xlab = "DNA & RNA binding factors",
                      ylab = "Aggregate directional factor contribution",
                      rotate = TRUE,
                      lab.size = 8,
                      main = "Directional factor model contributions"
    )+
    theme_pubr(base_size = plot.base.size, legend ="top")
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "shap_top_total_diretional_factor_contributions_",append, "_", train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), top.directional.effectors.plot)
  
  
  results <- list(class.contributions = class.contributions,
                  top.directional.effectors = top.directional.effectors,
                  aggregate.weighted.class.contributions.plot = aggregate.weighted.class.contributions.plot,
                  relative.aggregate.class.contributions.plot = relative.aggregate.class.contributions.plot,
                  aggregate.non.binding.feature.contrib = aggregate.non.binding.feature.contrib,
                  aggregate.non.binding.feature.contrib.plot = aggregate.non.binding.feature.contrib.plot,
                  factor.contrib = factor.contrib,
                  factor.contrib.plot = factor.contrib.plot,
                  avg.factor.contrib.plot = avg.factor.contrib.plot,
                  cum.factor.contrib = cum.factor.contrib,
                  cum.factor.contrib.plot = cum.factor.contrib.plot,
                  directional.factor.contrib = directional.factor.contrib,
                  directional.factor.contrib.plot = directional.factor.contrib.plot,
                  top.directional.effectors.plot = top.directional.effectors.plot,
                  percent.dna.sequence.contrib = percent.dna.sequence.contrib,
                  percent.dna.sequence.contrib.housekeeping.incl = percent.dna.sequence.contrib.housekeeping.incl)
  
  return(results)
  
}
###
# This function visualizes the top directional dynamic features 
###
visualize.individual.directional.dynamic.factor.effects <- function(top.directional.effectors, shap.scores, X, train.cline, matrix.type, feature.subspace, model.target, plot.base.size){
  
  # top.directional.effectors = cluster.effectors$all.cluster.factor.contribs[match(optimal.cluster.factor.directional.contrib$top.directional.binding.factors, cluster.effectors$all.cluster.factor.contribs$factor), ]
  # shap.scores
  # X = X
  # train.cline = train.cline
  # matrix.type = matrix.type
  # feature.subspace = feature.subspace
  # model.target = model.target
  # plot.base.size = plot.base.size
  
  # top.directional.effectors = optimal.cluster.factor.directional.contrib$top.directional.binding.factors
  # shap.scores
  # X = X
  # train.cline = train.cline
  # matrix.type = matrix.type
  # feature.subspace = feature.subspace
  # model.target = model.target
  # plot.base.size = plot.base.size
  
  top.factors <- top.directional.effectors
  shap_long <- shap.prep(shap_contrib = shap.scores, X_train = X)
  individual.directional.factor.contrib <- 
    lapply(top.factors, function(factor){
      #factor.feature.contrib <- factor.binding.contrib[ , grepl(factor, colnames(factor.binding.contrib))]
      shap_long.factor <- shap_long[grepl(factor, shap_long$variable),]
      shap_long.factor$variable <- factor
      shap_long.factor$rfvalue <- as.factor(shap_long.factor$rfvalue)
      shap_long.factor$stdfvalue <- as.factor(shap_long.factor$stdfvalue)
      binding.counts <- table(shap_long.factor$rfvalue)
      non.binding.max = max(shap_long.factor$value[shap_long.factor$rfvalue==0])+0.005
      binding.max = max(shap_long.factor$value[shap_long.factor$rfvalue==1])+0.005
      non.binding.count <- table(shap_long.factor$rfvalue)[1]
      binding.count <- table(shap_long.factor$rfvalue)[2]
      ylim.max <- max(shap_long.factor$value)
      ylim.min <- min(shap_long.factor$value)
      g <- shap.plot.dependence(data_long = shap_long.factor, x = factor, y = factor, smooth = F) +  
        ggtitle(factor)+
        xlab(paste0("Binding"))+
        ylab(paste0("Shap value"))+
        theme_pubr(base_size = 12)+
        theme(plot.margin=margin(t = 0.5, unit="cm"))+
        geom_text(x=1, y=non.binding.max+0.05, label=non.binding.count, color="blue", fontface="bold", size = 4)+
        geom_text(x=2, y=binding.max+0.05, label=binding.count, color="blue",fontface="bold", size = 4)+
        ylim(c(ylim.min, ylim.max+0.1))
      g
    })
  
  individual.directional.factor.contrib.plots <- 
    plot_grid(plotlist = individual.directional.factor.contrib,
              #labels = as.character(1:8), 
              align = "h", 
              nrow =2, 
              ncol = 5)
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "individual_diretional_factor_contributions_",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), 
                     individual.directional.factor.contrib.plots, 
                     nrow = 2, ncol = 4)
  
  return(individual.directional.factor.contrib.plots)
}
###
# This function visualizes the shap forces and shap force clusters 
###
visualize.shap.forces <- function(shap.scores,mean.shap.scores, X, train.cline, matrix.type, feature.subspace, model.target, plot.base.size){
  
  # Temp vars to test function
  # shap.scores = shap.scores
  # mean.shap.scores = mean.shap.scores
  # X = X
  # train.cline = train.cline
  # matrix.type = matrix.type
  # feature.subspace = feature.subspace
  # model.target = model.target
  # plot.base.size = plot.base.size
  
  force.plot_data <- shap.prep.stack.data(shap_contrib = shap.scores, 
                                          top_n = 5, 
                                          n_groups = 2, 
                                          cluster_method = "ward.D")
  force.plot <- 
    shap.plot.force_plot(force.plot_data, zoom_in_group = 2, y_parent_limit = c(-7,7))+
    ggtitle("Individual feature's impact on the traveling ratio of each transcript (top 5 features)" )+
    theme_pubr(base_size = plot.base.size)+
    guides(shape = guide_legend(override.aes = list(size = 9)))+
    guides(color = guide_legend(override.aes = list(size = 9)))+
    theme(legend.title = element_text(size = 9), 
          legend.text = element_text(size = 9),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "force_plot_",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), force.plot)
  
  # Shap force plot of clustered samples
  force.plot.cluster.plot <- 
    shap.plot.force_plot_bygroup(force.plot_data)+
    ggtitle("Clustered Shap Force Plot")+
    theme_pubr(base_size = plot.base.size)+
    rremove("legend")+
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))+
    scale_fill_jco()
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "cluster_force_plot_",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), force.plot.cluster.plot)
  
  return(list(force.plot_data = force.plot_data, 
              force.plot = force.plot,
              force.plot.cluster.plot = force.plot.cluster.plot))
}
###
#
###
visualize.shap.interactions <- function(model, X, all.optimal.factors){
  
  model <- minimal.model.results$select.feature.model$model
  X <- minimal.model.results$select.feature.model$X
  shap_long <- shap.prep(xgb_model = model, X_train = X)
  shap_int <- shap.prep.interaction(xgb_model = model, X_train = X)
  aggregate.shap_int <- as.data.frame(shap_int)
  abs.aggregate.shap_int <- colSums(abs(aggregate.shap_int))
  abs.aggregate.shap_int <- sort(abs.aggregate.shap_int, decreasing = T)
  shap.plot.dependence(data_long = shap_long,
                       data_int = shap_int,
                       x= "tx.len",
                       y = "chip.RBFOX2.five_prime", 
                       color_feature = "tx.len")
}
###
#
###
calculate.prior.knowledge.contributions <- function(shap.scores,mean.shap.scores, X, train.cline, matrix.type, feature.subspace, model.target, plot.base.size){
# 
#   shap.scores = shap.scores
#   mean.shap.scores = mean.shap.scores
#   X = X
#   train.cline = train.cline
#   matrix.type = matrix.type
#   feature.subspace = feature.subspace
#   model.target = model.target
#   plot.base.size = plot.base.size

  
  sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
  sub.factor.feature.spaces <- 
    build.sequence.specific.biologically.functional.feature.sub.spaces(stratify = sub.factor.feature.spaces)
  
  feature.types <- c("Chromatin", "Initiation", "Elongation", "7SK.Binding", "Elongation+7SK", "Termination","Splicing", "Processing")
  feature.types <- c(feature.types,  as.vector(outer(feature.types, c("ss", "nss"), paste, sep=":"))) 
  #feature.types <- rev(feature.types)
  
  process.contributions <- 
  lapply(feature.types, function(process){
    # Retrieve process specific factors
    process.factors <- sub.factor.feature.spaces[[process]][[train.cline]]
    # Check if there are any factors
    nfactors <- length(process.factors)
    if(nfactors>0){
      # Make a pattern out of factors to macth contributions matrix
      process.factors <- paste0(paste0(process.factors, "\\."), collapse = "|")
      # Retrieve mean contrib of each factor feature
      mean.contrib <- mean.shap.scores[grepl(process.factors, names(mean.shap.scores))]
      mean.contrib <- mean.contrib[mean.contrib>0]
      # Retrieve contributions
      contribs <- as.data.frame(shap.scores)
      contribs <- contribs[, grepl(process.factors, colnames(contribs))]
      if(length(mean.contrib)>0){
        # Calculate mean total contrib of the model factors
        mean.total.contrib <- mean(mean.contrib)
        # Subset contributions by valid factor features
        binding.signals <- colnames(contribs)[grepl("^chip|clip", colnames(contribs))]
        nfactors <-  gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*",
                          "", binding.signals) %>% unique() %>% length()
        total.process.contrib <- sum(abs(contribs))
        total.process.contrib <- 
        data.frame(contrib = total.process.contrib, 
                   nfactors = nfactors, 
                   mean.shap = mean.total.contrib)
        return(total.process.contrib)
      }
    }
    return(0)
  })
  names(process.contributions) <- feature.types
  process.contributions <- do.call(rbind, process.contributions)
  process.contributions$subspace <- rownames(process.contributions)
  process.contributions$subspace <- factor(process.contributions$subspace, levels = feature.types)
  
  # Transformations for ordering plot
  process.contributions$sequence.specific <- grepl(":ss|:nss", process.contributions$subspace)
  process.contributions$sequence.specific[!process.contributions$sequence.specific] <- "General"
  process.contributions$sequence.specific[process.contributions$sequence.specific == "TRUE"] <- 
    grepl(":ss", process.contributions$subspace[process.contributions$sequence.specific == "TRUE"])
  process.contributions$sequence.specific[process.contributions$sequence.specific=="TRUE"] <- "sequence.specific"
  process.contributions$sequence.specific[process.contributions$sequence.specific=="FALSE"] <- "non-sequence.specific"
  process.contributions$sequence.specific <- factor(process.contributions$sequence.specific)
  
  process.contributions$subspace <- gsub(":ss|:nss", "", process.contributions$subspace)
  process.contributions$subspace <- factor(  process.contributions$subspace, rev(feature.types))

  # Weight contributions by the number of factors
  #process.contributions$contrib <- (process.contributions$contrib*process.contributions$nfactors)/sum(process.contributions$nfactors)
  
  # Plot model results
  process.contributions.plot <- 
    ggbarplot(process.contributions, 
              x = "subspace", 
              y = "contrib",
              #fill = rgb(0.2,0.4,0.6,0.6),
              fill = "lightgrey",
              #fill = "mean.shap", 
              #label = paste0(paste0("Rsqrd. ", sub.model.performances$test.rsqrd), "\n#Factors ", sub.model.performances$nfactors, "\nContrib. ", format(round(sub.model.performances$mean.shap, 3), nsmall=3) ),
              label = paste0(process.contributions$nfactors),
              rotate = T,
              #lab.pos = "in",
              lab.col = "black",
              lab.font = "bold",
              lab.size = 4,
              lab.hjust = -0.5,
              lab.vjust = 0.5,
              position = position_dodge(0.9))+
    facet_grid(~sequence.specific)+
    xlab("Process")+
    ylab("Contribution")+
    theme_pubr(base_size = plot.base.size, legend = "right")+
    #scale_fill_material("grey")+
    #scale_fill_material("blue-grey")+
    #scale_fill_uchicago()+
    theme(panel.spacing = unit(2, "lines"))+
    scale_y_continuous(limits = c(0, (max(process.contributions$contrib)+1000)))
  
  # cowplot::save_plot(paste0(OUTPUT, PLOTS, "sequence_specific_go_term_model_performances",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"),
  #                    model.performance.plot)
  
  return(process.contributions.plot)
}
###
#
###
visualize.cluster.factors <- function(force.plot_data, shap.scores, X, Y){
  
  # Temp vars to test function
  # force.plot_data = shap.force.plots$force.plot_data
  # shap.scores = shap.scores
  # X = X
  # Y=Y
  
  sample.clustering <- data.frame(id = force.plot_data$ID, cluster = force.plot_data$group)
  sample.clustering <- sample.clustering[order(sample.clustering$id, decreasing = F),]
  model.target.clusters <- data.frame(target = Y, cluster = as.factor(sample.clustering$cluster))
  
  cluster.feature.contrib <- shap.scores
  cluster.feature.contrib$cluster[match(rownames(model.target.clusters), rownames(cluster.feature.contrib))] <- model.target.clusters$cluster
  cluster.feature.contrib <- split(cluster.feature.contrib, cluster.feature.contrib$cluster)
  
  cluster.factor.contrib <- 
    lapply(cluster.feature.contrib, function(cluster.contribs){
      #total.cluster.contrib <- sum(cluster.contribs)
      #cluster.contribs <- cluster.feature.contrib[[1]]
      # Retrieve samples from cluster
      cluster.samples <- rownames(cluster.contribs)
      # Retrieve binding features
      cluster.contribs <- cluster.contribs[, grepl("^chip|clip", colnames(cluster.contribs))]
      contributing.factors <-  gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*", "", colnames(cluster.contribs))
      # Retrieve unique binding factors
      contributing.factors <- unique(contributing.factors)
      # For each factor calculate its contribution
      cum.factor.contrib <- 
        lapply(contributing.factors, function(factor){
          factor.binding.patterns <- as.data.frame(X[match(cluster.samples, rownames(X)), grepl(factor, colnames(X))])
          colnames(factor.binding.patterns) <- colnames(X)[ grepl(factor, colnames(X))]
          factor.features <- colnames(factor.binding.patterns)
          factor.feature.contrib <- 
            lapply(factor.features, function(factor.binding.feature){
              bound.samples <- which(factor.binding.patterns[, factor.binding.feature]==1)
              sum(cluster.contribs[bound.samples, factor.binding.feature])
            })
          factor.feature.contrib <- sum(unlist(factor.feature.contrib))#/total.cluster.contrib
          data.frame(factor = factor, factor.contrib = factor.feature.contrib, direction = if (factor.feature.contrib<0) "decrease" else "increase" )
        })
      cum.factor.contrib <- do.call(rbind, cum.factor.contrib)
      #cum.factor.contrib <- cum.factor.contrib[cum.factor.contrib$factor.contrib != 0,]
      #cum.factor.contrib <- cum.factor.contrib[cum.factor.contrib$factor.contrib != 0,]$direction[which(cum.factor.contrib$factor.contrib == 0, arr.ind = T)] <- "none"
      cum.factor.contrib$direction <- as.factor(cum.factor.contrib$direction)
      cum.factor.contrib <- split(cum.factor.contrib, cum.factor.contrib$direction)
      cum.factor.contrib <- 
        lapply(cum.factor.contrib, function(contribs){contribs[order(abs(contribs$factor.contrib), decreasing = F),]})
      return(cum.factor.contrib)
    })
  
  dbp.factors <- 
    lapply(CHIPseq, function(peaks){
      unique(peaks$hgnc_symbol)
    }) %>% unlist %>% unique 
  
  # Retrieve all RNA binding factors
  rbp.factors <- 
    lapply(eCLIPseq, function(peaks){
      unique(peaks$hgnc_symbol)
    }) %>% unlist %>% unique 
  
  top.cluster.effectors <- 
    lapply(cluster.factor.contrib, function(cluster.factor.contributions){
      top.directional.effectors <- 
        lapply(cluster.factor.contributions, function(directional.effects){
          tail(directional.effects, 10)
        })
      top.directional.effectors <- do.call(rbind, top.directional.effectors)
      
      top.directional.effectors$factor.type <- "factor"
      top.directional.effectors$factor.type[top.directional.effectors$factor %in% dbp.factors] <- "DBP"
      top.directional.effectors$factor.type[top.directional.effectors$factor %in% rbp.factors] <- "RBP"
      top.directional.effectors$factor.type[top.directional.effectors$factor %in% rbp.factors &
                                              top.directional.effectors$factor %in% dbp.factors ] <- "DBP/RBP"
      
      top.directional.effectors.plot <- 
        ggpubr::ggbarplot(top.directional.effectors, 
                          x="factor",
                          y = "factor.contrib",
                          fill = "factor.type",
                          color = "black",            # Set bar border colors to white
                          palette = "uchicago",            # jco journal color palett. see ?ggpar
                          sort.val = "desc",          # Sort the value in descending order
                          legend.title = "Factor Type",
                          sort.by.groups = FALSE,     # Don't sort inside each group
                          x.text.angle = 90,          # Rotate vertically x axis texts
                          xlab = "DNA & RNA binding factors",
                          ylab = "Aggregate directional factor contribution",
                          rotate = TRUE,
                          lab.size = 10,
        )+
        theme_pubr(base_size = 9, legend ="top")
      
      return(top.directional.effectors.plot)
      #cowplot::save_plot(paste0(OUTPUT, PLOTS, "shap_top_total_diretional_factor_contributions_",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), top.directional.effectors.plot)
    })
  p <- 
    plot_grid(top.cluster.effectors[[1]],
              top.cluster.effectors[[2]], 
              ncol = 2,
              nrow = 1)
  
  # aggregate.cluster.contrib <- lapply(seq_along(cluster.factor.contrib), function(cluster){
  #   cluster.factor.contributions <- cluster.factor.contrib[[cluster]]
  #   all.cluster.factor.contribs <- rbind(cluster.factor.contributions$decrease, cluster.factor.contributions$increase)
  #   all.cluster.factor.contribs$cluster <- cluster
  #   all.cluster.factor.contribs$factor.contrib <- abs(all.cluster.factor.contribs$factor.contrib)
  #   all.cluster.factor.contribs <- all.cluster.factor.contribs[order(all.cluster.factor.contribs$factor.contrib,decreasing = T),]
  # })
  # 
  # top.cluster.factors <- lapply(aggregate.cluster.contrib, function(cluster.factors){head(cluster.factors$factor, n = 30)})
  # top.cluster.factors <- unique(unlist(top.cluster.factors))
  # aggregate.cluster.contrib <- do.call(rbind, aggregate.cluster.contrib)
  # aggregate.cluster.contrib <- aggregate.cluster.contrib[aggregate.cluster.contrib$factor %in% top.cluster.factors, ]
  # aggregate.cluster.contrib$cluster <- factor(aggregate.cluster.contrib$cluster)
  # 
  # ggplot(aggregate.cluster.contrib, aes(y=factor.contrib, x=factor)) +
  #   geom_bar(position="dodge", stat="identity")+
  #   scale_fill_viridis(discrete = T, end = 0.5)+
  #   coord_flip()+
  #   geom_text(aes(label=factor))+
  #   facet_grid(cols = vars(cluster))
  
  # aggregate.cluster.contrib <- lapply(seq_along(cluster.factor.contrib), function(cluster){
  #   cluster.factor.contributions <- cluster.factor.contrib[[cluster]]
  #   all.cluster.factor.contribs <- rbind(cluster.factor.contributions$decrease, cluster.factor.contributions$increase)
  #   all.cluster.factor.contribs$cluster <- cluster
  #   all.cluster.factor.contribs$factor.contrib <- abs(all.cluster.factor.contribs$factor.contrib)
  #   all.cluster.factor.contribs <- all.cluster.factor.contribs[order(all.cluster.factor.contribs$factor.contrib,decreasing = T),]
  # })
  # 
  # top.cluster.factors <- lapply(aggregate.cluster.contrib, function(cluster.factors){head(cluster.factors$factor, n = 30)})
  # top.cluster.factors <- unique(unlist(top.cluster.factors))
  # aggregate.cluster.contrib <- do.call(rbind, aggregate.cluster.contrib)
  # aggregate.cluster.contrib <- aggregate.cluster.contrib[aggregate.cluster.contrib$factor %in% top.cluster.factors, ]
  # aggregate.cluster.contrib$cluster <- factor(aggregate.cluster.contrib$cluster)
  # 
  # aggregate.cluster.contrib$factor <- factor(aggregate.cluster.contrib$factor, levels = rev(aggregate.cluster.contrib$factor[aggregate.cluster.contrib$cluster == 1])) 
  # aggregate.cluster.contrib$factor.contrib[aggregate.cluster.contrib$direction == "decrease"] <- aggregate.cluster.contrib$factor.contrib[aggregate.cluster.contrib$direction == "decrease"]*(-1)
  # 
  # ggplot(aggregate.cluster.contrib, aes(fill=cluster, y=factor.contrib, x=factor)) + 
  #   geom_bar(position="dodge", stat="identity")+
  #   scale_fill_viridis(discrete = T, end = 0.5)+
  #   coord_flip()+
  #   geom_text(aes(label=factor))
  
  
  ## Compare ranks of factors
  all.first.cluster.factor.contribs <- rbind(cluster.factor.contrib[[1]]$decrease, cluster.factor.contrib[[1]]$increase)
  all.first.cluster.factor.contribs$factor.contrib <- abs(all.first.cluster.factor.contribs$factor.contrib)
  all.first.cluster.factor.contribs <- all.first.cluster.factor.contribs[order(all.first.cluster.factor.contribs$factor.contrib, decreasing = T),]
  all.first.cluster.factor.contribs$pos <- 1:dim(all.first.cluster.factor.contribs)[1]
  
  all.second.cluster.factor.contribs <- rbind(cluster.factor.contrib[[2]]$decrease, cluster.factor.contrib[[2]]$increase)
  all.second.cluster.factor.contribs$factor.contrib <- abs(all.second.cluster.factor.contribs$factor.contrib)
  all.second.cluster.factor.contribs <- all.second.cluster.factor.contribs[order(all.second.cluster.factor.contribs$factor.contrib, decreasing = T),]
  all.second.cluster.factor.contribs$pos <- 1:dim(all.second.cluster.factor.contribs)[1]
  
  factor.pos.change <- 
    all.first.cluster.factor.contribs$pos - all.second.cluster.factor.contribs$pos[match(all.first.cluster.factor.contribs$factor,
                                                                                         all.second.cluster.factor.contribs$factor )]
  names(factor.pos.change) <- all.first.cluster.factor.contribs$factor
  factor.pos.change <- factor.pos.change[order(abs(factor.pos.change), decreasing = T)]
  all.first.cluster.factor.contribs$pos.change <- factor.pos.change[match(all.first.cluster.factor.contribs$factor, names(factor.pos.change))]
  all.first.cluster.factor.contribs$new.factor <- all.second.cluster.factor.contribs$factor
  all.first.cluster.factor.contribs$new.factor.pos.change <- factor.pos.change[match(all.first.cluster.factor.contribs$new.factor, names(factor.pos.change))]
  
  all.cluster.factor.contribs.plots <- 
    lapply(cluster.factor.contrib, function(cluster.factor.contributions){
      all.cluster.factor.contribs <- rbind(cluster.factor.contributions$decrease, cluster.factor.contributions$increase)
      all.cluster.factor.contribs$factor.contrib <- abs(all.cluster.factor.contribs$factor.contrib)
      all.cluster.factor.contribs <- all.cluster.factor.contribs[order(all.cluster.factor.contribs$factor.contrib, decreasing = T),]
      
      all.cluster.factor.contribs$factor.type <- "factor"
      all.cluster.factor.contribs$factor.type[all.cluster.factor.contribs$factor %in% dbp.factors] <- "DBP"
      all.cluster.factor.contribs$factor.type[all.cluster.factor.contribs$factor %in% rbp.factors] <- "RBP"
      all.cluster.factor.contribs$factor.type[all.cluster.factor.contribs$factor %in% rbp.factors &
                                                all.cluster.factor.contribs$factor %in% dbp.factors ] <- "DBP/RBP"
      
      known.7sk.binding.factors <-  build.pausing.associated.factor.sub.feature.spaces()
      known.7sk.binding.factors <- unique(unlist(known.7sk.binding.factors$Known))
      known.7sk.binding.factors.colors <- rev(ifelse(head(all.cluster.factor.contribs$factor, 30)%in% known.7sk.binding.factors, "red", "black"))
      
      all.cluster.factor.contribs.plot <- 
        ggpubr::ggbarplot(head(all.cluster.factor.contribs, 30), 
                          x="factor",
                          y = "factor.contrib",
                          fill = "factor.type",
                          color = "black",            # Set bar border colors to white
                          palette = "uchicago",            # jco journal color palett. see ?ggpar
                          sort.val = "asc",          # Sort the value in descending order
                          legend.title = "Factor Type",
                          sort.by.groups = FALSE,     # Don't sort inside each group
                          x.text.angle = 90,          # Rotate vertically x axis texts
                          lab.col = "black",
                          xlab = "DNA & RNA binding factors",
                          ylab = "Aggregate factor contribution",
                          rotate = TRUE,
                          lab.size = 10,
        )+
        theme_pubr(base_size = 9, legend ="top")+
        theme(axis.text.y = element_text(colour = known.7sk.binding.factors.colors))+
        scale_y_continuous(expand = expansion(mult = c(0, .1)))
      
      return(all.cluster.factor.contribs.plot)
    })
  
  all.cluster.factor.contribs.plots[[1]] <- all.cluster.factor.contribs.plots[[1]]+
    geom_text(aes(label=all.first.cluster.factor.contribs$pos.change[match(all.cluster.factor.contribs.plots[[1]]$data$factor, all.first.cluster.factor.contribs$factor)]), 
              position = position_dodge(width = 1), hjust = -0.5, size = 3)+
    ggtitle("Cluster 1 (Pause Release)")+
    theme( legend.position = "none", plot.title = element_text(size=10))
  
  all.cluster.factor.contribs.plots[[2]] <- all.cluster.factor.contribs.plots[[2]]+
    geom_text(aes(label=all.first.cluster.factor.contribs$pos.change[match(all.cluster.factor.contribs.plots[[2]]$data$factor, all.first.cluster.factor.contribs$factor)]), 
              position = position_dodge(width = 1), hjust = -1, size = 3)+
    ggtitle("Cluster 2 (Pause Enhance)")+
    theme( legend.position = "right", plot.title = element_text(size=10))
  
  
  undirectional.cluster.factor.contrib <- lapply(cluster.factor.contrib, function(cluster.factor.contributions){
    all.cluster.factor.contribs <- rbind(cluster.factor.contributions$decrease, cluster.factor.contributions$increase)
    all.cluster.factor.contribs$factor.contrib <- abs(all.cluster.factor.contribs$factor.contrib)
    all.cluster.factor.contribs <- all.cluster.factor.contribs[order(all.cluster.factor.contribs$factor.contrib, decreasing = T),]
    colnames(all.cluster.factor.contribs) <- c("factor", "contrib", "direction")
    all.cluster.factor.contribs
  })
  
  
  return(list(first.cluster.factors = top.cluster.effectors[[1]],
              second.cluster.factors = top.cluster.effectors[[2]], 
              top.cluster.effectors = top.cluster.effectors,
              #all.first.cluster.factors = all.cluster.factor.contribs[[1]], 
              #all.second.cluster.factors = all.cluster.factor.contribs[[2]],
              all.cluster.factor.contribs.plots = all.cluster.factor.contribs.plots, 
              undirectional.cluster.factor.contrib = undirectional.cluster.factor.contrib))
}
###
# This function visualizes the distribution of the target of the shap clusters
###
visualize.shap.cluster.target.distributions <- function(force.plot_data,shap.scores, X, Y, train.cline, matrix.type, feature.subspace, model.target, plot.base.size, optimal.factors = NULL){
  
  # # Temp vars to test function
  # force.plot_data = shap.force.plots$force.plot_data
  # shap.scores
  # X=X
  # train.cline = train.cline
  # matrix.type = matrix.type
  # feature.subspace = feature.subspace
  # model.target = model.target
  # plot.base.size = plot.base.size
  # optimal.factors = optimal.factors.plot$optimal.factors
  
  # append = ""
  # if(!is.null(optimal.factors)){
  #   non.binding.signals <- colnames(shap.scores)[!grepl("^chip|clip", colnames(shap.scores))]
  #   optimal.factors.pattern <- paste0(optimal.factors, collapse = "|")
  #   optimal.factor.features <- colnames(shap.scores)[grepl(optimal.factors.pattern, colnames(shap.scores))]
  #   sub.features <- c(non.binding.signals, optimal.factor.features)
  #   force.plot_data <- shap.prep.stack.data(shap_contrib = shap.scores[, sub.features], 
  #                                           top_n = 5, 
  #                                           n_groups = 2, 
  #                                           cluster_method = "ward.D")
  #   append <- "influential_factors"
  # }
  
  ### Average target per group as violinplot
  sample.clustering <- data.frame(id = force.plot_data$ID, cluster = force.plot_data$group)
  sample.clustering <- sample.clustering[order(sample.clustering$id, decreasing = F),]
  model.target.clusters <- data.frame(target = Y, cluster = as.factor(sample.clustering$cluster))
  comparison <- list( c("1", "2"))
  cluster.target.dist.plot <- 
    ggpubr::ggviolin(model.target.clusters, 
                     x = "cluster", 
                     y = "target", 
                     fill = "cluster",
                     palette =pal_uchicago("dark")(2),
                     add = "boxplot", 
                     add.params = list(fill = "white"))+
    ggpubr::stat_compare_means(comparisons = comparison, label = "p.signif")+
    ggpubr::stat_compare_means(label.y = (max(Y)+mean(Y)/2))+
    theme_pubr(base_size = plot.base.size)
  #ggtitle(paste0("Cluster target distributions ",  "(", train.cline, " / ", model.target.type, ")"))
  
  #cowplot::save_plot(paste0(OUTPUT, PLOTS, "cluster_target_distribution_", append, "_",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), cluster.target.dist.plot)
  return(cluster.target.dist.plot)
}
###
# This function visualizes deconvoluted target based on shap force clusters
###
visualize.target.deconvolution <- function(force.plot_data, Y, train.cline, matrix.type, feature.subspace, model.target, plot.base.size, optimal.factors = NULL, shap.scores = NULL, append = ""){
  
  # Temp vars to test function
  # force.plot_data = shap.force.plots$force.plot_data
  # Y = Y
  # train.cline = train.cline
  # matrix.type = matrix.type
  # feature.subspace = feature.subspace
  # model.target = model.target
  # plot.base.size = plot.base.size
  # optimal.factors = NULL
  
  # if(!is.null(optimal.factors)){
  #   non.binding.signals <- colnames(shap.scores)[!grepl("^chip|clip", colnames(shap.scores))]
  #   optimal.factors.pattern <- paste0(optimal.factors, collapse = "|")
  #   optimal.factor.features <- colnames(shap.scores)[grepl(optimal.factors.pattern, colnames(shap.scores))]
  #   sub.features <- c(non.binding.signals, optimal.factor.features)
  #   force.plot_data <- shap.prep.stack.data(shap_contrib = shap.scores[, sub.features],
  #                                           top_n = 5,
  #                                           n_groups = 2,
  #                                           cluster_method = "ward.D")
  # }
  
  ### Average target per group as violinplot
  sample.clustering <- data.frame(id = force.plot_data$ID, cluster = force.plot_data$group)
  sample.clustering <- sample.clustering[order(sample.clustering$id, decreasing = F),]
  model.target.clusters <- data.frame(target = Y, cluster = as.factor(sample.clustering$cluster))
  
  ### Plot deconvoluted signals of target
  conv.target.plot <-
    ggpubr::gghistogram(data.frame(target=Y),
                        x = "target",
                        add = "mean",
                        color = rgb(0.2,0.4,0.6,0.6),
                        binwidth = 0.05)+
    xlab(model.target.type)+
    ggtitle(paste0("Target distribution (", train.cline, " / ", model.target.type, ")"))+
    theme_pubr(base_size = plot.base.size)
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "conv_target__",append, train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"),
                     conv.target.plot)
  
  w.test.res <- 
  wilcox.test(model.target.clusters$target[model.target.clusters$cluster==1], 
              model.target.clusters$target[model.target.clusters$cluster==2])
  w.test.p <- ifelse(w.test.res$p.value==0, 2.2e-16, w.test.res$p.value)
  w.test.p <- signif(w.test.p, digits = 3)
  deconv.target.plot <-
    ggpubr::gghistogram(model.target.clusters,
                        x = "target",
                        add = "mean",
                        rug = TRUE,
                        color = "cluster",
                        fill = "cluster",
                        palette = pal_uchicago("dark")(2),
                        binwidth = 0.05)+
    xlab(model.target.type)+
    ggtitle(paste0("Deconvoluted target distribution (", train.cline, " / ", model.target.type, ")"))+
    theme_pubr(base_size = plot.base.size)+
    annotate("text", x = 12, y = 50, label = paste0("Wilcoxon, p=",w.test.p), fontface = "bold", color = "red", size = 6)
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "deconv_target__",append, train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"),
                     deconv.target.plot)
  
  ### Compare SHAP based deconvolution with gaussian mixture modeling approach
  bic <- mclustBIC(Y)
  clust <- Mclust(Y, x = bic)
  clusters <- clust$classification
  gaussian.mixture.model.clusters <- data.frame(target = Y, cluster = as.factor(clusters))
  
  gaussian.deconv.target.plot <-
    ggpubr::gghistogram(gaussian.mixture.model.clusters,
                        x = "target",
                        add = "mean",
                        rug = TRUE,
                        color = "cluster",
                        fill = "cluster",
                        palette = pal_uchicago("light")(2),
                        binwidth = 0.05)+
    xlab(model.target.type)+
    ggtitle(paste0("Gaussian Mixture Model (", train.cline, " / ", model.target.type, ")"))+
    theme_pubr(base_size = plot.base.size)
  
  deconv.approach.signals <- data.frame(shap = as.numeric(model.target.clusters$cluster),
                                        gaussian = as.numeric(gaussian.mixture.model.clusters$cluster))
  cor(deconv.approach.signals$shap, deconv.approach.signals$gaussian)
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "gaussian_mixture_model_target_deconv_",append, train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"),
                     gaussian.deconv.target.plot)
  
  
  return(list(conv.target.plot = conv.target.plot,
              deconv.target.plot = deconv.target.plot,
              gaussian.deconv.target.plot = gaussian.deconv.target.plot))
  
}
###
# This function calculates the minimal number of factors that have the highest shap contribtion
###
retrieve.optimal.number.of.factors <- function(cum.factor.contrib, train.cline, matrix.type, feature.subspace, model.target, plot.base.size, append = "", topN = 0){
  
  # Temp vars to test function
  # cum.factor.contrib = dynamic.features.plots$cum.factor.contrib
  # train.cline = train.cline
  # matrix.type = matrix.type
  # feature.subspace = feature.subspace
  # model.target = model.target
  # plot.base.size = plot.base.size
  
  
  ### Gene set enrichment analysis of top factors
  # Identify top factors (minimize number of factors, while maximizing aggregate model contributions)
  factor.contribution.distribution <- cum.factor.contrib
  factor.contribution.distribution$cumsum <- cumsum(cum.factor.contrib$contrib)
  factor.contribution.distribution$nfactors <- 1:dim(factor.contribution.distribution)[1]
  total.contrib <- sum(factor.contribution.distribution$contrib)
  factor.contribution.distribution$cumloss <- total.contrib-factor.contribution.distribution$cumsum
  
  cum.sum <- data.frame(x = factor.contribution.distribution$nfactors, y = factor.contribution.distribution$cumsum)
  cum.loss <- data.frame(x = factor.contribution.distribution$nfactors, y = factor.contribution.distribution$cumloss)
  optimum <- curve_intersect(cum.sum, cum.loss, empirical = TRUE, domain = NULL)
  optimum.num.factors <- ceiling(optimum$x)
  percent.contribution <- sum(factor.contribution.distribution$contrib[1:optimum.num.factors])/total.contrib
  
  colors <- c("cumsum" = "black", "cumloss" = "red")
  optimum.factor.plot <- 
    ggplot(factor.contribution.distribution, aes(x=nfactors)) + 
    geom_line(aes(y = cumsum, color = "cumsum")) + 
    geom_line(aes(y = cumloss, color= "cumloss"))+
    geom_vline(aes(xintercept = optimum$x), linetype="dashed", color = "darkred")+
    xlab("Number of factors")+
    #ggtitle("Optimal number of factors")+
    theme_pubr(base_size = 12)+
    theme(legend.position = c(0.9, 0.5))+
    labs(x = "Number of factors",
         y = "Factor contribution",
         color = "Legend") +
    scale_color_manual(values = colors)
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "opt_num_factors_",append, "_", train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), optimum.factor.plot)
  optimal.factors <- factor.contribution.distribution$factor[1:optimum.num.factors]
  if(topN != 0){
    optimal.factors <- factor.contribution.distribution$factor[1:topN]
  }
  return(list(optimal.factors = optimal.factors,
              percent.contribution = percent.contribution,
              optimum.factor.plot = optimum.factor.plot))
}
###
# This function visualizes the model contributions of the most influential factors
###
visualize.optimal.factor.contributions <- function(force.plot_data, shap.scores, optimal.factors, cluster.membership, X = X,train.cline, matrix.type, feature.subspace,model.target ){
  
  # sample.clustering <- data.frame(id = force.plot_data$ID, cluster = force.plot_data$group)
  # sample.clustering <- sample.clustering[order(sample.clustering$id, decreasing = F),]
  # model.target.clusters <- data.frame(target = Y, cluster = as.factor(sample.clustering$cluster))
  # 
  # cluster.feature.contrib <- shap.scores
  # cluster.feature.contrib$cluster[match(rownames(model.target.clusters), rownames(cluster.feature.contrib))] <- model.target.clusters$cluster
  # cluster.feature.contrib <- split(cluster.feature.contrib, cluster.feature.contrib$cluster)

  # force.plot_data = shap.force.plots$force.plot_data
  # shap.scores = shap.scores
  # optimal.factors = all.optimal.factors
  # cluster.membership = all.optimal.factors %in%  cluster1.optimal.factors.plot$optimal.factors
  # X = X
  
  optimal.factor.pattern <- paste0(optimal.factors, collapse = "|")
  binding.feature.contribs <- as.data.frame(shap.scores)[ , grepl("^chip|clip", colnames(shap.scores))]
  binding.patterns <-  X[ , grepl("^chip|clip", colnames(X))]
  
  optimal.factor.feature.contrib <- binding.feature.contribs[, grepl(optimal.factor.pattern, colnames(binding.feature.contribs))]
  optimal.factor.binding.pattern <- binding.patterns[, grepl(optimal.factor.pattern, colnames(binding.patterns))]
  optimal.factor.feature.contrib <- 
    lapply(optimal.factors, function(factor){
      factor.feature.contrib <- optimal.factor.feature.contrib[, grepl(factor, colnames(optimal.factor.feature.contrib))]
      factor.binding.pattern <- optimal.factor.binding.pattern[, grepl(factor, colnames(optimal.factor.binding.pattern))]
      factor.feature.contrib <- sum(abs(factor.feature.contrib))#/sum(factor.binding.pattern == 1)
      #factor.feature.contrib <- sum((factor.feature.contrib[factor.binding.pattern == 0]))
      #factor.feature.contrib <- sum(factor.feature.contrib)
    })
  optimal.factor.feature.contrib <- as.data.frame(do.call(rbind, optimal.factor.feature.contrib))
  optimal.factor.feature.contrib$factor <- optimal.factors
  colnames(optimal.factor.feature.contrib) <- c("Contribution", "Factor")
  #optimal.factor.feature.contrib$cluster = factor(revalue(factor(as.numeric(cluster.membership)), c("1" = "1", "0" = "2")), levels = c(1, 2))
  
  # Color by factor type
  dbp.factors <- 
    lapply(CHIPseq, function(peaks){
      unique(peaks$hgnc_symbol)
    }) %>% unlist %>% unique 
  
  # Retrieve all RNA binding factors
  rbp.factors <- 
    lapply(eCLIPseq, function(peaks){
      unique(peaks$hgnc_symbol)
    }) %>% unlist %>% unique 
  
  optimal.factor.feature.contrib$factor.type <- "factor"
  optimal.factor.feature.contrib$factor.type[optimal.factor.feature.contrib$Factor %in% dbp.factors] <- "DBP"
  optimal.factor.feature.contrib$factor.type[optimal.factor.feature.contrib$Factor %in% rbp.factors] <- "RBP"
  optimal.factor.feature.contrib$factor.type[optimal.factor.feature.contrib$Factor %in% rbp.factors &
                                               optimal.factor.feature.contrib$Factor %in% dbp.factors ] <- "DBP/RBP"
  
  sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
  sub.factor.feature.spaces <- 
    build.sequence.specific.biologically.functional.feature.sub.spaces(stratify = sub.factor.feature.spaces)
  known.elongation.factors <- unique(unlist(sub.factor.feature.spaces$`Elongation`))
  known.elongation.factors.colors <- (ifelse(optimal.factor.feature.contrib$Factor %in% known.elongation.factors, "red", "black"))
  names(known.elongation.factors.colors) <- optimal.factor.feature.contrib$Factor
  known.elongation.factors.colors <- known.elongation.factors.colors[ optimal.factor.feature.contrib$Factor[order(optimal.factor.feature.contrib$Contribution, decreasing = F)]]
  # 
  # all.cluster.factor.contribs.plot <- 
  #   ggpubr::ggbarplot(head(all.cluster.factor.contribs, 30), 
  #                     x="factor",
  #                     y = "factor.contrib",
  #                     fill = "factor.type",
  #                     color = "black",            # Set bar border colors to white
  #                     palette = "uchicago",            # jco journal color palett. see ?ggpar
  #                     sort.val = "asc",          # Sort the value in descending order
  #                     legend.title = "Factor Type",
  #                     sort.by.groups = FALSE,     # Don't sort inside each group
  #                     x.text.angle = 90,          # Rotate vertically x axis texts
  #                     lab.col = "black",
  #                     xlab = "DNA & RNA binding factors",
  #                     ylab = "Aggregate factor contribution",
  #                     rotate = TRUE,
  #                     lab.size = 10,
  #   )+
  #   theme_pubr(base_size = 9, legend ="top")+
  # 
  #   scale_y_continuous(expand = expansion(mult = c(0, .1)))
  
  
  optimal.factor.contrib.plot <- 
    ggpubr::ggbarplot(optimal.factor.feature.contrib, 
                      x="Factor",
                      y = "Contribution",
                      fill = "factor.type",            # Set bar border colors to white
                      palette = "uchicago",            # jco journal color palett. see ?ggpar
                      sort.val = "asc",          # Sort the value in descending order
                      sort.by.groups = FALSE,     # Don't sort inside each group
                      #x.text.angle = 90,          # Rotate vertically x axis texts
                      xlab = "Most influential factors",
                      ylab = "Total factor contribution",
                      rotate = TRUE,
                      ggtheme = theme_pubr(base_size = 10),
                      lab.col = "black",
                      color = "black"
                      # main = "Factor model contributions"
    )+
    theme(axis.text.y = element_text(colour = known.elongation.factors.colors))
  #cowplot::save_plot(paste0(OUTPUT, PLOTS, "shap_top_total_influential_factor_contributions_",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), optimal.factor.contrib.plot)
  
  optimal.factor.feature.contrib <- optimal.factor.feature.contrib[order(optimal.factor.feature.contrib$Contribution, decreasing = F), ]
  #top.directional.binding.factors <- c(head(optimal.factor.feature.contrib$Factor, n = 5), tail(optimal.factor.feature.contrib$Factor, n = 5))
  return(list(optimal.factor.contrib.plot = optimal.factor.contrib.plot))
#   return(list(optimal.factor.contrib.plot = optimal.factor.contrib.plot,
#               top.directional.binding.factors = top.directional.binding.factors))
}
###
# This function performs a gene set enrichment analysis
###
perform.gene.set.enrichment.analysis <- function(universe, node.size = 5){
  
  # # Check if a precalculated version exists
  # target.file <- paste(OUTPUT, "optimal.factor.gsea.RDS", sep = "")
  # if(!(CARGS$new) & file.exists(target.file)){
  #   optimal.factor.gsea.res <- readRDS(target.file)
  #   return(optimal.factor.gsea.res)
  # }
  
  GOdata <- new("topGOdata",
                description = "Model Feature GSEA", 
                ontology = "BP", 
                allGenes = universe, 
                nodeSize = node.size, 
                annot = annFUN.org,
                mapping="org.Hs.eg.db",
                ID="SYMBOL")
  
  resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  allRes <- GenTable(GOdata, 
                     classicFisher = resultFisher,
                     orderBy = "resultFisher", 
                     ranksOf = "classicFisher")
  allRes$adj.p <- p.adjust(allRes$classicFisher, method = "fdr")
  optimal.factor.gsea.res <- allRes[, c("GO.ID","Term", "classicFisher", "adj.p")]
  colnames(optimal.factor.gsea.res) <- c("GO.ID","Term","p","adj.p")
  saveRDS(optimal.factor.gsea.res, paste0(OUTPUT, "optimal.factor.gsea.RDS"))
  
  
  optimal.factor.gsea.res$p <- as.numeric(optimal.factor.gsea.res$p)
  optimal.factor.gsea.res$adj.p <- as.numeric(optimal.factor.gsea.res$adj.p)
  optimal.factor.gsea.res$p.rescaled <- -log10(optimal.factor.gsea.res$p)
  optimal.factor.gsea.res$adj.p.rescaled <- -log10(optimal.factor.gsea.res$adj.p)
  #optimal.factor.gsea.res <- optimal.factor.gsea.res[order(optimal.factor.gsea.res$p.rescaled, decreasing = F),]
  #optimal.factor.gsea.res$adj.p <- format(round(optimal.factor.gsea.res$adj.p, 3), nsmall = 3)
  optimal.factor.gsea.res <- optimal.factor.gsea.res[order(optimal.factor.gsea.res$p, decreasing = F),]
  
  gsea.plot <- 
    ggbarplot(optimal.factor.gsea.res, 
              x = "Term", 
              y = "p.rescaled",
              fill = "p.rescaled", 
              #label = paste0(paste0("p=", optimal.factor.gsea.res$p), " / ", paste0("adj.p=", optimal.factor.gsea.res$adj.p) ),
              rotate = T,
              lab.pos = "in",
              lab.col = "white",
              lab.size = 4,
              lab.hjust = 1.1,
              lab.vjust = 1,
              position = position_dodge(0.9))+
    xlab("GO Term")+
    ylab("-log10(p-value)")+
    scale_fill_gsea(reverse = F)+
    theme_pubr(base_size = plot.base.size, legend = "right")
  #geom_hline(yintercept = min(optimal.factor.gsea.res$adj.p.rescaled), colour = "red")
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "optimal_factor_gsea_",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"),
                     gsea.plot)
  
  return(list(optimal.factor.gsea = optimal.factor.gsea.res,
              gsea.plot = gsea.plot))
}
###
# This function visualizes the correlation of the traveling ratio and gene expression targets
###
visualize.target.correlations <- function(target.model.results, train.cline, matrix.type, feature.subspace, model.target, plot.base.size){
  
  targets <- data.frame(traveling.ratio = target.model.results$target_traveling.ratio$Y,
                        transcript.expression = target.model.results$target_transcript.expression$Y)
  
  # tryCatch({dev.off()}, error = function(e) {return(NULL)})
  # heatscatter(targets$traveling.ratio,
  #             targets$transcript.expression,
  #             xlab = "Traveling Ratio",
  #             ylab = "Transcript Expression",
  #             colpal="bl2gr2rd",
  #             main="Traveling Ratio vs. Transcript Expression",
  #             cor=T, 
  #             add.contour = T, 
  #             nrcol = 100,
  #             alpha = 25,
  #             method = "pearson",
  #             nlevels = 15)
  # abline(v = mean(targets$traveling.ratio), col="red", lwd=2, lty=2)
  # abline(h = mean(targets$transcript.expression), col="red", lwd=2, lty=2)
  # text(x=12, y=3, labels = paste("rho = ", format(round(cor(targets$traveling.ratio, targets$transcript.expression), 2))), font = 2, cex = 1.5)
  # p <- recordPlot()
  # return(p)
  
  tr.vs.exp <- 
    ggscatter(targets, 
              x = "traveling.ratio", 
              y = "transcript.expression",
              #color = "resid",
              add ="reg.line",
              conf.int = TRUE,
              add.params = list(color = "red",
                                fill = "lightgray",
                                size = 0.3),
              #palette = pal_uchicago("dark")(9)[9],
              alpha = 1/10
    ) +
    theme_pubr(base_size = plot.base.size)+
    #ggtitle("Traveling Ratio vs. Transcript Expression")+
    xlab("Traveling Ratio")+
    ylab("Transcript Expression")+
    ggpubr::stat_cor(cor.coef.name = c("rho"))+
    gradient_color(c("black", pal_uchicago("light")(9)[1]))
  return(tr.vs.exp)
}
###
# This function compares the feature contributions to a reverse target models feature contributions
###
visualize.reverse.model.feature.contributions <- function(optimal.factors, forward.model.directional.factor.contrib, model.results, X, train.cline, matrix.type, feature.subspace, model.target, plot.base.size){
  
  ### Inverse correlation of feature contributions with contributions from a reverse target model
  reverse.target <- if(model.target=="target_traveling.ratio") "target_transcript.expression" else "target_traveling.ratio"
  # Build the reverse target
  reverse.model.type <- gsub("target_", "", reverse.target)
  # Retrieve model results from the reverse target model
  reverse.target.model.results <- model.results$models[[train.cline]][[matrix.type]][[feature.subspace]][[reverse.target]]
  # Rretrieve shap values of reverse target model
  reverse.target.model.shap.scores <- as.data.frame(reverse.target.model.results$shap.values$shap_score)
  
  ### Feature class (DBP, RBP etc.) ranking with absolute shap values
  # Retrieve indices of binary binding features of factor bindings
  binding.signals <- colnames(reverse.target.model.shap.scores)[grepl("^chip|clip", colnames(reverse.target.model.shap.scores))]
  factor.binding.contrib <- reverse.target.model.shap.scores[ , binding.signals]
  # Make shap values absolute to make contributions undirectional
  absolute.feature.contrib <- as.data.frame(apply(factor.binding.contrib, 2, abs))
  # Filter for non-zero contributions
  aggregate.feature.contrib <- colSums(absolute.feature.contrib)
  aggregate.feature.contrib <- aggregate.feature.contrib[aggregate.feature.contrib>0]
  
  factor.contrib <- sort(aggregate.feature.contrib, decreasing = T)
  factor.contrib <- data.frame(contrib = factor.contrib, factor.binding = names(factor.contrib))
  
  contributing.factors <- factor.contrib$factor.binding
  contributing.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*", 
                               "", contributing.factors)
  contributing.factors <- unique(contributing.factors)
  
  directional.factor.contrib <- 
    lapply(contributing.factors, function(factor){
      factor.binding.patterns <- X[, grepl(factor, colnames(X))]
      factor.features <- colnames(X)[grepl(factor, colnames(X))]
      factor.feature.contrib <- 
        lapply(factor.features, function(factor.binding.feature){
          bound.samples <- which(factor.binding.patterns[, factor.binding.feature]==1)
          sum(factor.binding.contrib[bound.samples, factor.binding.feature])
        })
      factor.feature.contrib <- sum(unlist(factor.feature.contrib))
      data.frame(factor = factor, factor.contrib = factor.feature.contrib, direction = if (factor.feature.contrib<0) "decrease" else "increase" )
    })
  directional.factor.contrib <- do.call(rbind, directional.factor.contrib)
  directional.factor.contrib <- directional.factor.contrib[order(abs(directional.factor.contrib$factor.contrib), decreasing = F),]
  directional.factor.contrib$direction <- as.factor(directional.factor.contrib$direction)
  
  reverse.model.directional.factor.contrib <- directional.factor.contrib
  forward.model.directional.factor.contrib <- do.call(rbind, forward.model.directional.factor.contrib)
  matching.factors <- intersect(reverse.model.directional.factor.contrib$factor,
                                forward.model.directional.factor.contrib$factor)
  
  reverse.model.directional.factor.contrib <- reverse.model.directional.factor.contrib[match(matching.factors, reverse.model.directional.factor.contrib$factor), ]
  forward.model.directional.factor.contrib <- forward.model.directional.factor.contrib[match(matching.factors, forward.model.directional.factor.contrib$factor), ]
  
  rev.forw.model.factor.contrib.cor <- 
    cor(forward.model.directional.factor.contrib$factor.contrib,
        reverse.model.directional.factor.contrib$factor.contrib)
  
  df <- data.frame(forward.model.features = forward.model.directional.factor.contrib$factor.contrib, 
                   reverse.model.features = reverse.model.directional.factor.contrib$factor.contrib,
                   factor = matching.factors, 
                   model.contribution = as.factor(revalue(as.character(matching.factors %in% optimal.factors), c("TRUE" = "high", "FALSE" = "low"))))
  
  inverse.model.feature.correlations <- 
    ggpubr::ggscatter(df, 
                      x = "forward.model.features", 
                      y = "reverse.model.features",
                      color = "model.contribution",
                      palette = "uchicago",
                      add = "reg.line",
                      add.params = list(color = "blue",
                                        fill =  "lightgray"),
                      rug = T,
                      alpha = 1/2, 
                      label = "factor",
                      conf.int = T,
                      label.select = optimal.factors,
                      repel = TRUE,
                      font.label = c(12, "bold"),
    )+
    ggpubr::stat_cor(cor.coef.name = c("rho"), 
                     method = "pearson",
                     label.x=(max(forward.model.directional.factor.contrib$factor.contrib)-abs(min(forward.model.directional.factor.contrib$factor.contrib)))/2)+
    xlab(paste0("Forward model factor contributions \n"," (", model.target.type, " model)"))+
    ylab(paste0("Reverse model factor contributions \n"," (", reverse.model.type, " model)"))+
    theme_pubr( base_size = 12)
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "inverse_model_feature_correlations_",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"),
                     inverse.model.feature.correlations)
  
  return(list(feature.correlations = inverse.model.feature.correlations, 
              reverse.mdoel = reverse.target.model.results))
}
###
# This function compares the correlations of simulated and real knockdown expression profiles
###
visualize.knockdown.correlations <- function(optimal.factors, model.results){
  # Identify factors with KO expression data
  optimal.ko.factors <- optimal.factors[optimal.factors %in% names(KO.GENE.QUANT[[train.cline]])]
  ## Retrieve predictions 
  shap.forces <- model.results$models[[train.cline]][[matrix.type]][[feature.subspace]][["target_transcript.expression"]]$shap.values$shap_score
  shap.bias <- as.numeric(model.results$models[[train.cline]][[matrix.type]][[feature.subspace]][["target_transcript.expression"]]$shap.values$BIAS0)
  shap.forces <- as.data.frame(shap.forces)
  # Subset by only binding signals
  shap.forces <- shap.forces[, colnames(shap.forces)[grepl("^chip|clip", colnames(shap.forces))]]
  # Retrieve unperturbed expression profile
  expression.target <- model.results$models[[train.cline]][[matrix.type]][[feature.subspace]][["target_transcript.expression"]]$Y
  
  expression.profile.comparisons <- 
    lapply(optimal.ko.factors, function(factor){
      # Retrieve knockdown profile for optimal factor
      factor.ko.profile <- KO.GENE.QUANT[[train.cline]][[factor]]
      # Remove factor's contributions
      ko.shap.forces <- shap.forces[ , which(!grepl(factor, colnames(shap.forces)))]
      model.predictions <- rowSums(ko.shap.forces)+shap.bias
      # Retrieve common samples
      common.samples <- intersect(rownames(expression.target), factor.ko.profile$transcript_id)
      # Subset data sets by common samples
      common.factor.ko.profile <- factor.ko.profile[match(common.samples, factor.ko.profile$transcript_id), ]
      common.expression.target <- expression.target[match(common.samples, rownames(expression.target)), ]
      common.model.predictions <- model.predictions[match(common.samples, names(model.predictions))]
      # Aggregate expression profiles
      expression.profiles <- 
        data.frame(Exp = common.expression.target, 
                   Sim.KO = common.model.predictions,
                   KO = common.factor.ko.profile$log10_FPKM)
      
      # Calculate pairwise correlation of expression profiles
      expression.correlation.profiles <- cor(expression.profiles)
      profile.pairs <- combn(c("Exp", "Sim.KO", "KO"), 2)
      expression.correlation.profiles <- 
        apply(profile.pairs, 2, function(profile.pair){
          expression.correlation.profiles[profile.pair[1], profile.pair[2]]
        })
      names(expression.correlation.profiles) <- apply(profile.pairs, 2, paste0, collapse = " / ")
      expression.correlation.profiles <- as.data.frame(expression.correlation.profiles)
      expression.correlation.profiles$profile.pair <- rownames(expression.correlation.profiles)
      expression.correlation.profiles$factor <- factor
      expression.correlation.profiles
    })
  expression.profile.comparisons <- do.call(rbind, expression.profile.comparisons)
  expression.profile.comparisons$expression.correlation.profiles <- format(round(expression.profile.comparisons$expression.correlation.profiles, 2), nsmall = 2)
  expression.profile.comparisons$expression.correlation.profiles <- as.numeric(expression.profile.comparisons$expression.correlation.profiles)
  expression.profile.comparisons <- expression.profile.comparisons[order(expression.profile.comparisons$expression.correlation.profiles, decreasing = T),]
  expression.profile.comparisons$factor <- factor(expression.profile.comparisons$factor, levels=optimal.ko.factors)
  
  knockdown.correlation.plot <- 
    ggbarplot(expression.profile.comparisons, 
              x = "profile.pair", 
              y = "expression.correlation.profiles",
              fill = "profile.pair", 
              color = "profile.pair", 
              label = TRUE,
              palette = "uchicago",
              position = position_dodge(0.9))+
    facet_grid(~factor)+
    rremove("x.text")+
    rremove("x.ticks")+
    xlab("Profile Pair")+
    ylab("Profile Correlation")
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "knockout_profiles__",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"),
                     knockdown.correlation.plot)
  
  return(knockdown.correlation.plot)
  
}
###
# This function visualizes the model performances
###
visualize.go.term.model.results <- function(model.results, train.cline, matrix.type, model.target, plot.base.size, random.model.results){
  
  # Temp vars to test function
  # model.results = model.results
  # train.cline = train.cline
  # matrix.type = matrix.type
  # model.target = model.target
  # plot.base.size = plot.base.size
  
  # Retrieve model performances
  model.performances <- model.results$performances
  # Remove models with no factor binding features
  model.performances <- model.performances[model.performances$nfactors > 0,]
  
  # Subset by go term associated model performances
  sub.model.performances <- model.performances[model.performances$cline==train.cline &
                                                 model.performances$modeltype==matrix.type &
                                                 model.performances$target==model.target, ]
  # Transformations for ordering plot
  sub.model.performances$sequence.specific <- grepl(":ss|:nss", sub.model.performances$subspace)
  sub.model.performances$sequence.specific[!sub.model.performances$sequence.specific] <- "General"
  sub.model.performances$sequence.specific[sub.model.performances$sequence.specific == "TRUE"] <- 
    grepl(":ss", sub.model.performances$subspace[sub.model.performances$sequence.specific == "TRUE"])
  sub.model.performances$sequence.specific[sub.model.performances$sequence.specific=="TRUE"] <- "sequence.specific"
  sub.model.performances$sequence.specific[sub.model.performances$sequence.specific=="FALSE"] <- "non-sequence.specific"
  sub.model.performances$sequence.specific <- factor(sub.model.performances$sequence.specific)
  
  # Transformations for ordering plot
  # sub.model.performances$sequence.specific <- factor(grepl(":ss", sub.model.performances$subspace))
  # sub.model.performances$sequence.specific <- revalue(sub.model.performances$sequence.specific, c("TRUE"="sequence.specific", 
  #                                                                                                 "FALSE"="non-sequence.specific"))
  feature.types <- c("Chromatin", "Initiation",  "Elongation","7SK.Binding", "Elongation+7SK", "Termination","Splicing", "Processing",  "All", "Random")
  feature.types <- rev(feature.types)
  sub.model.performances$subspace <- gsub(":ss|:nss", "", sub.model.performances$subspace)
  sub.model.performances <- sub.model.performances[sub.model.performances$subspace %in% feature.types,]
  
  ## Append random model results
  random.results <- data.frame(matrix(vector(), ncol = dim(random.model.results)[2]))
  colnames(random.results) <- colnames(random.model.results)
  random.results <- rbind(random.results, random.model.results[1, ])
  random.results$train.rsqrd <- mean(random.model.results$train.rsqrd)
  random.results$test.rsqrd <- mean(random.model.results$test.rsqrd)
  random.results$nfeatures <- mean(random.model.results$nfeatures)
  random.results$nfactors <- mean(random.model.results$nfactors)
  random.results$train.power <- mean(random.model.results$train.power)
  random.results$test.power <- mean(random.model.results$test.power)
  random.results$mean.shap <- mean(random.model.results$mean.shap)
  random.results$sequence.specific <- "General"
  random.results$subspace <- factor("Random")
  
  target.cols <- c("subspace", "train.rsqrd", "test.rsqrd", "nfeatures", "nfactors", "sequence.specific", "mean.shap")
  sub.model.performances <- rbind(sub.model.performances[,target.cols], random.results[, target.cols])
  
  # Reverse order of models
  sub.model.performances$subspace <- factor(sub.model.performances$subspace, levels = feature.types)
  
  
  #sub.model.performances$test.rsqrd <- format(round(sub.model.performances$test.rsqrd, 3), nsmall=3)
  sub.model.performances$test.rsqrd <- as.numeric(sub.model.performances$test.rsqrd)
  #sub.model.performances <- sub.model.performances[order(sub.model.performances$mean.shap, sub.model.performances$test.rsqrd, decreasing = F), ]
  sub.model.performances <- sub.model.performances[order(sub.model.performances$test.rsqrd, decreasing = F), ]
  #sub.model.performances$mean.shap <- format(round(sub.model.performances$mean.shap, 3), nsmall=3)
  
  # Plot model results
  model.performance.plot <- 
    ggbarplot(sub.model.performances, 
              x = "subspace", 
              y = "test.rsqrd",
              #fill = rgb(0.2,0.4,0.6,0.8), 
              fill = "#69b3a2",
              #label = paste0(paste0("Rsqrd. ", sub.model.performances$test.rsqrd), "\n#Factors ", sub.model.performances$nfactors, "\nContrib. ", format(round(sub.model.performances$mean.shap, 3), nsmall=3) ),
              label = paste0(sub.model.performances$nfactors),
              rotate = T,
              #lab.pos = "in",
              lab.col = "black",
              lab.font = "bold",
              lab.size = 4,
              lab.hjust = -0.5,
              lab.vjust = 0.5,
              position = position_dodge(0.9))+
    facet_grid(~sequence.specific)+
    xlab("Model")+
    ylab("Test Data R-squared")+
    theme_pubr(base_size = plot.base.size, legend = "right")+
    #scale_fill_material("grey")+
    scale_fill_material("blue-grey")+
    #scale_fill_uchicago()+
    theme(panel.spacing = unit(2, "lines"))+
    scale_y_continuous(limits = c(0, 0.8), breaks = seq(0,0.8, 0.1))
  
  # cowplot::save_plot(paste0(OUTPUT, PLOTS, "sequence_specific_go_term_model_performances",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"),
  #                    model.performance.plot)
  
  return(list(model.performance.plot = model.performance.plot,
              sub.model.performances = sub.model.performances))
}
###
# This function visualizes the model performances
###
visualize.custom.model.results <- function(model.results, train.cline, matrix.type, model.target, plot.base.size){
  
  # model.results = model.results
  # train.cline = train.cline
  # matrix.type = matrix.type
  # model.target = model.target
  # plot.base.size = plot.base.size
  
  # Retrieve model performances
  model.performances <- model.results$performances
  # Remove models with no factor binding features
  model.performances <- model.performances[model.performances$nfactors>0,]
  # Define go term models
  feature.types <- c("chromatin", "initiation",  "elongation", "termination", "translation", "splicing", "processing", "pos_regulation", "neg_regulation", "sdom")
  feature.types <- paste0(feature.types, collapse ="|")
  # Subset by non go term associated model performances
  sub.model.performances <- model.performances[model.performances$cline==train.cline &
                                                 model.performances$modeltype==matrix.type &
                                                 model.performances$target==model.target &
                                                 !grepl(feature.types, model.performances$subspace), ]
  
  sub.model.performances <- sub.model.performances[!grepl("^RN7SK", sub.model.performances$subspace),]
  
  # Transformations for ordering plot
  sub.model.performances$sequence.specific <- grepl(":ss|:nss", sub.model.performances$subspace)
  sub.model.performances$sequence.specific[!sub.model.performances$sequence.specific] <- "General"
  sub.model.performances$sequence.specific[sub.model.performances$sequence.specific == "TRUE"] <- 
    grepl(":ss", sub.model.performances$subspace[sub.model.performances$sequence.specific == "TRUE"])
  sub.model.performances$sequence.specific[sub.model.performances$sequence.specific=="TRUE"] <- "sequence.specific"
  sub.model.performances$sequence.specific[sub.model.performances$sequence.specific=="FALSE"] <- "non-sequence.specific"
  sub.model.performances$sequence.specific <- factor(sub.model.performances$sequence.specific)
  sub.model.performances$subspace <- gsub(":ss|:nss", "", sub.model.performances$subspace)
  
  #sub.model.performances$test.rsqrd <- format(round(sub.model.performances$test.rsqrd, 3), nsmall=3)
  sub.model.performances$test.rsqrd <- as.numeric(sub.model.performances$test.rsqrd )
  #sub.model.performances <- sub.model.performances[order( sub.model.performances$test.rsqrd,sub.model.performances$nfactors, decreasing = T), ]
  sub.model.performances <- sub.model.performances[order(sub.model.performances$mean.shap, sub.model.performances$test.rsqrd, decreasing = F), ]
  #sub.model.performances$mean.shap <- format(round(sub.model.performances$mean.shap, 6), nsmall=6)
  
  
  # temp.general <- colnames(model.matrices$K562$individual.model.matrices$`all.RN7SK.binding`)
  # temp.ss <- colnames(model.matrices$K562$individual.model.matrices$`all.RN7SK.binding:ss`)
  # temp.nss <- colnames(model.matrices$K562$individual.model.matrices$`all.RN7SK.binding:nss`)
  # 
  # temp.general <- temp.general[grepl("^chip|clip", temp.general)]
  # temp.ss <- temp.ss[grepl("^chip|clip", temp.ss)]
  # temp.nss <- temp.nss[grepl("^chip|clip", temp.nss)]
  # 
  # temp.general <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*", 
  #                      "", temp.general)
  # temp.ss <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*", 
  #                      "", temp.ss)
  # temp.nss <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*", 
  #                      "", temp.nss)
  # temp.general <- unique(temp.general)
  # temp.ss <- unique(temp.ss)
  # temp.nss <- unique(temp.nss)
  
  # Plot model results
  model.performance.plot <- 
    ggbarplot(sub.model.performances, 
              x = "subspace", 
              y = "test.rsqrd",
              fill = "mean.shap", 
              #label = paste0(paste0("Rsqrd. ", sub.model.performances$test.rsqrd), "\n#Factors ", sub.model.performances$nfactors, "\nContrib. ", format(round(sub.model.performances$mean.shap, 3), nsmall=3) ),
              label = paste0(sub.model.performances$nfactors),
              rotate = T,
              #lab.pos = "in",
              lab.col = "black",
              lab.font = "bold",
              lab.size = 4,
              lab.hjust = -0.5,
              lab.vjust = 0.5,
              position = position_dodge(0.9))+
    facet_grid(~sequence.specific)+
    xlab("Model")+
    ylab("Test Data R-squared")+
    theme_pubr(base_size = plot.base.size, legend = "right")+
    scale_fill_material("grey")+
    theme(panel.spacing = unit(2, "lines"))+
    scale_y_continuous(limits = c(0, 0.8), breaks = seq(0,0.8, 0.1))
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "sequence_specific_model_performances",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"),
                     model.performance.plot)
  
  return(list(model.performance.plot = model.performance.plot,
              sub.model.performances = sub.model.performances))
}
### 
# This function performs gene set enrichment analyses of feature sub spaces
###
feature.sub.space.gsea <- function(model.results, train.cline, matrix.type, feature.subspace, model.target){
  
  # Check if a precalculated version exists
  # target.file <- paste0(OUTPUT, "feature.sub.space.gsea.RDS")
  # if(!(CARGS$new) & file.exists(target.file)){
  #   feature.subpace.gsea.results <- readRDS( target.file)
  #   return(feature.subpace.gsea.results)
  # }
  
  # X TEMP vars to test function
  # model.target="target_traveling.ratio"
  # matrix.type <- "individual.model.matrices"
  # feature.subspace <- "All"
  # model.target.type <- gsub("target_", "", model.target)
  # train.cline <- "K562"
  
  ## Retrieve all available factors from all assays
  clip.targets <- lapply(eCLIPseq, function(peaks){
    as.character(unique(peaks$hgnc_symbol))
  })
  chip.targets <- lapply(CHIPseq, function(peaks){
    as.character(unique(peaks$hgnc_symbol))
  })
  all.available.factors <- unique(c(unlist(clip.targets), unlist(chip.targets)))
  
  # Define feature sets to perform gsea on
  go.term.models <- c("chromatin", "initiation",  "elongation", "termination", "translation", "splicing", "processing", "pos_regulation", "neg_regulation")
  go.term.models <- paste0(go.term.models, collapse = "|")
  relevant.feature.subspaces <- unique(model.results$performances$subspace[!grepl(go.term.models, model.results$performances$subspace)])
  
  feature.subpace.gsea.results <- 
    lapply(relevant.feature.subspaces, function(feature.subspace){
      # Load relevant model results
      select.model.results <- model.results$models[[train.cline]][[matrix.type]][[feature.subspace]][[model.target]]
      
      # Extract relevant variables from model results
      X <- select.model.results$X
      Y <- select.model.results$Y
      X.test = select.model.results$X.test
      Y.test <- select.model.results$Y.test
      train.preds <- select.model.results$train.preds
      test.preds <- select.model.results$test.preds
      shap.values <- select.model.results$shap.values
      samples <- select.model.results$samples
      model.target.type <- gsub("target_", "", model.target)
      # Identify test cell line
      test.cline <- if(train.cline=="K562") "HepG2" else "K562"
      
      mean.shap.values <- shap.values$mean_shap_score
      mean.shap.values <- mean.shap.values[mean.shap.values>0]
      contributing.factors <- names(mean.shap.values)[grepl("^chip|clip", names(mean.shap.values))]
      
      # Perform gene set enrichment analysis with top scoring factors
      contributing.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*", 
                                   "", contributing.factors) %>% unique()
      
      if(length(contributing.factors) == 0){return()}
      
      universe <- factor(as.integer(all.available.factors %in% contributing.factors))
      names(universe) <- all.available.factors
      GOdata <- new("topGOdata",
                    description = "Model Feature GSEA", 
                    ontology = "BP", 
                    allGenes = universe, 
                    nodeSize = 50, 
                    annot = annFUN.org,
                    mapping="org.Hs.eg.db",
                    ID="SYMBOL")
      
      resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
      allRes <- GenTable(GOdata, 
                         classicFisher = resultFisher,
                         orderBy = "resultFisher", 
                         ranksOf = "classicFisher")
      #allRes[allRes$classicFisher < 0.05, ]
      allRes$adj.p <- p.adjust(allRes$classicFisher, method = "fdr")
      allRes <- allRes[allRes$adj.p < 0.05,]
      if(dim(allRes)[1] == 0){return()}
      allRes$Subspace <- feature.subspace
      # Retrieve top signfificant GO term associations
      allRes <- allRes[ , c("Subspace", "GO.ID", "Term", "classicFisher", "adj.p", "Annotated", "Significant", "Expected" )]
      colnames(allRes) <- c("Subspace", "GO.ID", "Term", "p", "adj.p", "Annotated", "Significant", "Expected" )
      allRes[1, ]
    })
  feature.subpace.gsea.results <- do.call(rbind, feature.subpace.gsea.results)
  
  # Format p-values to numeric type
  feature.subpace.gsea.results <- feature.subpace.gsea.results[!grepl("sdom", feature.subpace.gsea.results$Subspace),]
  feature.subpace.gsea.results$p <- as.numeric(feature.subpace.gsea.results$p)
  feature.subpace.gsea.results$adj.p <- -log10(as.numeric(feature.subpace.gsea.results$adj.p))
  feature.subpace.gsea.results <- feature.subpace.gsea.results[order(feature.subpace.gsea.results$adj.p, decreasing = F),]
  # Reformat adj.pavlues notation
  #x feature.subpace.gsea.results$adj.p <- format(feature.subpace.gsea.results$adj.p, digits = 2, scientific = T)
  # Sort in creasing order of pvalue
  feature.subpace.gsea.results <- feature.subpace.gsea.results[order(feature.subpace.gsea.results$Term, feature.subpace.gsea.results$adj.p, decreasing = F), ]
  # Format p values to character for markdown document
  #feature.subpace.gsea.results$p <- as.character(feature.subpace.gsea.results$p)
  #feature.subpace.gsea.results$adj.p <- as.character(feature.subpace.gsea.results$adj.p)
  saveRDS(feature.subpace.gsea.results, paste0(OUTPUT, "feature.sub.space.gsea.RDS"))
  
  
  
  gsea.plot <- 
    ggbarplot(feature.subpace.gsea.results, 
              x = "Subspace", 
              y = "adj.p",
              fill = "adj.p", 
              label = paste0(feature.subpace.gsea.results$Term),
              rotate = T,
              lab.pos = "out",
              lab.col = "black",
              lab.size = 3.9,
              lab.hjust = 1,
              lab.vjust = -4,
              position = position_dodge(0.9))+
    xlab("Subpace's best GO Term Enrichment")+
    ylab("-log10(adj. p-value)")+
    scale_fill_gsea(reverse = TRUE)+
    theme_pubr(base_size = plot.base.size, legend = "right")
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "custom_model_gsea_",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"),
                     gsea.plot)
  
  
  # # Transformations for ordering plot
  # feature.subpace.gsea.results <- feature.subpace.gsea.results[feature.subpace.gsea.results$adj.p<0.05,]
  # feature.subpace.gsea.results$sequence.specific <- grepl(":ss|:nss", feature.subpace.gsea.results$Subspace)
  # feature.subpace.gsea.results$sequence.specific[!feature.subpace.gsea.results$sequence.specific] <- "General"
  # feature.subpace.gsea.results$sequence.specific[feature.subpace.gsea.results$sequence.specific == "TRUE"] <-
  #   grepl(":ss", feature.subpace.gsea.results$Subspace[feature.subpace.gsea.results$sequence.specific == "TRUE"])
  # feature.subpace.gsea.results$sequence.specific[feature.subpace.gsea.results$sequence.specific=="TRUE"] <- "sequence.specific"
  # feature.subpace.gsea.results$sequence.specific[feature.subpace.gsea.results$sequence.specific=="FALSE"] <- "non-sequence.specific"
  # feature.subpace.gsea.results$sequence.specific <- factor(feature.subpace.gsea.results$sequence.specific)
  # feature.subpace.gsea.results$Subspace <- gsub(":ss|:nss", "", feature.subpace.gsea.results$Subspace)
  # 
  # feature.subpace.gsea.results.splitted <- split(feature.subpace.gsea.results, feature.subpace.gsea.results$sequence.specific)
  # gsea.heatmaps <- 
  # lapply(names(feature.subpace.gsea.results.splitted), function(enrichment.type){
  #   enrichments <- feature.subpace.gsea.results.splitted[[enrichment.type]]
  #   gsea.heatmap <- matrix(0,
  #                          nrow = length(unique(enrichments$Term)),
  #                          ncol = length(unique(enrichments$Subspace)))
  #   rownames(gsea.heatmap) <- unique(enrichments$Term)
  #   colnames(gsea.heatmap) <- unique(enrichments$Subspace)
  # 
  #   for(i in 1:dim(enrichments)[1]){
  #     subspace <- enrichments[i, "Subspace"]
  #     term <- enrichments[i, "Term"]
  #     p <- enrichments[i, "adj.p"]
  #     gsea.heatmap[term, subspace] <- -log10(p)
  #   }
  #   gsea.heatmap.plot <- 
  #   ggheatmap(gsea.heatmap, scale = "none", Colv = F, Rowv = F)
  # 
  # })
  # names(gsea.heatmaps) <- names(feature.subpace.gsea.results.splitted)
  # 
  # cowplot::save_plot(paste0(OUTPUT, PLOTS, "optimal_factor_gsea_",train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"),
  #                    gsea.plot)
  
  return(list(feature.subpace.gsea.results = feature.subpace.gsea.results, 
              feature.subpace.gsea.plot = gsea.plot))
}
### 
# This function trains a model based on all factors involved in RNA processing
###
train.minimal.model <- function(all.optimal.factors, train.cline, matrix.type, feature.subspace){
  
  # # Build a model with splicing and RN7SK factors
  # sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
  # sub.factor.feature.spaces <- 
  #   build.sequence.specific.biologically.functional.feature.sub.spaces(stratify = sub.factor.feature.spaces)
  # #sub.factor.feature.spaces <- append(sub.factor.feature.spaces, sequence.specific.biologically.functional.feature.sub.spaces)
  # 
  # select.features <- sub.factor.feature.spaces$chromatin
  # select.features <- combine.lists(select.features, sub.factor.feature.spaces$initiation)
  # select.features <- combine.lists(select.features, sub.factor.feature.spaces$elongation)
  # select.features <- combine.lists(select.features, sub.factor.feature.spaces$termination)
  # select.features <- combine.lists(select.features, sub.factor.feature.spaces$splicing)
  # select.features <- combine.lists(select.features, sub.factor.feature.spaces$translation)
  # select.features <- combine.lists(select.features, sub.factor.feature.spaces$pausing)
  # select.features <- combine.lists(select.features, sub.factor.feature.spaces$processing)
  
  select.features <- all.optimal.factors
  
  
  #select.features <- combine.lists(select.features, sub.factor.feature.spaces$termination)
  base.model.data <- model.matrices[[train.cline]][[matrix.type]][[feature.subspace]]
  #base.model.data[which(is.na(base.model.data), arr.ind = T)] <- 0
  #indx <- apply(base.model.data, 2, function(x) any(is.nan(x) | is.infinite(x)))
  #which(indx)
  
  nonbinding.features <- colnames(base.model.data)[!grepl("^chip|clip", colnames(base.model.data))]
  binding.features <- colnames(base.model.data)[grepl("^chip|clip", colnames(base.model.data))]
  binding.features <- binding.features[grepl(paste0(select.features, collapse = "|"), binding.features)]
  base.model.data <- base.model.data[, c(nonbinding.features, binding.features)]
  select.feature.model <- apply.xgboost(train.model.data = base.model.data, 
                                        test.model.data = NULL, 
                                        model.target = "target_traveling.ratio")
  
  ## Model performance
  # Extract relevant variables from model results
  model <- select.feature.model$model
  X <- select.feature.model$X
  Y <- select.feature.model$Y
  X.test = select.feature.model$X.test
  Y.test <- select.feature.model$Y.test
  train.preds <- select.feature.model$train.preds
  test.preds <- select.feature.model$test.preds
  shap.values <- select.feature.model$shap.values
  samples <- select.feature.model$samples
  model.target.type <- gsub("target_", "", model.target)
  shap.scores <- as.data.frame(select.feature.model$shap.values$shap_score)
  rownames(shap.scores) <- samples
  
  ## Evaluate model performance
  # Calculate root mean squared error of test data
  test.residuals <- Y.test - test.preds
  test.rmse <- sqrt(mean(test.residuals^2))
  ## Calculate the r-squared of the test data 
  # Calculate the mean of the test data observations
  mean.y.test <- mean(Y.test)
  # Calculate the total sum of squares
  tss =  sum((Y.test - mean.y.test)^2)
  # Calculate residual sum of squares
  rss =  sum(test.residuals^2)
  # Calculate R-squared
  rsq.test  =  1 - (rss/tss)
  
  # Observed vs predicted values with residuals
  obs.v.pred = data.frame(obs = Y.test, pred = test.preds, resid = abs(test.residuals))
  # Plot observed vs predicted targets without residual errors
  cons.obs.v.pred.plot <- 
    ggscatter(obs.v.pred, 
              x = "obs", 
              y = "pred",
              color = "resid",
              add ="reg.line",
              conf.int = TRUE,
              add.params = list(color = "red",
                                fill = "lightgray",
                                size = 0.3),
              #palette = pal_uchicago("dark")(9)[9],
              alpha = 1/10
    ) +
    theme_pubr(base_size = plot.base.size)+
    #ggtitle(paste0("Observed vs. Predicted traveling ratio (K562)"))+
    xlab("Observed traveling ratio")+
    ylab("Predicted traveling ratio")+
    ggpubr::stat_cor(cor.coef.name = c("rho"))+
    gradient_color(c("black", pal_uchicago("light")(9)[1]))+
    theme( legend.position = "right", plot.title = element_text(size=10))
  
  return(list(minimal.model.perf = cons.obs.v.pred.plot, 
              select.feature.model = select.feature.model))
  
  # cowplot::save_plot(paste0(OUTPUT, PLOTS,"obs_vs_pred_target_","rna_processing", ".pdf"), 
  #                    plot = obs.v.pred.plot)
  
  
  
  
  #length(select.features$K562)
  
  ## Model performance
  # Extract relevant variables from model results
  X <- select.feature.model$X
  Y <- select.feature.model$Y
  X.test = select.feature.model$X.test
  Y.test <- select.feature.model$Y.test
  train.preds <- select.feature.model$train.preds
  test.preds <- select.feature.model$test.preds
  shap.values <- select.feature.model$shap.values
  samples <- select.feature.model$samples
  model.target.type <- gsub("target_", "", model.target)
  shap.scores <- as.data.frame(select.feature.model$shap.values$shap_score)
  rownames(shap.scores) <- samples
  
  cons.dynamic.features.plots <- 
    visualize.dynamic.features(shap.scores = shap.scores, 
                               X = select.feature.model$X, 
                               train.cline = train.cline, 
                               matrix.type = matrix.type, 
                               feature.subspace = "consensus", 
                               model.target = model.target,
                               plot.base.size = 10, 
                               append = "rna_processing")
  
  optimal.factors.plot <- 
    retrieve.optimal.number.of.factors(cum.factor.contrib = cons.dynamic.features.plots$cum.factor.contrib,
                                       train.cline = train.cline,
                                       matrix.type = matrix.type,
                                       feature.subspace = "rna_processing",
                                       model.target = model.target,
                                       plot.base.size = plot.base.size,
                                       append = "rna_processing")
  optimal.factors <- optimal.factors.plot$optimal.factors
  
  # Visualize features
  optimal.factor.pattern <- paste0(optimal.factors, collapse = "|")
  binding.feature.contribs <- shap.scores[ , grepl("^chip|clip", colnames(shap.scores))]
  optimal.factor.feature.contrib <- binding.feature.contribs[, grepl(optimal.factor.pattern, colnames(binding.feature.contribs))]
  optimal.factor.feature.contrib <- 
    lapply(optimal.factors, function(factor){
      factor.feature.contrib <- optimal.factor.feature.contrib[, grepl(factor, colnames(optimal.factor.feature.contrib))]
      factor.feature.contrib <- sum(abs(factor.feature.contrib))
    })
  optimal.factor.feature.contrib <- as.data.frame(do.call(rbind, optimal.factor.feature.contrib))
  optimal.factor.feature.contrib$factor <- optimal.factors
  colnames(optimal.factor.feature.contrib) <- c("Contribution", "Factor")
  
  cons.optimal.factor.contrib.plot <- 
    ggpubr::ggbarplot(optimal.factor.feature.contrib, 
                      x="Factor",
                      y = "Contribution",
                      fill = "grey",            # Set bar border colors to white
                      palette = "uchicago",            # jco journal color palett. see ?ggpar
                      sort.val = "asc",          # Sort the value in descending order
                      sort.by.groups = FALSE,     # Don't sort inside each group
                      #x.text.angle = 90,          # Rotate vertically x axis texts
                      xlab = "Most influential factors",
                      ylab = "Total factor contribution",
                      rotate = TRUE,
                      ggtheme = theme_pubr(base_size = 10),
                      # main = "Factor model contributions"
    )
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "shap_top_total_influential_factor_contributions_", "rna_processing", "_", train.cline, "_", gsub("\\.", "_", matrix.type), "_" , gsub("\\.", "_", feature.subspace), "_", gsub("\\.", "_", model.target), ".pdf"), optimal.factor.contrib.plot)
  
  
  
  # binding.features <- names(shap.values$mean_shap_score)[grepl("^chip|clip", names(shap.values$mean_shap_score))]
  # all.factors <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*",
  #                     "", binding.features) %>% unique()
  # universe <- factor(as.integer(all.factors %in% optimal.factors))
  # names(universe) <- all.factors
  # cons.optimal.factor.gsea <- perform.gene.set.enrichment.analysis(universe, node.size = 2)
  
  
  ## Evaluate model performance
  # Calculate root mean squared error of test data
  test.residuals <- Y.test - test.preds
  test.rmse <- sqrt(mean(test.residuals^2))
  ## Calculate the r-squared of the test data 
  # Calculate the mean of the test data observations
  mean.y.test <- mean(Y.test)
  # Calculate the total sum of squares
  tss =  sum((Y.test - mean.y.test)^2)
  # Calculate residual sum of squares
  rss =  sum(test.residuals^2)
  # Calculate R-squared
  rsq.test  =  1 - (rss/tss)
  
  # Observed vs predicted values with residuals
  obs.v.pred = data.frame(obs = Y.test, pred = test.preds, resid = abs(test.residuals))
  # Plot observed vs predicted targets without residual errors
  cons.obs.v.pred.plot <- 
    ggscatter(obs.v.pred, 
              x = "obs", 
              y = "pred",
              color = "resid",
              add ="reg.line",
              conf.int = TRUE,
              add.params = list(color = "red",
                                fill = "lightgray",
                                size = 0.3),
              #palette = pal_uchicago("dark")(9)[9],
              alpha = 1/10
    ) +
    theme_pubr(base_size = plot.base.size)+
    ggtitle(paste0("Observed vs. Predicted traveling ratio (K562)"))+
    xlab("Observed traveling ratio")+
    ylab("Predicted traveling ratio")+
    ggpubr::stat_cor(cor.coef.name = c("rho"))+
    gradient_color(c("black", pal_uchicago("light")(9)[1]))
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS,"obs_vs_pred_target_","rna_processing", ".pdf"), 
                     plot = obs.v.pred.plot)
  
  ## Traveling ratio deconvolution
  deconv.target <- 
    visualize.target.deconvolution(force.plot_data = NULL, 
                                   Y = Y, 
                                   train.cline = train.cline, 
                                   matrix.type = matrix.type, 
                                   feature.subspace = "rna_processing", 
                                   model.target = model.target, 
                                   plot.base.size = plot.base.size, 
                                   optimal.factors = optimal.factors, 
                                   shap.scores = shap.scores, 
                                   append = "rna_processing")
  
  
  
  consensus.model.results <- 
    plot_grid(cons.obs.v.pred.plot,
              cons.dynamic.features.plots$top.directional.effectors.plot,
              cons.optimal.factor.contrib.plot,
              deconv.target$deconv.target.plot,
              labels = "AUTO",
              align = "h",
              ncol = 2,
              nrow = 2,
              axis = "b")
  
  shap.force.plots <-
    visualize.shap.forces(shap.scores = shap.scores, 
                          mean.shap.scores = mean.shap.scores,
                          X = X,
                          train.cline = train.cline,
                          matrix.type = matrix.type,
                          feature.subspace = feature.subspace,
                          model.target = model.target,
                          plot.base.size = plot.base.size)
  
  cluster.effectors <- visualize.cluster.factors(force.plot_data = shap.force.plots$force.plot_data,
                                                 shap.scores = shap.scores, 
                                                 X = X, 
                                                 Y=Y)
  
  
  
  
  
  deconvoluted.target.plots <- 
    visualize.target.deconvolution(force.plot_data = shap.force.plots$force.plot_data, 
                                   Y = Y,
                                   train.cline = train.cline,
                                   matrix.type = matrix.type,
                                   feature.subspace = feature.subspace,
                                   model.target = model.target,
                                   plot.base.size = plot.base.size, 
                                   optimal.factors = NULL)
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS,"paper_figure_rna_processing_model.pdf"), 
                     plot = consensus.model.results,
                     ncol = 2,
                     nrow = 2)
  
  return(list(cum.contrib = cons.dynamic.features.plots$cum.factor.contrib.plot,
              directional.contrib = cons.dynamic.features.plots$top.directional.effectors.plot,
              influential.factor.contrib = cons.optimal.factor.contrib.plot,
              obs.v.pred = cons.obs.v.pred.plot,
              deconv.target = deconv.target$deconv.target.plot,
  ))
}
###
# This function prepares the figures for the paper
###
prepare.paper.figures <- function(model.performance, 
                                  sync.model.performance,
                                  cell.type.specific.sample.model.performances, 
                                  shap.summary.plot, 
                                  static.features.plots, 
                                  dynamic.features.plots, 
                                  individual.directional.factor.contrib.plots, 
                                  shap.force.plots, 
                                  cluster.target.distribution.plots, 
                                  deconvoluted.target.plots, 
                                  optimal.factors.plot, 
                                  optimal.factor.gsea, 
                                  target.correlations, 
                                  reverse.model.feature.correlations, 
                                  knockdown.correlations.plot, 
                                  go.term.model.results.plot, 
                                  custom.model.results.plot, 
                                  train.cline){
  
  train.cline <-  "HepG2"
  
  #####
  #  Model performances
  #####
  
  indiv.full.model.perf <- model.performance$obs.v.pred.plot+
    #ggtitle("Observed vs. Predicted Traveling Ratio\nof independent test data set (K562)")+
    xlab("Observed Traveling Ratio")+
    ylab("Predicted Traveling Ratio")+
    theme( legend.position = "right", plot.title = element_text(size=10))
  
  sync.full.model.perf <- sync.model.performance$obs.v.pred.plot+
    #ggtitle("Observed vs. Predicted Traveling Ratio\nof independent cross cell type data set (HepG2)")+
    xlab("Observed Traveling Ratio")+
    ylab("Predicted Traveling Ratio")+
    theme( legend.position = "right", plot.title = element_text(size=10))
  
  cell.type.specific.sample <- cell.type.specific.sample.model.performances$train.v.test.target+
    #ggtitle("Cell type specific samples")+
    xlab("Traveling Ratio (K562 cell line)")+
    ylab("Traveling Ratio (HepG2 cell line)")+
    theme( plot.title = element_text(size=10))
  
  cell.type.specific.sample.prediction.performances <- cell.type.specific.sample.model.performances$cline.specific.sample.predictions.plot+
    #ggtitle("Cell type specific sample prediction performances")+
    xlab("log2(Obs_K562 / Obs_HepG2)")+
    ylab("log2(Pred_K562 / Pred_HepG2)")+
    theme( plot.title = element_text(size=10))
  
  model.perf.grid <-  
    plot_grid(indiv.full.model.perf,
              sync.full.model.perf,
              cell.type.specific.sample,
              cell.type.specific.sample.prediction.performances,
              align = "vh",
              labels="AUTO",
              ncol = 2,
              nrow = 2,
              axis = "b")
  
  cowplot::save_plot(paste0("plots/paper_figure_model_performances_", train.cline,".png"),
                     plot = model.perf.grid,
                     ncol = 2,
                     nrow = 2, 
                     dpi = 600,
                     base_height = 4.71)
  
  
  #####
  # Model array
  #####
  
  go.term.models <- go.term.model.results.plot$model.performance.plot+
    theme_pubr(base_size = 12)+
    theme( legend.position = "right")+
    xlab("")
  cowplot::save_plot(paste0("plots/paper_figure_model_array_", train.cline,".png"),
                     plot = go.term.models,
                     ncol = 2,
                     nrow = 2, 
                     dpi = 600,
                     base_height = 4.71)
  
  
  
  #####
  # Final model
  #####
  #####
  # Traveling Ratio Deconvolution
  #####
  forces.plot <- shap.force.plots$force.plot+
    xlab("Transcripts")+
    ylab("Feature Contribution")+
    ggtitle("")+
    scale_fill_jco()+
    theme(plot.margin=margin(r = 1, unit="cm"))
  
  
  force.cluster.plot <- shap.force.plots$force.plot.cluster.plot+
    ggtitle("")+
    xlab("Transcripts")+
    ylab("Feature Contribution")+
    scale_fill_jco()
  
  deconvoluted.target <- deconvoluted.target.plots$deconv.target.plot+
    ggtitle("")+
    xlab("Traveling Ratio")+
    ylab("Count")+
    theme(plot.margin=margin(t = 0, unit="cm"))
  
  aggregate.class.contributions.plot <- dynamic.features.plots$aggregate.class.contributions.plot+ggtitle("")
  
  feature.contrib.plot.grid <- 
    plot_grid(forces.plot,
              force.cluster.plot,
              deconvoluted.target,
              aggregate.class.contributions.plot,
              labels = "AUTO",
              align = "h",
              ncol = 2,
              nrow = 2,
              axis = "b", 
              rel_heights = c(1.3, 1.3, 0.9, 0.9))
  
  cowplot::save_plot(paste0("plots/paper_figure_feature_contrib_", train.cline,".png"),
                     plot = feature.contrib.plot.grid,
                     ncol = 2,
                     nrow = 2, 
                     dpi = 600)
  
  ###
  all.cluster.factor.contribs.plots <- cluster.effectors$all.cluster.factor.contribs.plots
  all.cluster.factor.contribs.plots.grid <-
    plot_grid(all.cluster.factor.contribs.plots[[1]],
              all.cluster.factor.contribs.plots[[2]],
              labels = "A",
              ncol = 2, 
              nrow = 1)
  
  contingency.test.results
  
  top.influential.grid <- 
    plot_grid(optimal.cluster.factor.directional.contrib$optimal.factor.contrib.plot,
              select.individual.directional.factor.contrib.plots,
              labels = c("B", "C"),
              align = "h",
              ncol = 2,
              nrow = 1,
              axis = "b", 
              rel_widths = c(0.8, 1.2))
  
  
  optimal.factor.target.deconvolution
  
  target.deconvolution.grid <- 
    plot_grid(optimal.factor.target.deconvolution,
              minimal.model.results,
              labels = c("D", "E"),
              align = "h",
              ncol = 2,
              nrow = 1,
              axis = "b")
  
  
  final.model.results.grid <- 
    plot_grid(all.cluster.factor.contribs.plots.grid, 
              top.influential.grid, 
              target.deconvolution.grid, 
              ncol = 1,
              nrow = 3,
              axis = "b")
  
  cowplot::save_plot(paste0("plots/paper_figure_factor_selection_", train.cline,".png"),
                     plot = final.model.results.grid,
                     ncol = 2,
                     nrow = 3, 
                     dpi = 600)
  
  
  ## SUPPLEMENTARY FIGURES
  
  
  #####
  # Feature effects
  #####
  feature.effect.plot <- shap.summary.plot+
    #ggtitle("Feature's effects on the traveling ratio (top 25 features)")+
    ggtitle("")+
    xlab("Feature")+
    theme( axis.text.y = element_text(size=10))
  
  static.feature.effect.plot <- static.features.plots+
    theme(plot.margin=margin(l=0.75,t=0.5, unit="cm"))
  
  top.directional.effectors <- dynamic.features.plots$top.directional.effectors.plot+
    #ggtitle("Greatest factor contributions (top 10 per direction)")+
    ggtitle("")+
    xlab("Factor")+
    ylab("Aggregate SHAP contribution")+
    theme( axis.text.y = element_text(size=10))
  
  indiv.directional.effectors <- individual.directional.factor.contrib.plots+
    theme(plot.margin=margin(t=0.5, unit="cm"))
  
  feature.effects.grid1 <- 
    plot_grid(feature.effect.plot,
              top.directional.effectors,
              labels = c("A", "B"),
              align = "h",
              ncol = 2,
              #nrow = 1,
              axis = "tb")
  
  feature.effects.grid2 <- 
    plot_grid(indiv.directional.effectors,
              static.feature.effect.plot,
              labels = c("C", "D"),
              align = "h",
              ncol = 2,
              #nrow = 1,
              hjust = c(-0.5, -1.5),
              axis = "tb", 
              rel_widths = c(0.6, 1))
  
  feature.effects.grid <- 
    plot_grid(feature.effects.grid1,
              feature.effects.grid2,
              #labels = "AUTO",
              align = "v",
              nrow = 2,
              axis = "b")
  
  cowplot::save_plot(paste0(OUTPUT, PLOTS, "paper_figure_feature_effects_", train.cline,".png"),
                     plot = feature.effects.grid,
                     ncol = 2,
                     nrow = 2, base_height = 5)
  
  
  
  
}

### 
# Create supplementary tables in xls file
###
create.supplementary.figures <- function(shap.summary.plot, static.features.plots){
  
  ## Inverse correlation of transcript traveling ratios with transcript expressions 
  inv.corr.tr.exp <- visualize.traveling.ratio.transcript.expression.distribution(plot.base.size)
  inv.corr.tr.exp.grid <- 
    plot_grid(inv.corr.tr.exp[[1]], 
              inv.corr.tr.exp[[2]],
              ncol = 1, 
              nrow = 2,
              align = "b",
              labels = "AUTO")
  cowplot::save_plot(paste0("plots/suppl_figure_inverse_corr_tr_exp.png"),
                     plot = inv.corr.tr.exp.grid,
                     ncol = 1,
                     nrow = 2, 
                     dpi = 600)
  
  ## Traveling ratio distribution
  traveling.ratio.distributions <- visualize.traveling.ratio(plot.base.size)
  traveling.ratio.distributions <- 
    plot_grid(traveling.ratio.distributions[[1]], 
              traveling.ratio.distributions[[2]],
              ncol = 1, 
              nrow = 2,
              align = "b",
              labels = "AUTO")
  cowplot::save_plot(paste0("plots/suppl_figure_traveling_ratios.png"),
                     plot = traveling.ratio.distributions,
                     ncol = 1,
                     nrow = 2, 
                     dpi = 600)
  
  sequence.feature.distributions <- visualize.features()
  cowplot::save_plot(paste0("plots/suppl_figure_sequence_feature_distribution_K562.png"),
                     plot = sequence.feature.distributions$K562,
                     ncol = 2,
                     nrow = 8, 
                     dpi = 600,  
                     base_width = 14, 
                     base_height = 5)
  
  cowplot::save_plot(paste0("plots/suppl_figure_sequence_feature_distribution_HepG2.png"),
                     plot = sequence.feature.distributions$HepG2,
                     ncol = 2,
                     nrow = 8, 
                     dpi = 600, 
                     base_width = 14, 
                     base_height = 5)
  
  ## SHAP contributions
  shap.summary.plot <- shap.summary.plot+
    theme( axis.text.y = element_text(size=12))
  
  cowplot::save_plot(paste0("plots/suppl_figure_shap_contributions.png"),
                     plot = shap.summary.plot,
                     ncol = 2,
                     nrow = 2, 
                     dpi = 600)
  
  
  cowplot::save_plot(paste0("plots/suppl_figure_shap_static_feature_distributions.png"),
                     plot = static.features.plots,
                     ncol = 2,
                     nrow = 8, 
                     dpi = 600,
                     base_width = 14, 
                     base_height = 5)
  
  
  cowplot::save_plot(paste0("plots/suppl_figure_shap_non_binding_feature_contributions.png"),
                     plot = dynamic.features.plots$aggregate.non.binding.feature.contrib.plot,
                     ncol = 1,
                     nrow = 1, 
                     dpi = 600)

  
  visualize.7SK.binders()
  
  
  
}
### 
# Create supplementary tables in xls file
###
create.supplementary.tables <- function(feature.vectors, model.matrices, model.results, random.model.results){
  
  ## CHIPseq Factors per cell line
  K562.chip.factors <- data.frame(Factor = unique(CHIPseq$K562$hgnc_symbol), 
                                  ID = CHIPseq$K562$ensembl_id[match(unique(CHIPseq$K562$hgnc_symbol), CHIPseq$K562$hgnc_symbol)])
  HepG2.chip.factors <- data.frame(Factor = unique(CHIPseq$HepG2$hgnc_symbol), 
                                  ID = CHIPseq$HepG2$ensembl_id[match(unique(CHIPseq$HepG2$hgnc_symbol), CHIPseq$HepG2$hgnc_symbol)])
  
  ## ENCODE Accession table
  encode.chipseq.K562.accession.table <- read.table(paste0(OUTPUT, "filtered.encode.chipseq.accessions.table.K562.txt"))
  colnames(encode.chipseq.K562.accession.table) <- c("Experiment", "Acession", "Target")
  encode.chipseq.HepG2.accession.table <- read.table(paste0(OUTPUT, "filtered.encode.chipseq.accessions.table.HepG2.txt"))
  colnames(encode.chipseq.HepG2.accession.table) <- c("Experiment", "Acession", "Target")
  # eCLIPseq factors
  K562.eclip.factors <- data.frame(Factor = unique(eCLIPseq$K562$hgnc_symbol), 
                                  ID = eCLIPseq$K562$ensembl_id[match(unique(eCLIPseq$K562$hgnc_symbol), eCLIPseq$K562$hgnc_symbol)])
  HepG2.eclip.factors <- data.frame(Factor = unique(eCLIPseq$HepG2$hgnc_symbol), 
                                   ID = eCLIPseq$HepG2$ensembl_id[match(unique(eCLIPseq$HepG2$hgnc_symbol), eCLIPseq$HepG2$hgnc_symbol)])
  
  encode.eclipseq.K562.accession.table <- read.table(paste0(OUTPUT, "filtered.encode.eclipseq.accessions.table.K562.txt"))
  colnames(encode.eclipseq.K562.accession.table) <- c("Experiment", "Acession", "Target")
  encode.eclipseq.HepG2.accession.table <- read.table(paste0(OUTPUT, "filtered.encode.eclipseq.accessions.table.HepG2.txt"))
  colnames(encode.eclipseq.HepG2.accession.table) <- c("Experiment", "Acession", "Target")
  
  K562.7SK.binding <- data.frame(Factor = unique(novel.rn7sk.binders$K562$hgnc_symbol), 
                                 Factor_ID = unique(novel.rn7sk.binders$K562$ensembl_id),
                                 RN7SK = unique(novel.rn7sk.binders$K562$RN7SK.tx.id))
  
  HepG2.7SK.binding <- data.frame(Factor = unique(novel.rn7sk.binders$HepG2$hgnc_symbol), 
                                 Factor_ID = unique(novel.rn7sk.binders$HepG2$ensembl_id),
                                 RN7SK = unique(novel.rn7sk.binders$HepG2$RN7SK.tx.id))
  
  ## Factor bindings on genomic regions
  pred.mat <- list(model.matrices$K562$individual.model.matrices$All,
                   model.matrices$HepG2$individual.model.matrices$All)
  factor.binding.dist <- 
    lapply(pred.mat, function(cline.features){
      binding.signals <- colnames(cline.features)[grepl("^chip|clip", colnames(cline.features))]
      factor.bindings <- cline.features[ , binding.signals]
      factor.binding.features <- colnames(factor.bindings)
      factors <- unique(gsub("chip\\.|clip\\.|\\.three_prime|\\.five_prime|\\.introns|\\.coding\\.exons|\\.Proximal.ncRNA.*", "", factor.binding.features))
      factor.stats <- 
        lapply(factors, function(factor){
          regions <- c("five_prime","introns","coding.exons","three_prime")
          region.bindings <- 
            lapply(regions, function(region){
              assays <- c("chip", "clip")
              assay.bindings <- 
                lapply(assays, function(assay){
                  target.factor.features <- paste0(assay, ".", factor, ".", region)
                  if(!target.factor.features %in% binding.signals){return(0)}
                  temp.binding.signals <- factor.bindings[, target.factor.features]
                  n <- table(temp.binding.signals)[["1"]]
                })
              names(assay.bindings) <- paste0(assays, ".", region)
              
              ncrna.assay.bindings <- 
                lapply(assays, function(assay){
                  target.factor.features <- paste0(assay, ".", factor, ".", region)
                  ncrna.target.factor.features <- paste0(target.factor.features, c(".Proximal.ncRNA.1", ".Proximal.ncRNA.2"))
                  ncrna.binding.counts <- 
                    lapply(ncrna.target.factor.features, function(ncrna.bindings){
                      if(!ncrna.bindings %in% binding.signals){return(0)}
                      temp.binding.signals <- factor.bindings[, ncrna.bindings]
                      n <- table(temp.binding.signals)[["1"]]
                    })
                  names(ncrna.binding.counts) <- paste0(assay, ".",region, ".", c("Proximal.ncRNA.1", "Proximal.ncRNA.2"))
                  ncrna.binding.counts
                })
              return(c(assay.bindings, unlist(ncrna.assay.bindings)))
            })
          region.bindings <- unlist(region.bindings)
          region.bindings <- data.frame(region.bindings)
          colnames(region.bindings) <- factor
          return(region.bindings)
        })
      factor.stats <- do.call(cbind, factor.stats)
    })
  names(factor.binding.dist) <- c("K562", "HepG2")
  
  factor.binding.dist <- lapply(factor.binding.dist, t)
  factor.binding.dist <- lapply(factor.binding.dist, as.data.frame)
  
  
  ## Factors per GO Term
  sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
  sub.factor.feature.spaces <- 
    build.sequence.specific.biologically.functional.feature.sub.spaces(stratify = sub.factor.feature.spaces)
  sub.factor.feature.spaces <- sub.factor.feature.spaces[c("Chromatin", "Initiation", "Elongation", "Termination", "Processing","Splicing", "Translation", "7SK.Binding", "Elongation+7SK")]
  factors.per.process <- lapply(CARGS$cell.line, function(cline){lapply(sub.factor.feature.spaces, function(x){x[[cline]]})})
  names(factors.per.process) <- CARGS$cell.line
  
  n <- max(lengths(factors.per.process$K562))
  factors.per.process$K562 <- lapply(factors.per.process$K562, function(x){length(x) <- n;x})
  n <- max(lengths(factors.per.process$HepG2))
  factors.per.process$HepG2 <- lapply(factors.per.process$HepG2, function(x){length(x) <- n;x})
  
  factors.per.process.K562 <- do.call(cbind, factors.per.process$K562)
  factors.per.process.HepG2 <- do.call(cbind, factors.per.process$HepG2)
  factors.per.process.K562[which(is.na(factors.per.process.K562))] <- ""
  factors.per.process.HepG2[which(is.na(factors.per.process.HepG2))] <- ""
  factors.per.process.K562 <- as.data.frame(factors.per.process.K562)
  factors.per.process.HepG2 <- as.data.frame(factors.per.process.HepG2)
  ## Sequence specificity
  cline.sequence.specific.factors <- 
    lapply(CARGS$cell.line, function(cline){
      sequence.specific.factors <-  data.frame(matrix(0, ncol = 4, nrow = length(unique(c(SEQ.SPEC$dbp.factors[[cline]], SEQ.SPEC$rbp.factors[[cline]])))))
      rownames(sequence.specific.factors) <- unique(c(SEQ.SPEC$dbp.factors[[cline]], SEQ.SPEC$rbp.factors[[cline]]))
      colnames(sequence.specific.factors) <- c("SS", "NSS", "DBP", "RBP")
      
      sequence.specific.factors[SEQ.SPEC$dbp.factors[[cline]], "DBP"] <- 1
      sequence.specific.factors[SEQ.SPEC$rbp.factors[[cline]], "RBP"] <- 1
      
      sequence.specific.factors[rownames(sequence.specific.factors) %in% SEQ.SPEC$sequence.specific, "SS"] <- 1
      sequence.specific.factors[rownames(sequence.specific.factors) %in% SEQ.SPEC$nonsequence.specific, "NSS"] <- 1
      return(sequence.specific.factors)
    })
  names(cline.sequence.specific.factors) <- CARGS$cell.line
  
  ## Factors per subspace
  dbp.factors <- 
    lapply(CHIPseq, function(peaks){
      unique(peaks$hgnc_symbol)
    }) %>% unlist %>% unique 
  
  # Retrieve all RNA binding factors
  rbp.factors <- 
    lapply(eCLIPseq, function(peaks){
      unique(peaks$hgnc_symbol)
    }) %>% unlist %>% unique 
  
  all.factors <- unique(c(dbp.factors, rbp.factors))
  subspace.factors <- data.frame(matrix(0, ncol = length(all.factors), nrow = length(names(model.matrices$K562$individual.model.matrices))))
  rownames(subspace.factors) <- names(model.matrices$K562$individual.model.matrices)
  colnames(subspace.factors) <- all.factors
  
  lapply(names(model.matrices), function(cline){
    lapply(names(model.matrices[[cline]]), function(modeltype){
      lapply(names(model.matrices[[cline]][[modeltype]]), function(subspace){
        features <- colnames(model.matrices[[cline]][[modeltype]][[subspace]])
        features <- features[grepl("chip|clip", features)]
        features <- gsub("^chip.|^clip.|.five_prime|.three_prime|.introns|.coding.exons|.Proximal.*", "", features)
        factors <- unique(features)
        subspace.factors[subspace, factors] <<- 1
        return()
      })
    })
  })
  subspace.factors <- subspace.factors[subspaces, ]
  subspace.factors <- as.data.frame(t(subspace.factors))
  
  ## 7SK associated factors
  sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
  rn7sk.factors <- unique(unlist(sub.factor.feature.spaces$All.7SK))

  rn7sk.factor.indicator <- data.frame(matrix(0, ncol = length(rn7sk.factors), nrow = 2))
  rownames(rn7sk.factor.indicator) <- CARGS$cell.line
  colnames(rn7sk.factor.indicator) <- rn7sk.factors

  lapply(names(sub.factor.feature.spaces$All.7SK), function(cline){
    lapply(sub.factor.feature.spaces$All.7SK[[cline]], function(factor){
      rn7sk.factor.indicator[cline, factor] <<- 1
    })
  })

  ## Established pausing factors
  sub.factor.feature.spaces <- build.pausing.associated.factor.sub.feature.spaces()
  sub.factor.feature.spaces <- 
    build.sequence.specific.biologically.functional.feature.sub.spaces(stratify = sub.factor.feature.spaces)
  known.pausing.factors <- unique(unlist(c(sub.factor.feature.spaces$`Elongation`$K562, sub.factor.feature.spaces$`Elongation`$HepG2)))
  known.pausing.factors <- unique(known.pausing.factors[known.pausing.factors %in% unique(c(unlist(dbp.factors), unlist(rbp.factors)))])
  known.pausing.factors <- data.frame(Factors = known.pausing.factors)
  
  ## Hyperparameters
  hyperparameters <- list(eta = 0.08,
                     max_depth = 7,
                     gamma = 0.5,
                     lambda = 0.1,
                     alpha = 0.01,
                     colsample_bytree = 0.7,
                     subsample = 0.7, 
                     min_child_weight = 50,
                     booster="gbtree")
  hyperparameters <- data.frame(Parameter = names(hyperparameters),
                           Value = unlist(hyperparameters))
  rownames(hyperparameters) <- NULL
  
  ## Model results
  subspaces <- c("Chromatin", 
                 "Initiation", 
                 "Elongation", 
                 "Termination", 
                 "Processing", 
                 "Splicing",
                 "Translation", 
                 "7SK.Binding",
                 "Elongation+7SK", 
                 "All")
  subspaces <- c(subspaces, paste0(subspaces, ":nss"), paste0(subspaces, ":ss"))
  subspaces <- sort(subspaces)
  select.model.results.table <- model.results$performances[model.results$performances$target == "target_traveling.ratio" &
                                                             model.results$performances$subspace %in% subspaces,]
  select.model.results.table <- select.model.results.table[order(select.model.results.table$cline, select.model.results.table$modeltype, select.model.results.table$subspace),]
  # remove target counts
  select.model.results.table$nfeatures <- select.model.results.table$nfeatures-2
  select.model.results.table <- select.model.results.table[, c("cline", "modeltype", "subspace", "train.rsqrd",  "test.rsqrd", "nfeatures", "nfactors", "mean.shap")]
  random.model.results <- random.model.results[, c("cline", "modeltype", "subspace", "train.rsqrd",  "test.rsqrd", "nfeatures", "nfactors", "mean.shap")]
  select.model.results.table <- rbind(select.model.results.table, random.model.results)
  rownames(select.model.results.table) <- 1:dim(select.model.results.table)[1]
  
  ## Write out full feature matrix for each cell line
  WriteXLS(list(K562.chip.factors,
                HepG2.chip.factors,
                encode.chipseq.K562.accession.table,
                encode.chipseq.HepG2.accession.table,
                K562.eclip.factors,
                HepG2.eclip.factors,
                encode.eclipseq.K562.accession.table,
                encode.eclipseq.HepG2.accession.table,
                K562.7SK.binding,
                HepG2.7SK.binding,
                factor.binding.dist$K562, 
                factor.binding.dist$HepG2,
                known.pausing.factors, 
                factors.per.process.K562,
                factors.per.process.HepG2,
                cline.sequence.specific.factors$K562,
                cline.sequence.specific.factors$HepG2, 
                subspace.factors,
                #model.matrices$K562$individual.model.matrices$All,
                #model.matrices$HepG2$individual.model.matrices$All,
                hyperparameters,
                select.model.results.table),
           ExcelFileName = "plots/Supplementary_Tables.xls",
           SheetNames = c("S1 K562 CHIPseq Factors",
                          "S2 HepG2 CHIPseq Factors",
                          "S3 K562 CHIPseq Accessions",
                          "S4 HepG2 CHIPseq Accessions",
                          "S5 K562 eCLIPseq Factors",
                          "S6 HepG2 eCLIPseq Factors",
                          "S7 K562 eCLIPseq Accessions",
                          "S8 HepG2 eCLIPseq Accessions",
                          "S9 K562 7SK Binding factors",
                          "S10 HepG2 7SK Binding Factors",
                          "S11 K562 Factor Bindings",
                          "S12 HepG2 Factor Bindings",
                          "S13 Known Pausing Factors",
                          "S14 K562 Factors per Process",
                          "S15 HepG2 Factors per Process",
                          "S16 K562 Sequence Specificity",
                          "S17 HepG2 Sequence Specificity", 
                          "S18 Subspace Factors Presence",
                          #"S19 Full K562 Model Matrix",
                          #"S20 Full HepG2 Model Matrix",
                          "S19 Hyperparameters",
                          "S20 All Model Results"),
           row.names = T,
           col.names = T,
           BoldHeaderRow = T, 
           verbose = T)
}
### 
# Create figures for concept drawing
###
crate.concept.figures <- function(model.performance, deconvoluted.target.plots){
  indiv.full.model.perf <- model.performance$obs.v.pred.plot+
    #ggtitle("Observed vs. Predicted Traveling Ratio\nof independent test data set (K562)")+
    xlab("Observed Traveling Ratio")+
    ylab("Predicted Traveling Ratio")+
    theme( legend.position = "none", plot.title = element_text(size=10))
  
  cowplot::save_plot(paste0("plots/concept_figure_model_performance.png"),
                     plot = indiv.full.model.perf,
                     ncol = 1,
                     nrow = 1, 
                     dpi = 600)
  
  
  traveling.ratio.dist <- deconvoluted.target.plots$conv.target.plot+ggtitle("")+xlab("Traveling Ratio")
  cowplot::save_plot(paste0("plots/concept_figure_traveling_ratio.png"),
                     plot = traveling.ratio.dist,
                     ncol = 1,
                     nrow = 1, 
                     dpi = 600)
  
  deconv.traveling.ratio <- deconvoluted.target.plots$deconv.target.plot+ggtitle("")+xlab("Deconvoluted Traveling Ratio")
  cowplot::save_plot(paste0("plots/concept_figure_deconv_target.png"),
                     plot = deconv.traveling.ratio,
                     ncol = 1,
                     nrow = 1, 
                     dpi = 600)
  
  
  cluster.plot <- shap.force.plots$force.plot.cluster.plot+ggtitle("Transcriptional States")+xlab("")
  cowplot::save_plot(paste0("plots/concept_figure_cluster_plot.png"),
                     plot = cluster.plot,
                     ncol = 1,
                     nrow = 1, 
                     dpi = 600)
  
  effectors.plot <- cluster.effectors$all.cluster.factor.contribs.plots[[1]]+ggtitle("Pause Modulators")
  cowplot::save_plot(paste0("plots/concept_figure_effectors.png"),
                     plot = effectors.plot,
                     ncol = 1,
                     nrow = 1, 
                     dpi = 600)
  
}
### 
# This function performs a fisher test
###
perform.fisher.test <- function(categories, counts){
  
  ## Check enrichment of pausing.related factors in full models selected features
  contingency <- matrix(0, nrow = 2, ncol = 2, dimnames = categories)
  contingency[categories[[1]][[1]], categories[[2]][[1]]] <- counts[1]
  contingency[categories[[1]][[1]], categories[[2]][[2]]] <- counts[2]
  contingency[categories[[1]][[2]], categories[[2]][[1]]] <- counts[3]
  contingency[categories[[1]][[2]], categories[[2]][[2]]] <- counts[4]
  fisher.test(contingency, alternative = "greater")
}
### TODO
# This function performs an enrichment test of proximal ncRNAs as
###
proximal.ncrna.enrichments <- function(){
  cline <- "K562"
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
                                   "linc.RNA.transcripts", "cline"), environment())
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
  
  
  
  nearest.lincRNA.metrics <- 
    lapply(nearest.lincRNA.metrics, function(ncrna.matrics){
      cbind(ncrna.matrics, index = 1:num.of.nearest.lincRNAs)
    })
  nearest.lincRNA.metrics <- do.call(rbind, nearest.lincRNA.metrics)
  
  q1 <- 
    names(q1) <- valid.coding.transcripts[[cline]]$transcript_id
  
  temp <- lapply(1:num.of.nearest.lincRNAs, function(index) {nearest.lincRNA.metrics[nearest.lincRNA.metrics$`1:num.of.nearest.lincRNAs`==index, "quant"]})
  
  library(corrplot)
  x <- do.call(cbind, temp)
  x <- apply(x, 2, rescale)
  x <- cbind(x, expression)
  corrplot(x, type="upper")
  
}
### 
# This function plots the distribution of transcript quantifications
###
visualize.transcript.quantifications <- function(){
  
  lapply(names(GENE.QUANT), function(cline){
    cline.quantifications <- GENE.QUANT[[cline]]
    p <- 
      ggpubr::gghistogram(cline.quantifications, 
                          x = "FPKM",
                          add = "mean", 
                          binwidth = 0.05)+
      xlab("log10(FPKM)")+
      ggtitle(paste0("Transcript quantification distribution (", cline, ")"))+
      theme_pubr(base_size = plot.base.size)
    save.plot(paste0("transcript_quantifications_", cline), p)
    return()
  })
  return()
}
### 
# This function plots the distribution of traveling.ratios
###
visualize.traveling.ratio <- function(plot.base.size){

  blue.color <- rgb(0.2,0.4,0.6,0.6)
  traveling.ratio.distributions <- 
    lapply(names(traveling.ratios), function(cline){
      cline.traveling.ratios <- traveling.ratios[[cline]]
      p <- 
        ggpubr::gghistogram(as.data.frame(cline.traveling.ratios), 
                            x = "traveling.ratio",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.2)+
        xlab("Pausing Index")+
        ggtitle(paste0("Pausing index distribution (", cline, ")"))+
        theme_pubr(base_size = plot.base.size)+
        scale_x_continuous(breaks = seq(-20, 20, by=5))
      return(p)
    })
  return(traveling.ratio.distributions)
  
}
### 
# This function plots the traveling ratios agains transcript expressions [DEPRECATED]
###
visualize.traveling.ratio.transcript.expression.distribution <- function(plot.base.size){
  
  p <- 
    lapply(CARGS$cell.line, function(cline){
      cline.traveling.ratios <- traveling.ratios[[cline]]
      cline.transcript.expressions <- GENE.QUANT[[cline]]
      common.samples <- intersect(cline.traveling.ratios$transcript_id, 
                                  cline.transcript.expressions$transcript_id)
      tr.v.exp <- data.frame(traveling.ratio = cline.traveling.ratios$traveling.ratio[match(common.samples, cline.traveling.ratios$transcript_id)],
                             transcript.expression = cline.transcript.expressions$FPKM[match(common.samples, cline.transcript.expressions$transcript_id)])
      
      p <- 
        ggscatter(tr.v.exp,
                  x = "traveling.ratio",
                  y = "transcript.expression",
                  palette = "uchicago",
                  #add = "reg.line",
                  add.params = list(color = "blue",
                                    fill =  "lightgray"),
                  conf.int = T,
                  alpha = 1/5,
                  rug = F) +
        stat_cor(cor.coef.name = c("rho"), 
                 method = "pearson")+
        xlab("Pausing Index")+
        ylab("Transcript Expression")+
        # stat_binhex(color = "white")+
        # scale_fill_gradient(low = "white", high = "black")+
        ggtitle(paste0("Pausing Indices vs. Transcript Expressions (", cline, ")"))+
        theme_pubr(base_size = plot.base.size)
      
      return(p)
    })
  return(p)
  
}
### 
# This function plots a Venn diagram of 7SK binders
###
visualize.7SK.binders <- function(){
  
  ## RN7SK binder overlaps between cell lines including pseudo 7Sk trascnript binders
  a1.clipped.prots.K562 <- as.character(unique(novel.rn7sk.binders$K562$ensembl_id))
  a2.clipped.prots.HepG2 <- as.character(unique(novel.rn7sk.binders$HepG2$ensembl_id))
  
  #pdf("")
  # Draw venn diagram
  pdf(paste(OUTPUT, PLOTS,"rn7sk_binders_overlap.pdf", sep = ""), width = 4, height = 4)
  venn.plot <- draw.pairwise.venn(
    area1 = length(a1.clipped.prots.K562),
    area2 = length(a2.clipped.prots.HepG2),
    cross.area =  length(intersect(a1.clipped.prots.K562, a2.clipped.prots.HepG2 )),
    category = c("K562" , "HepG2"),
    fill = viridis(2, end = 0.75),
    cat.cex = 1,
    scaled = F
  )
  dev.off()
  
  
}
### 
# This function plots the egnineered features
###
visualize.features <- function(){
  
  plot.base.size <- 18
  title.size <- 18
  blue.color <- rgb(0.2,0.4,0.6,0.6)
  p <- 
    lapply(CARGS$cell.line, function(cline){
      feature.vectors <- model.matrices[[cline]]$individual.model.matrices$All
      
      tx.housekeeping <- as.character(feature.vectors$housekeeping)
      tx.housekeeping <- table(tx.housekeeping)
      tx.housekeeping <- data.frame(housekeeping = gsub("chr", "", names(tx.housekeeping)), count = tx.housekeeping)
      p1 <- 
        ggpubr::ggbarplot(tx.housekeeping, 
                          x="housekeeping",
                          y = "count.Freq",
                          color = "black",            # Set bar border colors to white
                          fill = blue.color,
                          palette = "uchicago",            # jco journal color palett. see ?ggpar
                          sort.val = "asc",          # Sort the value in descending order
                          sort.by.groups = FALSE,     # Don't sort inside each group
                          x.text.angle = 90,          # Rotate vertically x axis texts
                          xlab = "Housekeeping",
                          ylab = "Count",
                          rotate = TRUE,
                          ggtheme = theme_pubr(plot.base.size),
                          main = "Number of Housekeeping Genes"
        )+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      tx.chr.loc <- as.character(feature.vectors$tx.chr.loc)
      tx.chr.loc <- table(tx.chr.loc)
      tx.chr.loc <- data.frame(chromosome = gsub("chr", "", names(tx.chr.loc)), count = tx.chr.loc)
      p2 <- 
        ggpubr::ggbarplot(tx.chr.loc, 
                          x="chromosome",
                          y = "count.Freq",
                          color = "black",            # Set bar border colors to white
                          fill = blue.color,
                          palette = "uchicago",            # jco journal color palett. see ?ggpar
                          sort.val = "asc",          # Sort the value in descending order
                          legend.title = "Impact on target",
                          sort.by.groups = FALSE,     # Don't sort inside each group
                          x.text.angle = 90,          # Rotate vertically x axis texts
                          xlab = "Chromosome",
                          ylab = "Number of coding transcripts",
                          rotate = TRUE,
                          ggtheme = theme_pubr(plot.base.size),
                          main = "Number of transcripts per chromosome"
        )+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      tx.strand <- as.character(feature.vectors$tx.strand)
      tx.strand <- table(tx.strand)
      names(tx.strand) <- c("+", "-")
      tx.strand <- data.frame(tx.strand = names(tx.strand), count = tx.strand)
      p3 <- 
        ggpubr::ggbarplot(tx.strand, 
                          x="tx.strand",
                          y = "count.Freq",
                          color = "black",            # Set bar border colors to white
                          fill = blue.color,
                          palette = "uchicago",            # jco journal color palett. see ?ggpar
                          sort.val = "asc",          # Sort the value in descending order
                          sort.by.groups = FALSE,     # Don't sort inside each group
                          x.text.angle = 90,          # Rotate vertically x axis texts
                          xlab = "Strand",
                          ylab = "Count",
                          rotate = TRUE,
                          ggtheme = theme_pubr(plot.base.size),
                          main = "Number of transcripts per strand"
        )+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      tx.loc <- feature.vectors$tx.loc
      p4 <- 
        ggpubr::gghistogram(data.frame(tx.loc = tx.loc), 
                            x = "tx.loc",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("Location (rescaled)")+
        ggtitle(paste0("Position on the linear genome"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      tx.len <- feature.vectors$tx.len
      p5 <- 
        ggpubr::gghistogram(data.frame(tx.len = tx.len), 
                            x = "tx.len",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("Length (rescaled)")+
        ggtitle(paste0("Transcrip length"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      tx.tss.width <- feature.vectors$tx.tss.width
      p6 <- 
        ggpubr::gghistogram(data.frame(tx.tss.width = tx.tss.width), 
                            x = "tx.tss.width",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("CAGE TSS cluster width (rescaled)")+
        ggtitle(paste0("CAGE TSS cluster width"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      tx.tss.at.cont <- feature.vectors$tx.tss.at.cont
      p7 <- 
        ggpubr::gghistogram(data.frame(tx.tss.at.cont = tx.tss.at.cont), 
                            x = "tx.tss.at.cont",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("CAGE TSS cluster AT frequency")+
        ggtitle(paste0("CAGE TSS cluster AT frequency"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      tx.gc.seq <- feature.vectors$tx.gc.seq
      p8 <- 
        ggpubr::gghistogram(data.frame(tx.gc.seq = tx.gc.seq), 
                            x = "tx.gc.seq",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("Transcript GC frequency")+
        ggtitle(paste0("Transcript GC frequency"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      tx.ex.count <- feature.vectors$tx.ex.count
      p9 <- 
        ggpubr::gghistogram(data.frame(tx.ex.count = tx.ex.count), 
                            x = "tx.ex.count",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("Transcript exon count (rescaled)")+
        ggtitle(paste0("Exon count"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      
      tx.ex.ratio <- feature.vectors$tx.ex.ratio
      p10 <- 
        ggpubr::gghistogram(data.frame(tx.ex.ratio = tx.ex.ratio), 
                            x = "tx.ex.ratio",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("Exon density (rescaled)")+
        ggtitle(paste0("Transcript exon density"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      
      tx.ex.width <- feature.vectors$tx.ex.width
      p11 <- 
        ggpubr::gghistogram(data.frame(tx.ex.width = tx.ex.width), 
                            x = "tx.ex.width",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("Average exon width (rescaled)")+
        ggtitle(paste0("Average exon width in transcript"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      tx.ex.seq <- feature.vectors$tx.ex.seq
      p12 <- 
        ggpubr::gghistogram(data.frame(tx.ex.seq = tx.ex.seq), 
                            x = "tx.ex.seq",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("Proportion of exonic sequence (rescaled)")+
        ggtitle(paste0("Proportion of exonic sequence per transcript"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      tx.cpg.island.dist <- feature.vectors$cpg.island.dist
      p13 <- 
        ggpubr::gghistogram(data.frame(tx.cpg.island.dist = tx.cpg.island.dist), 
                            x = "tx.cpg.island.dist",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("Distance to nearest CpG island from TSS (rescaled)")+
        ggtitle(paste0("Distance to nearest CpG Island"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      cpg.island.length <- feature.vectors$cpg.island.length
      p14 <- 
        ggpubr::gghistogram(data.frame(cpg.island.length = cpg.island.length), 
                            x = "cpg.island.length",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("Island length (rescaled)")+
        ggtitle(paste0("Length of nearest CpG Island"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      cpg.island.count <- feature.vectors$cpg.island.count
      p15 <- 
        ggpubr::gghistogram(data.frame(cpg.island.count = cpg.island.count), 
                            x = "cpg.island.count",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("CpG island counts (rescaled)")+
        ggtitle(paste0("Number of CpG islands"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      cpg.island.percent.cpg <- feature.vectors$cpg.island.percent.cpg
      p16 <- 
        ggpubr::gghistogram(data.frame(cpg.island.percent.cpg = cpg.island.percent.cpg), 
                            x = "cpg.island.percent.cpg",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("Percent CpG")+
        ggtitle(paste0("Percentage of CpG in the nearest island"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      cpg.island.percent.cg <- feature.vectors$cpg.island.percent.cg
      p17 <- 
        ggpubr::gghistogram(data.frame(cpg.island.percent.cg = cpg.island.percent.cg), 
                            x = "cpg.island.percent.cg",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("Percent CG")+
        ggtitle(paste0("Percentage of CG in the nearest island"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      cpg.island.exp.obs <- feature.vectors$cpg.island.exp.obs
      p18 <- 
        ggpubr::gghistogram(data.frame(cpg.island.exp.obs = cpg.island.exp.obs), 
                            x = "cpg.island.exp.obs",
                            add = "mean", 
                            fill = blue.color,
                            binwidth = 0.01)+
        xlab("Observed/Expected")+
        ggtitle(paste0("Ratio of observed(cpgNum) to expected(numC*numG/length) CpG in island"))+
        theme_pubr(base_size = plot.base.size)+
        theme(plot.title = element_text(size = title.size, face = "bold"))
      
      p <- 
        plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9,p10, p11, p12, p13, p14, p15, p16,p17, p18,
                  labels = "AUTO", 
                  label_size = 12, 
                  align = "h",
                  ncol = 2,
                  nrow = 9)
      cowplot::save_plot(paste0("./plots/supplementary_figure_feature_distributions_", cline, ".png"), 
                         plot = p, 
                         ncol = 2, 
                         nrow = 9,
                         base_width = 14, 
                         base_height = 5)
      
      return(p)
    })
  names(p) <- CARGS$cell.line
  return(p)
}

temp <- function(){
  
  
  # Calculate root mean squared error of training data
  train.residuals <- Y - train.preds
  train.rmse <- sqrt(mean(train.residuals^2))
  ## Calculate the r-squared of the test data 
  # Calculate the mean of the test data observations
  mean.y.train <- mean(Y)
  # Calculate the total sum of squares
  train.tss =  sum((Y - mean.y.train)^2)
  # Calculate residual sum of squares
  train.rss =  sum(train.residuals^2)
  # Calculate R-squared
  rsq.train  =  1 - (train.rss/train.tss)
  # Plot distribution of residuals
  residuals <- data.frame(residuals = train.residuals)
  
  ggplot(residuals, aes(x = residuals)) +
    geom_density(fill = 'grey50', color = 'white', alpha = 0.7) +
    theme_minimal()
  
  # Observed vs predicted values with residuals
  obs.v.pred = data.frame(obs = Y, pred = train.preds, resid = train.residuals)
  
  # Plot observed vs predicted targets without residual errors
  ggplot(obs.v.pred, aes(x = obs, y = pred)) +
    geom_point(alpha = 1/10) +
    theme_minimal()
  
  # Plot observed vs predicted targets with residual errors
  ggplot(obs.v.pred, aes(x = obs, y = pred, size = resid)) +
    geom_point(alpha = 1/10) +
    theme_minimal()
  
  
  ## Evaluate model performance
  # Calculate root mean squared error of test data
  test.residuals <- Y.test - test.preds
  test.rmse <- sqrt(mean(test.residuals^2))
  ## Calculate the r-squared of the test data 
  # Calculate the mean of the test data observations
  mean.y.test <- mean(Y.test)
  # Calculate the total sum of squares
  tss =  sum((Y.test - mean.y.test)^2)
  # Calculate residual sum of squares
  rss =  sum(test.residuals^2)
  # Calculate R-squared
  rsq.test  =  1 - (rss/tss)
  
  # Plot distribution of residuals
  residuals <- data.frame(residuals = test.residuals)
  ggplot(residuals, aes(x = residuals)) +
    geom_density(fill = 'grey50', color = 'white', alpha = 0.7) +
    theme_minimal()
  
  # Observed vs predicted values with residuals
  obs.v.pred = data.frame(obs = Y.test, pred = test.preds, resid = test.residuals)
  
  # Plot observed vs predicted targets without residual errors
  ggplot(obs.v.pred, aes(x = obs, y = pred)) +
    geom_point(alpha = 1/10) +
    theme_minimal()
  
  # Plot observed vs predicted targets with residual errors
  ggplot(obs.v.pred, aes(x = obs, y = pred, size = resid)) +
    geom_point(alpha = 1/10) +
    theme_minimal()
  
}
