###
# As opposed to the GRanges 'reduce' or 'findOverlaps' function, this function 
# allows relative overlapping and subsequent merging of overlapping ranges
##
reduceGRanges <- function(x, y, min=0.6, intersect = T, both = F) {
  #x <- sig.eClip.calls$rep1$AARS
  #y <- sig.eClip.calls$rep2$AARS
  
  hits <- findOverlaps(x, y)
  xhits <- x[queryHits(hits)]
  yhits <- y[subjectHits(hits)]
  # Fragments from both replicates must meet the "min" criterion
  if(both){
    overlap <- width(pintersect(xhits, yhits))
    frac1 <- overlap/ width(xhits)
    frac2 <- overlap/ width(yhits)
    frac1 <- frac1 >= min
    frac2 <- frac2 >= min
    merge <- frac1 & frac2
    # Only the smaller fragment must fulfill the min criterion
  }else{
    frac <- width(pintersect(xhits, yhits)) / pmin(width(xhits), width(yhits))
    merge <- frac >= min
  }
  # Return merged fragments
  if(!intersect){
    return(reduce(c(xhits[merge], yhits[merge])))
  }else{
    # Return only overlapping region of fragments
    return(pintersect(xhits, yhits)[merge])
  }
  
  # xhits[!merge], yhits[!merge],
  # x[-queryHits(hits)], y[-subjectHits(hits)])
}
###
# Removes a pattern in a character vector and optionally from names attribute
##
strip.version <- function(id, pattern = "\\..*", replacement = "", names = F){
  if(names & (!is.null(names(id)))){
    names(id)  <- gsub(pattern, replacement, names(id))
    return(gsub(pattern, replacement = replacement, id))
  }else{
    return(gsub(pattern, replacement = replacement, id))
  }
}
###
# Combine lists by list names
##
combine.lists <- function(L1, L2){
  element.names <- names(L1)
  new.list <- 
    lapply(names(L1), function(lname){
      unique(c(L1[[lname]], L2[[lname]]))
    })
  names(new.list) <- element.names
  new.list
}
###
# This function saves plots
##
save.plot <- function(dest, plot, ggsave = T){
  
  if(ggsave){
    ggsave(paste0(OUTPUT, PLOTS, dest,".pdf"), plot = plot, device = "pdf")
  }
  
  
}