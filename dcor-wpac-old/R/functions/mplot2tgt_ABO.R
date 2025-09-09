# Source required functions
source("R/functions/tgt.2.geno.table.R")
source("R/functions/haps.w.ns.R")

# CHANGED: Added library call to make the script's dependency on 'dplyr' explicit.
library(dplyr)

# Define the function to convert mplot output to a TGT-like genotype table
# CHANGED: The function signature was modified to make all key parameters required arguments.
# This makes the function more robust and self-contained.
mplot2tgt <- function(project, out.path, AB.min.het, AB.max.homo, min.read.depth) {
  
  why.removed <- list()
  
  # CHANGED: Replaced paste0() with file.path() for a more reliable way to build file paths
  # that works on any operating system.
  geno_file <- file.path(out.path, paste0(project, ".rds"))
  genos <- readRDS(geno_file)[,-1]
  
  genos <- mutate(genos, id.loc = paste0(id, "-", locus))
  genos$to.remove <- FALSE
  
  # Remove haplotypes that are NA
  genos$to.remove[which(is.na(genos$haplo))] <- TRUE
  why.removed$na <- length(which(genos$to.remove))
  
  # Remove genotypes that don't meet the minimum read depth criterion
  inds <- unique(genos$id)
  loci <- unique(genos$locus)
  num.locs <- length(loci)
  
  tot.depth <- data.frame(do.call(rbind, lapply(inds, function(i){
    df <- filter(genos, id == i) %>% filter(to.remove == FALSE)
    do.call(rbind, lapply(loci, function(l){
      haps <- filter(df, locus == l)
      data.frame(i, l, sum(haps$depth))
    }))
  })))
  names(tot.depth) <- c("id", "locus", "tot.depth")
  
  genos.to.drop <- filter(tot.depth, tot.depth < min.read.depth) %>% mutate(id.loc = paste0(id, "-", locus))
  genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE
  why.removed$minDepth <- length(which(genos$to.remove))
  print("Done filtering on read depth")
  
  # Remove genotypes for individuals whose major haplotype read depth is less than min.read.depth/2
  genos.to.drop <- filter(genos, depth < min.read.depth/2) %>% filter(rank == 1) %>% filter(to.remove == FALSE)
  genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE
  why.removed$majorAlleleDepth <- length(which(genos$to.remove))
  
  # Remove minor haplotypes (rank > 1) where the allelic balance suggests a homozygote
  genos$to.remove[which(genos$rank > 1 & genos$allele.balance < AB.max.homo)] <- TRUE
  why.removed$AB.max.hom <- length(which(genos$to.remove))
  
  # Remove genotypes where the minor allele frequency is too low to call a heterozygote
  genos.to.drop <- filter(genos, rank > 1) %>% filter(to.remove == FALSE) %>% 
    filter(allele.balance < AB.min.het) %>% select(id.loc)
  genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE
  why.removed$AB.min.het <- length(which(genos$to.remove))
  
  # Remove haplotypes with rank > 4
  genos$to.remove[which(genos$rank > 4)] <- TRUE
  why.removed$rank.gt.4 <- length(which(genos$to.remove))
  print("Done filtering on rank and allelic balance")
  
  # Create the TGT-like data structure
  tgt <- data.frame(do.call(rbind, lapply(unique(genos$id.loc), function(x){
    df <- filter(genos, id.loc == x)
    if (nrow(df) > 4) df <- df[1:4,]
    
    loc <- df$locus[1]
    ind <- df$id[1]
    haps <- df$haplo
    depth <- df$depth
    
    if(length(haps) < 4) {
      haps <- c(haps, rep(NA, (4-length(haps))))
      depth <- c(depth, rep(0, (4-length(depth))))
    }
    
    # CHANGED: The nested if/else statements were reformatted for better readability.
    # The logic is identical but easier to understand.
    if(df$to.remove[1]) {
      gt <- NA 
    } else {
      if(nrow(df) == 1 || df$to.remove[2]) {
        gt <- paste(haps[1], haps[1], sep= "/")
      } else {
        gt <- ifelse(haps[2] < haps[1], paste(haps[2], haps[1], sep= "/"), paste(haps[1], haps[2], sep= "/"))
      }
    }
    
    num.haps <- sum(!df$to.remove)
    res <- c(loc, ind, gt, haps, depth, num.haps)
    names(res) <- c("locus", "Indiv", "gt", "haplo.1", "haplo.2", "haplo.3", "haplo.4", 
                    "depth.1", "depth.2", "depth.3", "depth.4", "num.haps")
    return(res)
  })))
  
  # CHANGED: The four repetitive lines for type conversion were replaced with a single,
  # more efficient line that applies the change to all 'depth' columns at once.
  tgt[, grep("depth", names(tgt))] <- lapply(tgt[, grep("depth", names(tgt))], as.integer)
  
  # return(list(tgt = tgt, why.removed = why.removed)) # Optional detailed return
  return(tgt)
}