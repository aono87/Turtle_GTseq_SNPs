# ADDED: Explicitly load dplyr for grouping and mutations.
library(dplyr) 

#
# Function: mplot2tgt_ratio
# Description: Converts microhaplotype data to a Tidy Genotype Table (TGT)
#              using allelic ratio for filtering.
#
# Arguments:
#   project: The name of the project.
#   AR.min.het: The MINIMUM allelic ratio for the minor allele to call a heterozygote.
#               Genotypes below this threshold will be marked as missing.
#   AR.max.homo: The MAXIMUM allelic ratio for the minor allele to still call a homozygote.
#                Minor alleles below this are discarded.
#   min.read.depth: The minimum total read depth for a genotype call.
#
mplot2tgt_ratio <- function(project, AR.min.het = 0.30, AR.max.homo = 0.20, 
                            min.read.depth = min.read.depth) {
  
  why.removed <- list()
  genos <- readRDS(paste0(out.path, project, ".rds"))[,-1]
  genos <- mutate(genos, id.loc = paste0(id, "-", locus))
  genos$to.remove <- FALSE
  
  # --- ADDED: Calculate Allelic Ratio ---
  # First, calculate the total read depth for each individual at each locus.
  locus_depths <- genos %>%
    group_by(id.loc) %>%
    summarise(total_locus_depth = sum(depth, na.rm = TRUE), .groups = 'drop')
  
  # Now, join the total depth back to the main table and calculate the ratio for each allele.
  genos <- genos %>%
    left_join(locus_depths, by = "id.loc") %>%
    mutate(allelic_ratio = ifelse(total_locus_depth > 0, depth / total_locus_depth, 0))
  # --- END ADDED SECTION ---

  # remove haplotypes that are NA
  genos$to.remove[which(is.na(genos$haplo))] <- TRUE
  why.removed$na <- length(which(genos$to.remove))
  
  # remove genotypes that don't meet the minimum read depth criterion
  inds <- unique(genos$id)
  loci <- unique(genos$locus)
  num.locs <- length(loci)
  # This part of the original script is slightly inefficient; we can just use the 'locus_depths' table we already made.
  genos.to.drop <- filter(locus_depths, total_locus_depth < min.read.depth)
  genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE
  why.removed$minDepth <- length(which(genos$to.remove))
  print("Done filtering on read depth")
  
  # remove genotypes for individuals whose read depth of its major haplotype at a locus is less than min.read.depth/2
  genos.to.drop <- filter(genos, depth < min.read.depth/2) %>% filter(rank == 1) %>% filter(to.remove == FALSE)
  genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE
  why.removed$majorAlleleDepth <- length(which(genos$to.remove))
  
  # --- CHANGED: Filter based on Allelic Ratio ---
  # remove minor haplotypes where the ratio is low enough to call a homozygote
  genos$to.remove[which(genos$rank > 1 & genos$allelic_ratio < AR.max.homo)] <- TRUE
  why.removed$AR.max.hom <- length(which(genos$to.remove))
  
  # remove entire genotypes with a minor allele ratio not high enough to call a heterozygote
  genos.to.drop <- filter(genos, rank > 1) %>% filter(to.remove == FALSE) %>% 
    filter(allelic_ratio < AR.min.het) %>% select(id.loc)
  genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE
  why.removed$AR.min.het <- length(which(genos$to.remove))
  # --- END CHANGED SECTION ---
  
  # remove genotypes with rank > 4
  genos$to.remove[which(genos$rank > 4)] <- TRUE
  why.removed$rank.gt.4 <- length(which(genos$to.remove))
  print("Done filtering on rank and allelic ratio")
  
  # create tgt-like data structure
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
    #create gt field
    if(df$to.remove[1]) gt <- NA else{
      if(nrow(df) == 1) gt <- paste(haps[1],haps[1], sep= "/") else {
        if(df$to.remove[2]) gt <- paste(haps[1],haps[1], sep= "/") else {
          gt <- ifelse(haps[2] < haps[1], paste(haps[2],haps[1], sep= "/"), paste(haps[1],haps[2], sep= "/"))
        }
      }
    }
    num.haps <- sum(!df$to.remove) 
    res <- c(loc, ind, gt, haps, depth, num.haps)
    names(res) <- c("locus", "Indiv", "gt", "haplo.1", "haplo.2", "haplo.3", "haplo.4", 
                    "depth.1", "depth.2", "depth.3", "depth.4", "num.haps")
    return(res)
  })))
  
  tgt$depth.1 <- as.integer(tgt$depth.1)
  tgt$depth.2 <- as.integer(tgt$depth.2)
  tgt$depth.3 <- as.integer(tgt$depth.3)
  tgt$depth.4 <- as.integer(tgt$depth.4)
  return(tgt)
}