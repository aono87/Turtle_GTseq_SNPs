# --- QC Analysis Helper Function ---
# This function runs all QC steps and saves reports for a given dataset stage
# tgt_df=current tgt file
# stage_label=stage that you are in ("initial", "pass1", "final", etc)
# project=project name used for output file prefix
# total_loci=total number of loci at this particular stage ("num.locs.initial", etc)

run_qc_analysis <- function(tgt_df, stage_label, project, total_loci) {
  message(paste("\n--- Running QC Analysis for Stage:", stage_label, "---"))
  
  # 1. Questionable Haplotypes
  questionable.hap <- sapply(1:nrow(tgt_df), function(i) {
    (grepl("N", tgt_df$gt[i]) || grepl("X", tgt_df$gt[i]) || tgt_df$num.haps[i] > 2)
  })
  genos.to.check <- filter(tgt_df, questionable.hap)
  if (nrow(genos.to.check) > 0) {
    message(paste("  - Found", nrow(genos.to.check), "questionable genotypes."))
    write.csv(genos.to.check, file.path(results.raw.path, paste0(project, ".genos.to.check.", stage_label, ".csv")))
  } else {
    message("  - No questionable genotypes found.")
  }
  
  # 2. Replicate Mismatches
  LABIDs <- unique(tgt_df$Indiv) %>% substr(start = 1, stop = 8)
  replicates <- LABIDs[duplicated(LABIDs)]
  if (length(replicates) > 0) {
    mismatches.to.check <-do.call('rbind',lapply(replicates, function(r){
      rep.tgt <- tgt_df[grep(substr(r, start = 1, stop = 8), tgt_df$Indiv),]
      compare.replicates(rep.tgt)
    }))
    # KKM
    # It would be helpful to add here some code to remove duplicate rows from
    # mismatches.to.check and add a column with LABIDs (no a, b, c, etc). Here's
    # some untested code to do that:
    # mismatches.to.check <- mismatches.to.check |> 
    #   distinct() |> 
    #   mutate(LABID = substr(Indiv, start = 1, stop =  8))
    if (!is.null(mismatches.to.check) && nrow(mismatches.to.check) > 0) {
      message(paste("  - Found", nrow(mismatches.to.check), "mismatched genotypes between replicates."))
      write.csv(mismatches.to.check, file.path(results.raw.path, paste0(project, ".genotype.mismatches.", stage_label, ".csv")), row.names = FALSE)
    } else {
      message("  - No mismatched genotypes were found between replicates.")
    }
  } else {
    message("  - No replicate samples found to compare in this dataset.")
  }
  
  # 3. Locus and Individual Summaries
  loc.sum <- tgt_df %>% filter(!is.na(gt)) %>% group_by(locus) %>%
    summarise(inds.genoed = n_distinct(Indiv), .groups = 'drop') %>%
    mutate(prop.genoed = inds.genoed / n_distinct(tgt_df$Indiv))
  ind.sum <- tgt_df %>% filter(!is.na(gt)) %>% group_by(Indiv) %>%
    summarise(loci.genoed = n_distinct(locus), .groups = 'drop') %>%
    mutate(prop.genoed = loci.genoed / total_loci)
  
  write.csv(loc.sum, file.path(results.raw.path, paste0(project, ".locus.summary.", stage_label, ".csv")), row.names = FALSE)
  write.csv(ind.sum, file.path(results.raw.path, paste0(project, ".indiv.summary.", stage_label, ".csv")), row.names = FALSE)
  
  # 4. Histograms
  pdf(file.path(results.raw.path, paste0(project, ".qc.histograms.", stage_label, ".pdf")))
  hist(loc.sum$prop.genoed, breaks = 40, main = paste("Locus Completeness -", stage_label), xlab = "Proportion of Individuals Genotyped", col = "skyblue")
  hist(ind.sum$prop.genoed, breaks = 40, main = paste("Individual Completeness -", stage_label), xlab = "Proportion of Loci Genotyped", col = "salmon")
  dev.off()
  message(paste("  - QC reports saved with suffix '_", stage_label, ".csv/pdf'.", sep=""))
  message("  - QC reports saved for questionable genotypes ('genos.to.check'), locus summary, individual summary and a summary of mismatched genotypes for replicates (if present)")
}
