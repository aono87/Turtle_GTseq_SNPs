# --- QC Analysis Helper Function (Returns a List of All Summary Data) ---
# This function runs all QC steps, saves reports, and returns a comprehensive list
# containing all generated summary data frames.

run_qc_analysis <- function(tgt_df, stage_label, project, total_loci, out_path) {
  message(paste("\n--- Running QC Analysis for Stage:", stage_label, "---"))
  
  # NEW: Initialize a list to hold all summary data frames.
  qc_summaries <- list()
  
  # 1. Questionable Haplotypes
  questionable.hap <- (grepl("N", tgt_df$gt) | grepl("X", tgt_df$gt) | tgt_df$num.haps > 2)
  
  # UPDATED: Only create the 'genos.to.check' df if any questionable haplotypes are found.
  if (any(questionable.hap)) {
    genos.to.check <- dplyr::filter(tgt_df, questionable.hap)
    message(paste("  - Found", nrow(genos.to.check), "questionable genotypes."))
    write.csv(genos.to.check, file.path(out_path, paste0(project, ".genos.to.check.", stage_label, ".csv")))
    # NEW: Add the data frame to our summary list.
    qc_summaries$questionable_genotypes <- genos.to.check
  } else {
    message("  - No questionable genotypes found.")
  }
  
  # 2. Replicate Mismatches
  LABIDs <- unique(tgt_df$Indiv) %>% substr(start = 1, stop = 8)
  replicates <- LABIDs[duplicated(LABIDs)]
  if (length(replicates) > 0) {
    mismatches.to.check <- do.call('rbind', lapply(replicates, function(r) {
      rep.tgt <- tgt_df[grep(substr(r, start = 1, stop = 8), tgt_df$Indiv), ]
      compare.replicates(rep.tgt)
    }))
    
    if (!is.null(mismatches.to.check) && nrow(mismatches.to.check) > 0) {
      mismatches.to.check <- mismatches.to.check %>%
        dplyr::distinct() %>%
        dplyr::mutate(LABID = substr(Indiv, start = 1, stop = 8))
      
      message(paste("  - Found", nrow(mismatches.to.check), "mismatched genotypes between replicates."))
      write.csv(mismatches.to.check, file.path(out_path, paste0(project, ".genotype.mismatches.", stage_label, ".csv")), row.names = FALSE)
      # NEW: Add the data frame to our summary list.
      qc_summaries$replicate_mismatches <- mismatches.to.check
    } else {
      message("  - No mismatched genotypes were found between replicates.")
    }
  } else {
    message("  - No replicate samples found to compare in this dataset.")
  }
  
  # 3. Locus and Individual Summaries
  loc.sum <- tgt_df %>% 
    dplyr::filter(!is.na(gt)) %>% 
    dplyr::group_by(locus) %>%
    dplyr::summarise(inds.genoed = dplyr::n_distinct(Indiv), .groups = 'drop') %>%
    dplyr::mutate(prop.genoed = inds.genoed / dplyr::n_distinct(tgt_df$Indiv))
  
  ind.sum <- tgt_df %>% 
    dplyr::filter(!is.na(gt)) %>% 
    dplyr::group_by(Indiv) %>%
    dplyr::summarise(loci.genoed = dplyr::n_distinct(locus), .groups = 'drop') %>%
    dplyr::mutate(prop.genoed = loci.genoed / total_loci)
  
  write.csv(loc.sum, file.path(out_path, paste0(project, ".locus.summary.", stage_label, ".csv")), row.names = FALSE)
  write.csv(ind.sum, file.path(out_path, paste0(project, ".indiv.summary.", stage_label, ".csv")), row.names = FALSE)
  
  # NEW: Add these essential summaries to the list, maintaining original names for compatibility.
  qc_summaries$loc_summary <- loc.sum
  qc_summaries$ind_summary <- ind.sum
  
  # 4. Histograms
  pdf(file.path(out_path, paste0(project, ".qc.histograms.", stage_label, ".pdf")))
  hist(loc.sum$prop.genoed, breaks = 40, main = paste("Locus Completeness -", stage_label), xlab = "Proportion of Individuals Genotyped", col = "skyblue")
  hist(ind.sum$prop.genoed, breaks = 40, main = paste("Individual Completeness -", stage_label), xlab = "Proportion of Loci Genotyped", col = "salmon")
  dev.off()
  
  message(paste("  - QC reports saved with suffix '_", stage_label, ".csv/pdf'.", sep=""))
  message("  - QC reports saved for questionable genotypes, locus/individual summaries, and replicate mismatches (if present).")
  
  # Return the comprehensive list of all summary data frames.
  return(qc_summaries)
}