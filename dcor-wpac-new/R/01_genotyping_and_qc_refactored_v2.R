# 01_genotyping_and_qc_refactored_v2.R
#
# This script uses a multi-part, interactive approach to filter GTseq data.
#   - PART 1: Generates initial data and runs a full QC analysis.
#   - PART 2: Applies a first, loose filter based on user input, then
#     re-runs the entire QC analysis on the intermediate dataset.
#   - PART 3: Applies a second, stricter filter to create the final dataset,
#     then re-runs the entire QC analysis on the final dataset.
#   - PART 4: Saves all final outputs.
# WORKFLOW DESCRIPTION:
#   - For the typical user that has known replicate samples, run PART 1 and PART 2 with all replicates separate. 
#   - After PART 2, evaluate the replicates and merge SAM files as needed. Generate a new sample label file accounting for merged SAM files. 
#   - Re-run PART 1 and PART 2 with merged replicates. Evaluate questionable haplotypes at this point and decide how to address them. 
#   - Run PART 3 and PART 4 after replicates and questionable haplotypes have been accounted for. 

# --- LOAD LIBRARIES ----
library(tidyverse)
library(vcfR)
library(microhaplot)
library(vcftoolsR)
library(openxlsx) # Added for writing multi-sheet Excel files

# --- LOAD CUSTOM FUNCTIONS ----
source("R/functions/mplot2tgt_SR.R")
source("R/functions/Compare.replicates.R")
source("R/functions/QC.helper.function_CURRENT.R") # Should contain the latest run_qc_analysis
source("R/functions/tracking_functions.R")      # Source for all tracking functions
#source("R/functions/zpad_ids.R") # Only if necessary

# ----------------------------------------------------------------------
# --- MASTER CONFIGURATION ----
# ----------------------------------------------------------------------

# ❗ STEP 1: SET THE CURRENT ANALYSIS STAGE ❗
# Choose one: "wreps", "merged", or "final_nodups"
analysis_stage <- "merged"

# --- Project and File Naming ---
project_name <- "dcor_wpac_new"
run_label <- paste(project_name, analysis_stage, sep = "_")

message("🚀 Starting workflow for run: ", run_label)

# --- Input Paths ---
sam_path <- "data-raw/sam_files/"
vcf_path <- "data-raw/vcf/Dcor_DcPanel_205_maf.targetSNPs_012224.recode.vcf"
app_path <- "~/Documents/GitHub/Shiny/microhaplot/"

# --- Dynamic Label Path Selection ---
label_path <- switch(analysis_stage,
                     "wreps"        = "data-raw/metadata/dcor.wpac.labels_zpad.txt",
                     "merged"       = "data-raw/metadata/dcor.wpac.labels_zpad_merged.txt",
                     "final_nodups" = "data-raw/metadata/dcor.wpac.labels_zpad_merged_nodups.txt"
)

# --- Output Paths ---
base_results_r_path <- "results-R"
base_results_raw_path <- "results-raw"
results_r_path <- file.path(base_results_r_path, analysis_stage)
results_raw_path <- file.path(base_results_raw_path, analysis_stage)
dir.create(results_r_path, showWarnings = FALSE, recursive = TRUE)
dir.create(results_raw_path, showWarnings = FALSE, recursive = TRUE)

# --- Genotyping Parameters (for mplot2tgt) ---
min_read_depth <- 20
ab_min_het <- 3 / 7
ab_max_homo <- 2 / 8
ar_max_homo <- 0.2

# --- Progress Tracking Initialization ---
# Initialize the summary tracker
progress_tracker <- tibble::tibble(
  Stage = character(), Loci_Count = integer(), Loci_Removed = integer(), Manual_Loci_Removed = character(),
  Individuals_Count = integer(), Individuals_Removed = integer(), Manual_Individuals_Removed = character(), Notes = character()
)

# Initialize the detailed, individual-level tracker
all_individuals_df <- read.delim(label_path, header = FALSE, stringsAsFactors = FALSE)
individual_tracker <- tibble::tibble(Indiv = all_individuals_df[[2]], analysis_stage = analysis_stage)

# Initialize the detailed, locus-level tracker (will be populated after initial data load)
locus_tracker <- tibble::tibble()

# Initialize the manual genotype change tracker
manual_genotype_changes_log <- tibble::tibble()

# ----------------------------------------------------------------------
# --- PART 1: INITIAL DATA GENERATION & QC ----
# ----------------------------------------------------------------------
message("\n--- Starting PART 1: Initial Data Generation & QC ---")

# Step 1A: Genotype and clean data
haplo_read_tbl <- prepHaplotFiles(
  run.label = run_label, sam.path = sam_path, out.path = results_r_path, label.path = label_path,
  vcf.path = vcf_path, app.path = app_path, n.jobs = parallel::detectCores()
)
mplot_rds_file <- file.path(results_r_path, paste0(run_label, ".rds"))

tgt <- mplot2tgt(
  project = run_label,
  results.mplot.file = mplot_rds_file,
  AB.min.het = ab_min_het,
  AB.max.homo = ab_max_homo,      
  AR.max.homo = ar_max_homo,      
  min.read.depth = min_read_depth 
)

tgt <- tgt %>%
  group_by(locus) %>% filter(!all(is.na(gt))) %>% ungroup() %>%
  group_by(Indiv) %>% filter(!all(is.na(gt))) %>% ungroup()

# Populate the locus_tracker with all unique loci from the initial data
all_loci_df <- unique(tgt$locus)
locus_tracker <- tibble::tibble(locus = all_loci_df, analysis_stage = analysis_stage)

message("Initial data loaded: ", n_distinct(tgt$Indiv), " individuals and ", n_distinct(tgt$locus), " loci.")
write_rds(tgt, file = file.path(results_r_path, paste0(run_label, "_initial_tgt.rds")))

# Step 1B: Run initial QC and update trackers with baseline state
qc_results_initial <- run_qc_analysis(
  tgt_df = tgt,
  stage_label = "initial",
  project = run_label,
  total_loci = n_distinct(tgt$locus),
  out_path = results_raw_path
)

# Log initial state to summary tracker
progress_tracker <- bind_rows(progress_tracker,
                              tibble(Stage = "01_Initial_Load", Loci_Count = n_distinct(tgt$locus), Individuals_Count = n_distinct(tgt$Indiv)))

# Log initial state to detailed individual tracker
ind_sum_initial <- qc_results_initial$ind_summary
individual_tracker <- individual_tracker %>%
  left_join(ind_sum_initial, by = "Indiv") %>%
  rename(initial_loci_genoed = loci.genoed, initial_prop_genoed = prop.genoed) %>%
  mutate(initial_status = if_else(is.na(initial_loci_genoed), "No Data", "Kept")) %>%
  relocate(initial_status, .after = last_col())

# Log initial state to detailed locus tracker
loc_sum_initial <- qc_results_initial$loc_summary
locus_tracker <- locus_tracker %>%
  left_join(loc_sum_initial, by = "locus") %>%
  rename(initial_inds_genoed = inds.genoed, initial_prop_genoed = prop.genoed) %>%
  mutate(initial_status = if_else(is.na(initial_inds_genoed), "No Data", "Kept")) %>%
  relocate(initial_status, .after = last_col())


# Step 1C: Save initial tracker state to a single Excel file
message("Saving initial tracker states to Excel...")
initial_trackers_list <- list(
  "Progress_Summary" = progress_tracker,
  "Individual_Tracker" = individual_tracker,
  "Locus_Tracker" = locus_tracker
)
write.xlsx(initial_trackers_list, file = file.path(results_raw_path, paste0(run_label, "_tracker_01_initial.xlsx")), na.string = "NA")


message("✅ PART 1 COMPLETE.")
message("🛑 ACTION REQUIRED: Examine reports and set LOOSE thresholds for Part 2.")

# ----------------------------------------------------------------------
# --- PART 2: FIRST PASS FILTERING & QC ----
# ----------------------------------------------------------------------
message("\n--- Starting PART 2: First Pass (Loose) Filtering ---")

# Step 2A: Define LOOSE filtering thresholds
pass1.locus.threshold <- 0.10
pass1.indiv.threshold <- 0.10
pass1.manual.loci.to.remove <- c("locus019")
pass1.manual.inds.to.remove <- c("")

# Step 2B: Apply the first filter
loc_sum_initial <- qc_results_initial$loc_summary
ind_sum_initial <- qc_results_initial$ind_summary

loci_to_keep_pass1 <- loc_sum_initial %>% filter(prop.genoed >= pass1.locus.threshold, !locus %in% pass1.manual.loci.to.remove) %>% pull(locus)
inds_to_keep_pass1 <- ind_sum_initial %>% filter(prop.genoed >= pass1.indiv.threshold, !Indiv %in% pass1.manual.inds.to.remove) %>% pull(Indiv)

tgt_before_pass1 <- tgt
tgt <- tgt %>% filter(locus %in% loci_to_keep_pass1, Indiv %in% inds_to_keep_pass1)

# Step 2C: Update all trackers with a single consolidated function
updated_trackers_pass1 <- update_all_trackers(
  progress_tracker = progress_tracker,
  individual_tracker = individual_tracker,
  locus_tracker = locus_tracker,
  tgt_before = tgt_before_pass1,
  tgt_after = tgt,
  ind_summary_before = ind_sum_initial,
  loc_summary_before = loc_sum_initial,
  stage_prefix = "pass1",
  indiv_threshold = pass1.indiv.threshold,
  locus_threshold = pass1.locus.threshold,
  manual_inds_remove = pass1.manual.inds.to.remove,
  manual_loci_remove = pass1.manual.loci.to.remove,
  previous_ind_status_col = "initial_status",
  previous_locus_status_col = "initial_status"
)
progress_tracker <- updated_trackers_pass1$progress_tracker
individual_tracker <- updated_trackers_pass1$individual_tracker %>%
  relocate(pass1_status, .after = last_col())
locus_tracker <- updated_trackers_pass1$locus_tracker %>%
  relocate(pass1_status, .after = last_col())

# Save pass1 tracker state to a single Excel file
message("Saving pass1 tracker states to Excel...")
pass1_trackers_list <- list(
  "Progress_Summary" = progress_tracker,
  "Individual_Tracker" = individual_tracker,
  "Locus_Tracker" = locus_tracker
)
write.xlsx(pass1_trackers_list, file = file.path(results_raw_path, paste0(run_label, "_tracker_02_pass1.xlsx")), na.string = "NA")

# Step 2D: Save intermediate files and re-run QC
write_rds(tgt, file = file.path(results_r_path, paste0(run_label, "_pass1_tgt.rds")))
qc_results_pass1 <- run_qc_analysis(
  tgt_df = tgt,
  stage_label = "pass1",
  project = run_label,
  total_loci = n_distinct(tgt$locus),
  out_path = results_raw_path
)

message("✅ PART 2 COMPLETE.")
message("🛑 ACTION REQUIRED: Examine new reports and set FINAL thresholds for Part 3.")

# ----------------------------------------------------------------------
# --- PART 3: FINAL (STRICT) FILTERING & QC ----
# ----------------------------------------------------------------------
message("\n--- Starting PART 3: Final (Strict) Filtering ---")

# Step 3A: Define FINAL filtering thresholds
final.locus.threshold <- 0.50
final.indiv.threshold <- 0.60
final.manual.loci.to.remove <- c("locus016")
final.manual.inds.to.remove <- c("")

# Optional: Manually change genotypes
genos_to_change_path <- 'data-raw/genos_to_change.csv'
if (file.exists(genos_to_change_path)) {
  # Read the CSV, ensuring strings are not converted to factors and empty strings are read as NA.
  genos_to_change <- read.csv(genos_to_change_path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  message("Applying ", nrow(genos_to_change), " manual genotype changes.")
  
  # Check if the required 'gt' column exists for the new genotypes.
  if("gt" %in% names(genos_to_change)){
    # Rename 'gt' to 'new_gt' for clarity in the join.
    genos_to_change <- dplyr::rename(genos_to_change, new_gt = gt)
    
    # Log the changes before applying them.
    manual_genotype_changes_log <- tgt %>%
      dplyr::inner_join(genos_to_change, by = c("locus", "Indiv")) %>%
      # Filter for rows where the genotype is actually different, correctly handling NAs.
      dplyr::filter(gt != new_gt | (is.na(gt) & !is.na(new_gt)) | (!is.na(gt) & is.na(new_gt))) %>%
      # Select and rename columns for a clear log.
      dplyr::select(locus, Indiv, original_gt = gt, new_gt) %>%
      # Explicitly convert NA values to the string "NA" for clear reporting.
      dplyr::mutate(
        original_gt = ifelse(is.na(original_gt), "NA", as.character(original_gt)),
        new_gt = ifelse(is.na(new_gt), "NA", as.character(new_gt))
      )
    
    # Add a flag to the changes to identify which rows in 'tgt' need updating.
    genos_to_change_flagged <- genos_to_change %>% dplyr::mutate(..to_be_changed.. = TRUE)
    
    # Apply the changes to the main 'tgt' dataframe.
    tgt <- tgt %>%
      dplyr::left_join(genos_to_change_flagged, by = c("locus", "Indiv")) %>%
      # For rows that were flagged (i.e., they were in the change file), replace the old 'gt' with 'new_gt'.
      # This correctly handles changing genotypes TO NA.
      dplyr::mutate(gt = ifelse(is.na(..to_be_changed..), gt, new_gt)) %>%
      # Clean up by removing the temporary helper columns.
      dplyr::select(-new_gt, -..to_be_changed..)
    
  } else {
    warning("Manual genotype change file ('", genos_to_change_path, "') does not contain a 'gt' column. No changes will be logged or applied.")
  }
}

# Step 3B: Apply the final filter
loc_sum_pass1 <- qc_results_pass1$loc_summary
ind_sum_pass1 <- qc_results_pass1$ind_summary

#Create list of loci to remove
if (exists("final.manual.loci.to.remove") && any(final.manual.loci.to.remove != "")) {
  loci_to_keep_final <- loc_sum_pass1 %>%
    filter(prop.genoed >= final.locus.threshold, 
           !locus %in% final.manual.loci.to.remove) %>%
    pull(locus)
} else {
  loci_to_keep_final <- loc_sum_pass1 %>%
    filter(prop.genoed >= final.locus.threshold) %>%
    pull(locus)
}

tgt_final_loci <- filter(tgt, locus %in% loci_to_keep_final)

num_locs_final <- length(loci_to_keep_final)
ind_sum_final <- tgt_final_loci %>% filter(!is.na(gt)) %>% group_by(Indiv) %>%
  summarise(loci.genoed = n_distinct(locus), .groups = 'drop') %>%
  mutate(prop.genoed = loci.genoed / num_locs_final)

#create list of inds to remove
if (exists("final.manual.inds.to.remove") && any(final.manual.inds.to.remove != "")) {
  inds_to_keep_final <- ind_sum_final %>%
    filter(prop.genoed >= final.indiv.threshold, 
           !Indiv %in% final.manual.inds.to.remove) %>%
    pull(Indiv)
} else {
  inds_to_keep_final <- ind_sum_final %>%
    filter(prop.genoed >= final.indiv.threshold) %>%
    pull(Indiv)
}

tgt_before_final <- tgt
tgt_final <- filter(tgt_final_loci, Indiv %in% inds_to_keep_final)

# Step 3C: Update all trackers with a single consolidated function
updated_trackers_final <- update_all_trackers(
  progress_tracker = progress_tracker,
  individual_tracker = individual_tracker,
  locus_tracker = locus_tracker,
  tgt_before = tgt_before_final,
  tgt_after = tgt_final,
  ind_summary_before = ind_sum_final, # Using recalculated summary for this decision
  loc_summary_before = loc_sum_pass1, # Using summary from previous stage
  stage_prefix = "final",
  indiv_threshold = final.indiv.threshold,
  locus_threshold = final.locus.threshold,
  manual_inds_remove = final.manual.inds.to.remove,
  manual_loci_remove = final.manual.loci.to.remove,
  previous_ind_status_col = "pass1_status",
  previous_locus_status_col = "pass1_status"
)

progress_tracker <- updated_trackers_final$progress_tracker
individual_tracker <- updated_trackers_final$individual_tracker %>%
  relocate(final_status, .after = last_col())
locus_tracker <- updated_trackers_final$locus_tracker %>%
  relocate(final_status, .after = last_col())

# Save final filter tracker state to a single Excel file
message("Saving final filter tracker states to Excel...")
final_trackers_list <- list(
  "Progress_Summary" = progress_tracker,
  "Individual_Tracker" = individual_tracker,
  "Locus_Tracker" = locus_tracker
)
if (nrow(manual_genotype_changes_log) > 0) {
  final_trackers_list[["Manual_Genotype_Changes"]] <- manual_genotype_changes_log
}
write.xlsx(final_trackers_list, file = file.path(results_raw_path, paste0(run_label, "_tracker_03_final_filter.xlsx")), na.string = "NA")


# Step 3D: Re-run QC on the final dataset
run_qc_analysis(
  tgt_df = tgt_final,
  stage_label = "final",
  project = run_label,
  total_loci = n_distinct(tgt_final$locus),
  out_path = results_raw_path
)

# ----------------------------------------------------------------------
# --- PART 4: SAVE FINAL OUTPUTS ----
# ----------------------------------------------------------------------
message("\n--- Starting PART 4: Saving Final Outputs ---")

# Step 4A: Convert final data to geno.table format
geno.table <- tgt.2.geno.table(tgt_final)

# Step 4B: Save final R objects
save(geno.table, tgt_final, file = file.path(results_r_path, paste0(run_label, "_final_data.rda")))

# Step 4C: Save all trackers to a single final Excel file
message("Saving final tracking reports to a single Excel file...")
final_summary_list <- list(
  "Progress_Summary" = progress_tracker,
  "Individual_Tracker" = individual_tracker,
  "Locus_Tracker" = locus_tracker
)
if (nrow(manual_genotype_changes_log) > 0) {
  final_summary_list[["Manual_Genotype_Changes"]] <- manual_genotype_changes_log
}
write.xlsx(final_summary_list, file = file.path(results_raw_path, paste0(run_label, "_filtering_summary_report.xlsx")), na.string = "NA")


message("\n✅ WORKFLOW COMPLETE!")
message("Final dataset contains ", nrow(geno.table), " individuals and ", ncol(geno.table) - 1, " loci.")
message("All R objects and raw reports saved in their respective stage folders inside '", base_results_r_path, "' and '", base_results_raw_path, "'")

