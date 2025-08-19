# 01_genotyping_and_qc_refactored_v2.R
#
# Workflow for a Two-Step Interactive Filter approach to genotype and filter gtseq genotypes:
# This script uses a multi-part, interactive approach to filter data.
#   - PART 1: Generates initial data and runs a full QC analysis.
#   - PART 2: Applies a first, loose filter based on user input, then
#     re-runs the entire QC analysis on the intermediate dataset.
#   - PART 3: Applies a second, stricter filter to create the final dataset,
#     then re-runs the entire QC analysis on the final dataset.
#   - PART 4: Saves all final outputs.

# --- LOAD LIBRARIES ----
library(tidyverse)
library(vcfR)
library(microhaplot)
library(vcftoolsR)

# --- LOAD CUSTOM FUNCTIONS ----
source("R/functions/mplot2tgt_ABO.R")
source("R/functions/tgt.2.geno.table.R")
source("R/functions/Compare.replicates.R")
source("R/functions/haps.w.ns.R")
source("R/functions/QC.helper.function.R") # call run_qc_analysis
source("R/functions/zpad_ids.R") # Only if necessary

# ----------------------------------------------------------------------
# --- MASTER CONFIGURATION ----
# ----------------------------------------------------------------------

# ❗ STEP 1: SET THE CURRENT ANALYSIS STAGE ❗
# Choose one: "wreps", "merged", or "final_nodups"
analysis_stage <- "final_nodups"

# --- Project and File Naming ---
project_name <- "dcor_wpac"
run_label <- paste(project_name, analysis_stage, sep = "_") # This creates names like "dcor_wpac_wreps"

message("🚀 Starting workflow for run: ", run_label)

# --- Input Paths ---
sam_path <- "data-raw/sam_files/"
vcf_path <- "data-raw/vcf/Dcor_DcPanel_205_maf.targetSNPs_012224.recode.vcf"
app_path <- "~/Documents/GitHub/Shiny/microhaplot/" # Path to shinyhaplot app

# --- Dynamic Label Path Selection ---
# The script automatically selects the correct label file based on 'analysis_stage'
label_path <- switch(analysis_stage,
                     "wreps"        = "data-raw/metadata/dcor.wpac.labels_zpad.txt",
                     "merged"       = "data-raw/metadata/dcor.wpac.labels_zpad_merged.txt",
                     "final_nodups" = "data-raw/metadata/dcor.wpac.labels_zpad_merged_nodups.txt"
)

# --- Output Paths ---
# Outputs will now be saved in stage-specific subfolders (e.g., results-R/wreps/)
base_results_r_path <- "results-R"
base_results_raw_path <- "results-raw"

results_r_path <- file.path(base_results_r_path, analysis_stage)
results_raw_path <- file.path(base_results_raw_path, analysis_stage)

dir.create(results_r_path, showWarnings = FALSE, recursive = TRUE)
dir.create(results_raw_path, showWarnings = FALSE, recursive = TRUE)


# --- Genotyping Parameters (for mplot2tgt) ---
min_read_depth <- 20
ab_min_het <- 3 / 7 # Allele balance for heterozygotes
ab_max_homo <- 2 / 8 # Allele balance for homozygotes

# ----------------------------------------------------------------------
# --- PART 1: INITIAL DATA GENERATION & QC ----
# ----------------------------------------------------------------------
message("\n--- Starting PART 1: Initial Data Generation & QC ---")

# Step 1A: Run microhaplot
# This single block works for all stages now.
haplo_read_tbl <- prepHaplotFiles(
  run.label = run_label,
  sam.path = sam_path,
  out.path = results_r_path,
  label.path = label_path,
  vcf.path = vcf_path,
  app.path = app_path,
  n.jobs = parallel::detectCores()
)

# Convert to a tidy genotype table (tgt)
tgt <- mplot2tgt(
  project = run_label,
  out.path = results_r_path,
  AB.min.het = ab_min_het,
  AB.max.homo = ab_max_homo,
  min.read.depth = min_read_depth
)

# Remove loci and individuals with no genotype data
tgt <- tgt %>%
  group_by(locus) %>%
  filter(!all(is.na(gt))) %>%
  ungroup() %>%
  group_by(Indiv) %>%
  filter(!all(is.na(gt))) %>%
  ungroup()

message("Initial data loaded: ", n_distinct(tgt$Indiv), " individuals and ", n_distinct(tgt$locus), " loci.")

# Save the initial tgt object
write_rds(tgt, file = file.path(results_r_path, paste0(run_label, "_initial_tgt.rds")))

# Step 1B: Run initial QC analysis
num_locs_initial <- n_distinct(tgt$locus)
run_qc_analysis(tgt, "initial", num_locs_initial, project = run_label, out_path = results_raw_path)

message("✅ PART 1 COMPLETE.")
message("🛑 ACTION REQUIRED: Examine reports in '", results_raw_path, "' with the '_initial' suffix.")
message("   Then, set your LOOSE thresholds for Part 2 below.")

# ----------------------------------------------------------------------
# --- PART 2: FIRST PASS FILTERING & QC ----
# ----------------------------------------------------------------------
message("\n--- Starting PART 2: First Pass (Loose) Filtering ---")

# Step 2A: Define LOOSE filtering thresholds (USER INPUT REQUIRED)
pass1.locus.threshold <- 0.10 # Keep loci genotyped in >10% of individuals
pass1.indiv.threshold <- 0.10 # Keep individuals genotyped at >10% of loci
pass1.manual.loci.to.remove <- c("locus019") # e.g., c("locus019", "locus042")
pass1.manual.inds.to.remove <- c("")

# Step 2B: Apply the first filter
# Load initial summaries
loc_sum_initial <- read.csv(file.path(results_raw_path, paste0(run_label, ".locus.summary.initial.csv")))
ind_sum_initial <- read.csv(file.path(results_raw_path, paste0(run_label, ".indiv.summary.initial.csv")))

# Get lists of what to keep
loci_to_keep_pass1 <- loc_sum_initial %>%
  filter(prop.genoed >= pass1.locus.threshold, !locus %in% pass1.manual.loci.to.remove) %>%
  pull(locus)

inds_to_keep_pass1 <- ind_sum_initial %>%
  filter(prop.genoed >= pass1.indiv.threshold, !Indiv %in% pass1.manual.inds.to.remove) %>%
  pull(Indiv)

# Filter the tgt object
n_loci_before <- n_distinct(tgt$locus)
n_inds_before <- n_distinct(tgt$Indiv)

tgt <- tgt %>%
  filter(locus %in% loci_to_keep_pass1) %>%
  filter(Indiv %in% inds_to_keep_pass1)

message("   - First pass removed ", n_loci_before - n_distinct(tgt$locus), " loci and ",
        n_inds_before - n_distinct(tgt$Indiv), " individuals.")

# Save the pass1 tgt object
write_rds(tgt, file = file.path(results_r_path, paste0(run_label, "_pass1_tgt.rds")))

# Step 2C: Re-run QC analysis
num_locs_pass1 <- n_distinct(tgt$locus)
run_qc_analysis(tgt, "pass1", num_locs_pass1, project = run_label, out_path = results_raw_path)

message("✅ PART 2 COMPLETE.")
message("🛑 ACTION REQUIRED: If you have replicates, check for mismatches and merge them now and re-run from the top with analysis_stage = 'merged'.")
message("   Otherwise, examine the new reports and set FINAL thresholds for Part 3.")

# ----------------------------------------------------------------------
# --- PART 3: FINAL (STRICT) FILTERING & QC ----
# ----------------------------------------------------------------------
message("\n--- Starting PART 3: Final (Strict) Filtering ---")

# Step 3A: Define FINAL filtering thresholds (USER INPUT REQUIRED)
final.locus.threshold <- 0.50 # Keep loci genotyped in >50% of individuals
final.indiv.threshold <- 0.60 # Keep individuals genotyped at >60% of loci
final.manual.loci.to.remove <- c("locus016")
final.manual.inds.to.remove <- c("")

# Optional: Manually change genotypes from a CSV file
# CSV must have columns: locus, Indiv, gt
genos_to_change_path <- 'data-raw/genos_to_change.csv'
if (file.exists(genos_to_change_path)) {
  genos_to_change <- read.csv(genos_to_change_path)
  message("Applying ", nrow(genos_to_change), " manual genotype changes.")
  
  # A more efficient `dplyr` way to do this join-update
  tgt <- tgt %>%
    left_join(genos_to_change, by = c("locus", "Indiv")) %>%
    mutate(gt = ifelse(!is.na(gt.y), gt.y, gt.x)) %>%
    select(-gt.x, -gt.y)
}

# Step 3B: Apply final filters
# Load pass1 summaries
loc_sum_pass1 <- read.csv(file.path(results_raw_path, paste0(run_label, ".locus.summary.pass1.csv")))

# Filter Loci first
loci_to_keep_final <- loc_sum_pass1 %>%
  filter(prop.genoed >= final.locus.threshold, !locus %in% final.manual.loci.to.remove) %>%
  pull(locus)

tgt_final_loci <- filter(tgt, locus %in% loci_to_keep_final)

# Now, recalculate individual genotyping rates based on the *final locus set*
num_locs_final <- length(loci_to_keep_final)
ind_sum_final <- tgt_final_loci %>%
  filter(!is.na(gt)) %>%
  group_by(Indiv) %>%
  summarise(loci_genoed = n_distinct(locus), .groups = 'drop') %>%
  mutate(prop.genoed = loci_genoed / num_locs_final)

inds_to_keep_final <- ind_sum_final %>%
  filter(prop.genoed >= final.indiv.threshold, !Indiv %in% final.manual.inds.to.remove) %>%
  pull(Indiv)

# Final filter
tgt_final <- filter(tgt_final_loci, Indiv %in% inds_to_keep_final)

message("   - Final filter resulted in ", n_distinct(tgt_final$locus), " loci and ",
        n_distinct(tgt_final$Indiv), " individuals.")

# Step 3C: Re-run QC on the final dataset
run_qc_analysis(tgt_final, "final", n_distinct(tgt_final$locus), project = run_label, out_path = results_raw_path)

# ----------------------------------------------------------------------
# --- PART 4: SAVE FINAL OUTPUTS ----
# ----------------------------------------------------------------------
message("\n--- Starting PART 4: Saving Final Outputs ---")

# Step 4A: Convert final data to geno.table format
geno.table <- tgt.2.geno.table(tgt_final)

# Step 4B: Save final R objects
# Renamed tgt_final to keep it distinct from the intermediate 'tgt' for clarity in the final save
save(geno.table, tgt_final, file = file.path(results_r_path, paste0(run_label, "_final_data.rda")))

message("\n✅ WORKFLOW COMPLETE!")
message("Final dataset contains ", nrow(geno.table), " individuals and ", ncol(geno.table) - 1, " loci.")
message("All R objects and raw reports saved in their respective stage folders inside '", base_results_r_path, "' and '", base_results_raw_path, "'")


# ----------------------------------------------------------------------
# --- NEXT STEPS: ----
# ----------------------------------------------------------------------

# Convert the genotable to a gtypes file using the script 02_create_gtypes.R. 
# Check for duplicates using the script 03_check_duplicates.R
# Address duplicates and generate a final sample list for microhaplot: e.g. "dcor.wpac.labels_zpad_merged_nodups.txt"
# Re-run this script with the merged sam files, updated sample list, and change stage to: "final_nodups"