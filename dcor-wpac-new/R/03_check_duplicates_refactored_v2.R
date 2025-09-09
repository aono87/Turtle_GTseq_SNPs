# 03_CHECK_DUPLICATES.R
#
# Workflow Step 3 (Accessory Step):
# This script uses the final `gtypes` object to search for pairs of
# individuals with highly similar genotypes, which may indicate that they
# are duplicate samples from the same individual.
#
# !! RUN THIS SCRIPT BEFORE THE MAIN POPULATION GENETIC ANALYSES !!
#
# INPUTS:
#   - A `gtypes` object (`.rda` file) from `02_create_gtypes.R`.
#
# OUTPUTS:
#   - A CSV file listing pairs of individuals exceeding the similarity threshold.

# --- LOAD LIBRARIES ---
library(tidyverse)
library(strataG)

# ----------------------------------------------------------------------
# --- CONFIGURATION (UPDATED) ---
# ----------------------------------------------------------------------

# --- Project and File Naming ---
# These settings MUST match script 02.
project_name <- "dcor_wpac"
final_analysis_stage <- "final_nodups" #update with stage: "wreps", "merged" or "final_nodups"
min_reads <- 20

# --- Auto-Generated Paths (Do not change) ---
# Recreates the run label from previous scripts to find the correct files.
run_label <- paste(project_name, final_analysis_stage, sep = "_")

# Path to the gtypes .rda file created by `02_create_gtypes.R`.
gtypes.path <- file.path("data", paste0("gtypes_", run_label, "_minReads", min_reads, ".rda"))

# Path for the output CSV file.
results_raw_path <- file.path("results-raw", final_analysis_stage)
dir.create(results_raw_path, showWarnings = FALSE, recursive = TRUE)

# --- Analysis Parameters ---
# Set the genetic similarity threshold for identifying potential duplicates.
similarity.threshold <- 0.80

# -------------------------
# --- END CONFIGURATION ---
# -------------------------


# ====================================================================
# STEP 1: LOAD GTYPES DATA
# ====================================================================
message("Step 1: Loading gtypes object...")

# Check if the gtypes file exists before proceeding.
if (!file.exists(gtypes.path)) {
  stop(paste("Gtypes file not found at:", gtypes.path,
             "\nPlease ensure script 02 has been run successfully with matching configuration."))
}
load(gtypes.path)


# ====================================================================
# STEP 2: IDENTIFY DUPLICATE GENOTYPES
# ====================================================================
message(paste("Step 2: Searching for duplicate genotypes with a similarity threshold of", similarity.threshold, "..."))

# The `dupGenotypes` function compares all pairs of individuals.
duplicate.pairs <- dupGenotypes(g, num.shared = similarity.threshold)

# ====================================================================
# STEP 3: SAVE RESULTS AND PROVIDE NEXT STEPS (UPDATED)
# ====================================================================
if (is.null(duplicate.pairs)) {
  message("\n✅ No potential duplicate pairs found at the specified threshold.")
} else {
  # If not NULL, duplicates were found and duplicate.pairs is a data frame.
  message(paste("Found", nrow(duplicate.pairs), "pair(s) of potential duplicates."))
  
  # Define the output filename using the run_label for consistency.
  output.filename <- file.path(
    results_raw_path,
    paste0(run_label, ".duplicate.genotypes.", similarity.threshold * 100, "pct.csv")
  )
  
  # Write the results to a CSV file
  write.csv(duplicate.pairs, output.filename, row.names = FALSE)
  
  message("Results saved to: ", output.filename)
  print(duplicate.pairs)
  
  message("   🛑 ACTION REQUIRED:")
  message("   Review the CSV file of duplicate pairs. If these are unintended duplicates, you must address them.")
  message("   Once addressed (merged, removed, etc) re-run the entire workflow from script '01_genotyping_and_qc...' to ensure all filtering is accurate.")
  message("   Change stage to 'final_nodups' and generate a final sample list ('dcor.wpac.labels_zpad_merged_nodups.txt')")
}

message("\n✅ Duplicate check complete!")
