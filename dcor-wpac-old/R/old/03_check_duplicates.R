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

# ---------------------------
# --- CONFIGURATION ---
# ---------------------------

# --- Project and File Naming ---
# This 'project' variable MUST match the one from `02_create_gtypes.R`.
project <- "dcor_wpac_final"
min.reads <- 20 # Used for file naming to match previous scripts.

# --- Input File ---
# Path to the gtypes .rda file created by `02_create_gtypes.R`.
gtypes.path <- file.path("data", paste0("gtypes_", project, "_minReads", min.reads, ".rda"))

# --- Analysis Parameters ---
# Set the genetic similarity threshold for identifying potential duplicates.
# A value of 0.8 means two samples are flagged if they share 80% or more of their alleles.
similarity.threshold <- 0.80

# --- Output Paths ---
results.raw.path <- "results-raw/"

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
             "\nPlease ensure script 02 has been run successfully."))
}
load(gtypes.path)


# ====================================================================
# STEP 2: IDENTIFY DUPLICATE GENOTYPES
# ====================================================================
message(paste("Step 2: Searching for duplicate genotypes with a similarity threshold of", similarity.threshold, "..."))

# The `dupGenotypes` function compares all pairs of individuals.
duplicate.pairs <- dupGenotypes(g, num.shared = similarity.threshold)

# ====================================================================
# STEP 3: SAVE RESULTS AND PROVIDE NEXT STEPS
# ====================================================================
if (is.null(duplicate.pairs)) {
  message("No potential duplicate pairs found at the specified threshold.")
} else {
  # If not NULL, duplicates were found and duplicate.pairs is a data frame.
  message(paste("Found", nrow(duplicate.pairs), "pair(s) of potential duplicates."))
  
  # Define the output filename
  output.filename <- file.path(
    results.raw.path,
    paste0(project, ".duplicate.genotypes.", similarity.threshold * 100, "pct.csv")
  )
  
  # Write the results to a CSV file
  write.csv(duplicate.pairs, output.filename, row.names = FALSE)
  
  message("Results saved to: ", output.filename)
  print(duplicate.pairs)
  
  message("\n🛑 ACTION REQUIRED:")
  message("   Review the CSV file of duplicate pairs. Decide which individual from each pair to remove.")
  message("   You can either:")
  message("   1. Go back to your metadata file, remove the unwanted individuals, and re-run `02_create_gtypes.R`.")
  message("   2. Manually remove them from the `g` object in the next script (`04_population_genetics_analysis.R`).")
}

message("\n✅ Duplicate check complete!")
