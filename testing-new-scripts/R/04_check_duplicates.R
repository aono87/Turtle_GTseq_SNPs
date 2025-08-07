# 04_CHECK_DUPLICATES.R
#
# Accessory Step:
# This script uses the final `gtypes` object to search for pairs of
# individuals with highly similar genotypes, which may indicate that they
# are duplicate samples from the same individual.
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
# --- BEGIN CONFIGURATION ---
# ---------------------------

# --- Project and File Names ---
project.name <- "dcor.wpac.test"
min.reads <- 20 # Used for file naming to match original script

# --- Input Files ---
# Path to the gtypes .rda file created by `02_create_gtypes.R`
gtypes.path <- file.path("data/", paste0("gtypes_", project.name, "_minReads.", min.reads, ".rda"))

# --- Analysis Parameters ---
# Set the genetic similarity threshold for identifying potential duplicates.
# `dupGenotypes` uses the proportion of shared alleles. A value of 0.8 means
# two samples are flagged if they share 80% or more of their alleles.
# A value of 1.0 would only find perfect matches.
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
load(gtypes.path)


# ====================================================================
# STEP 2: IDENTIFY DUPLICATE GENOTYPES
# ====================================================================
message(paste("Step 2: Searching for duplicate genotypes with a similarity threshold of", similarity.threshold, "..."))

# The `dupGenotypes` function from strataG compares all pairs of individuals.
# It calculates the proportion of shared alleles between them.
# `num.shared` is the minimum proportion of shared alleles to be considered a duplicate.
duplicate.pairs <- dupGenotypes(g, num.shared = similarity.threshold)

# ====================================================================
# STEP 3: SAVE RESULTS
# ====================================================================
if (nrow(duplicate.pairs) > 0) {
  message(paste("Found", nrow(duplicate.pairs), "pair(s) of potential duplicates."))
  
  # Define the output filename
  output.filename <- file.path(
    results.raw.path,
    paste0(project.name, ".duplicate.genotypes.", similarity.threshold * 100, "pct.csv")
  )
  
  # Write the results to a CSV file
  write.csv(duplicate.pairs, output.filename, row.names = FALSE)
  
  message("Results saved to: ", output.filename)
  print(duplicate.pairs)
  
} else {
  message("No potential duplicate pairs found at the specified threshold.")
}

message("Duplicate check complete!")