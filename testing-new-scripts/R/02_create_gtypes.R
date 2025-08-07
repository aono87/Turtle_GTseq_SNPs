# 02_CREATE_GTYPES.R
#
# Workflow Step 2:
# This script loads the cleaned genotype table from Step 1, merges it with
# sample metadata, and creates a `gtypes` object from the `strataG` package.
# This object is the primary input for downstream population genetics analyses.
#
# INPUTS:
#   - The final `.rda` file from `01_genotyping_and_qc.R` (containing `geno.table`).
#   - A metadata file (Excel or CSV) with sample information and strata definitions.
#
# OUTPUTS:
#   - A `gtypes` object saved as an `.rda` file in the `data/` directory.

# --- LOAD LIBRARIES ---
library(tidyverse)
library(strataG)
library(readxl)

# ---------------------------
# --- BEGIN CONFIGURATION ---
# ---------------------------

# --- Project and File Names ---
project.name <- "dcor.wpac.test"
min.reads <- 20 # Used for file naming to match original script

# --- Input Files ---
# Path to the .rda file created by the `01_genotyping_and_qc.R` script
geno.data.path <- file.path("results-R/", paste0(project.name, ".final.geno.data.rda"))

# Path to your sample metadata file
# This file should contain a column that matches the individual IDs in your genotype data.
sample.info.path <- "data-raw/metadata/turtle.qa.qc.xlsx"
sample.info.sheet <- "SampleData" # Specify sheet name if using Excel

# --- Column Names ---
# The column in your metadata file that contains the individual IDs.
# This MUST match the IDs in the `geno.table` (e.g., "z0012526").
metadata.id.col <- "mplot.id"

# The column(s) in your metadata file to be used as strata.
# This will be the primary population definition for HWE/LD analyses.
strata.col.name <- "Stratum_ABO"

# -------------------------
# --- END CONFIGURATION ---
# -------------------------


# ====================================================================
# STEP 1: LOAD AND PREPARE DATA
# ====================================================================
message("Step 1: Loading genotype data and sample metadata...")

# Load the `geno.table` object from the previous script's output
load(geno.data.path)

# Read in sample metadata
sample.info <- read_xlsx(sample.info.path, sheet = sample.info.sheet)

# Check if required columns exist
if (!metadata.id.col %in% names(sample.info)) {
  stop(paste("Metadata ID column '", metadata.id.col, "' not found in sample info file."))
}
if (!strata.col.name %in% names(sample.info)) {
  stop(paste("Strata column '", strata.col.name, "' not found in sample info file."))
}

# ====================================================================
# STEP 2: MERGE GENOTYPES WITH METADATA
# ====================================================================
message("Step 2: Merging genotypes with metadata...")

# The `geno.table` has individuals as rows and loci as columns.
# The first column is named "Indiv".
# We will join the metadata to this table.

# Prepare the metadata by renaming the ID column for the join.
sample.info.renamed <- sample.info %>%
  rename(Indiv = all_of(metadata.id.col))

# The `geno.table` contains only the final, filtered set of individuals.
# An `inner_join` ensures we only keep metadata for these individuals.
df.merged <- inner_join(sample.info.renamed, geno.table, by = "Indiv")

# ====================================================================
# STEP 3: CREATE AND SAVE THE GTYPES OBJECT
# ====================================================================
message("Step 3: Creating the gtypes object...")

# The `df2gtypes` function requires a specific format.
# We need to specify which columns contain the ID, strata, and genotype data.
# Using column names is more robust than using column indices.

# Identify the column where the genetic loci start.
# We assume all columns after the metadata are loci.
first.locus.col <- which(!names(df.merged) %in% names(sample.info.renamed))[1]
message(paste("Locus data starts at column:", first.locus.col))

# Create the gtypes object
# We provide the full merged dataframe and tell the function where to find
# the individual IDs, the desired stratum, and where the locus data starts.
g <- df2gtypes(
  df.merged,
  ploidy = 2,
  id.col = "Indiv",
  strata.col = strata.col.name,
  loc.col = first.locus.col,
  sep = "/" # Specify the allele separator
)

# You can add more stratification schemes to the `gtypes` object if needed.
# For example, to add a scheme for 'Country':
# g <- addStrata(g, df.merged %>% select(Indiv, Country) %>% column_to_rownames("Indiv"))

print(g)
message("gtypes object created successfully.")

# Save the final `gtypes` object and the merged dataframe for inspection
output.filename <- file.path("data/", paste0("gtypes_", project.name, "_minReads.", min.reads, ".rda"))
save(g, df.merged, file = output.filename)

message("Workflow Step 2 complete!")
message("gtypes object saved to: ", output.filename)