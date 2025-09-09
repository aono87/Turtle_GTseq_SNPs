# 02_CREATE_GTYPES.R
#
# Workflow Step 2:
# This script loads the cleaned genotype table from Step 1, merges it with
# de-duplicated metadata, reshapes the genotype data, and creates a `gtypes`
# object with a named strata scheme.
#
# INPUTS:
#   - The final '.rda' file from Step 1 (containing the `geno.table` object).
#   - A metadata file (e.g., Excel or CSV) with sample information.
#
# OUTPUTS:
#   - A `gtypes` object saved as an '.rda' file in the 'data/' directory.

# --- LOAD LIBRARIES ---
library(tidyverse)
library(strataG)
library(swfscMisc)
library(readxl)

# ----------------------------------------------------------------------
# --- CONFIGURATION (UPDATED) ---
# ----------------------------------------------------------------------

# --- Project and File Naming ---
# These settings MUST match the final run of script 01.
project_name <- "dcor_wpac"
final_analysis_stage <- "final_nodups" # The stage of the data you want to load.
#"wreps", "merged", or "final_nodups"

# --- Auto-Generated Paths (Do not change) ---
# Recreates the run label from script 01 to find the correct file.
run_label <- paste(project_name, final_analysis_stage, sep = "_")
min_reads <- 20 # Assuming this is consistent with your project.

# Path to the .rda file created by script 01.
geno.data.path <- file.path("results-R", final_analysis_stage, paste0(run_label, "_final_data.rda"))

# --- Metadata Configuration (USER INPUT REQUIRED) ---
# Path to your sample metadata file.
sample.info.path <- "dcor.wpac.qa-qc.xlsx"
sample.info.sheet <- "SampleData"

# The column in your metadata that contains the individual IDs.
metadata.id.col <- "id"
# The column in your metadata to be used as the stratum.
#strata.col.name <- "Stratum_ABO"
strata.col.name <- "Stratum2_ABO" #for Fst tables

# -------------------------
# --- END CONFIGURATION ---
# -------------------------


# ====================================================================
# STEP 1: LOAD AND PREPARE DATA
# ====================================================================
message("Step 1: Loading genotype and metadata files...")

# Load the `geno.table` object from the previous script's output.
if (!file.exists(geno.data.path)) {
  stop(paste("Genotype data file not found at:", geno.data.path, "\nPlease check your configuration settings."))
}
load(geno.data.path)

# Read in sample metadata.
sample.info <- read_xlsx(sample.info.path, sheet = sample.info.sheet)

# Rename the ID column in the metadata to "Indiv" to facilitate merging.
sample.info.renamed <- sample.info %>%
  rename(Indiv = all_of(metadata.id.col))

#Check for and remove duplicate individuals from the metadata.
message("Checking for duplicate individuals in metadata...")
n_before <- nrow(sample.info.renamed)
sample.info.deduped <- sample.info.renamed %>%
  distinct(Indiv, .keep_all = TRUE) # Keeps the first unique entry for each ID
n_after <- nrow(sample.info.deduped)

if (n_before > n_after) {
  message(paste("NOTE:", n_before - n_after, "duplicate individual ID(s) found in metadata and were removed."))
}


# ====================================================================
# STEP 2: RESHAPE GENOTYPE DATA USING alleleSplit()
# ====================================================================
message("Step 2: Reshaping genotype data with alleleSplit()...")

# The df2gtypes function requires one column per allele.
# First, set the 'Indiv' column as row names for alleleSplit.
genos <- column_to_rownames(geno.table, var = "Indiv")

# Now, split the genotypes (e.g., "A/T") into two separate columns ("A", "T").
split.genos <- alleleSplit(genos, sep = "/") %>%
  data.frame() %>%
  rownames_to_column(var = "Indiv") # Convert row names back to an 'Indiv' column.


# ====================================================================
# STEP 3: CREATE THE GTYPES OBJECT
# ====================================================================
message("Step 3: Merging data and creating the gtypes object...")

# Join the de-duplicated metadata with the split genotype data.
df <- right_join(sample.info.deduped, split.genos, by = "Indiv")

# Prepare the data frame for the `schemes` argument.
# This defines the named stratification scheme(s).
df.schemes <- select(df, Indiv, all_of(strata.col.name)) %>%
  column_to_rownames(var = "Indiv")

# Identify the column where the locus data begins.
# This makes the script robust if you add/remove metadata columns.
first.locus.col <- which(!names(df) %in% names(sample.info.renamed))[1]

# Create the gtypes object.
g <- df2gtypes(
  df,
  ploidy = 2,
  id.col = "Indiv",
  strata.col = strata.col.name,
  loc.col = first.locus.col,
  schemes = df.schemes
)

print(g)


# ====================================================================
# STEP 4: SAVE FINAL OUTPUT (UPDATED)
# ====================================================================

# Save the final `gtypes` object and the merged dataframe.
# This filename will now match the input data stage (e.g., gtypes_dcor_wpac_final_nodups_...)
output.filename <- file.path("data", paste0("gtypes_", run_label, "_minReads", min_reads, "_Stratum2", ".rda"))
save(g, df, file = output.filename)

message("\n✅ Workflow Step 2 complete!")
message("   gtypes object saved to: ", output.filename)

