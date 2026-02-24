# 02_CREATE_GTYPES.R-REDUCED
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

geno.table #from final filtering step
# --- Metadata Configuration (USER INPUT REQUIRED) ---
# Path to your sample metadata file.
sample.info <- "qaqc.full-final-dataset.xlsx"
sample.info.sheet <- "sample_data"

# The column in your metadata that contains the individual IDs.
metadata.id.col <- "mplot.id"
# The column in your metadata to be used as the stratum.
strata.col.name <- "Stratum_ABO"
#strata.col.name <- "Stratum2_ABO" #for Fst tables

# -------------------------
# --- END CONFIGURATION ---
# -------------------------


# ====================================================================
# STEP 1: LOAD AND PREPARE DATA
# ====================================================================


# Read in sample metadata.
# Path to your sample metadata file.
sample.info <- "qaqc.full-final-dataset.xlsx" #Name of your sample info file
sample.info.sheet <- "sample_data" #Name of the sheet with the sample info data

# Rename the ID column in the metadata to "Indiv" to facilitate merging. This is specific to my file, just need to make sure it fits your file format
sample.info.renamed <- sample.info %>%
  rename(Indiv = all_of(metadata.id.col))

#Optional: Check for and remove duplicate individuals from the metadata.
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

#Update this with your filename

save(g, df, file = output.filename)

message("\n✅ Workflow Step 2 complete!")


