# 06_POPULATION_STRUCTURE_ANALYSIS.R
#
# Workflow Step 6:
# This script assesses population structure by calculating pairwise Fst and Chi-squared
# statistics between strata. It is designed to follow the data preparation and
# population genetics analyses in previous workflow steps.
#
# INPUTS:
#   - A `gtypes` object (`.rda` file) from `02_create_gtypes.R`.
#
# OUTPUTS:
#   - A scatter plot of missing genotypes vs. homozygosity.
#   - CSV files containing pairwise Fst and Chi-squared matrices.
#   - An RDA file saving key objects from the analysis.

# --- LOAD LIBRARIES ---
library(strataG)
library(genepop)
library(tidyverse)
library(ggplot2)

# ----------------------------------------------------------------------
# --- CONFIGURATION ---
# ----------------------------------------------------------------------

# --- Project and File Naming ---
# These settings MUST match your other scripts.
project_name <- "dcor_wpac"
final_analysis_stage <- "final_nodups" 
min_reads <- 20

# --- Auto-Generated Paths (Do not change) ---
run_label <- paste(project_name, final_analysis_stage, sep = "_")
#gtypes_path <- file.path("data", paste0("gtypes_", run_label, "_minReads", min_reads, ".rda"))
gtypes_path<- "data/gtypes_dcor_wpac_final_nodups_minReads20_Stratum2.rda"

# --- Output Paths ---
results_raw_path <- file.path("results-raw", final_analysis_stage)
results_r_path <- file.path("results-R", final_analysis_stage)
dir.create(results_raw_path, showWarnings = FALSE, recursive = TRUE)
dir.create(results_r_path, showWarnings = FALSE, recursive = TRUE)

# --- Analysis Parameters ---
# Stratification scheme to use for the analysis.
stratification_scheme <- "Stratum2_ABO" 
# Minimum number of individuals required per stratum for pairwise tests.
min_inds_for_pairwise <- 7
# Number of repetitions for p-value computation in pairwise tests.
num_reps_pvals <- 1000

# -------------------------
# --- END CONFIGURATION ---
# -------------------------


# ====================================================================
# STEP 1: LOAD DATA AND PREPARE STRATA
# ====================================================================
message("Step 1: Loading gtypes object and preparing strata...")

# Check if the gtypes file exists before proceeding.
if (!file.exists(gtypes_path)) {
  stop(paste("Gtypes file not found at:", gtypes_path,
             "\nPlease ensure script 02 has been run successfully with matching configuration."))
}
load(gtypes_path)
g

# Stratify the gtypes object based on the chosen scheme.
g_stratified <- stratify(g, scheme = stratification_scheme, drop = TRUE)
strata_to_drop <- c("In-water, CA", "In-water, Indonesia")
g_filtered <- g_stratified[!getStrata(g_stratified) %in% strata_to_drop, , ]
g_filtered

message(paste("Data stratified by:", stratification_scheme))


# ====================================================================
# STEP 2: INDIVIDUAL SUMMARY AND VISUALIZATION
# ====================================================================
message("\nStep 2: Summarizing individual data and creating diagnostic plots...")

ind_summary <- summarizeInds(g_filtered)

# Create a scatter plot of missing genotypes vs. homozygosity.
homozygosity_plot <- ggplot(ind_summary, aes(x = num.loci.missing.genotypes, y = pct.loci.homozygous)) +
  geom_point(alpha = 0.7, color = "#2c7fb8") +
  labs(
    title = "Individual Homozygosity vs. Missing Genotypes",
    x = "Number of Loci with Missing Genotypes",
    y = "Proportion of Homozygous Loci"
  ) +
  theme_minimal()

print(homozygosity_plot)

# Save the plot.
plot_path <- file.path(results_raw_path, paste0(run_label, "_missing_vs_homozygosity.pdf"))
ggsave(plot_path, plot = homozygosity_plot, device = "pdf")
message(paste("Diagnostic plot saved to:", plot_path))


# ====================================================================
# STEP 3: PAIRWISE POPULATION STRUCTURE TESTS
# ====================================================================
message("\nStep 3: Performing pairwise population structure tests...")

# Filter for strata with the minimum required number of individuals.
g_filtered <- g_filtered[,,which(getNumInd(g_filtered, by.strata = TRUE)$num.ind >= min_inds_for_pairwise)]
num_strata <- length(getStrataNames(g_filtered))

if (num_strata < 2) {
  stop(paste("Fewer than 2 strata remain after filtering for", min_inds_for_pairwise, "individuals. Cannot perform pairwise tests."))
}
message(paste("Running pairwise tests on", num_strata, "strata with at least", min_inds_for_pairwise, "individuals."))

# Run the pairwise tests. This may take some time.
pws_struct <- pairwiseTest(g_filtered, nrep = num_reps_pvals)
pws_summary <- pairwiseSummary(pws_struct)

message("Pairwise tests complete. Summary:")
print(pws_summary)


# ====================================================================
# STEP 4: EXTRACT AND SAVE RESULTS
# ====================================================================
message("\nStep 4: Extracting and saving pairwise matrices...")

# Extract Chi-squared and Fst matrices.
# Note: For genepop output, the statistic estimate is in the lower triangle and p-value is in the upper.
chi2_mat <- pairwiseMatrix(pws_struct, stat = 'CHIsq')
fst_mat <- pairwiseMatrix(pws_struct, stat = 'Fst')

message("\nPairwise Chi-squared Matrix (estimate below, p-value above):")
print(chi2_mat)
message("\nPairwise Fst Matrix (estimate below, p-value above):")
print(fst_mat)

# --- Save results to CSV ---
chi2_path <- file.path(results_raw_path, paste0(run_label, "_pairwise_chi2.csv"))
write.csv(chi2_mat, file = chi2_path)
message(paste("Chi-squared matrix saved to:", chi2_path))

fst_path <- file.path(results_raw_path, paste0(run_label, "_pairwise_fst.csv"))
write.csv(fst_mat, file = fst_path)
message(paste("Fst matrix saved to:", fst_path))

# --- Save R objects ---
rda_path <- file.path(results_r_path, paste0(run_label, "_pop_struct_results.rda"))
save(g_stratified, ind_summary, pws_struct, file = rda_path)
message(paste("R objects saved to:", rda_path))

message("\n✅ Population Structure Analysis Complete!")