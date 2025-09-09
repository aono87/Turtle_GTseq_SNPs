# 06_MIXED_STOCK_ANALYSIS.R
#
# Workflow Step 6:
# This script conducts a mixed-stock analysis using the 'rubias' package. It
# performs self-assignment of reference samples based on the "Nesting Area"
# stratification and estimates the mixing proportions of a specified mixture sample.
#
# INPUTS:
#   - A `gtypes` object (`.rda` file) from `02_create_gtypes.R`.
#   - Individual information data (`ind.info`).
#
# OUTPUTS:
#   - CSV files for self-assignment matrices (individual-to-group).
#   - A density plot comparing the z-scores of assignments to a normal distribution.
#   - A density plot of the MCMC traces for top contributing groups.

# --- LOAD LIBRARIES ---
library(rubias)
library(strataG)
library(tidyverse)
library(ggplot2)

# ----------------------------------------------------------------------
# --- CONFIGURATION ---
# ----------------------------------------------------------------------

# --- Project and File Naming ---
# These settings MUST match the previous scripts.
project_name <- "dcor_wpac"
analysis_stage <- "final_nodups"
min_reads <- 20

# --- Auto-Generated Paths (Do not change) ---
run_label <- paste(project_name, analysis_stage, sep = "_")
gtypes_file <- paste0("gtypes_", run_label, "_minReads", min_reads, ".rda")
gtypes_path <- file.path("data", gtypes_file)

# --- Output Paths ---
results_raw_path <- file.path("results-raw", analysis_stage)
dir.create(results_raw_path, showWarnings = FALSE, recursive = TRUE)

# --- Analysis Parameters ---
# The single stratification scheme to be used for the analysis.
stratification_scheme <- "Stratum_ABO"
# The name of the stratum to be treated as a mixture sample (set to NULL if none).
mixture_population_name <- NULL
# Minimum number of individuals per stratum to be included in the analysis.
min_inds_per_stratum <- 3
# MCMC burn-in period for mixture proportion trace analysis.
mcmc_burn_in <- 200

# -------------------------
# --- END CONFIGURATION ---
# -------------------------


# ====================================================================
# STEP 1: LOAD DATA AND PREPARE FOR RUBIAS
# ====================================================================
message("Step 1: Loading gtypes object and preparing data for rubias...")

# Check if the gtypes file exists before proceeding.
if (!file.exists(gtypes_path)) {
  stop(paste("Gtypes file not found at:", gtypes_path,
             "\nPlease ensure script 02 has been run successfully with matching configuration."))
}
load(gtypes_path)
#data("ind.info") # Load individual metadata, assumed to exist.

# Stratify data and filter for minimum sample size.
g_stratified <- stratify(g, scheme = stratification_scheme, drop = TRUE)
g_stratified <- g_stratified[,,
                             filter(getNumInd(g_stratified, by.strata = TRUE), num.ind >= min_inds_per_stratum) |>
                               pull(stratum)]
message(paste("Data stratified by '", stratification_scheme, "' and filtered for strata with at least ", min_inds_per_stratum, " individuals.", sep=""))

# Format the gtypes data into the required format for rubias.
# This assumes the 'schemes' data frame in the gtypes object contains the stratification column.
rubias_data <- right_join(
  # --- Part 1: Prepare the Metadata ---
  g_stratified@schemes %>%
    select(id, all_of(stratification_scheme)) %>%
    mutate(
      # Classify samples based on the stratum name
      sample_type = ifelse(
        startsWith(.data[[stratification_scheme]], "In-water"), 
        'mixture', 
        'reference'
      ),
      
      # CORRECTED LOGIC: Assign repunit conditionally.
      # If the sample is a 'reference', use the stratum name.
      # If it's a 'mixture', assign NA.
      repunit = ifelse(sample_type == 'reference', .data[[stratification_scheme]], NA),
      
      # Collection is always the stratum name
      collection = .data[[stratification_scheme]]
    ),
  
  # --- Part 2: Restructure the Genetic Data ---
  g_stratified@data %>%
    mutate(allele_copy = paste(locus, c(1,2), sep = '_')) %>%
    select(-locus, -stratum) %>%
    pivot_wider(names_from = allele_copy, values_from = allele)
) %>%
  rename(indiv = id) %>%
  select(indiv, sample_type, repunit, collection, everything())

message("Data successfully formatted for rubias analysis.")


# ====================================================================
# STEP 2: SELF-ASSIGNMENT OF REFERENCE SAMPLES
# ====================================================================
message("\nStep 2: Performing self-assignment of reference samples...")

reference_samples <- filter(rubias_data, sample_type == 'reference')
mixture_samples <- filter(rubias_data, sample_type == 'mixture')
self_assignment_results <- self_assign(reference = reference_samples, gen_start_col = 6)

# Create a summary matrix of assignments.
assignment_summary <- self_assignment_results %>%
  group_by(indiv) %>%
  top_n(1, scaled_likelihood) %>%
  ungroup() %>%
  mutate(across(where(is.list), as.character))

assignment_matrix <- table(
  Observed = assignment_summary$collection,
  Inferred = assignment_summary$inferred_collection
)
message("Self-assignment summary matrix:")
print(assignment_matrix)

nesting_z <- assignment_summary %>% 
  group_by(indiv) %>%
  top_n(1, z_score) %>% #picks the best population (highest z score)
  ungroup()
normo <- tibble(z_score = rnorm(1e06)) #create a million numbers that are normally distributed
ggplot(nesting_z, aes(x = z_score)) + #plots distribution of z scores (blue) against normal distribution (black)
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black")
#good fit: assignments are statistically sound and model is working well
#poor fit: issues with analysis including "ghost" populations or poor baseline data

# Save the detailed and summary assignment results.

write.csv(assignment_matrix, file = file.path(results_raw_path, paste0(run_label, "_rubias_self_assignment_matrix.csv")))
write.csv(assignment_summary, file = file.path(results_raw_path, paste0(run_label, "_rubias_self_assignment_summary.csv")))
message("Self-assignment results saved.")


# ====================================================================
# STEP 3: ESTIMATE PROPORTIONS OF MIXTURE SAMPLES
# ====================================================================
if (nrow(mixture_samples) > 0) {
  message("\nStep 3: Estimating mixing proportions of the mixture samples...")
  
  mix_est <- infer_mixture(reference = reference_samples, #uses genetic info from reference samples to estimate the proportion that each of those populations contributes to unknown mixture samples
                           mixture = mixture_samples,
                           gen_start_col = 6)
  
  # Summarize mixing proportions by reporting unit.
  rep_mix_ests <- mix_est$mixing_proportions %>%
    group_by(mixture_collection, repunit) %>%
    summarise(mixing_proportion = sum(pi), .groups = 'drop')
  
  message("Estimated mixing proportions:")
  print(rep_mix_ests)
  write.csv(rep_mix_ests, file = file.path(results_raw_path, paste0(run_label, "_rubias_mixture_proportions.csv")), row.names = FALSE)
  
  # --- Summarize Individual Proportions ---
  # For each individual, find the source population with the highest probability.
  indiv_assignments <- mix_est$indiv_posteriors %>%
    group_by(indiv) %>%
    # Sort by probability (PofZ) and then take the top row for each individual
    arrange(desc(PofZ)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    # Select and rename the correct PofZ column
    select(indiv, mixture_collection, repunit, "PofZ") %>%
    rename(
      assigned_population = repunit,
      assignment_probability = PofZ
    )
  
  message("Most likely source population for each individual:")
  print(indiv_assignments)
  
  # Save the individual assignment summary to top repunit to a CSV file.
  write.csv(indiv_assignments, file = file.path(results_raw_path, paste0(run_label, "_rubias_individual_assignments.csv")), row.names = FALSE)
  message(paste("\nIndividual assignment summary to top repunit saved to:", file.path(results_raw_path, paste0(run_label, "_rubias_individual_assignments.csv"))))
  
  # For each individual, show probabilities for all source populations
  indiv_assignments_all <- mix_est$indiv_posteriors %>%
    # Select and rename the correct PofZ column
    select(indiv, mixture_collection, repunit, "PofZ") %>%
    rename(
      assigned_population = repunit,
      assignment_probability = PofZ
    ) %>%
    # Add this line to format the column
    mutate(assignment_probability = sprintf("%.4g", assignment_probability))
  
  message("Most likely source population for each individual:")
  print(indiv_assignments_all)
  
  # Save the individual assignment summary to a CSV file.
  write.csv(indiv_assignments_all, file = file.path(results_raw_path, paste0(run_label, "_rubias_individual_assignments_all.csv")), row.names = FALSE)
  message(paste("\nIndividual assignment summary with all repunits saved to:", file.path(results_raw_path, paste0(run_label, "_rubias_individual_assignments_all.csv"))))
  
  # --- Analyze MCMC Traces ---
  # assess the stability and reliabiilty of mixed-stock analysis results
  # visualizes output of MCMC for two most important source populations
  message("\n   - Analyzing MCMC traces for top reporting units...")
  top2_rep_units <- rep_mix_ests %>% arrange(desc(mixing_proportion)) %>% slice(1:2)
  
  trace_subset <- mix_est$mix_prop_traces %>%
    filter(sweep > mcmc_burn_in) %>% #discards burn-in
    group_by(sweep, repunit) %>% #groups data by each sweep (step) and repunit in simulation
    summarise(pi_sum = sum(pi), .groups = 'drop') %>%
    filter(repunit %in% top2_rep_units$repunit) #filter to keep only top two repunits
  
  #create density plot showing range and likelihood of estimated mixing proportion for two top source populations
  mcmc_trace_plot <- ggplot(trace_subset, aes(x = pi_sum, colour = repunit, fill = repunit)) +
    geom_density(alpha = 0.5) +
    labs(
      title = "Posterior Distributions of Mixture Proportions",
      subtitle = paste("After discarding first", mcmc_burn_in, "sweeps as burn-in."),
      x = "Mixture Proportion (pi)",
      y = "Density"
    ) +
    theme_minimal()
  
  mcmc_plot_path <- file.path(results_raw_path, paste0(run_label, "_rubias_mcmc_trace_density.pdf"))
  ggsave(mcmc_plot_path, plot = mcmc_trace_plot, device = "pdf")
  message(paste("   - MCMC trace density plot saved to:", mcmc_plot_path))
  
} else {
  message("\nStep 3: Skipping mixture proportion estimation as no mixture samples were found.")
}

mcmc_trace_plot #peak is to the right (left skew), suggests the best estimate for populatio's contribution is high, with high confidence
#genotypes of mixture are a good match for this source population
message("\n✅ Rubias Analysis Complete!")