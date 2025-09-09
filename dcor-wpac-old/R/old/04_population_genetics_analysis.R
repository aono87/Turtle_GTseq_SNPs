# 04_POPULATION_GENETICS_ANALYSIS.R
#
# Workflow Step 4:
# This script performs downstream population genetic analyses using the `gtypes`
# object created in Step 2. It calculates locus summary statistics, visualizes
# data distributions, calculates individual homozygosity, and tests for
# Hardy-Weinberg Equilibrium (HWE) and Linkage Disequilibrium (LD).
#
# INPUTS:
#   - A `gtypes` object (`.rda` file) from `02_create_gtypes.R`.
#
# OUTPUTS:
#   - CSV and RDA files for locus summary statistics.
#   - A histogram plot of homozygosity.
#   - CSV files summarizing HWE and LD results.
#   - RDA files containing the full result objects for HWE and LD.

# --- LOAD LIBRARIES ---
library(tidyverse)
library(strataG)
library(purrr)
library(tibble)

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
# Homozygosity threshold.
# The script will first generate a histogram to help you choose this value.
# You can run the script once, inspect the plot, then adjust this value and re-run.
high.homo.threshold <- 0.7


# LD analysis: Loci with genotyping proportion below this threshold *within a stratum* will be removed.
ld.locus.geno.threshold <- 0.5

# P-value adjustment method for multiple testing. "holm" is a standard choice.
p.adjust.method <- "holm"

# --- Output Paths ---
results.raw.path <- "results-raw/"
data.path <- "data/"

# -------------------------
# --- END CONFIGURATION ---
# -------------------------


# ====================================================================
# STEP 1: LOAD DATA AND DEFINE STRATA
# ====================================================================
message("Step 1: Loading gtypes object and preparing data...")

# Check if the gtypes file exists before proceeding.
if (!file.exists(gtypes.path)) {
  stop(paste("Gtypes file not found at:", gtypes.path,
             "\nPlease ensure script 02 has been run successfully."))
}
load(gtypes.path)

# Automatically get the list of all strata from the gtypes object.
strata.to.analyze <- getStrataNames(g)
message(paste("Automatically detected", length(strata.to.analyze), "strata for analysis:"))
message(paste(" -", strata.to.analyze, collapse = "\n"))

# Optional: Manually remove specific strata from the analysis if needed.
# Removing the In-water strata since their population of origin is unknown
# Create a list of gtypes for each of the strata to keep
strata.to.analyze <- strata.to.analyze[!strata.to.analyze %in% c("In-water, CA", "In-water, Indonesia")]
print(strata.to.analyze)

# Create a new gtypes object containing only the strata to be analyzed.
pop.g <- g[i = which(getStrata(g) %in% strata.to.analyze)]

# Create a list of gtypes for each of the strata to keep
all.strats.g<-strataSplit(g)
pop.strats.g.split <- all.strats.g[names(all.strats.g) %in% strata.to.analyze]


# ====================================================================
# STEP 2: LOCUS SUMMARY STATISTICS
# ====================================================================
message("\nStep 2: Calculating locus summary statistics...")

# Calculate overall summary statistics for all loci.
loc.sum <- summarizeLoci(g)

# Calculate summary statistics for loci within each defined stratum.
loc.sum.strata <- summarizeLoci(g, by.strata = TRUE)

# SNP counts per locus ---
locus.snp.counts <- g@data |>
  select(locus, allele) |>
  distinct() |>
  filter(!is.na(allele)) |>
  mutate(num.snps = nchar(allele)) |>
  select(locus, num.snps) |>
  distinct()

# --- Join SNP counts with summary tables ---
loc.sum <- left_join(loc.sum, locus.snp.counts, by = "locus")
loc.sum.strat <- left_join(loc.sum.strata, locus.snp.counts, by = "locus")

# --- Save Locus Summary results ---
write.csv(loc.sum, file = file.path(results.raw.path, paste0(project, ".loc.sum.csv")), row.names = FALSE)
write.csv(loc.sum.strata, file = file.path(results.raw.path, paste0(project, ".loc.sums.by.strat.csv")), row.names = FALSE)
save(loc.sum.strata, loc.sum, file = file.path(data.path, paste0(project, ".loc.sum.rda")))

message("Locus summary statistics calculated and saved.")


# ====================================================================
# STEP 3: HOMOZYGOSITY CHECK
# ====================================================================
message("\nStep 3: Calculating and visualizing individual homozygosity...")

# First, calculate summary stats for each individual.
ind.summary <- summarizeInds(g)

# Print a numerical summary of the homozygosity distribution to the console.
message("Summary of homozygosity across all individuals:")
print(summary(ind.summary$pct.loci.homozygous))

# Generate a histogram to visualize the distribution.
homo_hist <- ggplot(ind.summary, aes(x = pct.loci.homozygous)) +
  geom_histogram(
    binwidth = 0.025,
    fill = "#2c7fb8",
    color = "black",
    alpha = 0.8,
    boundary = high.homo.threshold # Align bin edge with the threshold
  ) +
  geom_vline(
    xintercept = high.homo.threshold,
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  annotate(
    "text", x = high.homo.threshold, y = Inf,
    label = paste("Threshold =", high.homo.threshold),
    hjust = -0.1, vjust = 2, color = "red", angle = 90
  ) +
  labs(
    title = "Distribution of Individual Homozygosity-All samples",
    subtitle = "Inspect this plot to inform your filtering threshold.",
    x = "Proportion of Homozygous Loci",
    y = "Number of Individuals"
  ) +
  theme_minimal()

# Display the plot.
print(homo_hist)

message("\nACTION: A histogram has been generated in the 'Plots' pane.")
message("        The red dashed line shows the current 'high.homo.threshold'.")
message("        If you need to adjust it, change the value in the CONFIGURATION section and re-run the script.")


# Now, identify individuals exceeding the currently set threshold.
high.homo.samps <- filter(ind.summary, pct.loci.homozygous > high.homo.threshold) %>%
  pull(id)

if (length(high.homo.samps) > 0) {
  message("\nWarning: Found ", length(high.homo.samps), " individual(s) exceeding the threshold:")
  print(high.homo.samps)
  # To remove these individuals from further analysis, uncomment the following line:
  # g <- g[-which(getIndNames(g) %in% high.homo.samps),]
  # message("High-homozygosity individuals have been removed.")
} else {
  message("\nNo individuals found above the current homozygosity threshold.")
}

# Save histogram to drive
pdf(file.path(results.raw.path, paste0(project, ".homozygosity.histograms.pdf")))
homo_hist
dev.off()

# ====================================================================
# STEP 4: HARDY-WEINBERG EQUILIBRIUM (HWE) ANALYSIS
# ====================================================================
message("\nStep 4: Calculating Hardy-Weinberg Equilibrium...")

# Use imap to loop over the list and its names at the same time.
# .x will be the gtypes object (the list element)
# .y will be the stratum name (the list element's name)
hwe.list <- imap(pop.strats.g.split, ~{
  hweTest(.x) %>%
    data.frame() %>%
    rownames_to_column(var = "locus") %>%
    # Use setNames to name the columns dynamically
    setNames(c("locus", .y))
  
})

# Combine the list of data frames into a single table
hwe.res <- reduce(hwe.list, full_join, by = "locus")

# --- Calculate significance counts (Unadjusted) ---
hwe.res.unadj <- hwe.res %>%
  rowwise() %>%
  mutate(
    num.sig.p = sum(c_across(all_of(strata.to.analyze)) < 0.05, na.rm = TRUE)
  ) %>%
  ungroup()

# --- Calculate adjusted p-values ---
hwe.p.adj <- hwe.res %>%
  select(-locus) %>%
  as.matrix() %>%
  apply(1, function(p_vals) p.adjust(p_vals, method = p.adjust.method, n = sum(!is.na(p_vals)))) %>%
  t() %>%
  as.data.frame()

# Combine adjusted p-values with locus names and calculate significance
hwe.res.adj <- bind_cols(locus = hwe.res$locus, hwe.p.adj) %>%
  rowwise() %>%
  mutate(num.sig.adj.p = sum(c_across(all_of(strata.to.analyze)) < 0.05, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(num.sig.adj.p))


# --- Save HWE results ---
write.csv(hwe.res.unadj, file.path(results.raw.path, paste0(project, ".hwe.unadjusted.csv")), row.names = FALSE)
write.csv(hwe.res.adj, file.path(results.raw.path, paste0(project, ".hwe.adjusted.csv")), row.names = FALSE)
save(hwe.res.unadj, hwe.res.adj, file = file.path(data.path, paste0(project, ".hwe.results.rda")))

message("HWE analysis complete. Results saved.")

# ====================================================================
# STEP 5: LINKAGE DISEQUILIBRIUM (LD) ANALYSIS
# ====================================================================
message("\nStep 5: Calculating Linkage Disequilibrium...")

pop.strats.g.split

# --- Pre-filter loci within each stratum for missing data ---
g.filtered.for.ld <- imap(pop.strats.g.split, ~{
  loc_summary <- summarizeLoci(.x)
  loci_to_remove <- loc_summary %>% filter(prop.genotyped < ld.locus.geno.threshold) %>% pull(locus)
  if (length(loci_to_remove) > 0) {
    message(paste("   - For stratum '", .y, "', removing", length(loci_to_remove), "loci with <", ld.locus.geno.threshold * 100, "% genotyping."))
    .x <- .x[, -which(getLociNames(.x) %in% loci_to_remove),]
  }
  return(.x)
})

# --- Run LD analysis ---
message("Running LDgenepop on all strata. This may take a while...")
ld.strata.list <- map(g.filtered.for.ld, LDgenepop)

#Ensuring there is only one p-value per locus1 and locus2 combination, per stratum
list_pair_counts <- map_dfr(ld.strata.list, ~{
  .x %>%
    drop_na(1, 2) %>% # <-- ADD THIS LINE to remove NAs from columns 1 & 2
    count(across(1:2), name = "times_repeated")
}, .id = "source_df_name") %>%
  filter(times_repeated > 1)
list_pair_counts

# --- Process LD results ---
ld.sig.res <- imap(ld.strata.list, ~{
  .x %>%
    select(Locus.1, Locus.2, p.value) %>%  # Select only the columns you need
    drop_na(Locus.1, Locus.2) %>%  # Remove rows with NAs in locus columns
    group_by(Locus.1, Locus.2) %>%     # Ensure each locus pair is unique within each data frame. Here, we take the minimum p-value for any duplicates.
    summarise(p.value = min(p.value, na.rm = TRUE), .groups = "drop") %>%
    rename(!!sym(paste0('p.val.', .y)) := p.value) # Rename the p-value column dynamically using the list element's name
}) %>% reduce(full_join, by = c("Locus.1", "Locus.2")) # Sequentially join all data frames in the list by the two Locus columns

#Check that there are no locus 1 and locus 2 combos repeated
all_pair_counts <- ld.sig.res %>%
  count(across(1:2), name = "times_repeated") %>% # Count pairs
  filter(times_repeated > 1) # Keep only pairs that appear more than once
all_pair_counts

# Calculate unadjusted significance
ld.sig.res <- ld.sig.res %>%
  rowwise() %>%
  mutate(num.sig.p = sum(c_across(starts_with("p.val.")) < 0.05, na.rm = TRUE)) %>%
  ungroup()

# Calculate adjusted p-values
p.val.cols <- ld.sig.res %>% select(starts_with("p.val."))
ld.p.adj <- apply(p.val.cols, 1, function(p) p.adjust(p, method = p.adjust.method, n = sum(!is.na(p)))) %>%
  t() %>%
  as.data.frame()
colnames(ld.p.adj) <- paste0(names(p.val.cols), "_adj")

# Combine and summarize adjusted results
ld.res.adj <- ld.sig.res %>%
  select(Locus.1, Locus.2) %>%
  bind_cols(ld.p.adj) %>%
  rowwise() %>%
  mutate(num.sig.adj.p = sum(c_across(ends_with("_adj")) < 0.05, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(num.sig.adj.p))

# --- Save LD results ---
write.csv(ld.sig.res, file.path(results.raw.path, paste0(project, ".ld.unadjusted.csv")), row.names = FALSE)
write.csv(ld.res.adj, file.path(results.raw.path, paste0(project, ".ld.adjusted.csv")), row.names = FALSE)
save(ld.sig.res, ld.res.adj, ld.strata.list, file = file.path(data.path, paste0(project, ".ld.results.rda")))

message("LD analysis complete. Results saved.")
message("\n✅ Workflow Complete!")