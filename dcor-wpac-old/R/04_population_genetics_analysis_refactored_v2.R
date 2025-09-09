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

# ----------------------------------------------------------------------
# --- CONFIGURATION (UPDATED) ---
# ----------------------------------------------------------------------

# --- Project and File Naming ---
# These settings MUST match the previous scripts.
project_name <- "dcor_wpac"
final_analysis_stage <- "final_nodups"
min_reads <- 20

# --- Auto-Generated Paths (Do not change) ---
# Recreates the run label from previous scripts to find/save files correctly.
run_label <- paste(project_name, final_analysis_stage, sep = "_")

# Path to the gtypes .rda file created by `02_create_gtypes.R`.
gtypes.path <- file.path("data", paste0("gtypes_", run_label, "_minReads", min_reads, ".rda"))

# --- Output Paths ---
results_raw_path <- file.path("results-raw", final_analysis_stage)
data_path <- "data/"
dir.create(results_raw_path, showWarnings = FALSE, recursive = TRUE)

# --- Analysis Parameters ---
high.homo.threshold <- 0.7
ld.locus.geno.threshold <- 0.5
p.adjust.method <- "holm"

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
             "\nPlease ensure script 02 has been run successfully with matching configuration."))
}
load(gtypes.path)

# Automatically get the list of all strata from the gtypes object.
strata.to.analyze <- getStrataNames(g)
message(paste("Automatically detected", length(strata.to.analyze), "strata for analysis:"))
message(paste(" -", strata.to.analyze, collapse = "\n"))

# Optional: Manually remove specific strata
strata.to.analyze <- strata.to.analyze[!strata.to.analyze %in% c("In-water, CA", "In-water, Indonesia")]
print(strata.to.analyze)

# Create a new gtypes object containing only the strata to be analyzed.
pop.g <- g[i = which(getStrata(g) %in% strata.to.analyze)]

# Create a list of gtypes for each of the strata to keep
all.strats.g <- strataSplit(g)
pop.strats.g.split <- all.strats.g[names(all.strats.g) %in% strata.to.analyze]


# ====================================================================
# STEP 2: SUMMARY STATISTICS
# ====================================================================
message("\nStep 2: Calculating locus summary statistics...")

loc.sum <- summarizeLoci(g)
loc.sum.strata <- summarizeLoci(g, by.strata = TRUE)

locus.snp.counts <- g@data |>
  select(locus, allele) |>
  distinct() |>
  filter(!is.na(allele)) |>
  mutate(num.snps = nchar(allele)) |>
  select(locus, num.snps) |>
  distinct()

loc.sum <- left_join(loc.sum, locus.snp.counts, by = "locus")
loc.sum.strata <- left_join(loc.sum.strata, locus.snp.counts, by = "locus")

# --- Save Locus Summary results (filenames updated) ---
write.csv(loc.sum, file = file.path(results_raw_path, paste0(run_label, ".loc.sum.csv")), row.names = FALSE)
write.csv(loc.sum.strata, file = file.path(results_raw_path, paste0(run_label, ".loc.sums.by.strat.csv")), row.names = FALSE)
save(loc.sum.strata, loc.sum, file = file.path(data_path, paste0(run_label, ".loc.sum.rda")))

message("Locus summary statistics calculated and saved.")

message("\nStep 2: Calculating individual summary statistics...")
ind.sum <- summarizeInds(g)

# --- Save Individual Summary results (filenames updated) ---
write.csv(ind.sum, file = file.path(results_raw_path, paste0(run_label, ".ind.sum.csv")), row.names = FALSE)
save(ind.sum, file = file.path(data_path, paste0(run_label, ".ind.sum.rda")))
message("Individual summary statistics calculated and saved.")
# ====================================================================
# STEP 3: HOMOZYGOSITY CHECK
# ====================================================================
message("\nStep 3: Calculating and visualizing individual homozygosity...")

ind.summary <- summarizeInds(g)
message("Summary of homozygosity across all individuals:")
print(summary(ind.summary$pct.loci.homozygous))

homo_hist <- ggplot(ind.summary, aes(x = pct.loci.homozygous)) +
  geom_histogram(
    binwidth = 0.025, fill = "#2c7fb8", color = "black", alpha = 0.8,
    boundary = high.homo.threshold
  ) +
  geom_vline(
    xintercept = high.homo.threshold, color = "red", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = high.homo.threshold, y = Inf,
    label = paste("Threshold =", high.homo.threshold),
    hjust = -0.1, vjust = 2, color = "red", angle = 90
  ) +
  labs(
    title = "Distribution of Individual Homozygosity-All samples",
    subtitle = "Inspect this plot to inform your filtering threshold.",
    x = "Proportion of Homozygous Loci", y = "Number of Individuals"
  ) +
  theme_minimal()
print(homo_hist)

message("\nACTION: A histogram has been generated. Adjust 'high.homo.threshold' if needed and re-run.")

high.homo.samps <- filter(ind.summary, pct.loci.homozygous > high.homo.threshold) %>% pull(id)

if (length(high.homo.samps) > 0) {
  message("\nWarning: Found ", length(high.homo.samps), " individual(s) exceeding the threshold:")
  print(high.homo.samps)
} else {
  message("\nNo individuals found above the current homozygosity threshold.")
}

# --- Save histogram (filename updated) ---
pdf(file.path(results_raw_path, paste0(run_label, ".homozygosity.histograms.pdf")))
print(homo_hist)
dev.off()

# ====================================================================
# STEP 4: HARDY-WEINBERG EQUILIBRIUM (HWE) ANALYSIS
# ====================================================================
message("\nStep 4: Calculating Hardy-Weinberg Equilibrium...")

hwe.list <- imap(pop.strats.g.split, ~{
  hweTest(.x) %>%
    data.frame() %>%
    rownames_to_column(var = "locus") %>%
    setNames(c("locus", .y))
})
hwe.res <- reduce(hwe.list, full_join, by = "locus")

hwe.res.unadj <- hwe.res %>%
  rowwise() %>%
  mutate(num.sig.p = sum(c_across(all_of(strata.to.analyze)) < 0.05, na.rm = TRUE)) %>%
  ungroup()

hwe.p.adj <- hwe.res %>%
  select(-locus) %>%
  as.matrix() %>%
  apply(1, function(p_vals) p.adjust(p_vals, method = p.adjust.method, n = sum(!is.na(p_vals)))) %>%
  t() %>%
  as.data.frame()

hwe.res.adj <- bind_cols(locus = hwe.res$locus, hwe.p.adj) %>%
  rowwise() %>%
  mutate(num.sig.adj.p = sum(c_across(all_of(strata.to.analyze)) < 0.05, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(num.sig.adj.p))

# --- Save HWE results (filenames updated) ---
write.csv(hwe.res.unadj, file.path(results_raw_path, paste0(run_label, ".hwe.unadjusted.csv")), row.names = FALSE)
write.csv(hwe.res.adj, file.path(results_raw_path, paste0(run_label, ".hwe.adjusted.csv")), row.names = FALSE)
save(hwe.res.unadj, hwe.res.adj, file = file.path(data_path, paste0(run_label, ".hwe.results.rda")))

message("HWE analysis complete. Results saved.")

# ====================================================================
# STEP 5: LINKAGE DISEQUILIBRIUM (LD) ANALYSIS
# ====================================================================
message("\nStep 5: Calculating Linkage Disequilibrium...")

g.filtered.for.ld <- imap(pop.strats.g.split, ~{
  loc_summary <- summarizeLoci(.x)
  loci_to_remove <- loc_summary %>% filter(prop.genotyped < ld.locus.geno.threshold) %>% pull(locus)
  if (length(loci_to_remove) > 0) {
    message(paste("   - For stratum '", .y, "', removing", length(loci_to_remove), "loci with <", ld.locus.geno.threshold * 100, "% genotyping."))
    .x <- .x[, -which(getLociNames(.x) %in% loci_to_remove),]
  }
  return(.x)
})

message("Running LDgenepop on all strata. This may take a while...")
ld.strata.list <- map(g.filtered.for.ld, LDgenepop)

ld.sig.res <- imap(ld.strata.list, ~{
  .x %>%
    select(Locus.1, Locus.2, p.value) %>%
    drop_na(Locus.1, Locus.2) %>%
    group_by(Locus.1, Locus.2) %>%
    summarise(p.value = min(p.value, na.rm = TRUE), .groups = "drop") %>%
    rename(!!sym(paste0('p.val.', .y)) := p.value)
}) %>% reduce(full_join, by = c("Locus.1", "Locus.2"))

ld.sig.res <- ld.sig.res %>%
  rowwise() %>%
  mutate(num.sig.p = sum(c_across(starts_with("p.val.")) < 0.05, na.rm = TRUE)) %>%
  ungroup()

p.val.cols <- ld.sig.res %>% select(starts_with("p.val."))
ld.p.adj <- apply(p.val.cols, 1, function(p) p.adjust(p, method = p.adjust.method, n = sum(!is.na(p)))) %>%
  t() %>%
  as.data.frame()
colnames(ld.p.adj) <- paste0(names(p.val.cols), "_adj")

ld.res.adj <- ld.sig.res %>%
  select(Locus.1, Locus.2) %>%
  bind_cols(ld.p.adj) %>%
  rowwise() %>%
  mutate(num.sig.adj.p = sum(c_across(ends_with("_adj")) < 0.05, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(num.sig.adj.p))

# --- Save LD results (filenames updated) ---
write.csv(ld.sig.res, file.path(results_raw_path, paste0(run_label, ".ld.unadjusted.csv")), row.names = FALSE)
write.csv(ld.res.adj, file.path(results_raw_path, paste0(run_label, ".ld.adjusted.csv")), row.names = FALSE)
save(ld.sig.res, ld.res.adj, ld.strata.list, file = file.path(data_path, paste0(run_label, ".ld.results.rda")))

message("LD analysis complete. Results saved.")
message("\n✅ Workflow Complete!")
