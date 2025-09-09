# 05_DAPC_ANALYSIS.R
#
# Workflow Step 5:
# This script performs Discriminant Analysis of Principal Components (DAPC) to
# investigate population structure based on the "Nesting Area" stratification.
# It includes cross-validation to determine the optimal number of principal
# components (PCs), visualizes the results, and assigns individuals to
# populations.
#
# INPUTS:
#   - A `gtypes` object (`.rda` file) from `02_create_gtypes.R`.
#
# OUTPUTS:
#   - JPEG files for DAPC scatter plots and assignment tables.
#   - JPEG file showing the assignment of samples to reference populations.
#   - An RDA file saving the leave-one-out cross-validation results.

# --- LOAD LIBRARIES ---
library(adegenet)
library(strataG)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(ggrepel)
# The custom function script below is required for the leave-one-out validation step.
# source('R/functions/DAPC.fit.and.predict.R')

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
results_r_path <- file.path("results-R", analysis_stage)
dir.create(results_raw_path, showWarnings = FALSE, recursive = TRUE)
dir.create(results_r_path, showWarnings = FALSE, recursive = TRUE)

# --- Analysis Parameters ---
# The single stratification scheme to be used for the analysis.
stratification_scheme <- "Stratum_ABO"

# Gtypes file
g

# Minimum number of individuals per stratum for DAPC.
#min_inds_dapc <- 9

# Cross-validation parameters.
xval_n_pca_max <- 200
xval_training_set <- 0.9
xval_reps <- 30

# DAPC parameters for prediction (if using the custom predictAllIDsDAPC function).
loov_n_da <- 8
loov_n_pca <- 100

# -------------------------
# --- END CONFIGURATION ---
# -------------------------


# ====================================================================
# STEP 1: LOAD AND PREPARE DATA
# ====================================================================
message("Step 1: Loading gtypes object and preparing data for DAPC...")

# Check if the gtypes file exists before proceeding.
if (!file.exists(gtypes_path)) {
  stop(paste("Gtypes file not found at:", gtypes_path,
             "\nPlease ensure script 02 has been run successfully with matching configuration."))
}
load(gtypes_path)

# Stratify data and filter for minimum sample size.
g_stratified <- stratify(g, scheme = stratification_scheme, drop = TRUE)
#g_stratified <- g_stratified[,,
                             filter(getNumInd(g_stratified, by.strata = TRUE), num.ind >= min_inds_dapc) |>
                               pull(stratum)]
#message(paste("Data stratified by '", stratification_scheme, "' and filtered for strata with at least ", min_inds_dapc, " individuals.", sep=""))

# Convert to genind format for adegenet.
genind_all <- gtypes2genind(g_stratified)

# Separate reference and mixture populations.
#if (!is.null(mixture_population_name) && mixture_population_name %in% pop(genind_all)) {
#  message(paste("Separating mixture population:", mixture_population_name))
#  genind_ref <- genind_all[-which(pop(genind_all) == mixture_population_name)]
#  genind_mix <- genind_all[which(pop(genind_all) == mixture_population_name)]
#} else {
#  message("No mixture population specified or found; all samples will be used for DAPC training.")
#  genind_ref <- genind_all
#  genind_mix <- NULL
#}


# ====================================================================
# STEP 2: CROSS-VALIDATION TO DETERMINE OPTIMAL PCs
# ====================================================================
message("\nStep 2: Running cross-validation to find the optimal number of PCs...")

mat_ref <- tab(genind_all, NA.method = 'mean')
grp_ref <- pop(genind_all)

# Run the cross-validation. This can be time-consuming.
xval_results <- xvalDapc(mat_ref, grp_ref, n.pca.max = xval_n_pca_max,
                         training.set = xval_training_set, result = 'groupMean',
                         n.rep = xval_reps, xval.plot = TRUE)

dapc_results <- xval_results$DAPC
message(paste("   - Optimal number of PCs identified:", dapc_results$n.pca))
message(paste("   - Number of discriminant axes retained:", dapc_results$n.da))


# ====================================================================
# STEP 3: VISUALIZE AND SUMMARIZE DAPC RESULTS
# ====================================================================
message("\nStep 3: Generating plots and summarizing DAPC results...")

# Define a color palette for plotting.
num_groups <- nlevels(pop(genind_all))
plot_colors <- hcl.colors(num_groups, palette = "Zissou 1")

# --- Visualise DAPC Scatter Plot ---
scatter.dapc(dapc_results, col = plot_colors, scree.da = FALSE, legend = TRUE, posi.leg = "topright")

message("\nStep 3: Generating advanced plot with non-overlapping labels...")

# --- 1. Prepare Data for ggplot ---
# Extract individual coordinates and add population info
ind_coords_df <- as.data.frame(dapc_results$ind.coord)
ind_coords_df$population <- dapc_results$grp

population_counts <- ind_coords_df %>% 
  count(population)
print("Number of individuals per population:")
print(population_counts)

# Extract centroid coordinates
cent_coords_df <- as.data.frame(dapc_results$grp.coord)
cent_coords_df$population <- rownames(cent_coords_df)

# Calculate percentage of variance explained by the axes
percent_explained <- (dapc_results$eig[1:2] / sum(dapc_results$eig)) * 100

# --- 2. Create the ggplot object ---
dapc_plot <- ggplot() +
  # Add individual points
  geom_point(
    data = ind_coords_df,
    aes(x = LD1, y = LD2, color = population),
    size = 2, alpha = 0.6
  ) +
  # Add ellipses to group points by population
  stat_ellipse(
    data = ind_coords_df,
    aes(x = LD1, y = LD2, color = population, fill = population),
    geom = "polygon", alpha = 0.1, linetype = "dashed"
  ) +
  # Add the non-overlapping population labels at the centroids
  geom_label_repel(
    data = cent_coords_df,
    aes(x = LD1, y = LD2, label = population, fill = population),
    color = "white", fontface = "bold", point.padding = 1, box.padding = 0.5
  ) +
  # CHANGED: Apply the custom manual color scheme
  scale_color_manual(values = plot_colors, guide = "none") +
  scale_fill_manual(values = plot_colors, guide = "none") +
  # Add informative labels and a clean theme
  labs(
    title = "Discriminant Analysis of Principal Components (DAPC)",
    x = paste0("LD1 (", round(percent_explained[1], 1), "%)"),
    y = paste0("LD2 (", round(percent_explained[2], 1), "%)")
  ) +
  theme_classic()
dapc_plot

# --- Save DAPC Scatter Plot ---
scatter_path <- file.path(results_raw_path, paste0(run_label, "_DAPC_scatter.jpg"))
jpeg(filename = scatter_path, width = 800, height = 800)
scatter.dapc(dapc_results, col = plot_colors, scree.da = FALSE, legend = TRUE, posi.leg = "topright")
dev.off()
message(paste("   - DAPC scatter plot saved to:", scatter_path))

# --- Save DAPC Assignment Table ---
assign_path <- file.path(results_raw_path, paste0(run_label, "_DAPC_assignment_table.jpg"))
jpeg(filename = assign_path, width = 800, height = 800)
table.value(table(dapc_results$assign, grp_ref), 
            col.lab = levels(grp_ref))

dev.off()
message(paste("   - Assignment table plot saved to:", assign_path))

# --- Calculate and report assignment accuracy ---
assignment_accuracy <- mean(as.character(dapc_results$assign) == as.character(grp_ref))
message(paste("   - Overall model assignment accuracy for reference populations:", round(assignment_accuracy, 4)))
# - Overall model assignment accuracy for reference populations: 0.6528

message("\n✅ DAPC Analysis Complete!")

##Compoplot
lab<-pop(genind_all)
par(mar=c(12,4,5,1), xpd=TRUE)
compoplot(dapc_results, lab=lab, las=2, cex.names=0.7, cleg=.6, posi=list(x=0,y=1.2))

write.csv(dapc_results$posterior, "./dapc_results_posterior.csv")
