# 01_GENOTYPING_AND_QC.R
#
# Workflow Step 1:
# This script uses the `microhaplot` package to call genotypes from GT-seq data
# and then performs a rigorous quality control (QA/QC) process to filter
# out low-quality loci and individuals.
#
# INPUTS:
#   - A directory of SAM files for each individual.
#   - A VCF file containing the target SNP loci.
#   - A sample labels file for `microhaplot`.
#   - Custom functions located in `R/functions/`.
#
# OUTPUTS:
#   - A cleaned `geno.table` with high-quality genotypes.
#   - Summary CSV files for locus and individual filtering steps.
#   - Intermediate R objects (`.rds`, `.rda`) saved to `results-R/`.

# --- LOAD LIBRARIES ---
#Instructions to install libraries in readme.md file
library(tidyverse)
library(vcfR)
library(microhaplot)
library(vcftoolsR)

# --- LOAD CUSTOM FUNCTIONS ---
# Ensure the path to your custom functions is correct.
source("R/functions/mplot2tgt.R")
source("R/functions/tgt.2.geno.table.R")
# source("R/functions/Compare.replicates.R") # If needed

# ---------------------------
# --- BEGIN CONFIGURATION ---
# ---------------------------

# --- Project and File Paths ---
project <- "testing-new-scripts"
run.label <- "testing-new-scripts"

# Input paths
sam.path   <- "data-raw/sam_files/"
label.path <- "data-raw/metadata/w-leatherback-inds-zpad.txt"
vcf.path   <- "data-raw/vcf/Dcor_DcPanel_205_maf.targetSNPs_012224.recode.vcf"

# Output paths
out.path <- "results-R/" #need to use the "out.path" name for the
results.r.path <- "results-R/"
results.raw.path <- "results-raw/"
shiny.app.path <- "~/Documents/GitHub/Shiny/microhaplot/" # Path to shinyhaplot app

# --- Genotyping Parameters (for mplot2tgt) ---
# See `microhaplot` documentation for details on these parameters. 
min.read.depth <- 20
AB.min.het     <- 3/7   # Allele balance: min ratio for heterozygous call
AB.max.homo    <- 2/8   # Allele balance: max ratio for homozygous call
# min.AR.het   <- 3/10  # Not used in original script, but available
# max.AR.homo  <- 2/10  # Not used in original script, but available

# --- QA/QC Filtering Thresholds ---

##These will need to be configured for each project after seeing outputs of genotyping success
# Initial loose filter to remove worst offenders
initial.locus.missing.threshold <- 0.90 # Remove loci missing in >90% of individuals
initial.indiv.missing.threshold <- 0.90 # Remove individuals with >90% missing loci

# Final strict filter for analysis-ready dataset
final.locus.missing.threshold <- 0.50 # Keep loci genotyped in at least 50% of individuals
final.indiv.genotyped.threshold <- 0.60 # Keep individuals with at least 60% of loci genotyped

# -------------------------
# --- END CONFIGURATION ---
# -------------------------


# ====================================================================
# STEP 1: RUN MICROHAPLOT TO CALL GENOTYPES
# This step calls `microhaplot` to create the initial haplotype table.
# It can be skipped if you have already generated the .rds file.
# ====================================================================

  haplo.read.tbl <- prepHaplotFiles(
    run.label = run.label,
    sam.path = sam.path,
    out.path = results.r.path,
    label.path = label.path,
    vcf.path = vcf.path,
    app.path = shiny.app.path,
    n.jobs = parallel::detectCores() # Use all available cores
  )
  #You may get a warning that some commands used by microhaplot are deprecated
  #Run shiny app if you want to visualize the data. 
  # runShinyHaplot(shiny.app.path)

# ====================================================================
# STEP 2: CONVERT MICROHAPLOT OUTPUT & INITIAL QA/QC
# Converts `microhaplot` output and performs a first pass of QC.
# Specifically looking at questionable haplotypes (contain N, X or >2 haplotypes) and looking for mismatches between known duplicates individuals
# ====================================================================

# Convert microhaplot RDS file to a "tgt" (tidy genotype table) file
tgt.initial <- mplot2tgt(
  project = project,
  AB.min.het = AB.min.het,
  AB.max.homo = AB.max.homo,
  min.read.depth = min.read.depth
)

# Get total number of loci from the VCF file for calculating missing data (want the number of unique CHROM elements since we may have multiple SNPs per locus)
vcf <- read.vcfR(vcf.path)
num.locs <- length(unique(getCHROM(vcf)))
min.locs.per.ind<-num.locs*0.6

# Identify and visualize questionable genotypes (containing N, X, or >2 haplotypes)
questionable.hap <- sapply(1:nrow(tgt.initial), function(i) {
  (grepl("N", tgt.initial$gt[i]) || grepl("X", tgt.initial$gt[i]) || tgt.initial$num.haps[i] > 2)
})
genos.to.check <- filter(tgt.initial, questionable.hap)
#tgt.filtered <- filter(tgt.initial, !questionable.hap) this line can be used to filter out the questionable haplotypes. we will wait to do this for now. 
table(genos.to.check$locus) #visualise questionable genotypes

if (nrow(genos.to.check) > 0) {
  message(paste("Found ", nrow(genos.to.check), "questionable genotypes."))
  saveRDS(genos.to.check, file = file.path(results.r.path, paste0(project, ".genos.to.check.initial.rds")))
  write.csv(genos.to.check, file = file.path(results.raw.path, paste0(project, "genos.to.check.initial.csv")))
}

#Identify and compare genotypes of known replicates
LABIDs <- unique(tgt.initial$Indiv) %>% substr(start = 1, stop = 8)
replicates <- LABIDs[duplicated(LABIDs)]
mismatches.to.check <-do.call('rbind',lapply(replicates, function(r){
  rep.tgt <- tgt[grep(substr(r, start = 1, stop = 8), tgt$Indiv),]
  mismatches <- compare.replicates(rep.tgt)
}))

if (!is.null(mismatches.to.check) && is.data.frame(mismatches.to.check)) {
  if (nrow(mismatches.to.check) > 0) {
    message(paste("Found", nrow(mismatches.to.check), "mismatched genotypes."))
    saveRDS(mismatches.to.check, file = file.path(results.r.path, paste0(project, ".genotype.mismatches.initial.rds")))
    write.csv(mismatches.to.check, file = file.path(results.raw.path, paste0(project, ".genotype.mismatches.initial.csv")), row.names = FALSE)
  } else {
    message("No mismatched genotypes were found between replicates.")
  }
} else {
  # This handles the case where no replicates were found to begin with
  message("No mismatched genotypes were found between replicates.")
}

############STOP HERE#################
#Open and examine the saved files for genos.to.check and mismatches.to.check
#Save them to the qa/qc file
#We will not remove any loci/individuals yet due to these two checks since we are
#unsure whether they could be impacted by bad quality loci or individuals
#Proceed to the first filtering step (>90% missing loci/individuals)
#Will recalculate and evaluate genos.to.check and mismatches.to.check after this. 

# ====================================================================
# STEP 3: ITERATIVE FILTERING OF LOCI AND INDIVIDUALS
# This section applies filtering thresholds to create the final dataset.
# The thresholds were set in the configuration section at the top. 
# ====================================================================

# --- First Pass Filter (Remove worst offenders) ---

# Summarize and create a list of loci that don't pass first threshold
loc.sum.initial <- tgt.initial %>%
  filter(!is.na(gt)) %>%
  group_by(locus) %>%
  summarise(inds.genoed = n_distinct(Indiv), .groups = 'drop') %>%
  mutate(prop.genoed = inds.genoed / n_distinct(tgt.initial$Indiv))

loci.to.remove.initial <- loc.sum.initial %>%
  filter(prop.genoed < (1 - initial.locus.missing.threshold)) %>%
  pull(locus) %>%
  append("locus019") ##Can manually include in this filter any obviously bad loci from the genos.to.check step before


# Summarize and create a list of individuals that don't pass first threshold
#It looks like this step does not filters out individuals that have fewer genotypes than "min.genos.per.ind"
#This variable was not even included in the new script. 
ind.sum.initial <- tgt.initial %>%
  filter(!is.na(gt)) %>%
  group_by(Indiv) %>%
  summarise(loci.genoed = n_distinct(locus), .groups = 'drop') %>%
  mutate(prop.genoed = loci.genoed / num.locs)

inds.to.remove.initial <- ind.sum.initial %>%
  filter(prop.genoed < (1 - initial.indiv.missing.threshold)) %>%
  pull(Indiv)

# Apply the first filter
tgt.pass1 <- tgt.initial %>%
  filter(!locus %in% loci.to.remove.initial) %>%
  filter(!Indiv %in% inds.to.remove.initial)

message(paste("  - Initial filter removed", length(loci.to.remove.initial), "loci and", length(inds.to.remove.initial), "individuals."))

#######STOP HERE########
#Go back and recalculate the questionable and replicated genotypes now to see if
#Filtering out the worst loci and individuals cleared them up. 

# Identify and visualize questionable genotypes (containing N, X, or >2 haplotypes)
questionable.hap.pass1 <- sapply(1:nrow(tgt.pass1), function(i) {
  (grepl("N", tgt.pass1$gt[i]) || grepl("X", tgt.pass1$gt[i]) || tgt.pass1$num.haps[i] > 2)
})
genos.to.check.pass1 <- filter(tgt.pass1, questionable.hap.pass1)
#tgt.filtered <- filter(tgt.initial, !questionable.hap) this line can be used to filter out the questionable haplotypes. we will wait to do this for now. 
table(genos.to.check.pass1$locus) #visualise questionable genotypes

if (nrow(genos.to.check) > 0) {
  message(paste("Found ", nrow(genos.to.check), "questionable genotypes."))
  saveRDS(genos.to.check, file = file.path(results.r.path, paste0(project, ".genos.to.check.initial.rds")))
  write.csv(genos.to.check, file = file.path(results.raw.path, paste0(project, "genos.to.check.initial.csv")))
}

#Identify and compare genotypes of known replicates
LABIDs <- unique(tgt.initial$Indiv) %>% substr(start = 1, stop = 8)
replicates <- LABIDs[duplicated(LABIDs)]
mismatches.to.check <-do.call('rbind',lapply(replicates, function(r){
  rep.tgt <- tgt[grep(substr(r, start = 1, stop = 8), tgt$Indiv),]
  mismatches <- compare.replicates(rep.tgt)
}))

if (!is.null(mismatches.to.check) && is.data.frame(mismatches.to.check)) {
  if (nrow(mismatches.to.check) > 0) {
    message(paste("Found", nrow(mismatches.to.check), "mismatched genotypes."))
    saveRDS(mismatches.to.check, file = file.path(results.r.path, paste0(project, ".genotype.mismatches.initial.rds")))
    write.csv(mismatches.to.check, file = file.path(results.raw.path, paste0(project, ".genotype.mismatches.initial.csv")), row.names = FALSE)
  } else {
    message("No mismatched genotypes were found between replicates.")
  }
} else {
  # This handles the case where no replicates were found to begin with
  message("No mismatched genotypes were found between replicates.")
}




# --- Second Pass Filter (Apply strict thresholds for final dataset) ---
message("  - Applying final strict filter...")

# Create helper function to summarize locus stats to keep code DRY
summarize_loci <- function(tgt_df) {
  tgt_df %>%
    filter(!is.na(gt)) %>%
    group_by(locus) %>%
    summarise(inds.genoed = n_distinct(Indiv), .groups = 'drop') %>%
    mutate(prop.genoed = inds.genoed / n_distinct(tgt_df$Indiv))
}

# 1. Filter LOCI based on `final.locus.missing.threshold`
loc.sum.pass1 <- summarize_loci(tgt.pass1)
loci.to.keep.final <- loc.sum.pass1 %>%
  filter(prop.genoed >= final.locus.missing.threshold) %>%
  pull(locus)
tgt.pass2 <- filter(tgt.pass1, locus %in% loci.to.keep.final)
message(paste("  - Kept", length(loci.to.keep.final), "loci genotyped in at least", final.locus.missing.threshold * 100, "% of individuals."))

# 2. Filter INDIVIDUALS based on `final.indiv.genotyped.threshold`
num.locs.final <- length(loci.to.keep.final)
ind.sum.pass2 <- tgt.pass2 %>%
  filter(!is.na(gt)) %>%
  group_by(Indiv) %>%
  summarise(loci.genoed = n_distinct(locus), .groups = 'drop') %>%
  mutate(prop.genoed = loci.genoed / num.locs.final)

inds.to.keep.final <- ind.sum.pass2 %>%
  filter(prop.genoed >= final.indiv.genotyped.threshold) %>%
  pull(Indiv)
tgt.final <- filter(tgt.pass2, Indiv %in% inds.to.keep.final)
message(paste("  - Kept", length(inds.to.keep.final), "individuals with at least", final.indiv.genotyped.threshold * 100, "% of loci genotyped."))

# Final check on locus stats after final individual filtering
loc.sum.final <- summarize_loci(tgt.final)

# ====================================================================
# STEP 4: SAVE FINAL OUTPUTS
# ====================================================================
message("Step 4: Saving final outputs...")

# Convert final `tgt` object to a wide "geno.table" format
geno.table <- tgt.2.geno.table(tgt.final)

# Write summary files
write.csv(ind.sum.pass2, file.path(results.raw.path, paste0(project, ".indiv.summary.final.csv")), row.names = FALSE)
write.csv(loc.sum.final, file.path(results.raw.path, paste0(project, ".locus.summary.final.csv")), row.names = FALSE)

# Save the final R objects for the next script
save(geno.table, tgt.final, loc.sum.final, file = file.path(results.r.path, paste0(project, ".final.geno.data.rda")))

message("Workflow Step 1 complete!")
message("Final dataset contains ", nrow(geno.table), " individuals and ", ncol(geno.table) - 1, " loci.")
message("Outputs saved in ", results.r.path, " and ", results.raw.path)
