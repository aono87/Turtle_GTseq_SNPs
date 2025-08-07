# 01_genotyping_and_qc_interactive.R
# 
# Workflow for a Two-Step Interactive Filter approach to genotype and filter gtseq genotypes:
# This script uses a multi-part, interactive approach to filter data.
#   - PART 1: Generates initial data and runs a full QC analysis.
#   - PART 2: Applies a first, loose filter based on user input, then
#     re-runs the entire QC analysis on the intermediate dataset.
#   - PART 3: Applies a second, stricter filter to create the final dataset, 
#     then re-runs the entire QC analysis on the final dataset.
#   - PART 4: Saves all final outputs.
#
# Note: more filtering stpes can be added as needed. 
# Note: this is coded with messages for automation, but the interactive approach means this may not be helpful.
#
# INPUTS:
#   - A directory of SAM files, a VCF file, and a sample label file.
#   - Sample names must be z-padded in the format z0xxxxxx (8 digits per sample). Replicates should have an additional b or c after the 8 digit name. (e.g. sample 134108 and 134108b must be z0134108 and z0134108b).
#   - A function to convert sample names to z-padded names is provided if needed: R/functions/zpad_ids.R. 
#   - Function is meant to be run on a formatted label file (filename, sample id, population; no headers)
#
# OUTPUTS:
#   - A cleaned, final `geno.table` with high-quality genotypes.
#   - QC reports (summaries, plots) for the initial, intermediate, and final datasets.
#
# WORKFLOW DESCRIPTION:
#   - For the typical user that has known replicate samples, run PART 1 and PART 2 with all replicates separate. 
#   - After PART 2, evaluate the replicates and merge SAM files as needed. Generate a new sample label file accounting for merged SAM files. 
#   - Re-run PART 1 and PART 2 with merged replicates. Evaluate questionable haplotypes at this point and decide how to address them. 
#   - Run PART 3 and PART 4 after replicates and questionable haplotypes have been accounted for. 

# --- LOAD LIBRARIES ---
library(tidyverse)
library(vcfR)
library(microhaplot)
library(vcftoolsR)

# --- LOAD CUSTOM FUNCTIONS ---
source("R/functions/mplot2tgt_ABO.R")
source("R/functions/tgt.2.geno.table.R")
source("R/functions/Compare.replicates.R")
source("R/functions/haps.w.ns.R")
source("R/functions/QC.helper.function.R") #call run_qc_analysis
source("R/functions/zpad_ids.R") #Only if necessary to convert sample names to required format. 

# ---------------------------
# --- CONFIGURATION ---
# ---------------------------

# --- Project and File Paths ---
#project.wreps <- "dcor_wpac_wreps"
project.merged <- "dcor_wpac_merged" #need to switch project names when the replicates are removed. 
sam.path   <- "data-raw/sam_files/"
#label.path <- "data-raw/metadata/dcor.wpac.labels.csv" #If z-padding required, run Step 0 in PART 1 before setting this path
#label.path.wreps <- "data-raw/metadata/dcor.wpac.labels_zpad.txt" #need to switch label files once the z-padded sample names are updated. 
label.path.merged <- "data-raw/metadata/dcor.wpac.labels_zpad_merged.txt"
vcf.path   <- "data-raw/vcf/Dcor_DcPanel_205_maf.targetSNPs_012224.recode.vcf"
results.r.path <- "results-R/"
results.raw.path <- "results-raw/"
app.path <- "~/Documents/GitHub/Shiny/microhaplot/" # Path to shinyhaplot app.


# --- Genotyping Parameters (for mplot2tgt) ---
min.read.depth <- 20
AB.min.het     <- 3/7 #minor allele must be present in >30% of reads to be called a heterozygote
AB.max.homo    <- 2/8 #minor allele must be present in <20% of reads to be called a homozygote

# ----------------------------------------------------------------------
# --- PART 1: INITIAL DATA GENERATION & QC                           ---
# ----------------------------------------------------------------------
# Step 0: Convert sample IDs to correct z-padded format
#zpad_ids(input_file_path = label.path)

# Step 1A: Run microhaplot
#haplo.read.tbl <- prepHaplotFiles(
  run.label = project.wreps,
  sam.path = sam.path,
  out.path = results.r.path,
  label.path = label.path.wreps,
  vcf.path = vcf.path,
  app.path = app.path,
  n.jobs = parallel::detectCores() # Use all available cores
)

#Uncomment and run this when sam files of replicates have been merged
haplo.read.tbl <- prepHaplotFiles(
  run.label = project.merged,
  sam.path = sam.path,
  out.path = results.r.path,
  label.path = label.path.merged,
  vcf.path = vcf.path,
  app.path = app.path,
  n.jobs = parallel::detectCores() # Use all available cores
)

length(unique(haplo.read.tbl$id)) #175
length(unique(haplo.read.tbl$locus)) #185

#Convert microhaplot output into a tgt file (tidy genotype table)
tgt.initial.wreps <- mplot2tgt(project = project.wreps, out.path=results.r.path, AB.min.het = AB.min.het, AB.max.homo = AB.max.homo, min.read.depth = min.read.depth)
tgt.initial.merged <- mplot2tgt(project = project.merged, out.path=results.r.path, AB.min.het = AB.min.het, AB.max.homo = AB.max.homo, min.read.depth = min.read.depth)

length(unique(tgt.initial.merged$Indiv)) #175
length(unique(tgt.initial.merged$locus)) #185

#Remove any loci and/or individuals that have an NA across all genotypes
tgt.initial.merged <- tgt.initial.merged %>%
  group_by(locus) %>%
  filter(!all(is.na(gt))) %>%
  ungroup()

tgt.initial.merged<- tgt.initial.merged %>%
  group_by(Indiv) %>%
  filter(!all(is.na(gt))) %>%
  ungroup()

length(unique(tgt.initial.merged$locus)) #183
length(unique(tgt.initial.merged$Indiv)) #166

#write initial tgt file to an r file in case it's needed later
write_rds(tgt.initial.wreps, file=paste0("results-R/tgt.initial.", project.wreps, ".rds"))
write_rds(tgt.initial.merged, file=paste0("results-R/tgt.initial.", project.merged, ".rds")) #only use when replicates have been merged

# Step 1B: Run initial QC analysis on the raw data
num.locs.initial.wreps<- tgt.initial.wreps %>% filter(!is.na(gt)) %>% distinct(locus) %>% summarise(x = n()) %>% pull(x)
num.locs.initial.merged<- tgt.initial.merged %>% filter(!is.na(gt)) %>% distinct(locus) %>% summarise(x = n()) %>% pull(x)

run_qc_analysis(tgt.initial.wreps, "initial", num.locs.initial.wreps, project=project.wreps)
run_qc_analysis(tgt.initial.merged, "initial", num.locs.initial.merged, project=project.merged)

# ✅ PART 1 COMPLETE.
#🛑 ACTION REQUIRED FOR PART 2: 
#1. Examine all reports in '", results.raw.path, "' with the '_initial' suffix.
#2. Decide on LOOSE thresholds for a first filtering pass to remove the worst offenders.
#3. Enter your chosen thresholds in the configuration block for Part 2.
#--------------------------------------------------------------------------


# ----------------------------------------------------------------------
# --- PART 2: FIRST PASS FILTERING & QC                              ---
# ----------------------------------------------------------------------

# Step 2A: Define LOOSE filtering thresholds (USER INPUT REQUIRED)
# Goal: Remove clear outliers. E.g., loci genotyped in <10% of individuals
# or individuals with <10% of loci genotyped.
pass1.locus.threshold <- 0.10  # <-- SET your initial LOCUS threshold (proportion genotyped), only keep loci that were genotyped in >10% of individuals
pass1.indiv.threshold <- 0.10  # <-- SET your initial INDIVIDUAL threshold (proportion genotyped), only keep individuals that were genotyped at >10% of loci
pass1.manual.loci.to.remove <- c("locus019") # <-- Add any obviously bad or weird loci, e.g. c("locus019")
pass1.manual.inds.to.remove <- c("") # <-- Add any obviously bad or weird individuals. 

# Step 2B: Apply the first filter "pass1"
#Set the project (with replicates or when replicates are merged)
#project<-project.wreps
project<-project.merged #need to change to this once samples are merged

# Filter loci
loc.sum.initial <- read.csv(file.path(results.raw.path, paste0(project, ".locus.summary.initial.csv")))
loci.to.keep.pass1 <- loc.sum.initial %>%
  filter(prop.genoed >= pass1.locus.threshold) %>%
  filter(!locus %in% pass1.manual.loci.to.remove) %>% #Only include this line if there are loci to remove manually
  pull(locus)
#tgt.pass1.loci.wreps <- filter(tgt.initial.wreps, locus %in% loci.to.keep.pass1)
tgt.pass1.loci.merged <- filter(tgt.initial.merged, locus %in% loci.to.keep.pass1) #use this when samples have been merged

# Filter individuals
ind.sum.initial <- read.csv(file.path(results.raw.path, paste0(project, ".indiv.summary.initial.csv")))
inds.to.keep.pass1 <- ind.sum.initial %>%
  filter(prop.genoed >= pass1.indiv.threshold) %>%
#  filter(!Indiv %in% pass1.manual.inds.to.remove) %>% #Only include this line if there are individuals to remove manually
  pull(Indiv)
tgt.pass1.merged <- filter(tgt.pass1.loci.merged, Indiv %in% inds.to.keep.pass1)

#message(paste("  - First pass removed", length(unique(tgt.initial.wreps$locus)) - length(loci.to.keep.pass1), "loci and",
              length(unique(tgt.initial.wreps$Indiv)) - length(inds.to.keep.pass1), "individuals."))
message(paste("  - First pass removed", length(unique(tgt.initial.merged$locus)) - length(loci.to.keep.pass1), "loci and",
              length(unique(tgt.initial.merged$Indiv)) - length(inds.to.keep.pass1), "individuals."))

#write pass1 tgt file to an r file in case it's needed later
write_rds(tgt.pass1.merged, file=paste0("results-R/tgt.pass1.", project, ".rds"))

# Step 2C: Re-run QC analysis on the intermediate 'pass1' dataset
num.locs.pass1.merged <- length(unique(tgt.pass1.merged$locus))
run_qc_analysis(tgt.pass1.merged, "pass1", num.locs.pass1.merged, project=project)

#✅ PART 2 COMPLETE.
#--------------------------------------------------------------------------
#🛑 ACTION REQUIRED BEFORE PART 3 IF REPLICATE SAMPLES PRESENT: 
#1. Before proceeding, evaluate all replicate samples, merge desired SAM files, generate new sample label file
#2. Switch commented out lines in CONFIGURATION for 'project' and 'label.path'
#3. Re-run PART 1 and PART 2 using merged SAM files and updated sample label file

#🛑 ACTION REQUIRED BEFORE PART 3 (NO REPLICATES PRESENT): 
#1. Examine all NEW reports in '", results.raw.path, "' with the '_pass1' suffix.
#2. Address any questionable haplotypes present in the 'genos.to.check' file. If needed, open SAM files in Geneious to evaluate and decide how to move forward
#3. If individuals or loci need to be removed due to questionable haplotypes, add them to the 'final.manual.loc/inds.to.remove' lists
#4. If genotypes need to be manually changed, use the 'genos_to_change' script
#5. Decide on the FINAL, STRICTER thresholds for your analysis-ready dataset.
#6. Enter these final thresholds in the configuration block for Part 3.
#--------------------------------------------------------------------------


# ----------------------------------------------------------------------
# --- PART 3: FINAL (STRICT) FILTERING                               ---
# ----------------------------------------------------------------------

#Histograms of data
# Step 3A: Define FINAL filtering thresholds (USER INPUT REQUIRED)
final.locus.threshold <- 0.50  # <-- SET your final LOCUS threshold. Only keep loci genotyped at MORE than this percent of individuals.
final.indiv.threshold <- 0.60  # <-- SET your final INDIVIDUAL threshold. Only keep individuals genotyped at MORE than this percent of loci.
final.manual.loci.to.remove <- c("") # <-- Add any other loci to remove
final.manual.inds.to.remove <- c("") # <-- Add any other individuals to remove

# Step 3A-optional: Manually change genotypes of questionable haplotypes from a csv file
#CSV file will need to be made by user in the format: locus, Indiv, gt
#gt should reflect the new genotype
genos_to_change <- read.csv('data-raw/genos_to_change.csv')
for (i in 1:nrow(genos_to_change)){
  idx <- which(tgt$locus == genos_to_change$locus[i] & tgt$Indiv == genos_to_change$Indiv[i])
  tgt$gt[idx] <- genos_to_change$gt[i]
}

# Step 3B: Apply final filters to the 'pass1' dataset
# Filter loci
loc.sum.pass1 <- read.csv(file.path(results.raw.path, paste0(project, ".locus.summary.pass1.csv")))
loci.to.keep.final <- loc.sum.pass1 %>%
  filter(prop.genoed >= final.locus.threshold) %>%
  #filter(!locus %in% final.manual.loci.to.remove) %>% #only include this line if there are loci to remove manually
  pull(locus)
tgt.final.loci <- filter(tgt.pass1, locus %in% loci.to.keep.final)

# Filter individuals based on the new set of final loci
num.locs.final <- length(loci.to.keep.final)
ind.sum.pass1_final <- tgt.final.loci %>% filter(!is.na(gt)) %>% group_by(Indiv) %>%
  summarise(loci.genoed = n_distinct(locus), .groups = 'drop') %>%
  mutate(prop.genoed = loci.genoed / num.locs.final)
inds.to.keep.final <- ind.sum.pass1_final %>%
  filter(prop.genoed >= final.indiv.threshold) %>%
  #filter(!Indiv %in% final.manual.inds.to.remove) %>% #only include this line if there are individiuals to remove manually
  pull(Indiv)
tgt.final <- filter(tgt.final.loci, Indiv %in% inds.to.keep.final)

message(paste("  - Final filter resulted in", length(loci.to.keep.final), "loci and", length(inds.to.keep.final), "individuals."))

# Step 3C: Re-run QC analysis on the 'final' dataset
num.locs.final <- length(unique(tgt.final$locus))
run_qc_analysis(tgt.final, "final", num.locs.final)

# ----------------------------------------------------------------------
# --- PART 4: SAVE FINAL OUTPUTS                                     ---
# ----------------------------------------------------------------------
message("\n--- PART 4: SAVING FINAL OUTPUTS ---")

# Step 4A: Convert final data to geno.table format
geno.table <- tgt.2.geno.table(tgt.final)

# Step 4B: Save final R objects and summary tables
loc.sum.final <- tgt.final %>% filter(!is.na(gt)) %>% group_by(locus) %>%
  summarise(inds.genoed = n_distinct(Indiv), .groups = 'drop') %>%
  mutate(prop.genoed = inds.genoed / n_distinct(tgt.final$Indiv))
ind.sum.final <- ind.sum.pass1_final %>% filter(Indiv %in% inds.to.keep.final)

write.csv(loc.sum.final, file.path(results.raw.path, paste0(project, ".locus.summary.final.csv")), row.names = FALSE)
write.csv(ind.sum.final, file.path(results.raw.path, paste0(project, ".indiv.summary.final.csv")), row.names = FALSE)
save(geno.table, tgt.final, file = file.path(results.r.path, paste0(project, ".final.geno.data.rda")))

message("\n✅ WORKFLOW COMPLETE!")
message("Final dataset contains ", nrow(geno.table), " individuals and ", ncol(geno.table) - 1, " loci.")
message("All outputs saved in '", results.r.path, "' and '", results.raw.path, "'")
