# GT-seq Analysis Workflow for Population Genetics

This repository contains a standardized and interactive workflow for processing and analyzing SNP data generated via Genotyping-in-Thousands by sequencing (GT-seq). The workflow is designed to be run in stages, guiding the user from raw sequencing data through genotype calling, a multi-step quality control process, and finally to downstream population genetic analyses.

This workflow is **interactive**. After most steps, you will be required to examine the output files (summary tables and plots) to make informed decisions about filtering parameters for the next step.

---

## 🚀 Overall Workflow Philosophy

The core of this pipeline is an iterative filtering process. You don't just filter the data once. You will refine it through several stages to ensure high-quality genotypes for final analysis. The `analysis_stage` variable in `01_genotyping_and_qc_interactive.R` is crucial and controls this process.

The typical progression is as follows:

1. **`wreps` Stage**: Run the initial analysis with all your samples, including technical replicates.
2. **Duplicate Check**: After creating a preliminary `gtypes` object, you'll check for duplicate genotypes (i.e., the same individual sampled multiple times).
3. **Merge & Refine**: Based on the duplicate check and replicate comparisons, you will manually address these samples (e.g., merge, remove) and create a new, definitive sample list.
4. **`final_nodups` Stage**: You will re-run the entire workflow from the beginning using your cleaned sample list to produce the final, high-quality dataset for population genetic analysis.

[Image of a bioinformatics data analysis flowchart]

---

## 📋 Prerequisites

### R Packages

You will need the following R packages. Install them by running this code in your R console:

```r
# For installing packages from GitHub
# install.packages("devtools")

# CRAN Packages
install.packages(c("tidyverse", "vcfR", "readxl", "strataG", "swfscMisc", "purrr", "tibble"))

# GitHub Packages
devtools::install_github("ngthomas/microhaplot", build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
devtools::install_github("kmartien/vcftoolsR") # vcftools wrapper functions
```

### External Software

* **shinyhaplot**: The `microhaplot` package requires a Shiny app component for visualizing genotypes. Please follow the setup instructions on the `microhaplot` GitHub page.

### Recommended Directory Structure

For the scripts to run correctly, you should organize your project as follows:

```
your_project/
├── data/
│   └── (Generated gtypes and .rda files will be stored here)
├── data-raw/
│   ├── sam_files/          # Your raw SAM files
│   ├── vcf/                # Your VCF file of target SNPs
│   └── metadata/
│       └── sample_metadata.xlsx
│       └── sample_labels.txt
├── R/
│   ├── 01_genotyping_and_qc_interactive.R
│   ├── 02_create_gtypes.R
│   ├── 03_check_duplicates.R
│   ├── 04_population_genetics_analysis.R
│   └── functions/
│       ├── mplot2tgt_ABO.R
│       ├── tgt.2.geno.table.R
│       └── ... (other custom functions)
├── results-raw/
│   └── (Generated .csv files and plots will be saved here)
└── results-R/
    └── (Intermediate R objects will be saved here)
```

## ⚙️ The Workflow: Step-by-Step Guide

Run the scripts in the specified order. Before running each script, **carefully review and update the `--- CONFIGURATION ---` section at the top of the file.**

### Step 1: Genotype Calling and Interactive Quality Control

**Script:** `R/01_genotyping_and_qc_interactive.R`

**Goal:** This script is the engine of the workflow. It calls genotypes from your raw SAM files and guides you through a rigorous, multi-part filtering process to remove low-quality individuals and loci.

**Key Operations:**

* **Part 1: Initial Data Generation**:
  * Uses `microhaplot` to process SAM files based on your VCF of target SNPs.
  * Converts the output into a "tidy genotype table" (`tgt`) using custom functions.
  * Runs a preliminary QC analysis (`run_qc_analysis`) to generate initial summary reports and histograms on locus and individual completeness.

* **Part 2: First Pass (Loose) Filtering**:
  * Applies a very lenient filter to remove loci and individuals with extremely high rates of missing data.
  * Re-runs the QC analysis on this loosely filtered dataset to provide a clearer picture for the next step.

* **Part 3: Final (Strict) Filtering**:
  * Applies a second, more stringent set of filters based on your defined thresholds.
  * Crucially, it recalculates individual genotyping rates *after* the final set of loci has been selected, ensuring an accurate final filtering step.
  * Runs the QC analysis one last time on the final, cleaned data.

👉 **User Actions:**

1. Open `01_genotyping_and_qc_interactive.R`.
2. In the **MASTER CONFIGURATION** section, set the `analysis_stage` (start with `"wreps"`).
3. Update all file paths and genotyping parameters.
4. Run the script through **Part 1**.
5. **Stop and Examine**: Go to the `results-raw/wreps/` folder. Look at the CSV files and histograms with the `_initial` suffix. Use these to decide on your loose thresholds.
6. Set the **LOOSE filtering thresholds** in **Part 2** and run it.
7. **Stop and Examine Again**: Review the `_pass1` reports. Use these to decide on your final, stricter thresholds.
8. Set the **FINAL filtering thresholds** in **Part 3** and run the rest of the script.

**Output:** The script saves a final, cleaned `geno.table` and `tgt_final` object in an `.rda` file inside the `results-R/[analysis_stage]/` folder. All QC reports are saved in `results-raw/[analysis_stage]/`.

### Step 2: Create a `gtypes` Object for Analysis

**Script:** `R/02_create_gtypes.R`

**Goal:** To combine your clean genotype data with your sample metadata and format it into a `gtypes` object, which is required for the `strataG` analysis package.

**Key Operations:**

1. **Load Data**: Loads the final `geno.table` produced by Step 1.
2. **Load Metadata**: Reads sample information (e.g., population, location, year) from your specified Excel or CSV file.
3. **Merge & Reshape**: Joins the genotype and metadata tables and reshapes the genotypes into a one-allele-per-column format using `alleleSplit()`.
4. **Create `gtypes` Object**: Converts the merged data frame into a `gtypes` object, which efficiently stores genotypes, individual IDs, and population strata.

👉 **User Actions:**

1. Open `02_create_gtypes.R`.
2. Update the **CONFIGURATION** section to match the `project_name` and `final_analysis_stage` from Step 1.
3. Provide the correct path to your metadata file (`sample.info.path`) and specify the relevant column names for individual IDs and population strata.

**Output:** A `gtypes` object saved as an `.rda` file in the `data/` directory.

### Step 3: Check for Duplicate Genotypes

**Script:** `R/03_check_duplicates.R`

**Goal:** To identify pairs of samples with highly similar genotypes, which may indicate that they are unintentional duplicate samples of the same individual. **This step is critical for ensuring data independence.**

**Key Operations:**

1. **Load `gtypes` Object**: Loads the file created in Step 2.
2. **Compare Genotypes**: Uses the `strataG::dupGenotypes` function to compare all pairs of individuals.
3. **Flag Duplicates**: Identifies and lists pairs that exceed a genetic similarity threshold you define (e.g., 80% shared alleles).

👉 **User Actions:**

1. Open `03_check_duplicates.R`.
2. Update the **CONFIGURATION** to match the previous scripts. Set your desired `similarity.threshold`.
3. Run the script and examine the output CSV in `results-raw/`.
4. **🛑 ACTION REQUIRED**: If duplicates are found, you must:
   * Investigate the pairs and decide which samples to keep, merge, or remove.
   * Create a **new sample list `.txt` file** that reflects these changes.
   * **Return to Step 1** and re-run the entire workflow from the beginning, this time setting `analysis_stage <- "final_nodups"` and providing the path to your new, cleaned sample list.

### Step 4: Population Genetics Analysis

**Script:** `R/04_population_genetics_analysis.R`

**Goal:** To perform the final population genetic analyses on your clean, finalized dataset.

**Key Operations:**

1. **Individual Summaries**: Calculates per-individual homozygosity and generates a histogram to help identify outlier individuals (e.g., due to inbreeding or contamination).
2. **Hardy-Weinberg Equilibrium (HWE)**: Tests for deviations from HWE for each locus within each population, providing both raw and adjusted p-values.
3. **Linkage Disequilibrium (LD)**: Tests for LD between all pairs of loci within each population to check for non-random association of alleles.

👉 **User Actions:**

1. Open `04_population_genetics_analysis.R`.
2. Update the **CONFIGURATION** to match the `final_nodups` run from the previous steps.
3. Set your analysis parameters, such as the `high.homo.threshold` for flagging outlier individuals and the `ld.locus.geno.threshold` for filtering loci before the LD test.
4. Run the script.

**Output:** This script saves multiple `.csv` files (for HWE and LD results) in the `results-raw/final_nodups/` folder and `.rda` files containing the full result objects in the `data/` folder.

---

✅ Once you have completed all steps, your analysis is finished! The output files in the `results-raw` directory contain the key findings of your study.