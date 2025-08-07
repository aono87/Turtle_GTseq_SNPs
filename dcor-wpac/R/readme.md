---
title: "GT-seq Analysis Workflow for Population Genetics"
author: "Your Name"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: united
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# GT-seq Analysis Workflow for Population Genetics

This repository contains a standardized workflow for processing and analyzing SNP data generated via Genotyping-in-Thousands by sequencing (GT-seq). The workflow takes raw sequencing data through several stages: genotype calling with `microhaplot`, quality control (QA/QC), and downstream population genetic analyses, including tests for Hardy-Weinberg Equilibrium (HWE) and Linkage Disequilibrium (LD).

## Prerequisites

### R Packages

You will need the following R packages. You can install them by running the code chunk below in your R console. 

```{r install-packages, eval=FALSE}
# For installing packages from GitHub
# install.packages("devtools")

# CRAN Packages
install.packages(c("tidyverse", "vcfR", "readxl", "strataG"))

# GitHub Packages
devtools::install_github("ngthomas/microhaplot", build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
devtools::install_github("kmartien/vcftoolsR") # vcftools wrapper functions
```

### External Software & Functions

1.  **shinyhaplot**: The `microhaplot` package requires a Shiny app component. Follow the instructions on the `microhaplot` Github page for setting it up: https://github.com/ngthomas/microhaplot.
2.  **Custom Functions**: This workflow assumes you have a directory (e.g., `R/functions/`) containing custom helper functions like `mplot2tgt.R` and `tgt.2.geno.table.R`. Most of these functions are available from Karen Martiens github repository: https://github.com/kmartien/Mnov.gtseq.analysis/tree/main/R/functions.

### Recommended Directory Structure

To ensure the scripts run smoothly, we recommend the following directory structure:

```
your_project/
├── bash/ 
├── data/
│   └── (Microhaplot generated files will be stored here)
├── data-raw/
│   ├── sam_files/
│   ├── vcf/
│   └── metadata/
│       └── sample_metadata.xlsx
├── R/
│   ├── 01_genotyping_and_qc_interactive.R
│   ├── 02_create_gtypes.R
│   ├── 03_population_genetics_analysis.R
│   └── functions/
│       ├── mplot2tgt.R
│       └── ... (other custom functions)
├── results-raw/
│   └── (Generated .csv files will be saved here)
└── results-R/
    └── (Generated .rds and .rda files will be saved here)
```

Prior to starting, create a qa/qc file to store all of the relevant information including metadata, notes about filtering steps and samples/loci removed, and all csv outputs from the scripts. You can color code the tabs by filtering step if desired. 
## The Workflow

The workflow is divided into three main scripts. Run them in the specified order. Before running each script, **carefully review and update the `--- CONFIGURATION ---` section at the top of the file.**

### Step 1: Genotype Calling and Quality Control

**Script:** `R/01_genotyping_and_qc_interactive.R`

This script is the first step in the pipeline. It uses the `microhaplot` package to call genotypes from your raw SAM files and then performs a rigorous and iterative quality control process to filter out low-quality loci and individuals.

**Key Operations:**

1.  **Genotype Calling**: Runs `microhaplot::prepHaplotFiles` to process SAM files based on a VCF file of target SNPs, creating a `haplo.read.tbl`.
2.  **Genotype Conversion**: Converts the `microhaplot` output to a "tgt" ("tidy genotype table") table using custom functions.
3.  **Initial QC**: Identifies and flags questionable genotypes (e.g., containing 'N's, 'X's, or more than two haplotypes), looks for genotyping mismatches between known replicates, calculates locus and individual summaries using a helper function ('QC.helper.function.R')
4.  **Iterative Filtering**:
    * Performs an initial filtering pass to remove individuals and loci with very high rates of missing data (customizable by user).
    * Re-calculates summary statistics using helper function.
    * Applies a second, more stringent filtering based on user-defined thresholds (e.g., keep individuals with >60% of loci genotyped and loci genotyped in >50% of individuals).
    * Re-calculates summary statistics using helper function
5.  **Output**: Saves the final, cleaned `geno.table`, a summary of filtered loci, and the final `tgt` object.

### Step 2: Create a `gtypes` Object

**Script:** `R/02_create_gtypes.R`

This script prepares your data for population genetic analysis. It takes the clean genotype table from Step 1 and combines it with your sample metadata (e.g., population, sampling location, year). It then formats this combined data into a `gtypes` object, which is the required format for the `strataG` package.

**Key Operations:**

1.  **Load Data**: Loads the `geno.table` produced in Step 1.
2.  **Load Metadata**: Reads your sample information from an Excel or CSV file.
3.  **Merge Data**: Joins the genotype data with the sample metadata.
4.  **Create `gtypes` Object**: Converts the merged data frame into a `gtypes` object, embedding the genotype data along with desired stratification schemes (e.g., by population).
5.  **Output**: Saves the final `gtypes` object as an `.rda` file for use in the next step.

### Step 3: Population Genetics Analysis

**Script:** `R/03_population_genetics_analysis.R`

This is the final analytical script. Using the `gtypes` object from Step 2, it performs several standard population genetic analyses.

**Key Operations:**

1.  **Individual Summaries**: Calculates per-individual homozygosity and filters out any outlier individuals with unusually high homozygosity.
2.  **Hardy-Weinberg Equilibrium (HWE)**:
    * Tests for deviations from HWE for each locus within each defined population.
    * Calculates both unadjusted p-values and p-values adjusted for multiple comparisons using the Holm-Bonferroni method.
    * Summarizes the number of populations in which each locus significantly deviates from HWE.
3.  **Linkage Disequilibrium (LD)**:
    * First, filters out loci with a high proportion of missing data within each population.
    * Tests for LD between all pairs of loci within each population.
    * Calculates and summarizes unadjusted and Holm-Bonferroni adjusted p-values.
4.  **Output**: Saves `.csv` files containing the HWE and LD results (both unadjusted and adjusted) and `.rda` files with the full result objects.

## Accessory Step: Check for Duplicate Genotypes

**Script:** `R/04_check_duplicates.R`

This is an optional but recommended script to run after creating your final `gtypes` object in Step 2. It identifies pairs of samples that have highly similar genotypes, which could indicate accidental duplicate sampling of the same individual.

**Key Operations:**

1.  **Load Data**: Loads the final `gtypes` object.
2.  **Compare Genotypes**: Uses the `strataG::dupGenotypes` function to compare all pairs of individuals.
3.  **Flag Duplicates**: Identifies pairs that exceed a user-defined genetic similarity threshold (e.g., 80% shared alleles).
4.  **Output**: Saves a `.csv` file listing the potential duplicate pairs for further investigation.
