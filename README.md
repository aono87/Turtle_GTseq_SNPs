# Turtle_GTseq_SNPs
Extracting genotypes of SNPs from gtseq panel and conducting popgen analyses


Notes from meeting with Karen Martien on 17 April 2025
-Uses microhaplot to go from final gtseq snp loci to microhaplotype files
-Converts microhaplot into gtypes format (strataG)
-gtypes to genind for DAPC
-gtypes for popstruct analyses

Genotyping Workflow and Steps

Step 1: Set up and run microhaplot in R
  Required files:
    VCF file - This contains the locus information
    SAM files - Directory with the .sam files
    Sample labels file - File formatted to be read by microhaplot (filename sampleID  pop)
      Note: pop typically set to NA at this stage
  Directory layout:
    These required files are stored in data-raw/ 
  Required scripts:  Dcor.microhaplot.analysis.R
    Script stored in R/
  Required packages: vcfR, microhaplot
  Summary: Script will install and set up shiny and microhaplot, paths to required files and output locations will be assigned, microhaplot will be run. 
    Final file will be "haplo.read.tbl"
    Final file will also be saved in outpath (Shiny/microhaplot/*.rds)
  
Step 2: Call genotypes and QA/QC
  Required files: two .rds files saved in Shiny/microhaplot
  Directory layout: Store r files in results-R/, Store csv files in results-raw
  Required scripts: Dcor.call.final.genotypes.R
  Required packages: vcfR, vcftoolsR (Karen's wrapper), tidyverse, dplyr
  Required functions: mplot2tgt.R, Compare.replicates.R, tgt.2.geno.table.R, haps.w.ns.R
  Summary: This script takes the outpout of microhaplot and converts it to a tgt file for QA/QC
    If duplicates are present, this script can compare the genotypes and look for mismatches
    Identify questionable genotypes (N, X or more that 2 haplotypes)
    Generate locus and individual summaries to estimate missing data
    Filter out bad loci/individuals
    Re-run script when iteratively removing bad ind/loci until you reach a point where you're happy to move forward with analysis
    Document all output files in a qa.qc excel file
  Filtering steps:
    In the first filter, I removed any loci that were missing in >90% of individuals and any individual that was missing >90% of all genotypes.
      In doing this, the worst genotypes are removed and you will likely be left with a reasonable dataset
    After filtering, re run the questionable genotypes and mismatch samples and decide how to proceed. 
      In my case, there were no questionable genotypes after the bad loci/samples were removed
    Mismatches were opened in geneious and while the SNP looks like it likely is a heterozygote, the number of minor alleles was not high enough in one replicate or in the combined sam files to meet the threshold
      MAF threshold: MA must be present in >30% of reads to be considered a heterozygote, MA can't be present in >20% of reads to be a homozygote
      My replicates fell in the 20-30% range
      Ultimately, I was happy that the replicates were good and consistent, so I merged the sam files and re-ran the script
      In doing so, the mismatches were just changed to NAs and removed from the dataset in the loc/ind sum steps. 
  Final dataset:
    After filtering bad individuals and bad loci, the final dataset should be only individuals with >60% genotyped loci, and loci that were genotyped in >50% of individuals. 
  Extra step: finding duplicate samples
    There was one replicate that appeared to be a different sample (different genotypes across all loci)
    
    
#Ideas of what to include in a "helpful code/function"

reading and writing R data
writing csv files with and without pasted project name
    
    
    