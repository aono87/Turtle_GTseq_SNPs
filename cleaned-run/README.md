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
    Required files are stored in data-raw/ 
    Script stored in R/
  Required scripts:  
  	test.turtle.gtseq.genotyping.R
  Required R packages: 
  	vcfR, microhaplot
  Summary: Script will install and set up shiny and microhaplot, paths to required files and output locations will be assigned, microhaplot will be run. 
    Final file will be "haplo.read.tbl"
    Final file will also be saved in outpath (Shiny/microhaplot/*.rds)
  
Step 2: Call genotypes and QA/QC
  Required files: 
  	Two microhaplot .rds files 
  Directory layout: 
  	.rds files stored in Shiny/microhaplot
  	Script stored in R/
  	R files stored in results-R/
  	CSV files sotred in results-raw/
  Required scripts: 
  	test.turtle.gtseq.genotyping.R
  Required packages: 
  	vcfR, vcftoolsR (Karen's wrapper), tidyverse, dplyr
  Required functions: 
  	mplot2tgt.R, Compare.replicates.R, tgt.2.geno.table.R, haps.w.ns.R
  Summary: This script takes the outpout of microhaplot and converts it to a tgt file for QA/QC
    If duplicates are present, this script can compare the genotypes and look for mismatches
    Identify questionable genotypes (N, X or more that 2 haplotypes)
    Generate locus and individual summaries to estimate missing data
    Filter out bad loci/individuals
    Re-run script when iteratively removing bad ind/loci until you reach a point where you're happy to move forward with analysis
    Document all output files in a qa.qc excel file
  Filtering steps:
    In the first filter, remove any loci that are missing in >90% of individuals and any individual that was missing >90% of all genotypes.
      In doing this, the worst genotypes are removed and you will likely be left with a reasonable dataset
    After filtering, re-run the questionable genotypes and mismatch samples code and decide how to proceed. 
      In my case, there were no questionable genotypes after the bad loci/samples were removed
    Mismatches were opened in geneious and while the SNP looks like it likely is a heterozygote, the number of minor alleles (MA) was not high enough in one replicate or in the combined sam files to meet the threshold
      MAF threshold: MA must be present in >30% of reads to be considered a heterozygote, MA can't be present in >20% of reads to be a homozygote
      My replicates fell in the 20-30% range
      Ultimately, I was happy that the replicates were good and consistent, so I merged the sam files and re-ran the script
      In doing so, the mismatches were just changed to NAs and removed from the dataset in the loc/ind sum steps. 
  Final dataset:
    After filtering bad individuals and bad loci, the final dataset should be only individuals with >60% genotyped loci, and loci that were genotyped in >50% of individuals. 

Step 3: Create a gtypes file for extra QA/QC and futher analyses
  Required files:
  	Sample info file for with desired strata information
  	A file called "geno.table" which was generated in step 2 from the tgt file
  Directory layout:
  	Sample info read in from QA/QC file stored in results-raw/
  	gtypes file saved in data/
  Required scripts:
  	test.turtle.create.gtypes.R
  Required packages:
  	tidyverse, strataG, dplyr, swfscMisc, readxl
  Required functions:
  Summary:
    This script combines the sample info and desired stratification schemes to the final genotype data. The combined data are used to make a "gtypes" file (from R package strataG), which stores the strata info per sample alongside the genotypes. 
    
Step 4: Summarise final genotype data including homozygosity, Hardy-Weinberg Equilibrium and Linkage Disequilibrium
  Required files:
  	gtypes file from Step 3
  Directory layout:
  	R files stored in results-R/
  	CSV files sotred in results-raw/
  Required scripts:
  	test.turtle.indiv.summaries.HWE.LD.R
  Required packages:
  	tidyverse, strataG, dplyr
  Required functions:
  Summary: This script generates summaries of loci including homozygosity. It also filters out loci that are homozygous in >70% of individuals. This script calculates the likelihood that loci are out of Hardy-Weinberg Equilibrium (HWE) both across all samples and within individual populations. Finally, this script calculates the likelihood that a pair of loci are linked (Linkage Disequilibrium, LD) both across all samples and within individual populations. 
    
    
    