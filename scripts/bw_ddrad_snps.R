#######################Global analysis of Z. cav and M. den using one snp per locus##############

# Libraries ---------------------------------------------------------------
############Load required libraries
library(devtools)
library(adegenet)
library(ape)
library(ggplot2)
library(vcfR)
library(hierfstat)
library(pegas)
library(ade4)
library(SNPRelate)
library(qvalue)
library(rgl)
library(dartR)
library(genetics)
library(Rcpp)
library(tidyverse)
library(tess3r)
library(maps)
library(raster)  
library(marmap)
library(LEA)
library(mapplots)
library(RColorBrewer)
library(fields)
library(radiator)
library(rworldmap)
library(strataG)
library(swfscMisc)
library("PopGenome")
library(plyr)
library(cowplot)
library(Demerelate)
library(readxl)
library(ggpattern)
library(Rcpp)
library(RcppArmadillo)
library(assigner)
library(ggtree)
library(poppr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(wesanderson)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggtree)
library(treeio)
library("ggplotgui")

devtools::install_github("bcm-uga/TESS3_encho_sen")
library(Rcpp)
library(tess3r)
# Data Files --------------------------------------------------------------
#######Datafiles
mden43.gl
pop(mden43.gl)
indNames(mden43.gl)

mden43.pop.gl
pop(mden43.pop.gl)
indNames(mden43.gl)

mden43.gl
zczc125.gl

##Import new vcf into R for use in adegenet
#convert vcf file to vcfR file
zczc123.vcfR<-read.vcfR(("~/Dropbox/Phd/Bioinformatics/bw_ddrad_global/SNP_files/no_populations/one_snp_per_locus/zczc_snp_files/zczc123.snps.vcf"), verbose = TRUE)
#convert vcfR files to genlight files
zczc123.gl<-vcfR2genlight(zczc123.vcfR)
#this just shows that it worked and counts missing data, 
zczc123.gl #123 genotypes,  30,479 binary SNPs, size: 6.1 Mb 96228 (2.57 %) missing data
indNames(zczc123.gl)
summary(NA.posi(zczc123.gl))

pop(zczc123.gl)<-as.factor(c("Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Atlantic",
                             "Atlantic",
                             "Mediterranean",
                             "Atlantic",
                             "Mediterranean",
                             "Atlantic",
                             "Indo-Pacific",
                             "Mediterranean",
                             "Mediterranean",
                             "Indo-Pacific",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Atlantic",
                             "Mediterranean",
                             "Indo-Pacific",
                             "Atlantic",
                             "Indo-Pacific",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Mediterranean",
                             "Mediterranean",
                             "Indo-Pacific",
                             "Mediterranean",
                             "Mediterranean",
                             "Mediterranean",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Atlantic",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Mediterranean",
                             "Mediterranean",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Atlantic",
                             "Atlantic",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Mediterranean",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific",
                             "Indo-Pacific"))

zczc123.pop.gl<-zczc123.gl
pop(zczc123.pop.gl)<-as.factor(c("Med_West",
                                 "Med_West",
                                 "Med_West",
                                 "Med_West",
                                 "Med_West",
                                 "Med_West",
                                 "Med_West",
                                 "Med_Corfu",
                                 "Med_West",
                                 "Med_West",
                                 "Med_Corfu",
                                 "Med_West",
                                 "Med_West",
                                 "Med_East",
                                 "Indo_Cent",
                                 "Indo_Cent",
                                 "Atl_NE",
                                 "Atl_CanIs",
                                 "Med_East",
                                 "Atl_NE",
                                 "Med_East",
                                 "Atl_NE",
                                 "Indo_NE",
                                 "Med_East",
                                 "Med_West",
                                 "Indo_NE",
                                 "Med_East",
                                 "Med_Corfu",
                                 "Med_Corfu",
                                 "Med_East",
                                 "Med_East",
                                 "Med_East",
                                 "Atl_CanIs",
                                 "Atl_CanIs",
                                 "Atl_CanIs",
                                 "Atl_CanIs",
                                 "Atl_CanIs",
                                 "Atl_NE",
                                 "Atl_CanIs",
                                 "Atl_CanIs",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_CanIs",
                                 "Atl_NE",
                                 "Atl_NCarib",
                                 "Atl_NE",
                                 "Atl_CanIs",
                                 "Atl_Sp",
                                 "Atl_NE",
                                 "Atl_NCarib",
                                 "Atl_Sp",
                                 "Atl_CanIs",
                                 "Atl_CanIs",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_CanIs",
                                 "Atl_NCarib",
                                 "Atl_CanIs",
                                 "Indo_NE",
                                 "Indo_NE",
                                 "Atl_NCarib",
                                 "Med_East",
                                 "Indo_Sou",
                                 "Atl_NCarib",
                                 "Indo_NE",
                                 "Atl_NCarib",
                                 "Atl_NE",
                                 "Atl_CanIs",
                                 "Med_West",
                                 "Med_West",
                                 "Indo_Sou",
                                 "Med_West",
                                 "Med_East",
                                 "Med_West",
                                 "Atl_CanIs",
                                 "Atl_SCarib",
                                 "Atl_SCarib",
                                 "Atl_NCarib",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_CanIs",
                                 "Atl_NCarib",
                                 "Atl_NCarib",
                                 "Atl_NCarib",
                                 "Atl_NCarib",
                                 "Atl_NE",
                                 "Indo_NE",
                                 "Indo_NE",
                                 "Med_West",
                                 "Med_West",
                                 "Indo_NE",
                                 "Indo_NE",
                                 "Atl_SCarib",
                                 "Atl_SCarib",
                                 "Indo_NE",
                                 "Indo_NE",
                                 "Indo_Cent",
                                 "Indo_NE",
                                 "Indo_NE",
                                 "Indo_NE",
                                 "Indo_NE",
                                 "Indo_NE",
                                 "Indo_Cent",
                                 "Indo_NE",
                                 "Indo_NE",
                                 "Indo_Cent",
                                 "Indo_Sou",
                                 "Indo_Sou",
                                 "Indo_Sou",
                                 "Indo_Sou",
                                 "Indo_Sou",
                                 "Indo_Sou",
                                 "Indo_Sou",
                                 "Med_West",
                                 "Indo_NE",
                                 "Indo_Sou",
                                 "Indo_Sou",
                                 "Indo_Sou"))
pop(zczc123.pop.gl)

##NEW ZCAV ANALYSIS WITH TWO DUPLICATES REMOVED FROM MADEIRA
#Duplicate individuals:  
#Zcav20181_L2-19 and Zcav20181_L7-12
#Zcav20181_L2-12 andZcav20182_L7-20

#CALCULATING MISSINGNESS PER INDIVIDUAL. This counts the number of locations where NAs are in the data per individual and summarizes in a table. 
summary(NA.posi(zczc129.gl))

#Zcav20181_L2-19=329 and Zcav20181_L7-12=306 --> Keep Zcav20181_L7-12
#Zcav20181_L2-12=367 and Zcav20182_L7-20=1675 --> Keep Zcav20181_L2-12

#Remove worse duplicates: Zcav20181_L2-19 and Zcav20182_L7-20

#Interestingly, removing just the two individuals from the adegenet gl file gives different snp number sthan if you remove the individuals from the vcf file in vcftools and then re-import in as a gl. am going to use the re-imported file and make new adegenet files. 
#This is best to do with the DartR because when it removes the indivdiuals, 
#it also re-calculates call rates and genotypes and will remove any loci that become monomorphic. 
##Remove these indivdiuals from the GenLight file
#gl.make.recode.ind(zczc129.gl, outfile="new_ind_assignments.csv", outpath = "~/Dropbox/Phd/Bioinformatics/bw_ddrad/")
##Open .csv file in excel and in second column, replace sample name that needs to be removed with "Delete". Save file. 
#zczc127.gl<-gl.recode.ind(zczc129.gl, ind.recode="~/Dropbox/Phd/Bioinformatics/bw_ddrad/new_ind_assignments.csv")
##Check that individuals were removed
#zczc127.gl
#indNames(zczc127.gl)
##Re-calculate locus metadata with bad individuals removed
#zczc127.gl<-gl.recalc.metrics(zczc127.gl)
#zczc129.gl #129 individuals, 25059 snps
#zczc127.gl #127 individuals, 25057 snps
##Calculate how many SNPs per individual now
#summary(NA.posi(zczc127.gl))

##Import new vcf into R for use in adegenet
#convert vcf file to vcfR file
zczc127.2.vcfR<-read.vcfR(("~/Dropbox/Phd/Bioinformatics/bw_ddrad/SNP_files/no_populations/one_snp_per_locus/zczc_snp_files/zczc127.snps.vcf"), verbose = TRUE)
#convert vcfR files to genlight files
zczc127.2.gl<-vcfR2genlight(zczc127.2.vcfR)
#this just shows that it worked and counts missing data, 
zczc127.2.gl
indNames(zczc127.2.gl)
summary(NA.posi(zczc127.2.gl))
zczc127.gl<-zczc127.2.gl

##Adding population data to zczc127.gl
pop(zczc127.gl)<-as.factor(c("Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Indopacific", 
                             "Indopacific", 
                             "Atlantic", 
                             "Atlantic", 
                             "Mediterranean", 
                             "Atlantic", 
                             "Mediterranean", 
                             "Atlantic", 
                             "Indopacific", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Indopacific", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Indopacific", 
                             "Indopacific", 
                             "Atlantic", 
                             "Mediterranean", 
                             "Indopacific", 
                             "Atlantic", 
                             "Indopacific", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Indopacific", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Atlantic", 
                             "Indopacific", 
                             "Indopacific", 
                             "Atlantic", 
                             "Atlantic", 
                             "Mediterranean", 
                             "Mediterranean", 
                             "Indopacific", 
                             "Indopacific", 
                             "Atlantic", 
                             "Atlantic", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Mediterranean", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific", 
                             "Indopacific"))
popNames(zczc127.gl)
pop(zczc127.gl)
nPop(zczc127.gl)

zczc127.pop.gl<-zczc127.gl
pop(zczc127.pop.gl)<-as.factor(c("Med-Italy", 
                                 "Med-Italy", 
                                 "Med-Italy", 
                                 "Med-Italy", 
                                 "Med-Italy", 
                                 "Med-Italy", 
                                 "Med-France", 
                                 "Med-East", 
                                 "Med-Italy", 
                                 "Med-East", 
                                 "Med-Italy", 
                                 "Med-East", 
                                 "Med-France", 
                                 "Med-Italy", 
                                 "Med-Italy", 
                                 "Indo-Central", 
                                 "Indo-Central", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Med-East", 
                                 "Atl-East", 
                                 "Med-East", 
                                 "Atl-Madeira", 
                                 "Indo-East", 
                                 "Med-Italy", 
                                 "Med-Italy", 
                                 "Indo-East", 
                                 "Med-East", 
                                 "Med-East", 
                                 "Med-East", 
                                 "Med-East", 
                                 "Med-East", 
                                 "Med-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-Madeira", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-Other", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-Spain", 
                                 "Atl-East", 
                                 "Atl-Other", 
                                 "Atl-Spain", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-Other", 
                                 "Atl-Bahamas", 
                                 "Atl-East", 
                                 "Indo-East", 
                                 "Indo-East", 
                                 "Atl-Bahamas", 
                                 "Med-East", 
                                 "Indo-Spac", 
                                 "Atl-Bahamas", 
                                 "Indo-Mexico", 
                                 "Atl-Bahamas", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Med-Italy", 
                                 "Med-Italy", 
                                 "Indo-Spac", 
                                 "Med-Italy", 
                                 "Med-Italy", 
                                 "Med-Italy", 
                                 "Atl-East", 
                                 "Atl-Other", 
                                 "Atl-Other", 
                                 "Atl-Bahamas", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-East", 
                                 "Atl-Bahamas", 
                                 "Atl-Bahamas", 
                                 "Atl-Bahamas", 
                                 "Atl-Bahamas", 
                                 "Atl-East", 
                                 "Indo-East", 
                                 "Indo-East", 
                                 "Atl-Bahamas", 
                                 "Atl-Bahamas", 
                                 "Med-Italy", 
                                 "Med-Italy", 
                                 "Indo-East", 
                                 "Indo-East", 
                                 "Atl-Other", 
                                 "Atl-Other", 
                                 "Indo-East", 
                                 "Indo-East", 
                                 "Indo-Central", 
                                 "Indo-East", 
                                 "Indo-East", 
                                 "Indo-East", 
                                 "Indo-East", 
                                 "Indo-East", 
                                 "Indo-Central", 
                                 "Indo-Mexico", 
                                 "Indo-East", 
                                 "Indo-Central", 
                                 "Indo-South", 
                                 "Indo-South", 
                                 "Indo-Spac", 
                                 "Indo-Spac", 
                                 "Indo-South", 
                                 "Indo-South", 
                                 "Indo-South", 
                                 "Med-Italy", 
                                 "Indo-East", 
                                 "Indo-South", 
                                 "Indo-South", 
                                 "Indo-Spac"))

popNames(zczc127.pop.gl)
pop(zczc127.pop.gl)
nPop(zczc127.pop.gl)

##Creating new Zcav datafiles with the re-imported zczc127 data

zczc127.gi<-gl2gi(zczc127.gl)
zczc127.gt<-genind2gtypes(zczc127.gi)

zczc127.pop.gi<-gl2gi(zczc127.pop.gl)
zczc127.pop.gt<-genind2gtypes(zczc127.pop.gi)

##Found out that there were duplicates in Bahamas samples as well. Need to remove these samples from database:
#Zcav20182_L6-4 & Zcav20182_L6-9

#will use vcftools to remove the samples and then re-import back into R. 
##Import new vcf into R for use in adegenet
#convert vcf file to vcfR file
zczc125.vcfR<-read.vcfR(("~/Dropbox/Phd/Bioinformatics/bw_ddrad/SNP_files/no_populations/one_snp_per_locus/zczc_snp_files/zczc125.snps.vcf"), verbose = TRUE)
#convert vcfR files to genlight files
zczc125.gl<-vcfR2genlight(zczc125.vcfR)
#this just shows that it worked and counts missing data, 
zczc125.gl
indNames(zczc125.gl)
summary(NA.posi(zczc125.gl))
zczc127.gl<-zczc125.gl

##Adding population data to zczc127.gl
# pop(zczc125.gl)<-as.factor(c("Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Mediterranean", 
#                              "Atlantic", 
#                              "Mediterranean", 
#                              "Atlantic", 
#                              "Indopacific", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Indopacific", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Atlantic", 
#                              "Mediterranean", 
#                              "Indopacific", 
#                              "Atlantic", 
#                              "Indopacific", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Indopacific", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Mediterranean", 
#                              "Mediterranean", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Atlantic", 
#                              "Atlantic", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Mediterranean", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific", 
#                              "Indopacific"))
# popNames(zczc125.gl)
# pop(zczc125.gl)
# nPop(zczc125.gl)
# 
# zczc125.pop.gl<-zczc125.gl
# pop(zczc125.pop.gl)<-as.factor(c("Med-Italy", 
#                                  "Med-Italy", 
#                                  "Med-Italy", 
#                                  "Med-Italy", 
#                                  "Med-Italy", 
#                                  "Med-Italy", 
#                                  "Med-France", 
#                                  "Med-East", 
#                                  "Med-Italy", 
#                                  "Med-East", 
#                                  "Med-Italy", 
#                                  "Med-East", 
#                                  "Med-France", 
#                                  "Med-Italy", 
#                                  "Med-Italy", 
#                                  "Indo-Central", 
#                                  "Indo-Central", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Med-East", 
#                                  "Atl-East", 
#                                  "Med-East", 
#                                  "Atl-Madeira", 
#                                  "Indo-East", 
#                                  "Med-Italy", 
#                                  "Med-Italy", 
#                                  "Indo-East", 
#                                  "Med-East", 
#                                  "Med-East", 
#                                  "Med-East", 
#                                  "Med-East", 
#                                  "Med-East", 
#                                  "Med-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-Madeira", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-Other", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-Spain", 
#                                  "Atl-East", 
#                                  "Atl-Other", 
#                                  "Atl-Spain", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-Other", 
#                                  "Atl-Bahamas", 
#                                  "Atl-East", 
#                                  "Indo-East", 
#                                  "Indo-East", 
#                                  "Atl-Bahamas", 
#                                  "Med-East", 
#                                  "Indo-Spac", 
#                                  "Atl-Bahamas", 
#                                  "Indo-Mexico", 
#                                  "Atl-Bahamas", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Med-Italy", 
#                                  "Med-Italy", 
#                                  "Indo-Spac", 
#                                  "Med-Italy", 
#                                  "Med-Italy", 
#                                  "Med-Italy", 
#                                  "Atl-East", 
#                                  "Atl-Other", 
#                                  "Atl-Other", 
#                                  "Atl-Bahamas", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-East", 
#                                  "Atl-Bahamas", 
#                                  "Atl-Bahamas", 
#                                  "Atl-Bahamas", 
#                                  "Atl-Bahamas", 
#                                  "Atl-East", 
#                                  "Indo-East", 
#                                  "Indo-East", 
#                                  "Med-Italy", 
#                                  "Med-Italy", 
#                                  "Indo-East", 
#                                  "Indo-East", 
#                                  "Atl-Other", 
#                                  "Atl-Other", 
#                                  "Indo-East", 
#                                  "Indo-East", 
#                                  "Indo-Central", 
#                                  "Indo-East", 
#                                  "Indo-East", 
#                                  "Indo-East", 
#                                  "Indo-East", 
#                                  "Indo-East", 
#                                  "Indo-Central", 
#                                  "Indo-Mexico", 
#                                  "Indo-East", 
#                                  "Indo-Central", 
#                                  "Indo-South", 
#                                  "Indo-South", 
#                                  "Indo-Spac", 
#                                  "Indo-Spac", 
#                                  "Indo-South", 
#                                  "Indo-South", 
#                                  "Indo-South", 
#                                  "Med-Italy", 
#                                  "Indo-East", 
#                                  "Indo-South", 
#                                  "Indo-South", 
#                                  "Indo-Spac"))
# 
# popNames(zczc125.pop.gl)
# pop(zczc125.pop.gl)
# nPop(zczc125.pop.gl)


# Allele Frequency Spectra ------------------------------------------------
#####Allele Frequency Spectra. Distribution of counts of alternative (non-reference) alleles across all sites in the population. 

##AFS for all M. den. 
mden43.sum<-glSum(mden43.gl, alleleAsUnit = TRUE)
barplot(table(mden43.sum), xlab="Allele Counts", main="Distribution of ALT allele counts in total dataset")

##AFS for M. den, each ocean
mden43.seppop.gl<-seppop(mden43.gl, drop=TRUE)
mden43.seppop.gl$Atlantic
mden43.seppop.gl$Indopacific

#now must remove nonvarint positions within the population:
#How many alternative alleles at each locus?
n.alleles.atl<-colSums(as.matrix(mden43.seppop.gl$Atlantic))
n.alleles.pac<-colSums(as.matrix(mden43.seppop.gl$Indopacific))

#how many particular categories of alternative allele counts are in my population?
summary(as.factor(n.alleles.atl))
summary(as.factor(n.alleles.pac))

#remove the reference-only position AND any columns with NAs
mden43.seppop.atl.gl<-new("genlight", (as.matrix(mden43.seppop.gl$Atlantic))[,(colSums(as.matrix(mden43.seppop.gl$Atlantic))>0) & (colSums(is.na(as.matrix(mden43.seppop.gl$Atlantic))) == 0)])
mden43.seppop.atl.gl

mden43.seppop.pac.gl<-new("genlight", (as.matrix(mden43.seppop.gl$Indopacific))[,(colSums(as.matrix(mden43.seppop.gl$Indopacific))>0) & (colSums(is.na(as.matrix(mden43.seppop.gl$Indopacific))) == 0)])
mden43.seppop.pac.gl

#check that there are no zeros
summary(colSums(as.matrix(mden43.seppop.atl.gl)))
summary(colSums(as.matrix(mden43.seppop.pac.gl)))

#Plot Afs for each population
mden.atl.afs.sum<-glSum(mden43.seppop.atl.gl, alleleAsUnit = TRUE)
barplot(table(mden.atl.afs.sum), xlab="Allele Counts", main="Distribution of ALT allele counts in the Atlantic")

mden.pac.afs.sum<-glSum(mden43.seppop.pac.gl, alleleAsUnit = TRUE)
barplot(table(mden.pac.afs.sum), xlab="Allele Counts", main="Distribution of ALT allele counts in the Indopacific")

###########lapply function to do this over all populations:
# Now, plot AFS for all populations in a batch and save this into a pdf file (i.e., using lapply function
#                                                                             which goes over all elements of the list of genlights)
# #### plot AFS for all pops in a batch
# aa.genlight.sep <- seppop(aa.genlight, drop=TRUE) # separate genlight per
# population
# # remove the nonvariant positions AND columns with NA within that pop.
# aa.genlight.sep.2 <- lapply (aa.genlight.sep, function (pop)
# {new("genlight", (as.matrix(pop))[,(colSums(as.matrix(pop)) > 0)
#                                   & (colSums(is.na(as.matrix(pop))) == 0)])})
# ##add pop identity to list elements
# listnames<-names(aa.genlight.sep.2)
# for (i in seq(listnames)) {pop(aa.genlight.sep.2[[i]])<-
#   substr(indNames(aa.genlight.sep.2[[i]]),1,3)}
# # loop over each population in a list of populations and draw AFS into one
# fig
# pdf("AFS_all_barplot.pdf", width=5, height=5)
# par(mfrow=c(2,3),mar=c(2,2,2,0))
# mySum <- lapply (aa.genlight.sep.2, function (pop) {
#   barplot(table(glSum(pop, alleleAsUnit=T)), col="blue", space=0,
#           xlab="Allele counts",
#           main=paste(levels(pop(pop)),sum(table(glSum(pop, alleleAsUnit=T))),"SNPs",
#                      sep=" "))
# })
# dev.off()
# par(mfrow=c(1,1))




# Summary Statistics dartR------------------------------------------------------
#######dartR Summary statistics
library(dartR)
#data: 
mden43.gl
zczc129.gl
zczc127.gl
zczc125.gl

#splitting genlight files into ocean level files
mden43.gl<-seppop(mden43.gl)
mden43.atl.gl<-mden43.gl$Atlantic
mden43.pac.gl<-mden43.gl$Indopacific
mden43.atl.gl
mden43.pac.gl

zczc129.gl<-seppop(zczc129.gl)
zczc129.atl.gl<-zczc129.gl$Atlantic
zczc129.pac.gl<-zczc129.gl$Indopacific
zczc129.med.gl<-zczc129.gl$Mediterranean
zczc129.atl.gl
zczc129.med.gl
zczc129.pac.gl

mden43.pop.gl
mden43.pop.gl<-seppop(mden43.pop.gl)
mden.atl.bah.gl<-mden43.pop.gl$"Atl-Bahamas"
mden.atl.east.gl<-mden43.pop.gl$"Atl-East"
mden.atl.oth.gl<-mden43.pop.gl$"Atl-Other"
mden.pac.afr.gl<-mden43.pop.gl$"Indo-Africa"
mden.pac.haw.gl<-mden43.pop.gl$"Indo-Hawaii"
mden.pac.sou.gl<-mden43.pop.gl$`Indo-South`

zczc129.pop.gl
zczc129.pop.gl<-seppop(zczc129.pop.gl)
zczc.atl.bah.gl<-zczc129.pop.gl$`Atl-Bahamas`
zczc.atl.east.gl<-zczc129.pop.gl$`Atl-East`
zczc.atl.mad.gl<-zczc129.pop.gl$`Atl-Madeira`
zczc.atl.oth.gl<-zczc129.pop.gl$`Atl-Other`
zczc.atl.sp.gl<-zczc129.pop.gl$`Atl-Spain`
zczc.pac.cen.gl<-zczc129.pop.gl$`Indo-Central`
zczc.pac.east.gl<-zczc129.pop.gl$`Indo-East`
zczc.pac.mex.gl<-zczc129.pop.gl$`Indo-Mexico`
zczc.pac.sou.gl<-zczc129.pop.gl$`Indo-South`
zczc.pac.spac.gl<-zczc129.pop.gl$`Indo-Spac`
zczc.med.east.gl<-zczc129.pop.gl$`Med-East`
zczc.med.fra.gl<-zczc129.pop.gl$`Med-France`
zczc.med.it.gl<-zczc129.pop.gl$`Med-Italy`

#table of frequencies of alternate alleles per locus
table(glSum(mden43.gl))
table(glSum(zczc129.gl))
table(glSum(zczc127.gl))
table(glSum(zczc125.gl))

#table of frequencies of NAs per locus
table(glNA(mden43.gl))
table(glNA(zczc129.gl))

#call rate by locus
gl.report.callrate(mden43.gl)
gl.report.callrate(zczc129.gl)

#global stats using gl.basic.stats in dartR. 
mden43.stats<-gl.basic.stats(mden43.gl)
mden43.stats
mden43.stats$overall

zczc129.stats<-gl.basic.stats(zczc129.gl)
zczc129.stats$overall

zczc127.stats<-gl.basic.stats(zczc127.gl)
zczc127.stats$overall

zczc125.stats<-gl.basic.stats(zczc125.gl)
zczc125.stats

#ocean stats
#each row is a locus and each column is the population. the locus is the rowname. 
mden43.stats$Hs[1:3, 1:2]
mden43.stats$Ho[1:3, 1:2]
mden43.stats$ho[1:3, 1:2]
mden43.stats$n.ind.samp[1:3,1:2]
#each row is a locus and each column is a different metric. the locus ts the rowname. 
mden43.stats$perloc[1:3,1:3]
#one row with overall metrics for all samples. each column is a different metric. 
mden43.stats$overall

zczc129.stats$Hs[1:3, 1:3]
zczc129.stats$Ho[1:3,1:3]
zczc129.stats$ho[1:3,1:3]
zczc129.stats$n.ind.samp[1:3,1:3]
zczc129.stats$perloc[1:3,1:3]
zczc129.stats$overall

##########For some reason, the basic stats function does not like genlight files that have been through the seppop function. 

#population stats
mden43.pop.gl
zczc129.pop.gl
zczc127.pop.gl
zczc125.pop.gl


mden43.pop.stats<-gl.basic.stats(mden43.pop.gl)
zczc129.pop.stats<-gl.basic.stats(zczc129.pop.gl)
zczc127.pop.stats<-gl.basic.stats(zczc127.pop.gl)
zczc125.pop.stats<-gl.basic.stats(zczc125.pop.gl)
zczc125.pop.stats

mden43.pop.stats$Ho[1:3,1:7]
zczc129.pop.stats$Ho[1:3,1:13]
boxplot(mden43.pop.stats$Ho)
boxplot(mden43.pop.stats$Hs)
boxplot(mden43.pop.stats$ho)

boxplot(zczc129.pop.stats$ho)
boxplot(zczc129.stats$ho)

boxplot(mden43.pop.stats$ho)
boxplot(mden43.stats$ho)

zczc129.pop.ho<-zczc129.pop.stats$Ho
zczc129.pop.ho[zczc129.pop.ho == 0]<-NA
boxplot(zczc129.pop.ho)
colSums(zczc129.pop.ho)
colSums(zczc129.pop.ho, na.rm=TRUE)

###Loop to count how many non-zero values are in the HO file. 
df <- data.frame(snps=seq(1,13, by=1),pop=seq(1,13, by=1) )
for (i in 1:ncol(zczc129.pop.ho)){
  d <- data.frame(zczc129.pop.ho[,i])
  temp <- is.na(d)
  temp <- sum(temp==FALSE)
  df$snps[i] <- temp
  df$pop[i] <- colnames(zczc129.pop.ho)[i]
  
}
df

###Summary stats files
mden43.stats #All 43 M. den samples, organized into two ocean basins: Atlantic and Indopacific
mden43.pop.stats #All 43 M. den samples, organized into populations within ocean basins: Atl-Bahamas, Atl-East, Atl-Other, Indo-Africa, Indo-Hawaii, Indo-south, NA
zczc129.stats #All 129 Z. cav samples, organized into ocean basins: Atlantic, Indopacific, Mediterranean
zczc129.pop.stats #All 129 Z. cav samples, organized into populations within ocean basins:  Med-Italy, Atl-East, Med-France, Atl-Madeira, Atl-Other, Med-East, Atl-Spain, Indo-South, Indo-East, Indo-Spac, Atl-Bahamas, Indo-Central, Indo-Mexico

mden43.stats$overall
mden43.pop.stats$overall
zczc129.stats$overall
zczc129.pop.stats$overall
zczc127.stats$overall
zczc127.pop.stats$overall

summary(zczc127.stats$Ho)
summary(zczc127.stats$Hs)
summary(zczc127.stats$ho)
summary(zczc127.pop.stats$Ho)
summary(zczc127.pop.stats$Hs)
summary(zczc127.pop.stats$ho)

summary(mden43.stats$Ho)
summary(mden43.stats$Hs)
summary(mden43.stats$ho)
summary(mden43.pop.stats$Ho)
summary(mden43.pop.stats$Hs)
summary(mden43.pop.stats$ho)
summary(zczc129.stats$Ho)
summary(zczc129.stats$Hs)
summary(zczc129.stats$ho)
summary(zczc129.pop.stats$Ho)
summary(zczc129.pop.stats$Hs)
summary(zczc129.pop.stats$ho)

##AMOVA in dartR
mden43.amova<-gl.amova(mden43.gl, nperm=1000)
mden43.pop.amova<-gl.amova(mden43.pop.gl, nperm=1000)
zczc129.amova<-gl.amova(zczc129.gl, nperm=1000)
zczc129.pop.amova<-gl.amova(zczc129.pop.gl, nperm=1000)

summary(mden43.amova)



##Diversity indices in dartR
mden43.div<-gl.diversity(mden43.gl)
mden43.pop.div<-gl.diversity(mden43.pop.gl)
zczc129.div<-gl.diversity(zczc129.gl)
zczc129.pop.div<-gl.diversity(zczc129.pop.gl)
zczc127.div<-gl.diversity(zczc127.gl)

mden43.div
##Distance maps and heatmaps in dartR
mden43.dist<-gl.dist.pop(mden43.gl)
mden43.pop.dist<-gl.dist.pop(mden43.pop.gl) #Error in `.rowNamesDF<-`(x, value = value) : missing values in 'row.names' are not allowed
#This error may be because one individual didnt assign to a population and was assigned an NA. 
#will remove this as a population, make a new genlight file and try
#this worked!
mden42.pop.gl<-gl.drop.pop(mden43.pop.gl, pop.list=c("NA"))
mden42.pop.dist<-gl.dist.pop(mden42.pop.gl)
mden42.pop.dist
zczc129.dist<-gl.dist.pop(zczc129.gl)
zczc129.pop.dist<-gl.dist.pop(zczc129.pop.gl)
#heatmap code
gl.dist.heatmap(mden42.pop.dist, rank=TRUE)



# Fst CI dartR ------------------------------------------------------------
#calculating Fst values with 95% confidenc interval
mden43.boot.fst<-gl.fst.pop(mden43.gl, nboots=100, percent=95, nclusters=1)
mden43.boot.fst
# Lower bound CI limit  Upper bound CI limit    p-value             Fst
# 0.1264371             0.1133691               0.123526       0 0.1189114
# $Fsts
# Atlantic Indopacific
# Atlantic           NA          NA
# Indopacific 0.1189114          NA
# 
# $Pvalues
# Atlantic Indopacific
# Atlantic          NA          NA
# Indopacific        0          NA

mden43.pop.boot.fst<-gl.fst.pop(mden43.pop.gl, nboots=100, percent=95, nclusters=1)
# $Fsts
# Atl-East  Atl-Other Atl-Bahamas Indo-Hawaii  Indo-Africa   Indo-South NA
# Atl-East            NA         NA          NA          NA           NA           NA NA
# Atl-Other   0.01575297         NA          NA          NA           NA           NA NA
# Atl-Bahamas 0.04029976 0.01485074          NA          NA           NA           NA NA
# Indo-Hawaii 0.16212013 0.13287489  0.16127335          NA           NA           NA NA
# Indo-Africa 0.10472668 0.07008143  0.09975008  0.02681006           NA           NA NA
# Indo-South  0.16271209 0.13415492  0.16229989  0.01307980  0.009237563           NA NA
# NA          0.18999090 0.17061198  0.20102479  0.01005949 -0.016177580 0.0009779845 NA
# 
# $Pvalues
# Atl-East Atl-Other Atl-Bahamas Indo-Hawaii Indo-Africa Indo-South NA
# Atl-East          NA        NA          NA          NA          NA         NA NA
# Atl-Other          0        NA          NA          NA          NA         NA NA
# Atl-Bahamas        0         0          NA          NA          NA         NA NA
# Indo-Hawaii        0         0           0          NA          NA         NA NA
# Indo-Africa        0         0           0        0.00          NA         NA NA
# Indo-South         0         0           0        0.00        0.00         NA NA
# NA                 0         0           0        0.16        0.98       0.43 NA

# Population1 Population2 Lower bound CI limit Upper bound CI limit p-value           Fst
# 1     Atl-East   Atl-Other          0.012402064          0.018762808    0.00  0.0157529730
# 2     Atl-East Atl-Bahamas          0.036933317          0.043596201    0.00  0.0402997572
# 3     Atl-East Indo-Hawaii          0.154042688          0.167656353    0.00  0.1621201322
# 4     Atl-East Indo-Africa          0.100187640          0.109860827    0.00  0.1047266819
# 5     Atl-East  Indo-South          0.153208437          0.169678111    0.00  0.1627120867
# 6     Atl-East          NA          0.173216040          0.202998919    0.00  0.1899909001
# 7    Atl-Other Atl-Bahamas          0.010453935          0.018767517    0.00  0.0148507380
# 8    Atl-Other Indo-Hawaii          0.124662446          0.138942051    0.00  0.1328748873
# 9    Atl-Other Indo-Africa          0.064757166          0.075380172    0.00  0.0700814310
# 10   Atl-Other  Indo-South          0.124920657          0.140460825    0.00  0.1341549181
# 11   Atl-Other          NA          0.150023254          0.186100383    0.00  0.1706119806
# 12 Atl-Bahamas Indo-Hawaii          0.153429370          0.170219524    0.00  0.1612733478
# 13 Atl-Bahamas Indo-Africa          0.094380000          0.106380114    0.00  0.0997500759
# 14 Atl-Bahamas  Indo-South          0.152539832          0.169774842    0.00  0.1622998864
# 15 Atl-Bahamas          NA          0.184078892          0.214194852    0.00  0.2010247885
# 16 Indo-Hawaii Indo-Africa          0.023091940          0.031002517    0.00  0.0268100609
# 17 Indo-Hawaii  Indo-South          0.008171821          0.018034605    0.00  0.0130797977
# 18 Indo-Hawaii          NA         -0.009448912          0.027513390    0.16  0.0100594860
# 19 Indo-Africa  Indo-South          0.003550439          0.014004439    0.00  0.0092375633
# 20 Indo-Africa          NA         -0.033453984         -0.004014072    0.98 -0.0161775804
# 21  Indo-South          NA         -0.018314698          0.013746505    0.43  0.0009779845

zczc125.boot.fst<-gl.fst.pop(zczc125.gl, nboots=100, percent=95, nclusters=1)
# Population1	Population2	lower bound CI	upper bound CI	p-value	Fst
# Mediterranean	Indopacific	0.18579501	0.196399	0	0.19097749
# Mediterranean	Atlantic	0.1650606	0.1739923	0	0.16940166
# Indopacific	Atlantic	0.01679399	0.0186115	0	0.01783218

zczc125.pop1.boot.fst<-gl.fst.pop(zczc125.pop1.gl, nboots=100, percent=95, nclusters=1)
# $Fsts
# Med_A       Med_C     Med_B     Indo_A     Atl_AD     Atl_AC     Indo_C      Atl_AA    Atl_AE     Atl_AB Indo_B
# Med_A          NA          NA        NA         NA         NA         NA         NA          NA        NA         NA     NA
# Med_C  0.07267638          NA        NA         NA         NA         NA         NA          NA        NA         NA     NA
# Med_B  0.08884050 0.007056218        NA         NA         NA         NA         NA          NA        NA         NA     NA
# Indo_A 0.24317274 0.166374501 0.2698512         NA         NA         NA         NA          NA        NA         NA     NA
# Atl_AD 0.19403043 0.131194896 0.2085442 0.02886441         NA         NA         NA          NA        NA         NA     NA
# Atl_AC 0.20016537 0.134161394 0.2152729 0.03195148 0.00544131         NA         NA          NA        NA         NA     NA
# Indo_C 0.21190273 0.151797471 0.2273839 0.01775912 0.02431115 0.02957729         NA          NA        NA         NA     NA
# Atl_AA 0.21543057 0.142783897 0.2331287 0.03535070 0.01252524 0.01235090 0.03487647          NA        NA         NA     NA
# Atl_AE 0.36619480 0.327572194 0.4207762 0.21718246 0.14252337 0.14969321 0.17188012 0.156176759        NA         NA     NA
# Atl_AB 0.16172582 0.070547931 0.1693094 0.02794257 0.01004736 0.01235550 0.03325936 0.006053763 0.1704945         NA     NA
# Indo_B 0.20869151 0.141220666 0.2257712 0.01659260 0.01114737 0.01642387 0.01002424 0.021663719 0.1643889 0.01793933     NA
# 
# $Pvalues
# Med_A Med_C Med_B Indo_A Atl_AD Atl_AC Indo_C Atl_AA Atl_AE Atl_AB Indo_B
# Med_A     NA    NA    NA     NA     NA     NA     NA     NA     NA     NA     NA
# Med_C      0    NA    NA     NA     NA     NA     NA     NA     NA     NA     NA
# Med_B      0     0    NA     NA     NA     NA     NA     NA     NA     NA     NA
# Indo_A     0     0     0     NA     NA     NA     NA     NA     NA     NA     NA
# Atl_AD     0     0     0      0     NA     NA     NA     NA     NA     NA     NA
# Atl_AC     0     0     0      0      0     NA     NA     NA     NA     NA     NA
# Indo_C     0     0     0      0      0      0     NA     NA     NA     NA     NA
# Atl_AA     0     0     0      0      0      0      0     NA     NA     NA     NA
# Atl_AE     0     0     0      0      0      0      0      0     NA     NA     NA
# Atl_AB     0     0     0      0      0      0      0      0      0     NA     NA
# Indo_B     0     0     0      0      0      0      0      0      0      0     NA

zczc125.pop1.boot.1000.fst<-gl.fst.pop(zczc125.pop1.gl, nboots=1000, percent=95, nclusters=1)
zczc125.pop1.boot.1000.fst


zczc125.pop1.boot.fst$Bootstraps[,c(1:2,103:106)]
# Population1 Population2 Lower bound CI limit Upper bound CI limit p-value         Fst
# 1        Med_A       Med_C          0.066879070          0.076381189       0 0.072676384
# 2        Med_A       Med_B          0.084282395          0.093644848       0 0.088840500
# 3        Med_A      Indo_A          0.236406549          0.249384399       0 0.243172741
# 4        Med_A      Atl_AD          0.189659072          0.198459009       0 0.194030427
# 5        Med_A      Atl_AC          0.194548415          0.205259462       0 0.200165372
# 6        Med_A      Indo_C          0.207281787          0.217381947       0 0.211902730
# 7        Med_A      Atl_AA          0.210353573          0.220023290       0 0.215430571
# 8        Med_A      Atl_AE          0.356991936          0.374416539       0 0.366194796
# 9        Med_A      Atl_AB          0.155382741          0.166214526       0 0.161725820
# 10       Med_A      Indo_B          0.203998980          0.213456340       0 0.208691507
# 11       Med_C       Med_B          0.003328192          0.010490915       0 0.007056218
# 12       Med_C      Indo_A          0.160652316          0.173307340       0 0.166374501
# 13       Med_C      Atl_AD          0.126230163          0.136406857       0 0.131194896
# 14       Med_C      Atl_AC          0.128828619          0.139637556       0 0.134161394
# 15       Med_C      Indo_C          0.146133854          0.157932285       0 0.151797471
# 16       Med_C      Atl_AA          0.136839146          0.148222180       0 0.142783897
# 17       Med_C      Atl_AE          0.318464156          0.336000578       0 0.327572194
# 18       Med_C      Atl_AB          0.065514167          0.075361990       0 0.070547931
# 19       Med_C      Indo_B          0.136492333          0.147091900       0 0.141220666
# 20       Med_B      Indo_A          0.263035144          0.278667395       0 0.269851202
# 21       Med_B      Atl_AD          0.201627524          0.215142828       0 0.208544224
# 22       Med_B      Atl_AC          0.208947965          0.221397762       0 0.215272885
# 23       Med_B      Indo_C          0.220980860          0.234305093       0 0.227383922
# 24       Med_B      Atl_AA          0.226773492          0.239248626       0 0.233128742
# 25       Med_B      Atl_AE          0.411439837          0.427572595       0 0.420776161
# 26       Med_B      Atl_AB          0.163472856          0.175418214       0 0.169309355
# 27       Med_B      Indo_B          0.220076372          0.231457185       0 0.225771151
# 28      Indo_A      Atl_AD          0.025863908          0.031477521       0 0.028864413
# 29      Indo_A      Atl_AC          0.029181898          0.035116576       0 0.031951483
# 30      Indo_A      Indo_C          0.015877939          0.020186396       0 0.017759121
# 31      Indo_A      Atl_AA          0.032402594          0.037764120       0 0.035350704
# 32      Indo_A      Atl_AE          0.209276415          0.227358410       0 0.217182458
# 33      Indo_A      Atl_AB          0.024721131          0.031501489       0 0.027942568
# 34      Indo_A      Indo_B          0.013483688          0.019556753       0 0.016592605
# 35      Atl_AD      Atl_AC          0.004815766          0.006300853       0 0.005441310
# 36      Atl_AD      Indo_C          0.022916472          0.025578268       0 0.024311146
# 37      Atl_AD      Atl_AA          0.011231683          0.013804828       0 0.012525243
# 38      Atl_AD      Atl_AE          0.135348531          0.149893554       0 0.142523371
# 39      Atl_AD      Atl_AB          0.007994414          0.012422990       0 0.010047364
# 40      Atl_AD      Indo_B          0.009956554          0.012718405       0 0.011147368
# 41      Atl_AC      Indo_C          0.027748853          0.031058536       0 0.029577293
# 42      Atl_AC      Atl_AA          0.010888890          0.013585425       0 0.012350897
# 43      Atl_AC      Atl_AE          0.142672001          0.157387834       0 0.149693209
# 44      Atl_AC      Atl_AB          0.010028638          0.014507477       0 0.012355496
# 45      Atl_AC      Indo_B          0.015039545          0.017798357       0 0.016423868
# 46      Indo_C      Atl_AA          0.032791781          0.036835399       0 0.034876468
# 47      Indo_C      Atl_AE          0.166034836          0.179579899       0 0.171880123
# 48      Indo_C      Atl_AB          0.030664878          0.036406795       0 0.033259359
# 49      Indo_C      Indo_B          0.008612729          0.011240338       0 0.010024243
# 50      Atl_AA      Atl_AE          0.148220043          0.164603307       0 0.156176759
# 51      Atl_AA      Atl_AB          0.002828972          0.008443667       0 0.006053763
# 52      Atl_AA      Indo_B          0.019463654          0.023749567       0 0.021663719
# 53      Atl_AE      Atl_AB          0.159738917          0.178365950       0 0.170494461
# 54      Atl_AE      Indo_B          0.157284300          0.172211389       0 0.164388856
# 55      Atl_AB      Indo_B          0.015430528          0.020520458       0 0.017939332

zczc125.pop2.boot.fst<-gl.fst.pop(zczc125.pop2.gl, nboots=100, percent=95, nclusters=1)
zczc125.pop2.boot.fst$Bootstraps[,c(1:2,103:106)]
# Population1 Population2 Lower bound CI limit Upper bound CI limit p-value         Fst
# 1        Med_A       Med_C          0.069283594          0.076821913       0 0.072676384
# 2        Med_A       Med_B          0.084334345          0.093215525       0 0.088840500
# 3        Med_A      Indo_A          0.236760580          0.249240846       0 0.243172741
# 4        Med_A      Atl_BC          0.188368753          0.199225560       0 0.194030427
# 5        Med_A      Atl_BB          0.195242872          0.205163356       0 0.200165372
# 6        Med_A      Indo_C          0.206243653          0.217112966       0 0.211902730
# 7        Med_A      Atl_BA          0.179547931          0.190357880       0 0.185598626
# 8        Med_A      Atl_BD          0.355494096          0.372304879       0 0.366194796
# 9        Med_A      Indo_B          0.202364757          0.214695417       0 0.208691507
# 10       Med_C       Med_B          0.003386509          0.010745700       0 0.007056218
# 11       Med_C      Indo_A          0.159898371          0.171990166       0 0.166374501
# 12       Med_C      Atl_BC          0.126182588          0.135328768       0 0.131194896
# 13       Med_C      Atl_BB          0.129014437          0.138751395       0 0.134161394
# 14       Med_C      Indo_C          0.146080992          0.156025926       0 0.151797471
# 15       Med_C      Atl_BA          0.109619077          0.119048231       0 0.114863834
# 16       Med_C      Atl_BD          0.318494396          0.334280749       0 0.327572194
# 17       Med_C      Indo_B          0.136208892          0.146273756       0 0.141220666
# 18       Med_B      Indo_A          0.261627553          0.276465621       0 0.269851202
# 19       Med_B      Atl_BC          0.201338482          0.213191653       0 0.208544224
# 20       Med_B      Atl_BB          0.208641476          0.219959220       0 0.215272885
# 21       Med_B      Indo_C          0.221185057          0.232606080       0 0.227383922
# 22       Med_B      Atl_BA          0.190069352          0.201215645       0 0.196417339
# 23       Med_B      Atl_BD          0.409959524          0.427393376       0 0.420776161
# 24       Med_B      Indo_B          0.220259803          0.230772961       0 0.225771151
# 25      Indo_A      Atl_BC          0.025514935          0.031125987       0 0.028864413
# 26      Indo_A      Atl_BB          0.028679768          0.034154347       0 0.031951483
# 27      Indo_A      Indo_C          0.014819332          0.019829007       0 0.017759121
# 28      Indo_A      Atl_BA          0.027216642          0.032232365       0 0.030173927
# 29      Indo_A      Atl_BD          0.208681565          0.223380475       0 0.217182458
# 30      Indo_A      Indo_B          0.013464778          0.019086632       0 0.016592605
# 31      Atl_BC      Atl_BB          0.004592924          0.006533311       0 0.005441310
# 32      Atl_BC      Indo_C          0.022857261          0.025804659       0 0.024311146
# 33      Atl_BC      Atl_BA          0.008554202          0.010670200       0 0.009738826
# 34      Atl_BC      Atl_BD          0.135559084          0.147991602       0 0.142523371
# 35      Atl_BC      Indo_B          0.009780495          0.012276984       0 0.011147368
# 36      Atl_BB      Indo_C          0.027908689          0.031107786       0 0.029577293
# 37      Atl_BB      Atl_BA          0.009155002          0.011371042       0 0.010328401
# 38      Atl_BB      Atl_BD          0.142334777          0.155832196       0 0.149693209
# 39      Atl_BB      Indo_B          0.014900976          0.017883722       0 0.016423868
# 40      Indo_C      Atl_BA          0.029568271          0.033158613       0 0.031506636
# 41      Indo_C      Atl_BD          0.162845325          0.177920425       0 0.171880123
# 42      Indo_C      Indo_B          0.008747731          0.011353543       0 0.010024243
# 43      Atl_BA      Atl_BD          0.135822685          0.150862432       0 0.144353124
# 44      Atl_BA      Indo_B          0.016389842          0.019640059       0 0.018116823
# 45      Atl_BD      Indo_B          0.156785069          0.170840364       0 0.164388856

zczc125.pop2.boot.1000.fst<-gl.fst.pop(zczc125.pop2.gl, nboots=1000, percent=95, nclusters=1)
zczc125.pop2.boot.1000.fst$Bootstraps[,c(1:2,1003:1006)]
# Summary statistics - Poppr ----------------------------------------------
####pairwise population differentiation with poppr

library("apex")
library("adegenet")
library("pegas")
library("mmod")
library("poppr")

poppr(mden43.gi) #Produce a basic summary table for population genetic analyses.
# | Atlantic 
# | Indopacific 
# | Total 
#       Pop         N   MLG  eMLG   SE    H     G   lambda  E.5  Hexp    Ia     rbarD     File
# 1     Atlantic    28  28   15     0     3.33  28  0.964   1    0.114   46.7   0.00681   mden43.gi
# 2     Indopacific 15  15   15     0     2.71  15  0.933   1    0.140   10.9   0.00112   mden43.gi
# 3     Total       43  43   15     0     3.76  43  0.977   1    0.131   123.3  0.01077   mden43.gi

poppr(mden43.pop.gi)
# | Atl-Bahamas 
# | Atl-East 
# | Atl-Other 
# | Indo-Africa 
# | Indo-Hawaii 
# | Indo-South 
# | Total 
# Pop  N MLG eMLG       SE    H  G lambda E.5  Hexp     Ia    rbarD          File
# 1 Atl-Bahamas  7   7    7 0.00e+00 1.95  7  0.857   1 0.109 181.63 0.041052 mden43.pop.gi
# 2    Atl-East 16  16   10 0.00e+00 2.77 16  0.938   1 0.113  23.75 0.004091 mden43.pop.gi
# 3   Atl-Other  5   5    5 0.00e+00 1.61  5  0.800   1 0.113 119.31 0.029230 mden43.pop.gi
# 4 Indo-Africa  5   5    5 0.00e+00 1.61  5  0.800   1 0.140  26.15 0.004529 mden43.pop.gi
# 5 Indo-Hawaii  6   6    6 0.00e+00 1.79  6  0.833   1 0.136  14.32 0.002353 mden43.pop.gi
# 6  Indo-South  3   3    3 0.00e+00 1.10  3  0.667   1 0.136   1.42 0.000356 mden43.pop.gi
# 7          NA  1   1    1 0.00e+00 0.00  1  0.000 NaN 0.131     NA       NA mden43.pop.gi
# 8       Total 43  43   10 9.46e-07 3.76 43  0.977   1 0.131 123.26 0.010773 mden43.pop.gi

poppr(zczc129.gi)
# | Atlantic 
# | Indopacific 
# | Mediterranean 
# | Total 
# Pop   N MLG eMLG       SE    H   G lambda E.5   Hexp    Ia   rbarD       File
# 1      Atlantic  59  59   34 9.53e-07 4.08  59  0.983   1 0.1369  61.0 0.00312 zczc129.gi
# 2   Indopacific  36  36   34 3.94e-08 3.58  36  0.972   1 0.1326  41.5 0.00240 zczc129.gi
# 3 Mediterranean  34  34   34 0.00e+00 3.53  34  0.971   1 0.0992 198.5 0.02251 zczc129.gi
# 4         Total 129 129   34 1.61e-05 4.86 129  0.992   1 0.1364 133.9 0.00634 zczc129.gi

poppr(zczc129.pop.gi)
# | Med-Italy 
# | Atl-East 
# | Atl-Madeira 
# | Atl-Other 
# | Med-East 
# | Indo-South 
# | Indo-East 
# | Indo-Spac 
# | Atl-Bahamas 
# | Indo-Central 
# | Total 
# Pop   N MLG eMLG       SE     H   G lambda E.5   Hexp      Ia     rbarD           File
# 1     Med-Italy  20  20   10 0.00e+00 2.996  20  0.950   1 0.0975   47.26  6.30e-03 zczc129.pop.gi
# 2      Atl-East  35  35   10 8.89e-07 3.555  35  0.971   1 0.1356   48.97  2.71e-03 zczc129.pop.gi
# 3    Med-France   2   2    2 0.00e+00 0.693   2  0.500   1 0.0950      NA        NA zczc129.pop.gi
# 4   Atl-Madeira   4   4    4 0.00e+00 1.386   4  0.750   1 0.1148 4082.06  8.13e-01 zczc129.pop.gi
# 5     Atl-Other   7   7    7 0.00e+00 1.946   7  0.857   1 0.1358   70.53  6.35e-03 zczc129.pop.gi
# 6      Med-East  12  12   10 0.00e+00 2.485  12  0.917   1 0.0938  571.66  7.50e-02 zczc129.pop.gi
# 7     Atl-Spain   2   2    2 0.00e+00 0.693   2  0.500   1 0.0864      NA        NA zczc129.pop.gi
# 8    Indo-South   9   9    9 0.00e+00 2.197   9  0.889   1 0.1323    3.28  2.71e-04 zczc129.pop.gi
# 9     Indo-East  17  17   10 2.31e-07 2.833  17  0.941   1 0.1316   41.36  2.92e-03 zczc129.pop.gi
# 10    Indo-Spac   3   3    3 0.00e+00 1.099   3  0.667   1 0.1274   -0.35 -5.39e-05 zczc129.pop.gi
# 11  Atl-Bahamas  11  11   10 0.00e+00 2.398  11  0.909   1 0.1334  292.41  2.42e-02 zczc129.pop.gi
# 12 Indo-Central   5   5    5 0.00e+00 1.609   5  0.800   1 0.1257   40.60  4.80e-03 zczc129.pop.gi
# 13  Indo-Mexico   2   2    2 0.00e+00 0.693   2  0.500   1 0.1244      NA        NA zczc129.pop.gi
# 14        Total 129 129   10 0.00e+00 4.860 129  0.992   1 0.1364  133.86  6.34e-03 zczc129.pop.gi

poppr(zczc127.gi)
#| Atlantic 
#| Indopacific 
#| Mediterranean 
#| Total 
#            Pop   N   MLG   eMLG       SE    H   G  lambda E.5  Hexp    Ia   rbarD       File
#1      Atlantic  57   57    34    1.12e-06 4.04  57  0.982   1  0.141  53.8 0.00280 zczc127.gi
#2   Indopacific  36   36    34    3.94e-08 3.58  36  0.972   1  0.136  41.5 0.00243 zczc127.gi
#3 Mediterranean  34   34    34    0.00e+00 3.53  34  0.971   1  0.102 198.6 0.02266 zczc127.gi
#4         Total 127   127   34    0.00e+00 4.84 127  0.992   1  0.140 134.7 0.00652 zczc127.gi

poppr(zczc127.pop.gi)
# Pop   N MLG eMLG       SE     H   G lambda E.5   Hexp      Ia    rbarD           File
# 1   Atl-Bahamas  11  11   10 0.00e+00 2.398  11  0.909   1 0.1376 292.134 2.42e-02 zczc127.pop.gi
# 2      Atl-East  35  35   10 8.89e-07 3.555  35  0.971   1 0.1396  49.244 2.77e-03 zczc127.pop.gi
# 3   Atl-Madeira   2   2    2 0.00e+00 0.693   2  0.500   1 0.1330      NA       NA zczc127.pop.gi
# 4     Atl-Other   7   7    7 0.00e+00 1.946   7  0.857   1 0.1399  71.843 6.50e-03 zczc127.pop.gi
# 5     Atl-Spain   2   2    2 0.00e+00 0.693   2  0.500   1 0.0890      NA       NA zczc127.pop.gi
# 6  Indo-Central   5   5    5 0.00e+00 1.609   5  0.800   1 0.1293  40.132 4.77e-03 zczc127.pop.gi
# 7     Indo-East  17  17   10 2.31e-07 2.833  17  0.941   1 0.1354  41.513 2.96e-03 zczc127.pop.gi
# 8   Indo-Mexico   2   2    2 0.00e+00 0.693   2  0.500   1 0.1280      NA       NA zczc127.pop.gi
# 9    Indo-South   7   7    7 0.00e+00 1.946   7  0.857   1 0.1360   3.429 3.18e-04 zczc127.pop.gi
# 10    Indo-Spac   5   5    5 0.00e+00 1.609   5  0.800   1 0.1342   0.861 9.46e-05 zczc127.pop.gi
# 11     Med-East  12  12   10 0.00e+00 2.485  12  0.917   1 0.0967 570.267 7.51e-02 zczc127.pop.gi
# 12   Med-France   2   2    2 0.00e+00 0.693   2  0.500   1 0.0981      NA       NA zczc127.pop.gi
# 13    Med-Italy  20  20   10 0.00e+00 2.996  20  0.950   1 0.1006  48.132 6.44e-03 zczc127.pop.gi
# 14        Total 127 127   10 0.00e+00 4.844 127  0.992   1 0.1403 134.739 6.52e-03 zczc127.pop.gi
# Exploring vcftools data -------------------------------------------------
######playing around with vcftools data
sum<-read.table("~/Dropbox/Phd/Bioinformatics/bw_ddrad/Results/R-summaries/r-working-summaries.csv", sep=",", header=TRUE)
head(sum)
sum<-sum[,1:15]
head(sum)
sum.zc<-subset(sum, sum$species=="Zcav")
sum.md<-subset(sum, sum$species=="Mden")
sum.zc
sum.md

hist(sum.zc$no.snps.ind)
hist(sum.md$no.snps.ind)

zc.no.snps<-(sum.zc$no.snps.ind)
md.no.snps<-(sum.md$no.snps.ind)
zc.no.snps

plot(sum.zc$no.snps.ind, sum.zc$f)
boxplot(sum.zc$f ~ sum.zc$population)
boxplot(sum.zc$f ~ sum.zc$ocean.basin)
dev.off()

plot(sum.md$no.snps.ind, sum.md$f)
boxplot(sum.md$f ~ sum.md$population)
boxplot(sum.md$f ~ sum.md$ocean.basin)

#for some reason the M. densirostris plots were still including mediteranean factors. have converted those columns to characters and they work now!
sum.md$ocean.basin
sum.md$population

sum.md$ocean.basin<-as.character(sum.md$ocean.basin)
sum.md$ocean.basin
sum.md$population<-as.character(sum.md$population)
sum.md$population

sum.zc$ocean.basin<-as.character(sum.zc$ocean.basin)
sum.zc$ocean.basin
sum.zc$population<-as.character(sum.zc$population)
sum.zc$population

# Inbreeding Coefficients -------------------------------------------------
###plotting inbreeding coefficient by ocean basin and population using ggplot
png("mden.fstats.pop.png",width=6, height=4, units='in', res=300 )
sum.md %>%
mutate(population=fct_relevel(population, c("Atl-Bah", "Atl-East", "Atl-Oth","Indo-Haw", "Indo-Afr", "Indo-Sou"))) %>%
ggplot(aes(x=population, y=f, fill=ocean.basin)) +
  scale_fill_manual(values = my_zissou)+
  geom_boxplot() + 
  theme_classic() +
  labs(fill="Ocean Basin") + 
  labs(x=NULL, y="Inbreeding Coefficient (F)") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()

png("mden.fstats.ocean.png",width=6, height=4, units='in', res=300 )
ggplot(sum.md, aes(x=sum.md$ocean.basin, y=sum.md$f, fill=sum.md$ocean.basin)) +
  scale_fill_manual(values = my_zissou)+
  geom_boxplot() + 
  theme_classic() +
  labs(fill="Ocean Basin") + 
  labs(x=NULL, y="Inbreeding Coefficient (F)") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()

png("zcav.fstats.pop.png",width=6, height=4, units='in', res=300 )
ggplot(sum.zc, aes(x=sum.zc$population, y=sum.zc$f, fill=sum.zc$ocean.basin)) +
  scale_fill_manual(values = my_zissou)+
  geom_boxplot() + 
  theme_classic() +
  labs(fill="Ocean Basin") + 
  labs(x=NULL, y="Inbreeding Coefficient (F)") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()
 
png("zcav.fstats.ocean.png",width=6, height=4, units='in', res=300 )
ggplot(sum.zc, aes(x=sum.zc$ocean.basin, y=sum.zc$f, fill=sum.zc$ocean.basin)) +
  scale_fill_manual(values = my_zissou)+
  geom_boxplot() + 
  theme_classic() +
  labs(fill="Ocean Basin") + 
  labs(x=NULL, y="Inbreeding Coefficient (F)") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()

sum.md
sum.zc

tapply(sum.zc.no.admix$f, sum.zc.no.admix$population, mean)
tapply(sum.zc$f, sum.zc$ocean.basin, mean)

boxplot(sum.zc.no.admix$f ~sum.zc.no.admix$population)
boxplot(sum.md.no.admix$f ~sum.md.no.admix$population)

boxplot(sum.zc.no.admix$f ~sum.zc.no.admix$ocean.basin)
boxplot(sum.md.no.admix$f ~sum.md.no.admix$ocean.basin)

tapply(sum.md.no.admix$f, sum.md.no.admix$population, mean)

tapply(sum.md.no.admix$f, sum.md.no.admix$ocean.basin, mean)
tapply(sum.zc.no.admix$f, sum.zc.no.admix$ocean.basin, mean)

###New inbreeding plots with admixed populations removed
sum.zc.no.admix
sum.md.no.admix

###plotting inbreeding coefficient by ocean basin and population using ggplot
png("mden.fstats.noadmix.pop.png",width=6, height=4, units='in', res=300 )
sum.md.no.admix %>%
  mutate(population=fct_relevel(population, c("Atl-Bah", "Atl-East", "Indo-Haw", "Indo-Afr", "Indo-Sou"))) %>%
  ggplot(aes(x=population, y=f, fill=ocean.basin)) +
  scale_fill_manual(values = my_zissou)+
  geom_boxplot() + 
  theme_classic() +
  labs(fill="Ocean Basin") + 
  labs(x=NULL, y="Inbreeding Coefficient (F)") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()

png("zcav.fstats.noadmix.pop.png",width=6, height=4, units='in', res=300 )
ggplot(sum.zc.no.admix, aes(x=sum.zc.no.admix$population, y=sum.zc.no.admix$f, fill=sum.zc.no.admix$ocean.basin)) +
  scale_fill_manual(values = my_zissou)+
  geom_boxplot() + 
  theme_classic() +
  labs(fill="Ocean Basin") + 
  labs(x=NULL, y="Inbreeding Coefficient (F)") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()

# Site Depth --------------------------------------------------------------
##makin histograms of site depth

#first made a dataframe with the mean site depth per ocean basin
md.dat<-ddply(sum.md, "ocean.basin", summarise, rating.mean=mean(mean.site.depth))
md.dat

#then plot histogram
png("mden.site-depth.ocean.png",width=6, height=4, units='in', res=300 )
ggplot(sum.md, aes(x=sum.md$mean.site.depth, fill=sum.md$ocean)) + 
  scale_fill_manual(values = wes_palette(name="Darjeeling1"))+
  geom_histogram(binwidth=15, alpha=0.5, position="identity") +
  geom_vline(data=md.dat, aes(xintercept=rating.mean, colour=ocean.basin, legend=NULL), linetype="dashed", size=1, show.legend=FALSE) +
  scale_color_manual(values=wes_palette(name="Darjeeling1"))+
  theme_classic() +
  labs(x="Mean site depth per locus, per individual", y="Frequency") +
  labs(fill="Ocean Basin") 
dev.off()

#making dataframe of averages by ocean basin
zc.dat<-ddply(sum.zc, "ocean.basin", summarise, rating.mean=mean(mean.site.depth))
zc.dat

#plotting zc histogram
png("zcav.site-depth.ocean.png",width=6, height=4, units='in', res=300 )
ggplot(sum.zc, aes(x=sum.zc$mean.site.depth, fill=sum.zc$ocean)) + 
  scale_fill_manual(values = wes_palette(name="Darjeeling1"))+
  geom_histogram(binwidth=15, alpha=0.5, position="identity") +
  geom_vline(data=md.dat, aes(xintercept=rating.mean, colour=ocean.basin, legend=NULL), linetype="dashed", size=1, show.legend=FALSE) +
  scale_color_manual(values=wes_palette(name="Darjeeling1"))+
  theme_classic() +
  labs(x="Mean site depth per locus, per individual", y="Frequency") +
  labs(fill="Ocean Basin") 
dev.off()

summary(sum.md$mean.site.depth)
mean(sum.md$mean.site.depth)
sd(sum.md$mean.site.depth)
describe(sum.md$mean.site.depth)

summary(sum.zc$mean.site.depth)
mean(sum.zc$mean.site.depth)
sd(sum.zc$mean.site.depth)
describe(sum.zc$mean.site.depth)



# Missingness -------------------------------------------------------------
##plotting missingness per indivdiaul by ocean basin
#making dataframe of averages by ocean basin
zc.miss<-ddply(sum.zc, "ocean.basin", summarise, rating.mean=mean(prop.missing.snps))
zc.miss

#plotting zc histogram
png("zcav.missing.ocean.png",width=6, height=4, units='in', res=300 )
ggplot(sum.zc, aes(x=sum.zc$prop.missing.snps, fill=sum.zc$ocean)) + 
  scale_fill_manual(values = wes_palette(name="Darjeeling1"))+
  geom_histogram(binwidth=.05,alpha=0.5, position="identity") +
  geom_vline(data=zc.miss, aes(xintercept=rating.mean, colour=ocean.basin), linetype="dashed", size=1, show.legend=FALSE) +
  scale_color_manual(values=wes_palette(name="Darjeeling1"))+
  theme_classic() +
  labs(x="Mean missingness per individual", y="Frequency") +
  labs(fill="Ocean Basin")
dev.off()

#making dataframe of averages by ocean basin
md.miss<-ddply(sum.md, "ocean.basin", summarise, rating.mean=mean(prop.missing.snps))
md.miss

#plotting md histogram
png("mden.missing.ocean.png",width=6, height=4, units='in', res=300 )
ggplot(sum.md, aes(x=sum.md$prop.missing.snps, fill=sum.md$ocean)) + 
  scale_fill_manual(values = wes_palette(name="Darjeeling1"))+
  geom_histogram(binwidth=.05,alpha=0.5, position="identity") +
  geom_vline(data=md.miss, aes(xintercept=rating.mean, colour=ocean.basin), linetype="dashed", size=1, show.legend=FALSE) +
  scale_color_manual(values=wes_palette(name="Darjeeling1"))+
  theme_classic() +
  labs(x="Mean missingness per individual", y="Frequency") +
  labs(fill="Ocean Basin")
dev.off()


# No. Snps ----------------------------------------------------------------
###boxplots of number of snps


library(wesanderson)
#title="Number of SNP loci per individual  by population, Z. cavirostris"
png("zcav.snp-no.pop.png",width=6, height=4, units='in', res=300 )
ggplot(sum.zc, aes(x=sum.zc$population, y=sum.zc$no.snps.ind, fill=sum.zc$ocean.basin)) +
  scale_fill_manual(values = my_zissou)+
  geom_boxplot() +
  theme_classic() +
  labs(fill="Ocean Basin") + 
  labs(x=NULL, y="Number of SNP loci") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()

png("zcav.snp-no.ocean.png",width=6, height=4, units='in', res=300 )
ggplot(sum.zc, aes(x=sum.zc$ocean.basin, y=sum.zc$no.snps.ind, fill=sum.zc$ocean.basin)) +
  scale_fill_manual(values = my_zissou)+
  geom_boxplot() + 
  theme(legend.position = "none") + 
  labs(x=NULL, y="Number of SNP loci")+
  labs(fill="Ocean Basin") +
  theme_classic()
dev.off()

#need to re-order here, hawaii needs to come before africa to match rest of plots
png("mden.snp-no.pop.png",width=6, height=4, units='in', res=300 )
sum.md %>%
mutate(population=fct_relevel(population, c("Atl-Bah", "Atl-East", "Atl-Oth","Indo-Haw", "Indo-Afr", "Indo-Sou"))) %>%
ggplot(aes(x=population, y=no.snps.ind, fill=ocean.basin)) +
  scale_fill_manual(values = my_zissou)+
  geom_boxplot() +
  theme_classic() +
  labs(fill="Ocean Basin") + 
  labs(x=NULL, y="Number of SNP loci") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()

png("mden.snp-no.ocean.png",width=6, height=4, units='in', res=300 )
ggplot(sum.md, aes(x=sum.md$ocean.basin, y=sum.md$no.snps.ind, fill=sum.md$ocean.basin)) +
  scale_fill_manual(values = my_zissou)+
  geom_boxplot() +
  theme_classic() +
  labs(fill="Ocean Basin") + 
  labs(x=NULL, y="Number of SNP loci") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()

summary(sum.md$no.snps.ind)
mean(sum.md$no.snps.ind)
sd(sum.md$no.snps.ind)
describe(sum.md$no.snps.ind)

summary(sum.zc$no.snps.ind)
mean(sum.zc$no.snps.ind)
sd(sum.zc$no.snps.ind)
mean(sum.zc$no.snps.ind)
max(sum.zc$no.snps.total)
describe(sum.zc$no.snps.ind)


#calculates the mean number of snps per individual by population or ocean basin
tapply(sum.zc$no.snps.ind, sum.zc$population, mean) 
tapply(sum.zc$no.snps.ind, sum.zc$ocean, mean)

tapply(sum.md$no.snps.ind, sum.md$population, range)
tapply(sum.md$no.snps.ind, sum.md$population, mean)

##Remove populations that were admixed
#Mexico
#Atl-Other
#Atl-Other

table(sum.zc$pop)
#remove Atl-Oth and Indo-Mex
table(sum.md$pop)
#remove Atl-Oth

sum.zc.no.admix<-subset(sum.zc, sum.zc$population!="Atl-Oth")
sum.zc.no.admix<-subset(sum.zc.no.admix, sum.zc.no.admix$population!="Indo-Mex")
table(sum.zc.no.admix$population)

sum.md.no.admix<-subset(sum.md, sum.md$population!="Atl-Oth")
table(sum.md.no.admix$population)

##New SNP counts with admixed populations removed
tapply(sum.zc.no.admix$no.snps.ind, sum.zc.no.admix$population, range)
tapply(sum.zc.no.admix$no.snps.ind, sum.zc.no.admix$population, mean)

tapply(sum.md.no.admix$no.snps.ind, sum.md.no.admix$population, range)
tapply(sum.md.no.admix$no.snps.ind, sum.md.no.admix$population, mean)

#new plots without admixed individuals
#need to re-order here, hawaii needs to come before africa to match rest of plots
#title="Number of SNP loci per individual  by population, Z. cavirostris"
png("zcav.snp-no.noadmix.pop.png",width=6, height=4, units='in', res=300 )
ggplot(sum.zc.no.admix, aes(x=sum.zc.no.admix$population, y=sum.zc.no.admix$no.snps.ind, fill=sum.zc.no.admix$ocean.basin)) +
  scale_fill_manual(values = my_zissou)+
  geom_boxplot() +
  theme_classic() +
  labs(fill="Ocean Basin") + 
  labs(x=NULL, y="Number of SNP loci") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()

png("mden.snp-no.noadmix.pop.png",width=6, height=4, units='in', res=300 )
sum.md.no.admix %>%
  mutate(population=fct_relevel(population, c("Atl-Bah", "Atl-East", "Indo-Haw", "Indo-Afr", "Indo-Sou"))) %>%
  ggplot(aes(x=population, y=no.snps.ind, fill=ocean.basin)) +
  scale_fill_manual(values = my_zissou)+
  geom_boxplot() +
  theme_classic() +
  labs(fill="Ocean Basin") + 
  labs(x=NULL, y="Number of SNP loci") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()


# Maps of sampling locations ----------------------------------------------
####making maps of data

map<-read.csv("~/Dropbox/Phd/Bioinformatics/bw_ddrad_global/Results/R-summaries/map-files.csv", header=TRUE, sep=",")
map
map.zc<-subset(map, map$species=="Zcav")
map.md<-subset(map, map$species=="Mden")
head(map.zc)
head(map.md)

is.factor(map$species)
map$species<-as.factor(map$species)

is.numeric(map$lat)
is.numeric(map$long)


#basic world map
wm<-map_data("world")
wm<-ggplot(wm, aes(long, lat, group=group)) +
  geom_polygon(fill="grey", colour="white") +
  coord_fixed(1.3) +
  theme_classic()
wm

png("ddrad.sample.map.png", width=6.75, height=4, units="in", res=300)
all.wm<-wm +
  geom_jitter(data=map, position=position_jitter(width = 1, height = 1),
              aes(x=long, y=lat, group=species, colour=species, shape=species), alpha=.6, size=2.5) +
  labs(colour=NULL, shape=NULL, x="Longitude", y="Latitude") +
  theme(legend.position=c(.93,0.6)) +
  scale_color_manual(values=wes_palette("Zissou1")[c(1,5,1,5)], breaks=c("Mden", "Zcav", "Blainville's", "Cuvier's"),labels=c("Blainville's-ddRAD", "Cuvier's-ddRAD", "Blainville's-mtDNA", "Cuvier's-mtDNA")) +
  scale_shape_manual(values=c(16,16, 8,8), breaks=c("Mden", "Zcav", "Blainville's", "Cuvier's"),labels=c("Blainville's-ddRAD", "Cuvier's-ddRAD", "Blainville's-mtDNA", "Cuvier's-mtDNA"))
all.wm
dev.off()

##Splitting maps by species
zc.map <- rbind(map[map$species=="Cuvier's",], map[map$species=="Zcav",])
png("./zc.map.png", width=6.75, height=4, units="in", res=500)
zc.wm<-wm +
  geom_jitter(data=zc.map, position=position_jitter(width = 1, height = 2),
              aes(x=long, y=lat, group=species, colour=species, shape=species), alpha=.9, size=3) +
  labs(colour=NULL, shape=NULL, x="Longitude", y="Latitude") +
  theme(legend.position=c(.14,0.2)) +
  scale_color_manual(values=paired.col, breaks=c("Zcav", "Cuvier's"),labels=c("Cuvier's-ddRAD", "Cuvier's-mtDNA")) +
  scale_shape_manual(values=c(16,17), breaks=c("Zcav","Cuvier's"),labels=c("Cuvier's-ddRAD", "Cuvier's-mtDNA"))
zc.wm
dev.off()

md.map <- rbind(map[map$species=="Blainville's",], map[map$species=="Mden",])
png("md.map.png", width=6.75, height=4, units="in", res=500)
md.wm<-wm +
  geom_jitter(data=md.map, position=position_jitter(width = 1, height = 2),
              aes(x=long, y=lat, group=species, colour=species, shape=species), alpha=.9, size=3) +
  labs(colour=NULL, shape=NULL, x="Longitude", y="Latitude") +
  theme(legend.position=c(.16,0.2)) +
  scale_color_manual(values=paired.col[c(3,4)], breaks=c("Mden", "Blainville's"),labels=c("Blainville's-ddRAD", "Blainville's-mtDNA")) +
  scale_shape_manual(values=c(16,17), breaks=c("Mden","Blainville's"),labels=c("Blainville's-ddRAD", "Blainville's-mtDNA"))
md.wm
dev.off()


##Splitting maps by marker type-using wes anderson palet
# rad.map<-subset(map, map$species==c("Mden", "Zcav"))
# png("./ddrad.map.png", width=6.75, height=4, units="in", res=500)
# rad.wm<-wm +
#   geom_jitter(data=rad.map, position=position_jitter(width = 1, height = 1),
#               aes(x=long, y=lat, group=species, colour=species, shape=species), alpha=.6, size=3) +
#   labs(colour=NULL, shape=NULL, x="Longitude", y="Latitude") +
#   theme(legend.position=c(.93,0.6)) +
#   scale_color_manual(values=wes_palette("Zissou1")[c(5,1)], breaks=c("Zcav", "Mden"),labels=c("Cuvier's-ddRAD", "Blainville's-ddRAD")) +
#   scale_shape_manual(values=c(16,17), breaks=c("Zcav","Mden"),labels=c("Cuvier's-ddRAD", "Blainville's-ddRAD")) 
# rad.wm
# dev.off()

#splitting maps by marker type using paired Color Brewer
display.brewer.pal(n=12, name='Paired')
png("./ddrad.map.png", width=6.75, height=4, units="in", res=500)
rad.wm<-wm +
  geom_jitter(data=rad.map, position=position_jitter(width = 1, height = 1),
              aes(x=long, y=lat, group=species, colour=species, shape=species), alpha=.9, size=3) +
  labs(colour=NULL, shape=NULL, x="Longitude", y="Latitude") +
  theme(legend.position=c(.14,0.2)) +
 scale_color_manual(values=paired.col, breaks=c("Zcav","Mden"),labels=c("Cuvier's-ddRAD", "Blainville's-ddRAD")) +
  scale_shape_manual(values=c(16,17), breaks=c("Zcav","Mden"),labels=c("Cuvier's-ddRAD", "Blainville's-ddRAD")) 
rad.wm
dev.off()

mito.map<-subset(map, map$species==c("Blainville's", "Cuvier's"))
png("./mito.map.png", width=6.75, height=4, units="in", res=500)
mito.wm<-wm +
  geom_jitter(data=mito.map, position=position_jitter(width = 1, height = 1),
              aes(x=long, y=lat, group=species, colour=species, shape=species), alpha=.9, size=3) +
  labs(colour=NULL, shape=NULL, x="Longitude", y="Latitude") +
  scale_color_manual(values=paired.col[c(3,4)], breaks=c("Cuvier's", "Blainville's"),labels=c("Cuvier's-mitogenomes", "Blainville's-mitogenomes")) +
  scale_shape_manual(values=c(16,17), breaks=c("Cuvier's","Blainville's"),labels=c("Cuvier's-mitogenomes", "Blainville's-mitogenomes")) +
  theme(legend.position=c(.16,.2), legend.background = element_rect(fill = NA) )
mito.wm
dev.off()

##Map of ALL ITABW samples
itabw.map<-read.csv("~/Dropbox/Phd/Bioinformatics/bw_ddrad/maps/ITABW_map.csv", header=TRUE, sep=",")
head(itabw.map)

png("all.ITABW.samples.map.png",width = 6, height = 4, units = 'in', res = 300 )
itabw.wm<-wm +
  geom_jitter(data=itabw.map, position=position_jitter(width = .5, height = .5),
              aes(x=as.numeric(long), y=as.numeric(lat), group=species, colour=species), alpha=.5, size=2) +
scale_colour_discrete(breaks=c("H. ampullatus", "Mbid", "Mden", "Meur", "Mmir", "Unknown", "Zcav"), labels=c("H. ampullatus", "M. bidens", "M. densirostris", "M. europaeus", "M. mirus", "Unknown Ziphiid", "Z. cavirostris")) +
  scale_color_manual(values=wes_palette(name="Darjeeling1", 7, type="continuous"))+
  labs(colour=NULL, x="Longitude", y="Latitude") +
  theme(legend.position=c(1.05,0.5)) +
  theme(legend.background = element_blank()) +
  theme(legend.text = element_text(size = 8)) +
  theme(plot.margin = margin(0.1, 2.1, 0.1, 0.1, "cm"))
itabw.wm
dev.off()

# #ddRAD Mden only map - samples in Red
# map.md<-subset(map, map$species=="Mden")
# head(map.md)
# 
# png("mden.allddRAD.map.png", width = 6, height = 4, units = 'in', res = 300)
# mdenwm<-wm + 
#   ##this is the information needed to add the data to the plot. 
#   geom_jitter(data=map.md, position=position_jitter(width=1, height=1),
#               aes(x=long, y=lat, group=species), colour="red", size=1) +
#   labs(colour=NULL, x="Longitude", y="Latitude") +
#   theme(legend.title = element_blank())
# ##this will show the new map
# mdenwm
# dev.off()

# #ddRAD Zcav only map-samples in Blue
# map.zc<-subset(map, map$species=="Zcav")
# head(map.zc)
# 
# 
# png("zcav.allddRAD.map.png", width = 6, height = 4, units = 'in', res = 300)
# zcavwm<-wm + 
#   ##this is the information needed to add the data to the plot. 
#   geom_jitter(data=map.zc, position=position_jitter(width=1, height=1),
#               aes(x=long, y=lat, group=species), colour="blue", size=1) +
#   labs(colour=NULL, x="Longitude", y="Latitude") +
#   theme(legend.title = element_blank())
# ##this will show the new map
# zcavwm
# dev.off()

##Blank world map
wm<-map_data("world")
wm<-ggplot(wm, aes(long, lat, group=group)) +
  geom_polygon(fill="grey", colour="white") +
  coord_fixed(1.3) +
  theme_classic()
wm


# mden43.map<-read.csv("~/Dropbox/Phd/Bioinformatics/bw_ddrad/SNP_files/no_populations/one_snp_per_locus/mden_results/mden43.coord.csv", sep=",", header=TRUE)
# head(mden43.map)
# 
# png("mden.allddRAD.map.png", width = 6, height = 4, units = 'in', res = 300)
# mdenwm<-wm + 
#   ##this is the information needed to add the data to the plot. 
#   geom_jitter(data=mden43.map, position=position_jitter(width=1, height=1),
#               aes(x=long, y=lat, group=ocean, col=ocean), size=1.5) +
#   labs(colour=NULL, x="Longitude", y="Latitude") +
#   theme(legend.title = element_blank())
# ##this will show the new map
# mdenwm
# dev.off()

# mdenwm<-wm + 
#   ##this is the information needed to add the data to the plot. 
#   geom_jitter(data=mden43.map, position=position_jitter(width=1, height=1),
#               aes(x=long, y=lat, group=pop1, col=pop1), size=1.5) +
#   labs(colour=NULL, x="Longitude", y="Latitude") +
#   theme(legend.title = element_blank())
# ##this will show the new map
# mdenwm


# Tajima’s D - Z. cavirostris -------------------------------------------------------------
######need to use vcf file format to do tajima's test. 
######Need to use original vcf file and remove the indivdiuals that didnt pass the glplots. 
#Mden individuals to remove: BWLib_L5-10, BWLib_L5-13, BWLib_L5-12, BWLib_L5-16, BWLib_L5-15, BWLib_L5-11
#Zcav individuals to remove: BWLib_L3-4, BWLib_L3-6, BWLib_L3-7
######looking at Tajima's D
#Atlantic Ocean Z. cav.
zczc.atl.D.1000<-read.table("./zczc129.atl.1000.Tajima.D", header=TRUE)
zczc.atl.D.5000<-read.table("./zczc129.atl.5000.Tajima.D", header=TRUE)
zczc.atl.D.10000<-read.table("./zczc129.atl.10000.Tajima.D", header=TRUE)
zczc.atl.D.100000<-read.table("./zczc129.atl.100000.Tajima.D", header=TRUE)

zc.atl.1<-ggplot(zczc.atl.D.1000, aes(x=zczc.atl.D.1000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc.atl.D.1000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc.atl.D.1000$TajimaD, na.rm=TRUE) +1), y=5000, 
           label=(paste("mean Tajima's D=", round(mean(zczc.atl.D.1000$TajimaD, na.rm=TRUE), 4))))

zc.atl.5<-ggplot(zczc.atl.D.5000, aes(x=zczc.atl.D.5000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc.atl.D.5000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc.atl.D.5000$TajimaD, na.rm=TRUE) +1), y=5000, 
           label=(paste("mean Tajima's D=", round(mean(zczc.atl.D.5000$TajimaD, na.rm=TRUE), 4))))

zc.atl.10<-ggplot(zczc.atl.D.10000, aes(x=zczc.atl.D.10000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc.atl.D.10000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc.atl.D.10000$TajimaD, na.rm=TRUE) +1), y=5000, 
           label=(paste("mean Tajima's D=", round(mean(zczc.atl.D.10000$TajimaD, na.rm=TRUE), 4))))

zc.atl.100<-ggplot(zczc.atl.D.100000, aes(x=zczc.atl.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc.atl.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc.atl.D.100000$TajimaD, na.rm=TRUE) +1), y=5000, 
           label=(paste("mean Tajima's D=", round(mean(zczc.atl.D.100000$TajimaD, na.rm=TRUE), 4))))
##combine all plots together
plot_grid(zc.atl.1, zc.atl.5, zc.atl.10, zc.atl.100)


#Indopacific Ocean Z. cav.
zczc.indo.D.1000<-read.table("./zczc129.indo.1000.Tajima.D", header=TRUE)
zczc.indo.D.5000<-read.table("./zczc129.indo.5000.Tajima.D", header=TRUE)
zczc.indo.D.10000<-read.table("./zczc129.indo.10000.Tajima.D", header=TRUE)
zczc.indo.D.100000<-read.table("./zczc129.indo.100000.Tajima.D", header=TRUE)

zc.indo.1<-ggplot(zczc.indo.D.1000, aes(x=zczc.indo.D.1000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc.indo.D.1000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc.indo.D.1000$TajimaD, na.rm=TRUE) +1), y=5000, 
           label=(paste("mean Tajima's D=", round(mean(zczc.indo.D.1000$TajimaD, na.rm=TRUE), 4))))

zc.indo.5<-ggplot(zczc.indo.D.5000, aes(x=zczc.indo.D.5000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc.indo.D.5000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc.indo.D.5000$TajimaD, na.rm=TRUE) +1), y=5000, 
           label=(paste("mean Tajima's D=", round(mean(zczc.indo.D.5000$TajimaD, na.rm=TRUE), 4))))

zc.indo.10<-ggplot(zczc.indo.D.10000, aes(x=zczc.indo.D.10000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc.indo.D.10000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc.indo.D.10000$TajimaD, na.rm=TRUE) +1), y=5000, 
           label=(paste("mean Tajima's D=", round(mean(zczc.indo.D.10000$TajimaD, na.rm=TRUE), 4))))

zc.indo.100<-ggplot(zczc.indo.D.100000, aes(x=zczc.indo.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc.indo.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc.indo.D.100000$TajimaD, na.rm=TRUE) +1), y=5000, 
           label=(paste("mean Tajima's D=", round(mean(zczc.indo.D.100000$TajimaD, na.rm=TRUE), 4))))
##combine all plots together
plot_grid(zc.indo.1, zc.indo.5, zc.indo.10, zc.indo.100)

#Mediterranean Sea Z. cav.
zczc.med.D.1000<-read.table("./zczc129.med.1000.Tajima.D", header=TRUE)
zczc.med.D.5000<-read.table("./zczc129.med.5000.Tajima.D", header=TRUE)
zczc.med.D.10000<-read.table("./zczc129.med.10000.Tajima.D", header=TRUE)
zczc.med.D.100000<-read.table("./zczc129.med.100000.Tajima.D", header=TRUE)

zc.med.1<-ggplot(zczc.med.D.1000, aes(x=zczc.med.D.1000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc.med.D.1000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc.med.D.1000$TajimaD, na.rm=TRUE) +1), y=2000, 
           label=(paste("mean Tajima's D=", round(mean(zczc.med.D.1000$TajimaD, na.rm=TRUE), 4))))

zc.med.5<-ggplot(zczc.med.D.5000, aes(x=zczc.med.D.5000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc.med.D.5000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc.med.D.5000$TajimaD, na.rm=TRUE) +1), y=2000, 
           label=(paste("mean Tajima's D=", round(mean(zczc.med.D.5000$TajimaD, na.rm=TRUE), 4))))

zc.med.10<-ggplot(zczc.med.D.10000, aes(x=zczc.med.D.10000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc.med.D.10000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc.med.D.10000$TajimaD, na.rm=TRUE) +1), y=2000, 
           label=(paste("mean Tajima's D=", round(mean(zczc.med.D.10000$TajimaD, na.rm=TRUE), 4))))

zc.med.100<-ggplot(zczc.med.D.100000, aes(x=zczc.med.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc.med.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc.med.D.100000$TajimaD, na.rm=TRUE) +1), y=2000, 
           label=(paste("mean Tajima's D=", round(mean(zczc.med.D.100000$TajimaD, na.rm=TRUE), 4))))
##combine all plots together
plot_grid(zc.med.1, zc.med.5, zc.med.10, zc.med.100)



###Tajima's D for Z. cav populations

#Atlantic - Bahamas
zcav.atl.bahamas.D.100000<-read.table("./zczc129.atl-bahamas.100000.Tajima.D", header=TRUE)

png("zcav.atl.bahamas.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.atl.bahamas<-ggplot(zcav.atl.bahamas.D.100000, aes(x=zcav.atl.bahamas.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.atl.bahamas.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.atl.bahamas.D.100000$TajimaD, na.rm=TRUE) +1.25), y=2500, 
           label=(paste("Atl-Bahamas mean Tajima's D=", round(mean(zcav.atl.bahamas.D.100000$TajimaD, na.rm=TRUE), 4)))) +
  annotate("text", x=(mean(zcav.atl.bahamas.D.100000$TajimaD, na.rm=TRUE) +1.25), y=2300, 
           label=(paste("n=11"))) +
  xlab("Tajima's D")
zc.atl.bahamas
dev.off()

#Atlantic-East
zcav.atl.east.D.100000<-read.table("./zczc129.atl-east.100000.Tajima.D", header=TRUE)

png("zcav.atl.east.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.atl.east<-ggplot(zcav.atl.east.D.100000, aes(x=zcav.atl.east.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.atl.east.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.atl.east.D.100000$TajimaD, na.rm=TRUE) +1.5), y=3500, 
           label=(paste("Atl-East mean Tajima's D=", round(mean(zcav.atl.east.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(zcav.atl.east.D.100000$TajimaD, na.rm=TRUE) +1.5), y=3200, 
           label=(paste("n=35")))+
  xlab("Tajima's D")
zc.atl.east
dev.off()

#Atlantic - Madeira
zcav.atl.madeira.D.100000<-read.table("./zczc129.atl-madeira.100000.Tajima.D", header=TRUE)

png("zcav.atl.madeira.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.atl.madeira<-ggplot(zcav.atl.madeira.D.100000, aes(x=zcav.atl.madeira.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.atl.madeira.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.atl.madeira.D.100000$TajimaD, na.rm=TRUE) +.5), y=2550, 
           label=(paste("Atl-Madeira mean Tajima's D=", round(mean(zcav.atl.madeira.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(zcav.atl.madeira.D.100000$TajimaD, na.rm=TRUE) +.5), y=2350, 
           label=(paste("n=4")))+
  xlab("Tajima's D")
zc.atl.madeira
dev.off()

#Atlantic - Other
zcav.atl.other.D.100000<-read.table("./zczc129.atl-other.100000.Tajima.D", header=TRUE)

png("zcav.atl.other.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.atl.other<-ggplot(zcav.atl.other.D.100000, aes(x=zcav.atl.other.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.atl.other.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.atl.other.D.100000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("Atl-Other mean Tajima's D=", round(mean(zcav.atl.other.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(zcav.atl.other.D.100000$TajimaD, na.rm=TRUE) +1), y=2300, 
           label=(paste("n=7")))+
  xlab("Tajima's D")
zc.atl.other
dev.off()

#Atlantic - Spain
zcav.atl.spain.D.100000<-read.table("./zczc129.atl-spain.100000.Tajima.D", header=TRUE)

png("zcav.atl.spain.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.atl.spain<-ggplot(zcav.atl.spain.D.100000, aes(x=zcav.atl.spain.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.atl.spain.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.atl.spain.D.100000$TajimaD, na.rm=TRUE) -1), y=2500, 
           label=(paste("Atl-Spain mean Tajima's D=", round(mean(zcav.atl.spain.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(zcav.atl.spain.D.100000$TajimaD, na.rm=TRUE) -1), y=2300, 
           label=(paste("n=2")))+
  xlab("Tajima's D")
zc.atl.spain
dev.off()

#Indopacific - central
zcav.indo.central.D.100000<-read.table("./zczc129.indo-central.100000.Tajima.D", header=TRUE)

png("zcav.indo.central.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.indo.central<-ggplot(zcav.indo.central.D.100000, aes(x=zcav.indo.central.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.indo.central.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.indo.central.D.100000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("Indo-Central mean Tajima's D=", round(mean(zcav.indo.central.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(zcav.indo.central.D.100000$TajimaD, na.rm=TRUE) +1), y=2300, 
           label=(paste("n=5")))+
  xlab("Tajima's D")
zc.indo.central
dev.off()

#Indopacific - east
zcav.indo.east.D.100000<-read.table("./zczc129.indo-east.100000.Tajima.D", header=TRUE)

png("zcav.indo.east.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.indo.east<-ggplot(zcav.indo.east.D.100000, aes(x=zcav.indo.east.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.indo.east.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.indo.east.D.100000$TajimaD, na.rm=TRUE) +1.5), y=2500, 
           label=(paste("Indo-East mean Tajima's D=", round(mean(zcav.indo.east.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(zcav.indo.east.D.100000$TajimaD, na.rm=TRUE) +1.5), y=2300, 
           label=(paste("n=17")))+
  xlab("Tajima's D")
zc.indo.east
dev.off()

#Indopacific - mexico
zcav.indo.mexico.D.100000<-read.table("./zczc129.indo-mexico.100000.Tajima.D", header=TRUE)

png("zcav.indo.mexico.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.indo.mexico<-ggplot(zcav.indo.mexico.D.100000, aes(x=zcav.indo.mexico.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.indo.mexico.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.indo.mexico.D.100000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("Indo-Mexico mean Tajima's D=", round(mean(zcav.indo.mexico.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(zcav.indo.central.D.100000$TajimaD, na.rm=TRUE) +1), y=2300, 
           label=(paste("n=2")))+
  xlab("Tajima's D")
zc.indo.mexico
dev.off()

#Indopacific - south
zcav.indo.south.D.100000<-read.table("./zczc129.indo-south.100000.Tajima.D", header=TRUE)

png("zcav.indo.south.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.indo.south<-ggplot(zcav.indo.south.D.100000, aes(x=zcav.indo.south.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.indo.south.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.indo.south.D.100000$TajimaD, na.rm=TRUE) +1.25), y=2500, 
           label=(paste("Indo-South mean Tajima's D=", round(mean(zcav.indo.south.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(zcav.indo.south.D.100000$TajimaD, na.rm=TRUE) +1.25), y=2300, 
           label=(paste("n=9")))+
  xlab("Tajima's D")
zc.indo.south
dev.off()

#Indopacific - spac
zcav.indo.spac.D.100000<-read.table("./zczc129.indo-spac.100000.Tajima.D", header=TRUE)

png("zcav.indo.spac.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.indo.spac<-ggplot(zcav.indo.spac.D.100000, aes(x=zcav.indo.spac.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.indo.spac.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.indo.spac.D.100000$TajimaD, na.rm=TRUE) +1.1), y=2500, 
           label=(paste("Indo=S. Pacific mean Tajima's D=", round(mean(zcav.indo.spac.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(zcav.indo.spac.D.100000$TajimaD, na.rm=TRUE) +1.1), y=2300, 
           label=(paste("n=3")))+
  xlab("Tajima's D")
zc.indo.spac
dev.off()

#Mediterranean - east
zcav.med.east.D.100000<-read.table("./zczc129.med-east.100000.Tajima.D", header=TRUE)

png("zcav.med.east.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.med.east<-ggplot(zcav.med.east.D.100000, aes(x=zcav.med.east.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.med.east.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.med.east.D.100000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("Med-East mean Tajima's D=", round(mean(zcav.med.east.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(zcav.med.east.D.100000$TajimaD, na.rm=TRUE) +1), y=2300, 
           label=(paste("n=12")))+
  xlab("Tajima's D")
zc.med.east
dev.off()

#Mediterranean - france
zcav.med.france.D.100000<-read.table("./zczc129.med-france.100000.Tajima.D", header=TRUE)

png("zcav.med.france.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.med.france<-ggplot(zcav.med.france.D.100000, aes(x=zcav.med.france.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.med.france.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.med.france.D.100000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("Med-France mean Tajima's D=", round(mean(zcav.med.france.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(zcav.med.france.D.100000$TajimaD, na.rm=TRUE) +1), y=2300, 
           label=(paste("n=2")))+
  xlab("Tajima's D")
zc.med.france
dev.off()

#Mediterranean - Italy
zcav.med.italy.D.100000<-read.table("./zczc129.med-italy.100000.Tajima.D", header=TRUE)

png("zcav.med.italy.D.png", width = 6, height = 4, units = 'in', res = 300)
zc.med.italy<-ggplot(zcav.med.italy.D.100000, aes(x=zcav.med.italy.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zcav.med.italy.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zcav.med.italy.D.100000$TajimaD, na.rm=TRUE) +1.1), y=2500, 
           label=(paste("Med-Italy mean Tajima's D=", round(mean(zcav.med.italy.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(zcav.med.italy.D.100000$TajimaD, na.rm=TRUE) +1.1), y=2300, 
           label=(paste("n=20")))+
  xlab("Tajima's D")
zc.med.italy
dev.off()


# Tajima’s D - M. densirostris --------------------------------------------
#Atlantic Ocean M. den.
mden.atl.D.1000<-read.table("./mden43.atl.1000.Tajima.D", header=TRUE)
mden.atl.D.5000<-read.table("./mden43.atl.5000.Tajima.D", header=TRUE)
mden.atl.D.10000<-read.table("./mden43.atl.10000.Tajima.D", header=TRUE)
mden.atl.D.100000<-read.table("./mden43.atl.100000.Tajima.D", header=TRUE)

md.atl.1<-ggplot(mden.atl.D.1000, aes(x=mden.atl.D.1000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.atl.D.1000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.atl.D.1000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("mean Tajima's D=", round(mean(mden.atl.D.1000$TajimaD, na.rm=TRUE), 4))))

md.atl.5<-ggplot(mden.atl.D.5000, aes(x=mden.atl.D.5000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.atl.D.5000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.atl.D.5000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("mean Tajima's D=", round(mean(mden.atl.D.5000$TajimaD, na.rm=TRUE), 4))))

md.atl.10<-ggplot(mden.atl.D.10000, aes(x=mden.atl.D.10000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.atl.D.10000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.atl.D.10000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("mean Tajima's D=", round(mean(mden.atl.D.10000$TajimaD, na.rm=TRUE), 4))))

md.atl.100<-ggplot(mden.atl.D.100000, aes(x=mden.atl.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.atl.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.atl.D.100000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("mean Tajima's D=", round(mean(mden.atl.D.100000$TajimaD, na.rm=TRUE), 4))))
##combine all plots together
plot_grid(md.atl.1, md.atl.5, md.atl.10, md.atl.100)

#Indopacific Ocean M. den.
mden.indo.D.1000<-read.table("./mden43.indo.1000.Tajima.D", header=TRUE)
mden.indo.D.5000<-read.table("./mden43.indo.5000.Tajima.D", header=TRUE)
mden.indo.D.10000<-read.table("./mden43.indo.10000.Tajima.D", header=TRUE)
mden.indo.D.100000<-read.table("./mden43.indo.100000.Tajima.D", header=TRUE)

md.indo.1<-ggplot(mden.indo.D.1000, aes(x=mden.indo.D.1000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.indo.D.1000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.indo.D.1000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("mean Tajima's D=", round(mean(mden.indo.D.1000$TajimaD, na.rm=TRUE), 4))))

md.indo.5<-ggplot(mden.indo.D.5000, aes(x=mden.indo.D.5000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.indo.D.5000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.indo.D.5000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("mean Tajima's D=", round(mean(mden.indo.D.5000$TajimaD, na.rm=TRUE), 4))))

md.indo.10<-ggplot(mden.indo.D.10000, aes(x=mden.indo.D.10000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.indo.D.10000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.indo.D.10000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("mean Tajima's D=", round(mean(mden.indo.D.10000$TajimaD, na.rm=TRUE), 4))))

md.indo.100<-ggplot(mden.indo.D.100000, aes(x=mden.indo.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.indo.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.indo.D.100000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("mean Tajima's D=", round(mean(mden.indo.D.100000$TajimaD, na.rm=TRUE), 4))))
##combine all plots together
plot_grid(md.indo.1, md.indo.5, md.indo.10, md.indo.100)

##Populations for Mden


# Atlantic - Bahamas
#Mediterranean - Italy
mden.atl.bahamas.D.100000<-read.table("./mden43.atl-bahamas.100000.Tajima.D", header=TRUE)

png("mden.atl.bahamas.D.png", width = 6, height = 4, units = 'in', res = 300)
md.atl.bahamas<-ggplot(mden.atl.bahamas.D.100000, aes(x=mden.atl.bahamas.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.atl.bahamas.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.atl.bahamas.D.100000$TajimaD, na.rm=TRUE) +1.15), y=2500, 
           label=(paste("Atl-Bahamas mean Tajima's D=", round(mean(mden.atl.bahamas.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(mden.atl.bahamas.D.100000$TajimaD, na.rm=TRUE) +1.15), y=2300, 
           label=(paste("n=7")))+
  xlab("Tajima's D")
md.atl.bahamas
dev.off()

# Atlantic - East
mden.atl.east.D.100000<-read.table("./mden43.atl-east.100000.Tajima.D", header=TRUE)

png("mden.atl.east.D.png", width = 6, height = 4, units = 'in', res = 300)
md.atl.east<-ggplot(mden.atl.east.D.100000, aes(x=mden.atl.east.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.atl.east.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.atl.east.D.100000$TajimaD, na.rm=TRUE) +1.5), y=2500, 
           label=(paste("Atl-East mean Tajima's D=", round(mean(mden.atl.east.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(mden.atl.east.D.100000$TajimaD, na.rm=TRUE) +1.5), y=2300, 
           label=(paste("n=16")))+
  xlab("Tajima's D")
md.atl.east
dev.off()

# Atlantic - Other
mden.atl.other.D.100000<-read.table("./mden43.atl-other.100000.Tajima.D", header=TRUE)

png("mden.atl.other.D.png", width = 6, height = 4, units = 'in', res = 300)
md.atl.other<-ggplot(mden.atl.other.D.100000, aes(x=mden.atl.other.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.atl.other.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.atl.other.D.100000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("Atl-Other mean Tajima's D=", round(mean(mden.atl.other.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(mden.atl.other.D.100000$TajimaD, na.rm=TRUE) +1), y=2300, 
           label=(paste("n=5")))+
  xlab("Tajima's D")
md.atl.other
dev.off()

# Indopacific - Africa
mden.indo.africa.D.100000<-read.table("./mden43.indo-africa.100000.Tajima.D", header=TRUE)

png("mden.indo.africa.D.png", width = 6, height = 4, units = 'in', res = 300)
md.indo.africa<-ggplot(mden.indo.africa.D.100000, aes(x=mden.indo.africa.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.indo.africa.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.indo.africa.D.100000$TajimaD, na.rm=TRUE) +1.1), y=2500, 
           label=(paste("Indo-S. Africa mean Tajima's D=", round(mean(mden.indo.africa.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(mden.indo.africa.D.100000$TajimaD, na.rm=TRUE) +1.1), y=2300, 
           label=(paste("n=5")))+
  xlab("Tajima's D")
md.indo.africa
dev.off()

# Indopacific - Hawaii
mden.indo.hawaii.D.100000<-read.table("./mden43.indo-hawaii.100000.Tajima.D", header=TRUE)

png("mden.indo.hawaii.D.png", width = 6, height = 4, units = 'in', res = 300)
md.indo.hawaii<-ggplot(mden.indo.hawaii.D.100000, aes(x=mden.indo.hawaii.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.indo.hawaii.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.indo.hawaii.D.100000$TajimaD, na.rm=TRUE) +1.2), y=2500, 
           label=(paste("Indo-Hawaii mean Tajima's D=", round(mean(mden.indo.hawaii.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(mden.indo.hawaii.D.100000$TajimaD, na.rm=TRUE) +1.2), y=2300, 
           label=(paste("n=6")))+
  xlab("Tajima's D")
md.indo.hawaii
dev.off()

# Indopacific - South
mden.indo.south.D.100000<-read.table("./mden43.indo-south.100000.Tajima.D", header=TRUE)

png("mden.indo.south.D.png", width = 6, height = 4, units = 'in', res = 300)
md.indo.south<-ggplot(mden.indo.south.D.100000, aes(x=mden.indo.south.D.100000$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(mden.indo.south.D.100000$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden.indo.south.D.100000$TajimaD, na.rm=TRUE) +1), y=2500, 
           label=(paste("Indo-South mean Tajima's D=", round(mean(mden.indo.south.D.100000$TajimaD, na.rm=TRUE), 4))))+
  annotate("text", x=(mean(mden.indo.south.D.100000$TajimaD, na.rm=TRUE) +1), y=2300, 
           label=(paste("n=3")))+
  xlab("Tajima's D")
md.indo.south
dev.off()


# Phylogenetic trees ------------------------------------------------------
###
##Using ggtree and tibble GOOD

library(tibble)
library(dbplyr)
library(ggplot2)
library(ggtree)

#convert tree to a tibble file
zczc125.tree.tib

#make another tibble file with ocean data
d<-tibble(label=zczc125.ind$ind, pop=zczc125.ind$ocean)
d

#join tibble file with ocean data
zczc125.tib<-full_join(zczc125.tree.tib, d, by='label')
zczc125.tib

#convert tibble to a treedata file
zczc125.tibtree<-as.treedata(zczc125.tib)
zczc125.tibtree

#plot new tibble tree file in ggtree
ggtree(zczc125.tibtree, layout="circular") +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=my_zissou) +
  theme(legend.position=c(0.2,0), legend.title = element_blank())


plot(zczc125.tree, typ="fan", show.tip=FALSE)
tiplabels(pch=20, col=my_zissou, cex=1.5)
legend('bottomleft', legend = c("ATL", "INDO-PAC", "MED"), bty="n", fill = cols, border = FALSE, cex =1)

##Using ggtree and tibble

library(tibble)
library(dplyr)
library(ggplot2)
library(ggtree)
library(treeio)
#convert tree to a tibble file
mden43.tree<-nj(dist(as.matrix(mden43.gl)))
mden43.tree.tib<-as.tibble(mden43.tree)
boot.phylo(mden43.tree)

#make another tibble file with ocean data
mden43.indnames<-as.data.frame(mden43.gl$ind.names)
mden43.oceans<-as.data.frame(mden43.gl$pop)
mden43.ind<-cbind(mden43.indnames, mden43.oceans)
colnames(mden43.ind)<-c("label","ocean")
mden43.ind
e<-tibble(label=mden43.ind$label, pop=mden43.ind$ocean)
e

#another tibble file with population data
mden43.pops<-as.data.frame(mden43.pop.gl$pop)
mden43.pop.ind<-cbind(mden43.indnames, mden43.pops)
colnames(mden43.pop.ind)<-c("label", "pop")
mden43.pop.ind
f<-tibble(label=mden43.pop.ind$label, pop=mden43.pop.ind$pop)
f

#join tibble file with ocean data
mden43.tib<-full_join(mden43.tree.tib, e, by='label')
mden43.tib

#join the tibble file with the population data
mden43.pop.tib<-full_join(mden43.tree.tib, f, by="label")
mden43.pop.tib

#convert tibble to a treedata file
mden43.tibtree<-as.treedata(mden43.tib)
mden43.tibtree

mden43.pop.tibtree<-as.treedata(mden43.pop.tib)
mden43.pop.tibtree

#plot new tibble tree file in ggtree
ggtree(mden43.tibtree, layout="circular") +
  geom_tippoint(aes(col=factor(pop)), size=3) +
  scale_colour_manual(values=my_zissou[c(2,1)]) +
  theme(legend.position=c(0.08,0.2), legend.title = element_blank())

ggtree(mden43.tibtree, layout="slanted") +
  geom_tippoint(aes(col=factor(pop)), size=2) +
  scale_colour_manual(values=my_zissou) +
  theme(legend.position=c(0.1,.9), legend.title = element_blank())

p<-ggtree(mden43.tibtree, branch.length = "none") +
  geom_tiplab()+
  geom_tippoint(aes(col=factor(pop)), size=2) +
  xlim(0, 15) +
  scale_colour_manual(values=paired.col[c(2,4)]) +
  theme(legend.position=c(0.1,.9), legend.title = element_blank())
p

  flip(p,45,47)
flip(p, 54, 68)

ggtree(mden43.pop.tibtree, branch.length = "none") +
  geom_tippoint(aes(col=factor(pop)), size=2) +
  theme(legend.position=c(0.1,.85), legend.title = element_blank())

####Actual Paper Trees
#convert tree to a tibble file
zczc125.tree<-bionjs(dist(as.matrix(zczc125.gl)))
zczc125.tree.tib<-as_tibble(zczc125.tree)


zczc125.tree.tib

zczc125.ind[,3]<-zczc125.pop1.gl$pop
zczc125.ind[,4]<-zczc125.pop2.gl$pop
colnames(zczc125.ind)<-c("ind", "ocean", "pop1", "pop2")
head(zczc125.ind)

zczc125.ind[,5]<-c("Italy-Ligurian1",
"Italy-Ligurian2",
"Italy-Ligurian3",
"Italy-Ligurian4",
"Italy-Ligurian5",
"Italy-Ligurian6",
"France-MED1",
"Corfu1",
"Italy-Ligurian7",
"Corfu2",
"Italy-Ligurian8",
"Corfu3",
"France-MED2",
"Italy-Ligurian9",
"Italy-Ionian1",
"USA-HI1",
"USA-HI2",
"Ireland4",
"Canary Islands, TF2",
"Crete2",
"Ireland5",
"Croatia",
"Madeira Islands2",
"USA-CA4",
"Italy-Ionian2",
"Italy-Ligurian13",
"USA-AK",
"Crete3",
"Corfu4",
"Corfu5",
"Greece",
"Israel1",
"Israel2",
"Canary Islands, FV1",
"Canary Islands, FV2",
"Canary Islands, EH1",
"Canary Islands, TF1",
"Canary Islands, EH2",
"France-ATL1",
"Canary Islands, EH3",
"Canary Islands, EH4",
"France-ATL2",
"France-ATL3",
"Ireland1",
"Ireland2",
"Madeira Islands1",
"Canary Islands, FV3",
"Scotland1",
"USA-FL",
"France-ATL4",
"Canary Islands, EH5",
"Spain1",
"France-ATL5",
"USA-NC",
"Spain2",
"Canary Islands, EH6",
"Canary Islands, EH7",
"Scotland2",
"Scotland3",
"Ireland3",
"Canary Islands, FV4",
"Puerto Rico1",
"Bahamas1",
"Canary Islands, EH8",
"USA-CA2",
"USA-CA3",
"Bahamas2",
"Crete1",
"Samoa",
"Bahamas3",
"Mexico1",
"Bahamas4",
"Scotland4",
"Canary Islands, EH9",
"Italy-Ligurian10",
"Italy-Ligurian15",
"Australia1",
"Italy-Ligurian16",
"Italy-Ionian3",
"Italy-Ligurian17",
"Canary Islands, TF3",
"Puerto Rico3",
"Puerto Rico4",
"Bahamas5",
"Scotland5",
"Scotland6",
"Scotland7",
"Canary Islands",
"Bahamas6",
"Bahamas7",
"Bahamas8",
"Bahamas9",
"Scotland8",
"USA-CA1",
"USA-OR1",
"Italy-Ligurian11",
"Italy-Ligurian14",
"Canada-West1",
"Canada-West2",
"Puerto Rico2",
"Virgin Islands",
"USA-CA6",
"USA-CA7",
"USA-HI3",
"USA-CA8",
"USA-CA9",
"USA-CA10",
"USA-OR2",
"USA-CA11",
"USA-HI4",
"Mexico2",
"USA-CA12",
"USA-Johnston Atoll",
"South Africa",
"New Zealand1",
"Philippines",
"Chile",
"New Zealand2",
"New Zealand3",
"New Zealand4",
"Italy-Ligurian12",
"USA-CA5",
"New Zealand5",
"New Zealand6",
"Australia2")

colnames(zczc125.ind)<-c("ind"  , "ocean" ,"pop1" , "pop2" , "region" )
zczc125.ind
#make another tibble file with ocean data
d<-tibble(label=zczc125.ind$ind, pop=zczc125.ind$ocean)
p1<-tibble(label=zczc125.ind$ind, pop=zczc125.ind$pop1)
p2<-tibble(label=zczc125.ind$ind, pop=zczc125.ind$pop2)
r<-tibble(label=zczc125.ind$ind, )

#make another tibble file with population data
e<-tibble(label=zczc125.ind$ind, pop=zczc125.ind$ocean)
ep1<-tibble(label=zczc125.ind$ind, pop=zczc125.ind$pop1)
ep2<-tibble(label=zczc125.ind$ind, pop=zczc125.ind$pop2)

#join tibble file with ocean data
zczc125.tib<-full_join(zczc125.tree.tib, d, by='label')
zczc125.pop1.tib<-full_join(zczc125.tree.tib, p1, by='label')
zczc125.pop2.tib<-full_join(zczc125.tree.tib, p2, by='label')

#convert tibble to a treedata file
zczc125.tibtree<-as.treedata(zczc125.tib)
zczc125.tibtree
zczc125.pop1.tibtree<-as.treedata(zczc125.pop1.tib)
zczc125.pop2.tibtree<-as.treedata(zczc125.pop2.tib)



#plot new tibble tree file in ggtree
ggtree(zczc125.tibtree, layout="circular") +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=my_zissou) +
  theme(legend.position=c(0.2,.1), legend.title = element_blank())

ggtree(zczc125.tibtree) +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=my_zissou) +
  theme(legend.position=c(.90,.2), legend.title = element_blank())

ggtree(zczc125.tibtree, branch.length = "none") +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  geom_tiplab()+
  xlim(0, 30) +
  ylim(100,150) +
  scale_colour_manual(values=my_zissou) +
  theme(legend.position=c(.1,.9), legend.title = element_blank())


ggtree(zczc125.pop1.tibtree, branch.length = "none") +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  theme(legend.position=c(.1,.8), legend.title=element_blank())


ggtree(zczc125.pop2.tibtree, branch.length="none") +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  theme(legend.position=c(.1,.8), legend.title=element_blank())


##Plotting trees differently so i can see which individuals are shoping up in different groups
install.packages("ggrepel")
library(ggrepel)
ggtree(zczc125.tibtree, branch.length = "none") +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=my_zissou) +
  geom_text(aes(label=zczc125.tib$label), position = position_nudge(x = 1), size=1) +
  theme(legend.position=c(0.2,.8), legend.title = element_blank())
ggtree(mden43.tibtree, branch.length = "none") +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=my_zissou) +
  geom_text(aes(label=mden43.tib$label), position = position_nudge(x = 1), size=1) +
  theme(legend.position=c(0.2,.8), legend.title = element_blank())


ggtree(zczc125.pop1.tibtree, branch.length = "none") +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  geom_text(aes(label=zczc125.tib$label), position = position_nudge(x = 1), size=1) +
  theme(legend.position=c(0.2,.8), legend.title = element_blank())

ggtree(zczc125.tibtree, layout="circular") +
  geom_tippoint(aes(col=factor(pop)), size=3) +
  scale_colour_manual(values=my_zissou) +
  geom_text(aes(label=zczc125.tib$label), position = position_nudge(x = 1), size=1) +
  theme(legend.position=c(0.2,.8), legend.title = element_blank())

ggtree(mden43.tibtree) +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=my_zissou) +
  geom_text(aes(label=mden43.tib$label), position = position_nudge(x = 1), size=2) +
  theme(legend.position=c(0.9,.8), legend.title = element_blank())

ggtree(mden43.pop.tibtree) +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=my_zissou) +
  geom_text(aes(label=mden43.pop.tib$label), position = position_nudge(x = .5), size=2) +
  theme(legend.position=c(0.2,.8), legend.title = element_blank())

##Trying to boostrap tree
zczc125.gl
mden43.gl

zczc125.dist<-dist(zczc125.gl)
mden43.dist<-dist(mden43.gl)

zczc125.boot.tree<-aboot(zczc125.gl, tree="upgma", distance=bitwise.dist, sample=100, showtree=F, cutoff=50, quiet=T)
zczc125.boot.tree<-aboot(zczc125.gi, tree="upgma", sample=100, showtree=F, cutoff=50, quiet=T)


cols <- brewer.pal(n = nPop(zczc125.gl), name = "Dark2")
plot.phylo(zczc125.boot.tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(zczc125.gl)])
nodelabels(zczc125.boot.tree$node.label, adj = c(2, -1), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend('topleft', legend = c("Atlantic","Indo-Pacific","Mediterranean"), fill = cols, border = FALSE, bty = "n", cex =1)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")
dev.off()

#My trees
zczc125.tibtree@phylo$tip.label<-zczc125.tipname.pop2
zczc125.tibtree@data
zczc125.pop1.tibtree
zczc125.pop2.tibtree

zczc125.tipname.pop1<-paste(zczc125.gl$pop, zczc125.pop1.gl$pop, zczc125.gl$ind.names, sep="_")
zczc125.tipname.pop2<-paste(zczc125.gl$pop, zczc125.pop2.gl$pop, zczc125.gl$ind.names, sep="_")

zczc125.region.tree<-zczc125.tibtree
zczc125.region.tree@phylo$tip.label<-zczc125.ind$region


write.tree(zczc125.tibtree@phylo, file="./zczc125.tibtree.pop2")
write.tree(zczc125.region.tree@phylo, file="./zczc125.region.tree")

mden43.tibtree
mden43.pop.tibtree

ggtree(zczc125.region.tree) +
  geom_tiplab()

# Snps per locus ----------------------------------------------------------

###Loading vcfs with all snps data
##
###Workout how many snps per locus in the vcf file (with no write_random_snps option (all_snps_per_locus))

##M. densirostris
#use vcftools either on cluster or laptop to remove the bad glPlot individuals. Need to also re-calculate allele frequencies and take out newly monomorphic loci using --maf function

vcftools --vcf mden49.snps.vcf --remove mden.ind.rem --maf 0.001 --recode --recode-INFO-all

#using terminal
grep -v "^#" mden43.snps.vcf | cut -f 1 | sort | uniq -c | sort -n > mden43_snps_per_locus.txt

##read in file to R
mden43.snps.locus<-read.table("./mden43_snps_per_locus.txt")
head(mden43.snps.locus)
tail(mden43.snps.locus)

#make a histogram of values
png("mden43_snps_per_locus.png", width = 6, height = 4, units = 'in', res = 300)
qplot(mden43.snps.locus$V1, geom="histogram", main="M. densirostris", xlab="Number of SNPs per locus", ylab="Frequency", binwidth=1, width=1) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept=mean(mden43.snps.locus$V1, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(mden43.snps.locus$V1, na.rm=TRUE) +15), y=2500, 
           label=(paste("mean no. SNPs per locus=", round(mean(mden43.snps.locus$V1, na.rm=TRUE), 4))))
dev.off()

##Z. cavirostris
#recoding vcf file removing bad snps
vcftools --vcf zczc132.snps.vcf --remove zczc.ind.rem --maf 0.001 --recode --recode-INFO-all

#using terminal
grep -v "^#" zczc129.snps.vcf | cut -f 1 | sort | uniq -c | sort -n > zczc129_snps_per_locus.txt

##read in file to R
zczc129.snps.locus<-read.table("./zczc129_snps_per_locus.txt")
head(zczc129.snps.locus)

#make a histogram of values
png("zczc129_snps_per_locus.png", width = 6, height = 4, units = 'in', res = 300)
qplot(zczc129.snps.locus$V1, geom="histogram", main="Z. cavirostris", xlab="Number of SNPs per locus", ylab="Frequency", binwidth=1, width=1) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept=mean(zczc129.snps.locus$V1, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc129.snps.locus$V1, na.rm=TRUE) +20), y=2500, 
           label=(paste("mean no. SNPs per locus=", round(mean(zczc129.snps.locus$V1, na.rm=TRUE), 4))))
dev.off()



# Heterozygosity ----------------------------------------------------------
###Looking at Heterozygosity in other ways since dartR results had lots of Zeros. 

##Heterozygosity with vcftools (per individual value)
##code for vcftools in Terminal
##using vcf files with only one snp per locus

#remove bad glPlot individuals, remove monomorphic loci
vcftools --vcf mden49.snps.vcf --remove mden.ind.rem --maf 0.001 --out mden43.snps --recode --recode-INFO-all
#calculate heterozygous sites
vcftools --vcf mden43.snps.vcf --het 

#load in table
mden43.het<-read.table("./mden43.het", header=TRUE, sep="")
#calculate observed and expected heterozygosities
Ho<-mden43.het$O.HOM./mden43.het$N_SITES
He<-mden43.het$E.HOM./mden43.het$N_SITES
#add as new columns in table
mden43.het<-cbind(mden43.het, Ho, He)
mden43.het
mden43.ho.hist<-hist(mden43.het$Ho)
mden43.he.hist<-hist(mden43.het$He)

#remove bad glPlot individuals, remove monomorphic loci
vcftools --vcf zczc132.snps.vcf --remove zczc.ind.rem --maf 0.001 --out zczc129.snps --recode --recode-INFO-all
#calculate heterozygous sites
vcftools --vcf zczc129.snps.vcf --het 

#load in table
zczc129.het<-read.table("./zczc129.het", header=TRUE, sep="")
#calculate observed and expected heterozygosities
Ho<-zczc129.het$O.HOM./zczc129.het$N_SITES
He<-zczc129.het$E.HOM./zczc129.het$N_SITES
#add as new columns in table
zczc129.het<-cbind(zczc129.het, Ho, He)
zczc129.het
zczc129.ho.hist<-hist(zczc129.het$Ho)
zczc129.he.hist<-hist(zczc129.het$He)

##Calculating EXPECTED Heterozygosity with adegenet 
mden43.gl
mden43.gp

mden43.pop.gl
mden43.pop.gp

zczc129.gl
zczc129.gp

zczc129.pop.gl
zczc129.pop.gp

Hs(mden43.gp)
#Atlantic Indopacific 
#0.1119443   0.1355113
Hs(zczc129.gp)
#Atlantic   Indopacific Mediterranean 
#0.13574301    0.13067124    0.09757143
Hs(mden43.pop.gp)
#Atl-Bahamas    Atl-East   Atl-Other Indo-Africa Indo-Hawaii  Indo-South          NA 
#0.1012754   0.1090775   0.1000685   0.1258402   0.1244822   0.1126243   0.0725622
Hs(zczc129.pop.gp)
#Med-Italy     Atl-East   Med-France  Atl-Madeira    Atl-Other     Med-East    Atl-Spain
#0.09486371   0.13354704   0.07293787   0.09994635   0.12566630   0.08940503   0.07395048 
#Indo-South    Indo-East    Indo-Spac  Atl-Bahamas Indo-Central  Indo-Mexico 
#0.12475436   0.12756310   0.10581805   0.12720595   0.11248742   0.09009238

mden43.gl.sum<-summary(mden43.gl)
mden43.gl.sum

#dartR population level calculations of OBSERVED heterozygosity
gl.report.heterozygosity(mden43.gl)
gl.report.heterozygosity(zczc129.gl)
gl.report.heterozygosity(mden43.pop.gl)
gl.report.heterozygosity(zczc129.pop.gl)

##STRATAG heterozygosity --> only calcualtes per locus for all samples together. Not very helpful
mean(exptdHet(mden43.gt))
mean(obsvdHet(mden43.gt))

#If you just call up the genotypes file though from STRATAG it will give the observed Heterozygosity per population. 
mden43.gt
# <<< gtypes created on 2019-02-12 10:53:48 >>>
#   
#   Contents: 43 samples, 13988 loci, 2 strata
# 
# Strata summary:
#   num.samples num.missing num.alleles prop.unique.alleles heterozygosity
# Atlantic             28   0.5331713    1.552116          0.07159708      0.1053599
# Indopacific          15   0.1607092    1.765156          0.18058336      0.1270102

mden43.pop.gt
# <<< gtypes created on 2019-02-12 12:00:59 >>>
#   
#   Contents: 43 samples, 13988 loci, 7 strata
# 
# Strata summary:
#   num.samples num.missing num.alleles prop.unique.alleles heterozygosity
# Atl-Bahamas           7  0.11409780    1.340149          0.06051616      0.1020844
# Atl-East             16  0.17586503    1.448813          0.05672719      0.1064709
# Atl-Other             5  0.24320846    1.316271          0.07888905      0.1063829
# Indo-Africa           5  0.06233915    1.443666          0.13311410      0.1247117
# Indo-Hawaii           6  0.04746926    1.461038                 NaN            NaN
# Indo-South            3  0.04310838    1.320060                 NaN            NaN
# 1 samples are unstratified

zczc129.gt
# <<< gtypes created on 2019-02-11 16:17:53 >>>
#   
#   Contents: 129 samples, 25059 loci, 3 strata
# 
# Strata summary:
#   num.samples num.missing num.alleles prop.unique.alleles heterozygosity
# Atlantic               59   1.4567221    1.896205          0.05089988     0.11848688
# Indopacific            36   0.8740971    1.780239          0.07458398     0.11597910
# Mediterranean          34   1.5721298    1.394868          0.04180135     0.08498699

zczc129.pop.gt
# <<< gtypes created on 2019-02-12 17:51:30 >>>
#   
#   Contents: 129 samples, 25059 loci, 13 strata
# 
# Strata summary:
#   num.samples num.missing num.alleles prop.unique.alleles heterozygosity
# Atl-Bahamas           11  0.17893771    1.522327          0.08144778     0.11805175
# Atl-East              35  0.95442755    1.813440          0.08208628     0.11903858
# Atl-Madeira            4  0.10682789    1.249292                 NaN            NaN
# Atl-Other              7  0.18855501    1.478750          0.11877968     0.11271162
# Atl-Spain              2  0.02797398    1.116565                 NaN            NaN
# Indo-Central           5  0.14002953    1.366335          0.09469652     0.10853918
# Indo-East             17  0.40863562    1.617822          0.09000758     0.11720942
# Indo-Mexico            2  0.12287003    1.204318                 NaN            NaN
# Indo-South             9  0.15503412    1.517579          0.11205555     0.11759576
# Indo-Spac              3  0.04752783    1.293906                 NaN            NaN
# Med-East              12  0.55872142    1.329143          0.05562872     0.08082447
# Med-France             2  0.04014526    1.165011                 NaN            NaN
# Med-Italy             20  0.97326310    1.320085          0.01783790     0.08672468
# Preparing SNAPP files ---------------------------------------------------
 
#using gl2snapp command
mden43.gl
mden43.pop.gl
zczc129.gl
zczc129.pop.gl

gl2snapp(mden43.gl, outfile="mden43.snapp.nex", outpath="./", v=5)
gl2snapp(mden43.pop.gl, outfile="mden43.pop.snapp.nex", outpath="./", v=5)
gl2snapp(zczc129.gl, outfile="zczc129.snapp.nex", outpath="./", v=5)
gl2snapp(zczc129.pop.gl, outfile="zczc129.pop.snapp.nex", outpath="./", v=5)



# Preparing Treemix files -------------------------------------------------

#use dartT to convert from genlight to treemix

source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(c("SNPRelate", "qvalue"))
install.packages("dartR")
library(dartR)

gl2treemix(mden43.gl, outfile = "mden43.treemix.gz", outpath = "./", v = 2)
gl2treemix(mden43.pop.gl, outfile="mden43.pop.treemix.gz", outpath="./", v=2)

source("plotting_functions.R")
plot_tree("stem")

##may  need to assign individuals as separate populations to make the tree... need >2 pops and it looks like this will be plotting the populations as a whole, not individuals. 

# Global Paper - Structure Plots ------------------------------------------

#Mden43

library(ggplot2)
library(ggthemes)
library(readxl)
mden43.q.data <- read_excel("Results/global_paper/structure-plots/q-matrix-data.xlsx", sheet = "mden43")
as.factor(mden43.q.data$q.num)
as.factor(zczc129.q.data$q.num)

ggplot(mden43.q.data, aes( x = reorder(filename, -sort), y = q.val, fill=q.num) ) + 
  geom_bar( stat = "identity", position = "stack" ) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Ancestry Coefficient") +
 theme(axis.title.x = element_blank(), axis.text.x  = element_text(angle=90)) +
  theme(legend.position="none") +
  facet_grid(~ocean, scales="free", space="free_x")
#with no sample labels
ggplot(mden43.q.data, aes( x = reorder(filename, -sort), y = q.val, fill=q.num) ) + 
  geom_bar( stat = "identity", position = "stack" , width=1) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Ancestry Coefficient") +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none") +
  facet_grid(~ocean, scales="free", space="free_x")



  ##best two plots for paper
# Install
install.packages("wesanderson")
# Load
library(wesanderson)

zc.facets.ocean<-c("Atlantic", "Indo-Pacific", "Mediterranean")
names(zc.facets.ocean)<-c("Atlantic", "Indopacific", "Mediterranean")
md.facets.ocean<-c("Atlantic", "Indo-Pacific")
names(md.facets.ocean)<-c("Atlantic", "Indopacific")

# 
# dpi=600    #pixels per square inch
# png("mden43.qplot.png", width=6*dpi, height=3*dpi, res=dpi)
# 
#   ggplot(mden43.q.data, aes( x = reorder(filename, -sort), y = q.val, fill=factor(q.num)) ) + 
#     geom_bar( stat = "identity", position = "stack" , width=1) +
#     theme_classic() +
#     scale_y_continuous(expand = c(0,0)) +
#     ylab("Ancestry Coefficient") +
#     theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
#     theme(legend.position="none") +
#     facet_grid(~ocean, scales="free", space="free_x", labeller=labeller(ocean=md.facets.ocean)) +
#     scale_fill_manual(values=my_zissou)
# dev.off()
#   

dpi=600    #pixels per square inch

# #Atlantic Only samples tess plot
# zc.facets.atl<-c("Atl-Bah", "Atl-East", "Atl-Mad", "Atl-Oth", "Atl-Spa")
# names(zc.facets.atl)<-c("Bahamas", "East", "Madeira", "Other", "Spain")
# 
# png("zczc129.atl.qplot.png", width=6*dpi, height=3*dpi, res=dpi)
# ggplot(zczc129.pop.atl.q, aes( x = filename, y = q.val, fill=factor(q.num))) + 
#   geom_bar( stat = "identity", position = "stack", width=1 ) +
#   theme_classic() +
#   scale_y_continuous(expand = c(0,0)) +
#   ylab("Ancestry Coefficient") +
#   theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
#   theme(legend.position="none") +
#   facet_grid(~pop, scales="free", space="free_x", labeller=labeller(pop=zc.facets.atl)) +
#   theme(strip.text.x=element_text(size=4))+
#   scale_fill_manual(values=my_zissou)
# dev.off()
# 
# #Indopacific Only samples tess plot
# zc.facets.indo<-c("Indo-Cen", "Indo-East", "Indo-Mex", "Indo-Sou", "Indo-Spac")
# names(zc.facets.indo)<-c("Central Pacific", "East", "Mexico", "South Indopacific","South Pacific" )
# 
# png("zczc129.indo.qplot.png", width=6*dpi, height=3*dpi, res=dpi)
# ggplot(zczc129.pop.indo.q, aes( x = filename, y = q.val, fill=factor(q.num))) + 
#   geom_bar( stat = "identity", position = "stack", width=1 ) +
#   theme_classic() +
#   scale_y_continuous(expand = c(0,0)) +
#   ylab("Ancestry Coefficient") +
#   theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
#   theme(legend.position="none") +
#   facet_grid(~pop, scales="free", space="free_x", labeller=labeller(pop=zc.facets.indo)) +
#   theme(strip.text.x=element_text(size=4))+
#   scale_fill_manual(values=wes_palette("Zissou1")[c(2,3,1,5)])
# dev.off()
# 
# #Mediterranean samples only tess plot
# zc.facets.med<-c("Med-East", "Med-Fra", "Med-Ita")
# names(zc.facets.med)<-c("East", "France", "Italy")
# png("zczc129.med.qplot.png", width=6*dpi, height=3*dpi, res=dpi)
# ggplot(zczc129.pop.med.q, aes( x = filename, y = q.val, fill=factor(q.num))) + 
#   geom_bar( stat = "identity", position = "stack", width=1 ) +
#   theme_classic() +
#   scale_y_continuous(expand = c(0,0)) +
#   ylab("Ancestry Coefficient") +
#   theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
#   theme(legend.position="none") +
#   facet_grid(~pop, scales="free", space="free_x", labeller=labeller(pop=zc.facets.med)) +
#   theme(strip.text.x=element_text(size=4))+
#   scale_fill_manual(values=my_zissou)
# dev.off()

####structure plots for zczc125 dataset
zczc125.ocean.q.data <- read_excel("Results/global_paper/structure-plots/q-matrix-data.xlsx", sheet = "zczc125.ocean")
zczc125.ocean.q.data

   #pixels per square inch
png("zczc125.ocean.qplot.png", width=4, height=2, unit="in", res=600)
ggplot(zczc125.ocean.q.data, aes( x = reorder(filename, -sort.ocean), y = q.val.ocean, fill=factor(q.num.ocean)) ) + 
  geom_bar( stat = "identity", position = "stack", width=1 ) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Ancestry Coefficient") +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none") +
  facet_grid(~ocean, scales="free", space="free_x", labeller=labeller(ocean=zc.facets.ocean)) +
  scale_fill_manual(values=paired.col[c(6,4,2)])
dev.off()

##zcav populations tess plots
zczc125.pop.q.data <- read_excel("Results/global_paper/structure-plots/q-matrix-data.xlsx", sheet = "zczc125.pops")
zczc125.pop.q.data

#atl
zc.facets1.atl<-as_labeller(c('Atl_AA'="N. Caribbean", 'Atl_AB'="S. Caribbean", 'Atl_AC'="Canary Islands", 'Atl_AD'="NE Atlantic", 'Atl_AE'="Spain"))
zc.facets2.atl<-as_labeller(c('Atl_BA'="Caribbean", 'Atl_BB'="Canary Islands", 'Atl_BC'="NE Atlantic", 'Atl_BD'="Spain"))

zczc125.atl.q.data<-subset(zczc125.pop.q.data, zczc125.pop.q.data$ocean=="Atlantic")
zczc125.atl.q.data

png("zczc125.atl1.qplot.png", width=6, height=3, unit="in", res=600)
ggplot(zczc125.atl.q.data, aes( x = reorder(filename, -sort.pop1), y = q.val.pop1, fill=factor(q.num.pop1))) + 
  geom_bar( stat = "identity", position = "stack", width=1 ) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Ancestry Coefficient") +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none") +
  facet_grid(.~pop1, scales="free", space="free_x", labeller=labeller(pop1=as_labeller(zc.facets1.atl))) +
  theme(strip.text.x = element_text(size = 5.5)) +
  scale_fill_manual(values=brewer.pal(n=4,name="Blues"))
dev.off()

png("zczc125.atl2.qplot.png", width=6, height=3, unit="in", res=600)
ggplot(zczc125.atl.q.data, aes( x = reorder(filename, -sort.pop1), y = q.val.pop1, fill=factor(q.num.pop1))) + 
  geom_bar( stat = "identity", position = "stack", width=1 ) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Ancestry Coefficient") +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none") +
  facet_grid(~pop2, scales="free", space="free_x", labeller=labeller(pop2=as_labeller(zc.facets2.atl))) +
  theme(strip.text.x = element_text(size = 5.5)) +
  scale_fill_manual(values=brewer.pal(n=4, name="Blues"))
dev.off()

#Indo
zczc125.indo.q.data<-subset(zczc125.pop.q.data, zczc125.pop.q.data$ocean=="Indopacific")
zc.facets.indo<-as_labeller(c("Indo_A"="Central", "Indo_B"="South", "Indo_C"="NE Indo-Pacific"))


png("zczc125.indo.qplot.png", width=6, height=3, unit="in", res=600)
ggplot(zczc125.indo.q.data, aes( x = reorder(filename, -sort.pop1), y = q.val.pop1, fill=factor(q.num.pop1))) + geom_bar( stat = "identity", position = "stack", width=1 ) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Ancestry Coefficient") +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none") +
  facet_grid(~pop1, scales="free", space="free_x", labeller=as_labeller(zc.facets.indo))+
  theme(strip.text.x = element_text(size = 7)) +
  scale_fill_manual(values=brewer.pal(n=3, name="Greens"))
dev.off()

#Med
zczc125.med.q.data<-subset(zczc125.pop.q.data, zczc125.pop.q.data$ocean=="Mediterranean")
zc.facets.med<-as_labeller(c("Med_A"="West", "Med_B"="East", "Med_C"="Corfu"))


png("zczc125.med.qplot.png", width=6, height=3, unit="in", res=600)
ggplot(zczc125.med.q.data, aes( x = reorder(filename, -sort.pop1), y = q.val.pop1, fill=factor(q.num.pop1))) + geom_bar( stat = "identity", position = "stack", width=1 ) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Ancestry Coefficient") +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none") +
  facet_grid(~pop1, scales="free", space="free_x", labeller=as_labeller(zc.facets.med)) +
  theme(strip.text.x = element_text(size = 7)) +
  scale_fill_manual(values=rev(brewer.pal(n=12, name="YlOrRd")))
dev.off()

# #Atlantic Only samples tess plot
# zc.facets.atl<-c("Atl-Bah", "Atl-East", "Atl-Mad", "Atl-Oth", "Atl-Spa")
# names(zc.facets.atl)<-c("Bahamas", "East", "Madeira", "Other", "Spain")
# 
# png("zczc129.atl.qplot.png", width=6*dpi, height=3*dpi, res=dpi)
# ggplot(zczc129.pop.atl.q, aes( x = filename, y = q.val, fill=factor(q.num))) + 
#   geom_bar( stat = "identity", position = "stack", width=1 ) +
#   theme_classic() +
#   scale_y_continuous(expand = c(0,0)) +
#   ylab("Ancestry Coefficient") +
#   theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
#   theme(legend.position="none") +
#   facet_grid(~pop, scales="free", space="free_x", labeller=labeller(pop=zc.facets.atl)) +
#   theme(strip.text.x=element_text(size=4))+
#   scale_fill_manual(values=my_zissou)
# dev.off()
# 
# #Indopacific Only samples tess plot
# zc.facets.indo<-c("Indo-Cen", "Indo-East", "Indo-Mex", "Indo-Sou", "Indo-Spac")
# names(zc.facets.indo)<-c("Central Pacific", "East", "Mexico", "South Indopacific","South Pacific" )
# 
# png("zczc129.indo.qplot.png", width=6*dpi, height=3*dpi, res=dpi)
# ggplot(zczc129.pop.indo.q, aes( x = filename, y = q.val, fill=factor(q.num))) + 
#   geom_bar( stat = "identity", position = "stack", width=1 ) +
#   theme_classic() +
#   scale_y_continuous(expand = c(0,0)) +
#   ylab("Ancestry Coefficient") +
#   theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
#   theme(legend.position="none") +
#   facet_grid(~pop, scales="free", space="free_x", labeller=labeller(pop=zc.facets.indo)) +
#   theme(strip.text.x=element_text(size=4))+
#   scale_fill_manual(values=wes_palette("Zissou1")[c(2,3,1,5)])
# dev.off()
# 
# #Mediterranean samples only tess plot
# zc.facets.med<-c("Med-East", "Med-Fra", "Med-Ita")
# names(zc.facets.med)<-c("East", "France", "Italy")
# png("zczc129.med.qplot.png", width=6*dpi, height=3*dpi, res=dpi)
# ggplot(zczc129.pop.med.q, aes( x = filename, y = q.val, fill=factor(q.num))) + 
#   geom_bar( stat = "identity", position = "stack", width=1 ) +
#   theme_classic() +
#   scale_y_continuous(expand = c(0,0)) +
#   ylab("Ancestry Coefficient") +
#   theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
#   theme(legend.position="none") +
#   facet_grid(~pop, scales="free", space="free_x", labeller=labeller(pop=zc.facets.med)) +
#   theme(strip.text.x=element_text(size=4))+
#   scale_fill_manual(values=my_zissou)
# dev.off()


##M. densirostris tess plots for populations within ocean basins

mden43.q.data <- read_excel("q-matrix-data.xlsx", sheet = "mden43")
mden43.q.data

mden43.q.atl.data<-subset(mden43.q.data, mden43.q.data$ocean=="Atlantic")
mden43.q.atl.data$pop

mden43.q.indo.data<-subset(mden43.q.data, mden43.q.data$ocean=="Indopacific")
mden43.q.indo.data$pop

#Mden Atlantic samples only tess plot
md.facet.atl<-as_labeller(c("Bahamas", "East", "Other"))
names(md.facet.atl)<-c("Bahamas", "East", "Other")
zc.facets.med<-as_labeller(c("Med_A"="West", "Med_B"="East", "Med_C"="Corfu"))
md.facets.ocean

png("mden43.ocean.qplot.png", width=4, height=2, unit="in", res=600)
ggplot(mden43.q.data, aes( x = reorder(filename, -sort), y = q.val.ocean, fill=factor(q.num.ocean)) ) + 
  geom_bar( stat = "identity", position = "stack", width=1 ) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Ancestry Coefficient") +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none") +
  facet_grid(~ocean, scales="free", space="free_x", labeller=labeller(ocean=zc.facets.ocean)) +
  scale_fill_manual(values=paired.col[c(4,2)])
dev.off()


png("mden43.atl.qplot.png", width=6*dpi, height=3*dpi, res=dpi)
ggplot(mden43.q.atl.data, aes( x = filename, y = q.val.pop, fill=factor(q.num.pop))) + 
  geom_bar( stat = "identity", position = "stack", width=1 ) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Ancestry Coefficient") +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none") +
  facet_grid(~pop, scales="free", space="free_x") +
  theme(strip.text.x=element_text(size=10))+
  scale_fill_manual(values=rev(brewer.pal(n=4, name="Blues")))
dev.off()

#Mden Indopacific samples only tess plot

png("mden43.indo.qplot.png", width=6*dpi, height=3*dpi, res=dpi)
ggplot(mden43.q.indo.data, aes( x = filename, y = q.val.pop, fill=factor(q.num.pop))) + 
  geom_bar( stat = "identity", position = "stack", width=1 ) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Ancestry Coefficient") +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none") +
  facet_grid(~pop, scales="free", space="free_x") +
  theme(strip.text.x=element_text(size=10))+
  scale_fill_manual(values=rev(brewer.pal(n=3, name="Greens")))
dev.off()



# Other Random Test Code --------------------------------------------------
####Correlation Plots with Fst Data
install.packages("corrplot")
library(corrplot)


mden43.fst.mat <- read_excel("Results/global_paper/structure-plots/q-matrix-data.xlsx", sheet = "Sheet5")
dim(mden43.fst.mat)
rownames(mden43.fst.mat)
#delete first blank column
mden43.fst.mat<-mden43.fst.mat[1:8, 2:9]
mden43.fst.mat

rownames(mden43.fst.mat)<-colnames(mden43.fst.mat)
rownames(mden43.fst.mat)
colnames(mden43.fst.mat)

mden43.fst.mat
as.numeric(mden43.fst.mat)

#Corrplots

#convert to matrix
mden43.fst.mat<-data.matrix(mden43.fst.mat)
mden43.fst.mat

library(reshape2)
melt_mden<-melt(mden43.fst.mat)
melt_mden

library(ggplot2)
#All the data together
ggplot(data=melt_mden, aes(x=Var1, y=Var2, fill=value)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + 
  geom_tile(color="white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Fst")

#only fst values
lower_mden43_fst<-function(mden43.fst.mat){
  mden43.fst.mat[upper.tri(mden43.fst.mat)]<-NA
  return(mden43.fst.mat)
}
low_mden43_fst<-lower_mden43_fst(mden43.fst.mat)
low_mden43_fst

melt_low_mden43<-melt(low_mden43_fst)
melt_low_mden43<-read.csv("Results/global_paper/fst/mden.csv", header=TRUE, sep=",", )
melt_low_mden43

##Only with Fst values
pal<-wes_palette("Zissou1", 100, type="continuous")
ggplot(data=melt_low_mden43, aes(melt_mden$Var1, melt_mden$Var2, fill=value)) +
  geom_tile(colour="white") +
  scale_fill_gradientn(colours=pal, limit=c(0,.20), space="Lab", name="Fst Value", na.value="grey") +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=45, vjust=1, size=12, hjust=1), axis.text.y=element_text(size=12)) + 
  coord_fixed()

library(gplots)
heatmap.2(low_mden43_fst,na.rm=TRUE, dendrogram = NULL)
heatmap(mden43.fst.mat)

library(superheat)
superheat(low_mden43_fst, heat.na.col="white", bottom.label.text.angle = 45, scale=TRUE, legend.breaks=c(.05, .1, .15))

indNames(mden43.gl)

# DAPC 

mden43.gl
mden43.pop.gl
zczc127.gl
zczc127.pop.gl


dapc.zczc127.pop<-dapc(zczc127.pop.gl, zczc127.pop.gl$pop) #120 PCs, 6 discriminants, all individuals within a population are clustered and only the group name is there, not individual dots
dapc.zczc127.pop
scatter(dapc.zczc127.pop)

dapc2.zczc127.pop<-dapc(zczc127.pop.gl, zczc127.pop.gl$pop) #50 PCs, 7 discriminants, Not good. All populations except Spain are clustered tightly!
dapc2.zczc127.pop
scatter(dapc2.zczc127.pop)

dapc3.zczc127.pop<-dapc(zczc127.pop.gl, zczc127.pop.gl$pop) #120 PCs, 13 discriminants, Not good. All populations except Spain are clustered tightly!
dapc3.zczc127.pop
scatter(dapc3.zczc127.pop)

dapc4.zczc127.pop<-dapc(zczc127.pop.gl, zczc127.pop.gl$pop) #10, 2 - Spain is still so far away and the rest are all clustered super close. Will try to remove this pop and re-do
dapc4.zczc127.pop
scatter(dapc4.zczc127.pop)

#removing spanish population


# Wes anderson colour pallets ------------------------------------------

library(wesanderson)
wes_palette("Zissou1")
wes_palette("Darjeeling1")

my_zissou<-wes_palette("Zissou1")[c(2,3,5,1,4)]
my_zissou2<-wes_palette("Zissou1")[c(5,3,1)]
my_zissou3<-wes_palette("Zissou1")[c(3,5)]
wes_palette("Zissou1")[c(2,3,5,1,4)]

# diveRsity- Fst and confidence intervales -----------------------------------------
#Confidence Intervals for Fst values based on Ocean basins
zczc125.diffCalc<-diffCalc(infile="./zczc125.genepop.txt", outfile="./zczc125.fststats.noloc.ind.txt", fst = TRUE, pairwise=TRUE, para=TRUE, bs_locus = TRUE, boots=100, bs_pairwise=TRUE)
tail(zczc125.diffCalc$std_stats)
names(zczc125.diffCalc)

#diversity indices
install.packages("xlsx")
install.packages("sendplot")
install.packages("doParallel")
install.packages("parallel")
install.packages("foreach")
install.packages("iterators")
install.packages("plotrix")
library(xlsx)
library(sendplot)
library(doParallel)
library(parallel)
library(foreach)
library(iterators)
library(plotrix)
mden43.ho.ci<-basicStats(infile="./mden43.genepop.txt", outfile="./mden43.out2", ho_ci=TRUE, ho_boots=100, ho_alpha=)
names(mden43.ho.ci)
# [1] "ar"        "lps"       "obs_het"   "exp_het"   "uexp_het" 
# [6] "ho"       "hwe_llr_p" "hwe_hom"   "hwe_het"   "main_tab"
head(mden43.ho.ci$ho)
tail(mden43.ho.ci$ho)
mden43.ho.ci[ho]
mden43.ho.ci[["ho"]]

mden_data<-readGenepop(infile="./m=den43.genepop.txt", bootstrap=F)
div_basic<-divBasic(infile="./zczc125_genepop.txt", outfile="./out", bootstraps=100)

#Confidence Intervals for Fst values based on populations 1
diffCalc(infile="./zczc125.pop1.genepop.txt", outfile="./zczc125.fststats.pop1.txt", fst = TRUE, pairwise=TRUE, bs_locus=TRUE, para=TRUE, boots=100, bs_pairwise=TRUE)

#Confidence Intervals for Fst values based on populations 2
diffCalc(infile="./zczc125.pop2.genepop.txt", outfile="./zczc125.fststats.pop2.txt", fst = TRUE, pairwise=TRUE, bs_locus=TRUE, para=TRUE, boots=100, bs_pairwise=TRUE)
  
#Emma's code that worked
dif.diversity<-diffCalc(infile = "/Users/ecar026/Dropbox/Bioinformatics/Cluster_UOA/Output_high_quality/withrep/SRWNoRep_populations.snps.genepop", outfile = "fst.calc.pairwise", fst = TRUE, pairwise = TRUE, 
                        bs_locus = TRUE, para = TRUE, boots = 999, bs_pairwise = TRUE)

# ###Dont know why, but this is now not working for Mden data. get error message below. will try in another package. 
# Error in apply(pwDLoc, 1, function(x) { : 
#     dim(X) must have a positive length
# 
# 
# #Confidence intervals for Mden Fst values based on ocean basin level
# mden43.dif.fst<-diffCalc(infile="./mden43.genepop.txt", outfile="mden43.fstats.ocean.txt", fst=TRUE, pairwise=TRUE, bs_locus=TRUE, para=TRUE, boots=100, bs_pairwise = TRUE)
# 
# #Confidence intervals for Mden Fst values based on population level
# diffCalc(infile="~/Dropbox/Phd/Bioinformatics/bw_ddrad/Analysis Files/genepop/mden43/mden43.pop.genepop.txt", outfile="./mden43.fstats.pop.txt", fst=TRUE, pairwise=TRUE, bs_locus=TRUE, para=TRUE, boots=100, bs_pairwise = TRUE)

mden43.basic<-divBasic(infile="./mden43/mden43.pop.genepop.txt", outfile=NULL, bootstraps=50)

install.packages("diveRsity")
library(diveRsity)
data("Test_data")

basicRes <- divBasic(infile = Test_data, outfile = NULL, gp = 2, 
                     bootstraps = 1000)

zcav.div<-divBasic(infile="./zczc125_genepop.txt", outfile=NULL, gp=3, bootstraps=100)
zcav.basic<-basicStats(infile="./zczc125_genepop.txt", outfile="zczc125_basic_stats", fis_ci=TRUE, fis_boots=100, fis_alpha=0.05)

mden_indo_div<-read.table("./mden_indo_div.txt", header=TRUE)
mden_indo_div[1:5, 1:5]
mden_indo_div<-t(mden_indo_div)
colnames(mden_indo_div)<-mden_indo_div[1,]
rownames(mden_indo_div)

mden_indo_div<-mden_indo_div[-1,]
mden_indo_div<-as.data.frame(mden_indo_div)
mden_indo_div$fis<-as.numeric(mden_indo_div$fis)
mden_indo_div$fis_lo<-as.numeric(mden_indo_div$fis_lo)
mden_indo_div$fis_hi<-as.numeric(mden_indo_div$fis_hi)




summary(mden_indo_div$fis)
su

# Demerlate ---------------------------------------------------------------
library(Demerelate)

mden43.gl
mden43.pop.gl
zczc125.gl
zczc125.pop1.gl
zczc125.pop2.gl


#Demerelate will only work if there are no missing data....
#Will Only work if you have NO NAs at all. So must filter the dataset first!
md.gl<-gl.compliance.check(mden43.gl)
md.gl<-gl.filter.callrate(md.gl, method="loc", threshold=1)
md.gl # 43 genotypes,  9,144 binary SNPs, size: 1.7 Mb, 0 (0 %) missing data
md.gl<-gl.filter.monomorphs(md.gl)
md.gl<-gl.recalc.metrics(md.gl)

md.pop.gl<-gl.compliance.check(mden43.pop.gl)
md.pop.gl<-gl.filter.callrate(md.pop.gl, method="loc", threshold=1)
md.pop.gl # 43 genotypes,  9,144 binary SNPs, size: 1.7 Mb, 0 (0 %) missing data
md.pop.gl<-gl.filter.monomorphs(md.pop.gl)
md.pop.gl<-gl.recalc.metrics(md.pop.gl)

zc.gl<-gl.compliance.check(zczc125.gl)
zc.gl<-gl.filter.callrate(zc.gl, method="loc", threshold=1)
zc.gl # 125 genotypes,  6,888 binary SNPs, size: 3.3 Mb, 0 (0 %) missing data
zc.gl<-gl.filter.monomorphs(zc.gl)
zc.gl<-gl.recalc.metrics(zc.gl)
zc.gl<-gl.compliance.check(zc.gl)
ploidy(zc.gl)<-2
zc.gl@ploidy

zc.pop1.gl<-gl.compliance.check(zczc125.pop1.gl)
zc.pop1.gl<-gl.filter.callrate(zc.pop1.gl, method="loc", threshold=1)
zc.pop1.gl # 125 genotypes,  6,888 binary SNPs, size: 3.3 Mb, 0 (0 %) missing data
zc.pop1.gl<-gl.filter.monomorphs(zc.pop1.gl)
zc.pop1.gl<-gl.recalc.metrics(zc.pop1.gl)
ploidy(zc.pop1.gl)<-2

zc.pop2.gl<-gl.compliance.check(zczc125.pop2.gl)
zc.pop2.gl<-gl.filter.callrate(zc.pop2.gl, method="loc", threshold=1)
zc.pop2.gl # 125 genotypes,  6,888 binary SNPs, size: 3.3 Mb, 0 (0 %) missing data
zc.pop2.gl<-gl.filter.monomorphs(zc.pop2.gl)
zc.pop2.gl<-gl.recalc.metrics(zc.pop2.gl)
ploidy(zc.pop2.gl)<-2

#example data format: 
#Sample-ID      Population locus.1.a locus.1.b locus.2.a locus.2.b ...
#Ind.Norway.01  Norway      001       002       001       002 ...

#convert data to correct format using dart

md.dem<-gl2demerelate(md.gl)
md.pop.dem<-gl2demerelate(md.pop.gl)
zc.dem<-gl2demerelate(zc.gl)
zc.pop1.dem<-gl2demerelate(zc.pop1.gl)
zc.pop2.dem<-gl2demerelate(zc.pop2.gl)

md.fis<-F.stat(md.dem, object=TRUE, iteration=10, directory.name="./md_fis/", out.name="md_fis")
md.pop.fis<-F.stat(md.pop.dem, object=TRUE, iteration=10, directory.name="./md_fis", out.name="md_pop_fis")

zc.fis<-F.stat(zc.dem, object=TRUE, iteration=10, directory.name="./zc_fis/", out.name="zc_fis")
zc.pop1.fis<-F.stat(zc.pop1.dem, object=TRUE, iteration=10, directory.name="./zc_fis", out.name="zc_pop1_fis")
zc.pop2.fis<-F.stat(zc.pop2.dem, object=TRUE, iteration=10, directory.name="./zc_fis", out.name="zc_pop2_fis")



mden43.dem.fis<-Fis.calc(mden43.dem, iteration=100,number.loci =13989, directory.name = "./", out.name="mden43.dem.fis")
mden43.deme<-F.stat(mden43.dem, object=TRUE,iteration=100, directory.name="./", out.name="mden43.deme")

mden43.seppop.gl
mden43.seppop.atl.gl<-mden43.seppop.gl$Atlantic
mden43.seppop.atl.gl
mden43.seppop.indo.gl<-mden43.seppop.gl$Indopacific
mden43.seppop.indo.gl

mden43.atl.dem<-gl2demerelate(mden43.seppop.atl.gl)
mden43.atl.dem
mden43.indo.dem<-gl2demerelate(mden43.seppop.indo.gl)
mden43.indo.dem

zczc129.dem.stats<-Demerelate(zczc129.dem, object=TRUE, ho=TRUE)
zczc129.dem.stats

zczc125.dem<-gl2demerelate(zczc125.gl)

F.stat(zczc125.dem, object=TRUE, iteration=100, directory.name = "./", out.name = "zczc125.fstat")
zczc125.fstat

#Atlantic Weighted mean ho value: 0.128
#Indopacific Weighted mean ho value: 0.116
#Mediterranean Weighted mean ho value: 0.132

F.stat(mden43.dem, object=TRUE, iteration=100, directory.name = "./", out.name = "mden43.fstat.100")
mden43.fstat

F.stat(mden43.atl.dem, object=TRUE, iteration=100, directory.name="./", out.name="mden43.atl.100")
F.stat(mden43.indo.dem, object=TRUE, iteration=10, directory.name = "./", out.name = "mden43.indo.10")

mden43.indo.ho<-ho(mden43.atl.dem, allele.column = 1)

mden43.indo.dem[1:3,1:5]

mden43.indo.dem[is.na(mden43.indo.dem)]<-0
mden43.atl.dem[is.na(mden43.atl.dem)]<-0
mden43.dem[is.na(mden43.dem)]<-0

zczc125.dem[1:5,1:5]

zczc125.dem[is.na(zczc125.dem)]<-0

F.stat(zczc125.dem, object=TRUE, iteration=100, directory.name="./", out.name="zczc125.100")


##bootstrap ho results for CI per ocean basin
#import files of wc ho values calculated by demerelate
zczc.indo.ho<-read.table("./zczc_indo_ho.txt", header=T)
zczc.atl.ho<-read.table("./zczc_atl_ho.txt", header=T)
zczc.med.ho<-read.table("./zczc_med_ho.txt", header=T)
#check files
head(zczc.indo.ho)
head(zczc.atl.ho)
head(zczc.med.ho)

#function to calculate mean from ho values
indo.ho.mean = function(x, indices){
  return(mean(x[indices]))
}
indo.ho.boot<-boot(zczc.indo.ho$ho_wc, indo.ho.mean, R=10000)
indo.ho.boot.ci<-boot.ci(indo.ho.boot, conf=0.95, type="basic")
indo.ho.boot.ci
#Intervals : 
#Level      Basic         
#95%   ( 0.3760,  0.3863 ) 

atl.ho.mean = function(x, indices){return(mean(x[indices]))}
atl.ho.boot<-boot(zczc.atl.ho$ho_wc, atl.ho.mean, R=10000)
atl.ho.boot.ci<-boot.ci(atl.ho.boot, conf=0.95, type="basic")
atl.ho.boot.ci
#Intervals : 
#  Level      Basic         
#95%   ( 0.3114,  0.3204 ) 

med.ho.mean = function(x, indices){return(mean(x[indices]))}
med.ho.boot<-boot(zczc.med.ho$ho_wc, med.ho.mean, R=10000)
med.ho.boot.ci<-boot.ci(med.ho.boot, conf=0.95, type="basic")
med.ho.boot.ci
# Level      Basic         
# 95%   ( 0.7023,  0.7124 ) 
# SNPRelate ---------------------------------------------------------------
library(SNPRelate)
library(dartR)
library(radiator)

zczc125.gds<-gl2gds(zczc125.gl)
snpgds


# assigner ----------------------------------------------------------------
library(assigner)
library(dplyr)
library(tidyr)
library(purrr)
library(parallel)
library(SNPRelate)

mden43.gl
zczc125.gl

mden.tidy<-tidy_genlight(mden43.gl)
zcav.tidy<-tidy_genlight(zczc125.gl)

mden.tidy
zcav.tidy

mden.tidy.fis<-fst_WC84()

# Global paper- diversity stats -------------------------------------------
##From my github!

####Codes for measuring Diversity Statistics: Ho, Hs, ho, F and Tajima's D
#Packages used:
#Hierfstat or DartR: Basic stats. Wrapper function in DartR does same calculation with a genlight file as Hierfstat does. Calculates Ho, Hs and ho using Nei 1987 method.
#Poppr: Various summary statistics and diversity indices by population
#Vcftools: Can calculate number of observed loci and the F inbreeding coefficient (method of moments) for each individual. Can also calculate site depth and Tajima's D for each population. 
#Genepop: Calculates ho (Weir and cockerham), 1-Qinter, 1-Qintra for each population (and locus)
#Adegenet: Can calulate Hs (Expected heterozygosity/Gene diversity)

####Hierfstat####
library(hierfstat)

#convert file from genind (Adegenet) to hierfstat format. This file has ocean basin level populations: Atlantic, Indopacific, Mediterranean
zczc125.hfstat<-genind2hierfstat(zczc125.gi)
zczc125.atl.hfstat<-as.data.frame(subset(zczc125.hfstat, zczc125.hfstat$pop=="Atlantic"))
zczc125.indo.hfstat<-as.data.frame(subset(zczc125.hfstat, zczc125.hfstat$pop=="Indopacific"))
zczc125.med.hfstat<-as.data.frame(subset(zczc125.hfstat, zczc125.hfstat$pop=="Mediterranean"))
zczc125.atl.hfstat[,1]<-factor(zczc125.atl.hfstat[,1])
zczc125.indo.hfstat[,1]<-factor(zczc125.indo.hfstat[,1])
zczc125.med.hfstat[,1]<-factor(zczc125.med.hfstat[,1])

#convert mden genind to hfstat with two ocean populations: Atlantic and Indopacific
mden43.hfstat<-genind2hierfstat(mden43.gi)
mden43.atl.hfstat<-as.data.frame(subset(mden43.hfstat, mden43.hfstat$pop=="Atlantic"))
mden43.indo.hfstat<-as.data.frame(subset(mden43.hfstat, mden43.hfstat$pop=="Indopacific"))
mden43.atl.hfstat[,1]<-factor(mden43.atl.hfstat[,1])
mden43.indo.hfstat[,1]<-factor(mden43.indo.hfstat[,1])
  
#calculate basic stats (Ho, Hs and ho using Nei 1987) for each of the ocean basins
zczc125.hfstat.basicstats<-basic.stats(zczc125.hfstat, diploid=TRUE)
mden43.hfstat.basicstats<-basic.stats(mden43.hfstat, diploid=TRUE)

#view results. This will give basic summary statistics for each measure for each of the three defined ocean basins in the genind file. 

summary(zczc125.hfstat.basicstats$Ho)
summary(zczc125.hfstat.basicstats$Hs)
summary(zczc125.hfstat.basicstats$ho)

summary(mden43.hfstat.basicstats$Ho)
summary(mden43.hfstat.basicstats$Hs)
summary(mden43.hfstat.basicstats$ho)

#Calculate boostraps around fis values
# zczc125.bootppfis<-boot.ppfis(zczc125.hfstat, nboot=100)
# zczc125.bootppfis$fis.ci
# # ll     hl
# # 1 0.1322 0.1397
# # 2 0.1225 0.1311
# # 3 0.1405 0.1510

#For the weir and cockerham fis calculations to work, need to use a genind file and need to use seppop to generate separate genind files by population
zczc125.sep.gi<-seppop(zczc125.gi)
zczc125.med.gi<-zczc125.sep.gi$Mediterranean
zczc125.atl.gi<-zczc125.sep.gi$Atlantic
zczc125.indo.gi<-zczc125.sep.gi$Indopacific

mden43.sep.gi<-seppop(mden43.gi)
mden43.atl.gi<-mden43.sep.gi$Atlantic
mden43.indo.gi<-mden43.sep.gi$Indopacific

#Caclulate fis point estimates:
zczc125.med.wc<-wc(zczc125.med.gi)
zczc125.atl.wc<-wc(zczc125.atl.gi)
zczc125.indo.wc<-wc(zczc125.indo.gi)
mden43.atl.wc <-wc(mden43.atl.gi)
mden43.indo.wc<-wc(mden43.indo.gi)

zczc125.med.wc  #$FIS [1] 0.1450542
zczc125.atl.wc  #$FIS [1] 0.1362893
zczc125.indo.wc #$FIS [1] 0.1269361
mden43.atl.wc   #$FIS [1] 0.07763681
mden43.indo.wc  #$FIS [1] 0.09777664

#Calculate fis 95% confidence intervals
zczc125.med.boot.fis<-boot.ppfis(zczc125.med.gi, nboot=1000)
zczc125.atl.boot.fis<-boot.ppfis(zczc125.atl.gi, nboot=1000)
zczc125.indo.boot.fis<-boot.ppfis(zczc125.indo.gi, nboot=1000)
mden43.atl.boot.fis<-boot.ppfis(mden43.atl.gi, nboots=1000)
mden43.indo.boot.fis<-boot.ppfis(mden43.indo.gi, nboots=1000)

zczc125.med.boot.fis  #0.1386 0.1516
zczc125.atl.boot.fis  #0.1321 0.1404
zczc125.indo.boot.fis #0.1223 0.1315
mden43.atl.boot.fis   #0.0695 0.0853
mden43.indo.boot.fis  #0.0899 0.105

mden43.pop.gi
zczc125.pop1.gl
zczc125.pop2.gl

mden43.pop.sep.gi<-seppop(mden43.pop.gi)
zczc125.pop1.sep.gi<-seppop(zczc125.pop1.gi)
zczc125.pop2.sep.gi<-seppop(zczc125.pop2.gi)

zc.atl.ncarib.wc<-wc(zczc125.pop1.sep.gi$Atl_AA)  #0.1313147
zc.atl.scarib.wc<-wc(zczc125.pop1.sep.gi$Atl_AB)  #0.2101498
zc.atl.carib.wc<- wc(zczc125.pop2.sep.gi$Atl_BA)  #0.1565372
zc.atl.canis.wc<- wc(zczc125.pop1.sep.gi$Atl_AC)  #0.1223991
zc.atl.east.wc<-  wc(zczc125.pop1.sep.gi$Atl_AD)  #0.1217424
zc.ind.cent.wc<-  wc(zczc125.pop1.sep.gi$Indo_A)  #0.1524863
zc.ind.sou.wc<-   wc(zczc125.pop1.sep.gi$Indo_B)  #0.1269939
zc.ind.ne.wc<-    wc(zczc125.pop1.sep.gi$Indo_C)  #0.1074694
zc.med.west.wc<-  wc(zczc125.pop1.sep.gi$Med_A)   #0.08550426
zc.med.east.wc<-  wc(zczc125.pop1.sep.gi$Med_B)   #0.09467633
zc.med.corf.wc<-  wc(zczc125.pop1.sep.gi$Med_C)   #0.2051271

md.atl.oth<-wc(mden43.pop.sep.gi$`Atl-Other`)   #0.06433392
md.atl.bah<-wc(mden43.pop.sep.gi$`Atl-Bahamas`) #0.07250989
md.atl.east<-wc(mden43.pop.sep.gi$`Atl-East`)   #0.0568955
md.ind.afr<-wc(mden43.pop.sep.gi$`Indo-Africa`) #0.1238679
md.ind.haw<-wc(mden43.pop.sep.gi$`Indo-Hawaii`) #0.06385982
md.ind.sou<-wc(mden43.pop.sep.gi$`Indo-South`)  #0.07913545

boot.ppfis(dat = zczc125.pop1.gi, nboot = 1000)

#$fis.ci
#ll      hl
#1   0.1241  0.1383
#2   0.1993  0.2207
#3   0.1162  0.1289
#4   0.1160  0.1269
#5  -1.0455 -1.0215
#6   0.1416  0.1640
#7   0.1199  0.1338
#8   0.1014  0.1132
#9   0.0775  0.0930
#10  0.0831  0.1057
#11  0.1930  0.2176

boot.ppfis(dat = zczc125.pop2.gi, nboot = 1000)

#$fis.ci
#ll      hl
#1   0.1505  0.1629
#2   0.1167  0.1289
#3   0.1162  0.1268
#4  -1.0454 -1.0223
#5   0.1421  0.1637
#6   0.1202  0.1341
#7   0.1016  0.1131
#8   0.0779  0.0933
#9   0.0839  0.1059
#10  0.1925  0.2182

boot.ppfis(dat = mden43.pop.gi, nboot = 1000)

#$fis.ci
#ll     hl
#1 0.0586 0.0871
#2 0.0478 0.0658
#3 0.0478 0.0800
#4 0.1111 0.1369
#5 0.0515 0.0764
#6 0.0596 0.0947

# zczc125.hfstat.sorted<-zczc125.hfstat[order(zczc125.hfstat$pop),]
# dim(zczc125.hfstat.sorted)
# zczc125.hfstat.sorted[,1]
# zczc125.med.hfstat[1:5,1:5]
# 
# zczc125.sorted.bootppfis<-boot.ppfis(zczc125.hfstat.sorted)
# zczc125.sorted.bootppfis$fis.ci
# zczc125.sorted.bootppfis$fis.ci
# # ll     hl
# # 1 0.1330 0.1402
# # 2 0.1217 0.1318
# # 3 0.1388 0.1526
# 
# zczc125.bootppfis$fis.ci
# # ll     hl
# # 1 0.1322 0.1397
# # 2 0.1225 0.1311
# # 3 0.1405 0.1510
# 
# 
# # zczc125.zero.hfstat<-zczc125.hfstat.sorted
# # zczc125.zero.hfstat[is.na(zczc125.zero.hfstat)]<-0
# # zczc125.zero.bootppfis<-boot.ppfis(zczc125.zero.hfstat)
# # zczc125.zero.bootppfis$fis.ci
# # ll     hl
# # 1 0.3350 0.3430
# # 2 0.3314 0.3416
# # 3 0.5135 0.5256
# 
# mden43.bootppfis<-boot.ppfis(mden43.hfstat, nboot=100)
# mden43.bootppfis$fis.ci
# # ll     hl
# # 1 0.0718 0.0848
# # 2 0.0892 0.1048
# 
# zczc125.pop1.hfstat<-genind2hierfstat(zczc125.pop1.gi)
# zczc125.pop2.hfstat<-genind2hierfstat(zczc125.pop2.gi)
# mden43.pop.hfstat<-genind2hierfstat(mden43.pop.gi)
# 
# dim(zczc125.pop1.hfstat)
# dim(zczc125.pop2.hfstat)
# dim(mden43.pop.hfstat)
# 
# zczc125.pop1.ppfis<-boot.ppfis(zczc125.pop1.hfstat, nboot=100)
# # ll      hl
# # 1   0.1255  0.1383
# # 2   0.2009  0.2214
# # 3   0.1154  0.1281
# # 4   0.1172  0.1265
# # 5  -1.0468 -1.0217
# # 6   0.1424  0.1628
# # 7   0.1214  0.1336
# # 8   0.1022  0.1120
# # 9   0.0778  0.0932
# # 10  0.0840  0.1042
# # 11  0.1946  0.2153
# zczc125.pop2.ppfis<-boot.ppfis(zczc125.pop2.hfstat, nboot=100)
# # ll      hl
# # 1   0.1510  0.1622
# # 2   0.1155  0.1289
# # 3   0.1159  0.1256
# # 4  -1.0443 -1.0241
# # 5   0.1415  0.1626
# # 6   0.1221  0.1325
# # 7   0.1018  0.1130
# # 8   0.0783  0.0946
# # 9   0.0836  0.1046
# # 10  0.1922  0.2173
# mden43.pop.ppfis<-boot.ppfis(mden43.pop.hfstat, nboot=100)
# # ll     hl
# # 1 0.0595 0.0865
# # 2 0.0495 0.0664
# # 3 0.0511 0.0802
# # 4 0.1114 0.1385
# # 5 0.0525 0.0756
# # 6 0.0574 0.0927
# # 7   -Inf   -Inf
# 
# zczc125.atl.hfstat[1:5,1:5]
# wc(zczc125.atl.hfstat,diploid=TRUE)

#calculating dA using genet.dist
#“Da” This is Nei's et al genetic distance (eqn 7), performing nearly as well as “Dch”
dim(zczc125.hfstat)
zczc125.hfstat[,1]
dim(mden43.hfstat)
mden43.hfstat[,1]

zczc125.da<-genet.dist(zczc125.hfstat, method="Da")
zczc125.da
#     1           2
# 2 0.009789475            
# 3 0.034973667   0.038426921
mden43.da<-genet.dist(mden43.hfstat, method="Da")
mden43.da
#   1
# 2 0.03394864

zczc125.pop1.da<-genet.dist(zczc125.pop1.hfstat, method="Da")
zczc125.pop1.da
#Atl_AA Atl_AB Atl_AC Atl_AD Atl_AE Indo_A Indo_B Indo_C Med_A Med_B Med_C
#     1          2          3          4          5          6          7          8          9         10
# 2  0.02388664                                                                                                   
# 3  0.01582804 0.02346080                                                                                        
# 4  0.01516022 0.02257742 0.01138916                                                                             
# 5  0.05900240 0.06296778 0.05742398 0.05569772                                                                  
# 6  0.02935771 0.03441095 0.02690840 0.02597298 0.06535603                                                       
# 7  0.01992740 0.02625256 0.01676673 0.01508555 0.05916122 0.02385625                                            
# 8  0.02068594 0.02714409 0.01728869 0.01572295 0.06015079 0.02143742 0.01369381                                 
# 9  0.04660761 0.03655587 0.04363308 0.04242428 0.07767316 0.05129282 0.04421051 0.04480601                      
# 10 0.05240648 0.03864387 0.04947844 0.04827245 0.08190408 0.05601666 0.04993511 0.05060307 0.01397572           
# 11 0.04348304 0.03369204 0.04114356 0.04025574 0.07608629 0.04927431 0.04246346 0.04329398 0.01890403 0.01322522

zczc125.pop2.da<-genet.dist(zczc125.pop2.hfstat, method="Da")
zczc125.pop2.da
#Atl_BA Atl_BB Atl_BC Atl_BD Indo_A Indo_B Indo_C Med_A Med_B Med_C
#     1          2          3          4          5          6          7          8          9
# 2  0.01372393                                                                                        
# 3  0.01286668 0.01138916                                                                             
# 4  0.05783274 0.05742398 0.05569772                                                                  
# 5  0.02769556 0.02690840 0.02597298 0.06535603                                                       
# 6  0.01784729 0.01676673 0.01508555 0.05916122 0.02385625                                            
# 7  0.01862719 0.01728869 0.01572295 0.06015079 0.02143742 0.01369381                                 
# 8  0.03962202 0.04363308 0.04242428 0.07767316 0.05129282 0.04421051 0.04480601                      
# 9  0.04439439 0.04947844 0.04827245 0.08190408 0.05601666 0.04993511 0.05060307 0.01397572           
# 10 0.03653539 0.04114356 0.04025574 0.07608629 0.04927431 0.04246346 0.04329398 0.01890403 0.01322522

mden43.pop.da<-genet.dist(mden43.pop.hfstat, method="Da")
mden43.pop.da
#Atl-Bahamas Atl-East Atl-Other Indo-Africa Indo-Hawaii Indo-South NA
#   1          2          3          4          5          6
# 2 0.01654475                                                       
# 3 0.01925939 0.01663157                                            
# 4 0.03887539 0.03569815 0.03852625                                 
# 5 0.04971283 0.04627049 0.04896594 0.03419427                      
# 6 0.05341589 0.05067143 0.05291720 0.03906672 0.03558537           
# 7 0.07395675 0.07121956 0.07301784 0.06133184 0.05858583 0.06340805

####dartR####
#producing same summary stats using dartR with a genlight file instead of genind converted to hierfstat file
library(dartR)

#run basic stats function (wrapper for hierfstats) for each genlight file (one ocean basin level and the two potential population level files)
zczc125.gl.basicstats<-gl.basic.stats(zczc125.gl)
zczc125.pop1.basicstats<-gl.basic.stats(zczc125.pop1.gl)
zczc125.pop2.basicstats<-gl.basic.stats(zczc125.pop2.gl)

mden43.gl.basicstats<-gl.basic.stats(mden43.gl)
mden43.pop.basicstats<-gl.basic.stats(mden43.pop.gl)

#Summary statistics for each metric for ocean basin pops
summary(zczc125.gl.basicstats$Ho)
summary(zczc125.gl.basicstats$Hs)
summary(zczc125.gl.basicstats$ho)

summary(mden43.gl.basicstats$Ho)
summary(mden43.gl.basicstats$Hs)
summary(mden43.gl.basicstats$ho)

#Summary statistics for each metric for population groups 1
summary(zczc125.pop1.basicstats$Ho)
summary(zczc125.pop1.basicstats$Hs)
summary(zczc125.pop1.basicstats$ho)

#Summary statistics for each metric for population gorups 2
summary(zczc125.pop2.basicstats$Ho)
summary(zczc125.pop2.basicstats$Hs)
summary(zczc125.pop2.basicstats$ho)

#Summary statistics for each metric for Mden populations
summary(mden43.pop.basicstats$Ho)
summary(mden43.pop.basicstats$Hs)
summary(mden43.pop.basicstats$ho)



####Poppr####
library(poppr)

#calculate summary statistics for ocean basins and populations using genind files (adegenet)
zczc125.poppr<-poppr(zczc125.gi)
zczc125.pop1.poppr<-poppr(zczc125.pop1.gi)
zczc125.pop2.poppr<-poppr(zczc125.pop2.gi)

mden43.poppr<-poppr(mden43.gi)
mden43.pop.poppr<-poppr(mden43.pop.gi)

####vcftools####
#vcftools is a linux program. first need to run scripts there to get the data and then can load a summary table into r for the analysis. 

#Tajima's D: once the tajima's D values are generated for each population in vcftools, export out of linux and import into R. Then can calculate summary statistics. Will need to do this for each population. 

zczc125.tajima.atl<-read.table("./ocean_basins/zczc125.atl/zczc125.atl.100000.Tajima.D", header=TRUE)
summary(zczc125.tajima.atl$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.6104 -0.8074 -0.5964 -0.2228  0.2050  2.4817 

mden43.tajima.atl<-read.table("./ocean_basins/mden43.atl/mden43.atl.100000.Tajima.D", header=TRUE)
summary(mden43.tajima.atl$TajimaD)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -1.68380 -0.88483 -0.31246 -0.07331  0.75051  2.31782      107 

summary(mden43.tajima.indo$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -1.8778 -1.1470 -0.7637 -0.4502  0.0057  2.1304      83 
mden43.tajima.indo<-read.table("./ocean_basins/mden43.indo/mden43.indo.100000.Tajima.D", header=TRUE)

summary(read.table("./mden43.atl-bahamas.100000.Tajima.D", header=TRUE)$TajimaD)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -1.67053 -1.15524 -0.10305  0.03746  0.84228  1.93304      115
summary(read.table("./mden43.atl-east.100000.Tajima.D", header=TRUE)$TajimaD)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -1.71706 -0.78295 -0.13836  0.01493  0.85334  2.15921      120 
summary(read.table("./mden43.atl-other.100000.Tajima.D", header=TRUE)$TajimaD)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -1.56222 -1.11173  0.01499 -0.01658  0.81980  1.84427      109 
summary(read.table("./mden43.indo-africa.100000.Tajima.D", header=TRUE)$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -1.5622 -1.1117 -0.7999 -0.2846  0.5259  1.8443     113 
summary(read.table("./mden43.indo-hawaii.100000.Tajima.D", header=TRUE)$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -1.7469 -1.1405 -0.1949 -0.2625  0.5406  1.8912     125 
summary(read.table("./mden43.indo-south.100000.Tajima.D", header=TRUE)$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -1.5347 -1.0548 -1.0548 -0.2173  0.3335  1.7937     112 

#Heterozygosity/F: like Tajima's D, this will first need to be calculated in vcftools and a resulting file has values for each individual. Each file has one row for each individual and the following columns: average site depth across all loci/individual, observed homozygous sites, expected homozygous sites, no sites with data, no sites across whole dataset, F. Export these files out of linux and upload into R for summary statistics. 

#upload vcftools summaries into R and subset dataset for just Z. cavirostris data
sum<-read.table("~/Dropbox/Phd/Bioinformatics/bw_ddrad/Results/R-summaries/r-working-summaries.csv", sep=",", header=TRUE)
sum.zc<-subset(sum, sum$species=="Zcav")
sum.md<-subset(sum, sum$species=="Mden")

#population variables need to be converted to factors and characters for summaries to work
sum.zc$ocean.basin<-as.factor(sum.zc$ocean.basin)
sum.zc$pop1<-as.character(sum.zc$pop1)
sum.zc$pop1<-as.factor(sum.zc$pop1)
sum.zc$pop2<-as.character(sum.zc$pop2)
sum.zc$pop2<-as.factor(sum.zc$pop2)

sum.md$ocean.basin<-as.factor(sum.md$ocean.basin)
sum.md$population<-as.character(sum.md$population)
sum.md$population<-as.factor(sum.md$population)



#use "Psych" library to generate tables of summary statistics of F inbreeding coefficient
library(psych)
describeBy(sum.zc$f, sum.zc$ocean.basin, mat=TRUE, digits=3)
describeBy(sum.zc$f, sum.zc$pop1, mat=TRUE, digits=3)
describeBy(sum.zc$f, sum.zc$pop2, mat=TRUE, digits=3)

describeBy(sum.md$f, sum.md$ocean.basin, mat=TRUE, digits=3)
# item        group1 vars  n  mean    sd median trimmed   mad    min   max range  skew kurtosis    se
# X11    1      Atlantic    1 28 0.197 0.042  0.195   0.196 0.035  0.119 0.290 0.171 0.267   -0.303 0.008
# X12    2  Indo-Pacific    1 14 0.028 0.042  0.019   0.024 0.029 -0.032 0.145 0.177 1.295    1.809 0.011

describeBy(sum.md$f, sum.md$population, mat=TRUE, digits=3)
# item   group1 vars  n  mean    sd median trimmed   mad    min   max range  skew kurtosis    se
# X11    1  Atl-Bah    1  7 0.221 0.015  0.219   0.221 0.009  0.196 0.246 0.051 0.047   -0.825 0.006
# X12    2 Atl-East    1 16 0.187 0.043  0.187   0.185 0.022  0.119 0.290 0.171 0.818    0.317 0.011
# X13    3  Atl-Oth    1  5 0.197 0.056  0.206   0.197 0.063  0.128 0.279 0.151 0.211   -1.653 0.025
# X14    4 Indo-Afr    1  5 0.044 0.064  0.031   0.044 0.023 -0.032 0.145 0.177 0.458   -1.361 0.029
# X15    5 Indo-Haw    1  6 0.019 0.027  0.008   0.019 0.017 -0.010 0.057 0.067 0.403   -1.849 0.011
# X16    6 Indo-Sou    1  3 0.021 0.019  0.012   0.021 0.005  0.009 0.043 0.034 0.373   -2.333 0.011

#Calculate observed heterozygosity: 1-(observed homozygous snps/no snps per individual genotyped) and add this column to your dataframe
obs.het.ind<-(1-sum.zc$obs.homo.snps/sum.zc$no.snps.ind)
sum.zc<-cbind(sum.zc, obs.het.ind)

mden43.obs.het.ind<-(1-sum.md$obs.homo.snps/sum.md$no.snps.ind)
sum.md<-cbind(sum.md, mden43.obs.het.ind)

#Calculate summary statistics for each population
describeBy(sum.zc$obs.het.ind, sum.zc$ocean.basin, mat=TRUE, digits=3)
describeBy(sum.zc$obs.het.ind, sum.zc$pop1, mat=TRUE, digits=3)
describeBy(sum.zc$obs.het.ind, sum.zc$pop2, mat=TRUE, digits=3)

describeBy(sum.md$mden43.obs.het.ind, sum.md$ocean.basin, mat=TRUE, digits=3)
# item        group1 vars  n  mean    sd median trimmed   mad   min   max range   skew kurtosis    se
# X11    1      Atlantic    1 28 0.104 0.007  0.105   0.104 0.004 0.081 0.115 0.034 -1.175    1.758 0.001
# X12    2  Indo-Pacific    1 14 0.126 0.005  0.126   0.127 0.005 0.112 0.135 0.023 -1.017    1.160 0.001
describeBy(sum.md$mden43.obs.het.ind, sum.md$population, mat=TRUE, digits=3)
# item   group1 vars  n  mean    sd median trimmed   mad   min   max range   skew kurtosis    se
# X11    1  Atl-Bah    1  7 0.101 0.003  0.100   0.101 0.003 0.098 0.105 0.007  0.354   -1.548 0.001
# X12    2 Atl-East    1 16 0.106 0.007  0.106   0.106 0.003 0.087 0.115 0.028 -1.343    1.772 0.002
# X13    3  Atl-Oth    1  5 0.101 0.012  0.103   0.101 0.004 0.081 0.114 0.033 -0.678   -1.258 0.005
# X14    4 Indo-Afr    1  5 0.124 0.008  0.124   0.124 0.003 0.112 0.135 0.023 -0.239   -1.453 0.004
# X15    5 Indo-Haw    1  6 0.128 0.004  0.129   0.128 0.002 0.123 0.131 0.008 -0.425   -1.919 0.001
# X16    6 Indo-Sou    1  3 0.127 0.002  0.126   0.127 0.002 0.125 0.129 0.005  0.264   -2.333 0.001

describeBy(sum.md$no.snps.ind, sum.md$ocean.basin, mat=TRUE, digits=3)
# item        group1 vars  n     mean      sd median  trimmed    mad   min   max range   skew kurtosis     se
# X11    1      Atlantic    1 28 13721.64 494.056  13898 13834.21 29.652 11539 13933  2394 -3.312   11.097 93.368
# X12    2  Indo-Pacific    1 14 13835.21 139.324  13889 13860.00 22.239 13460 13913   453 -1.832    1.755 37.236
describeBy(sum.md$no.snps.ind, sum.md$population, mat=TRUE, digits=3)
# item   group1 vars  n     mean       sd median  trimmed     mad   min   max range   skew kurtosis      se
# X11    1  Atl-Bah    1  7 13760.00  230.418  13886 13760.00  34.100 13280 13909   629 -1.177   -0.228  87.090
# X12    2 Atl-East    1 16 13834.25  279.475  13911 13901.86  17.791 12789 13933  1144 -3.252    9.234  69.869
# X13    3  Atl-Oth    1  5 13307.60 1003.740  13730 13307.60 277.246 11539 13917  2378 -0.997   -1.018 448.886
# X14    4 Indo-Afr    1  5 13813.60  197.893  13903 13813.60  14.826 13460 13913   453 -1.067   -0.927  88.501
# X15    5 Indo-Haw    1  6 13877.33   27.573  13886 13877.33  13.343 13826 13906    80 -0.863   -0.792  11.257
# X16    6 Indo-Sou    1  3 13787.00  189.660  13896 13787.00   1.483 13568 13897   329 -0.385   -2.333 109.500
describeBy(sum.md$mean.site.depth, sum.md$ocean.basin, mat=TRUE, digits=3)
# item        group1 vars  n   mean     sd median trimmed    mad    min     max   range  skew kurtosis    se
# X11    1      Atlantic    1 28 47.624 25.960 46.608  47.071 30.564  4.768 105.429 100.661 0.234   -0.820 4.906
# X12    2  Indo-Pacific    1 14 64.310 30.484 61.991  61.985 19.388 13.543 142.972 129.430 0.815    0.919 8.147
describeBy(sum.md$mean.site.depth, sum.md$population, mat=TRUE, digits=3)
# item   group1 vars  n   mean     sd median trimmed    mad    min     max   range   skew kurtosis     se
# X11    1  Atl-Bah    1  7 36.043 24.365 33.293  36.043 22.556 10.993  82.084  71.091  0.731   -0.929  9.209
# X12    2 Atl-East    1 16 58.644 24.001 57.516  58.927 24.660  7.903 105.429  97.526 -0.098   -0.495  6.000
# X13    3  Atl-Oth    1  5 28.574 18.004 23.852  28.574 28.294  4.768  49.587  44.819 -0.052   -1.924  8.052
# X14    4 Indo-Afr    1  5 62.425 29.472 67.160  62.425 26.567 13.543  85.788  72.246 -0.721   -1.338 13.180
# X15    5 Indo-Haw    1  6 55.617  7.137 54.176  55.617  7.802 45.559  64.091  18.532  0.001   -1.758  2.913
# X16    6 Indo-Sou    1  3 84.837 57.516 83.579  84.837 82.459 27.961 142.972 115.011  0.022   -2.333 33.207
describeBy(sum.md$prop.missing.snps, sum.md$ocean.basin, mat=TRUE, digits=3)
# item        group1 vars  n  mean    sd median trimmed   mad   min   max range  skew kurtosis    se
# X11    1      Atlantic    1 28 0.019 0.035  0.006   0.011 0.002 0.004 0.175 0.171 3.312   11.097 0.007
# X12    2  Indo-Pacific    1 14 0.011 0.010  0.007   0.009 0.002 0.005 0.038 0.032 1.832    1.755 0.003
describeBy(sum.md$prop.missing.snps, sum.md$population, mat=TRUE, digits=3)
# item   group1 vars  n  mean    sd median trimmed   mad   min   max range  skew kurtosis    se
# X11    1  Atl-Bah    1  7 0.016 0.016  0.007   0.016 0.002 0.006 0.051 0.045 1.177   -0.228 0.006
# X12    2 Atl-East    1 16 0.011 0.020  0.006   0.006 0.001 0.004 0.086 0.082 3.252    9.234 0.005
# X13    3  Atl-Oth    1  5 0.049 0.072  0.018   0.049 0.020 0.005 0.175 0.170 0.997   -1.018 0.032
# X14    4 Indo-Afr    1  5 0.012 0.014  0.006   0.012 0.001 0.005 0.038 0.032 1.067   -0.927 0.006
# X15    5 Indo-Haw    1  6 0.008 0.002  0.007   0.008 0.001 0.006 0.012 0.006 0.864   -0.792 0.001
# X16    6 Indo-Sou    1  3 0.014 0.014  0.007   0.014 0.000 0.007 0.030 0.024 0.385   -2.333 0.008


####Genepop####
#ho and Gene diversities: Calculates for each locus as well as population. 
library(genepop)

#First you will need to convert your vcf file into a genepop file using PGD spider. once that is done, rename the individuals using text wrangler so that the first characters are the population. 

#run genedivho function for each dataset. The resulting file will be saved in your working directory and at the very bottom will have values for ho, 1-Qinter and 1-Qintra gene diversities. 
zczc125.gp.ho<-genedivho("./zczc125_genepop.txt")
zczc125.pop1.gp.ho<-genedivho("./zczc125.pop1.genepop.txt")
zczc125.pop2.gp.ho<-genedivho("./zczc125.pop2.genepop.txt")

mden43.gp.ho<-genedivho("./mden43.genepop.txt")
mden43.pop.gp.ho<-genedivho("./mden43.pop.genepop.txt")

test_diff("./zczc125_genepop.txt")

zczc125.genedivFis<-genedivFis(inputFile = "./zczc125_genepop_oceans.txt")

####Adegenet####
#Calculate Hs for each population using genpop files (Adegenet). Can also use to calculate observed heterozygosity for each population

library(adegenet)

#convert genind files to genpop files
zczc125.gp<-genind2genpop(zczc125.gi)
zczc125.pop1.gp<-genind2genpop(zczc125.pop1.gi)
zczc125.pop2.gp<-genind2genpop(zczc125.pop2.gi)

mden43.gp
mden43.pop.gp

#Calculate Hs for each population
Hs(zczc125.gp)
# Atlantic   Indopacific Mediterranean 
# 0.13569395    0.13068689    0.09758311 
Hs(zczc125.pop1.gp)
# Atl_AA     Atl_AB     Atl_AC     Atl_AD     Atl_AE     Indo_A     Indo_B     Indo_C      Med_A      Med_B      Med_C 
# 0.12909228 0.12075041 0.13128402 0.13173343 0.07395933 0.11250089 0.12674538 0.12796312 0.09385036 0.08047348 0.09431465 
Hs(zczc125.pop2.gp)
# Atl_BA     Atl_BB     Atl_BC     Atl_BD     Indo_A     Indo_B     Indo_C      Med_A      Med_B      Med_C 
# 0.13202948 0.13128402 0.13173343 0.07395933 0.11250089 0.12674538 0.12796312 0.09385036 0.08047348 0.09431465 

Hs(mden43.gp)
# Atlantic Indopacific 
# 0.1119443   0.1355113 
Hs(mden43.pop.gp)
# Atl-Bahamas    Atl-East   Atl-Other Indo-Africa Indo-Hawaii  Indo-South          NA 
# 0.1012754   0.1090775   0.1000685   0.1258402   0.1244822   0.1126243   0.0725622 

#use seppop and lapply to calculate observed heterozygosity for each population
# lapply(seppop(zczc125.gi), function(e) mean(summary(e)$Hobs, na.rm = TRUE))
# lapply(seppop(zczc125.pop1.gi), function(e) mean(summary(e)$Hobs, na.rm = TRUE))
# lapply(seppop(zczc125.pop2.gi), function(e) mean(summary(e)$Hobs, na.rm = TRUE))
# lapply(seppop(mden43.gi), function(e) mean(summary(e)$Hobs, na.rm = TRUE))
# # $Atlantic
# # [1] 0.1053599
# # 
# # $Indopacific
# # [1] 0.1270102
# lapply(seppop(mden43.pop.gi), function(e) mean(summary(e)$Hobs, na.rm = TRUE))
# # $`Atl-Bahamas`
# # [1] 0.1020844
# # 
# # $`Atl-East`
# # [1] 0.1064709
# # 
# # $`Atl-Other`
# # [1] 0.1063829
# # 
# # $`Indo-Africa`
# # [1] 0.1247117
# # 
# # $`Indo-Hawaii`
# # [1] 0.1280213
# # 
# # $`Indo-South`
# # [1] 0.1274278
# # 
# # $`NA`
# # [1] 0.130557
# 

zczc125.gl
mden43.gl

mden.atl.gi
mden.pac.gi

test<-Hs.test(mden43.gi[pop="Atlantic"], mden43.gi[pop="Indopacific"], n.sim=100)
# Monte-Carlo test
# Call: Hs.test(x = mden43.gi[pop = "Atlantic"], y = mden43.gi[pop = "Indopacific"], 
#               n.sim = 100)
# 
# Observation: -0.02356708 
# 
# Based on 100 replicates
# Simulated p-value: 0.00990099 
# Alternative hypothesis: two-sided 
# 
# Std.Obs.y   Expectation      Variance 
# -3.372702e+00 -7.749360e-04  4.566827e-05 

test<-Hs.test(mden43.gi[pop="Atlantic"], mden43.gi[pop="Indopacific"], n.sim=100)


####Diversity index differences ----------------------------------------------------------

#testing for significance between ocean basins in Ho and ho (per locus calculations Nei)
library(PMCMRplus)

zczc125.ho
head(zczc125.ho)
mden43.ho

zczc125.ho<-zczc125.hfstat.basicstats$Ho
zczc125.hs<-zczc125.hfstat.basicstats$Hs
zczc125.ho<-zczc125.hfstat.basicstats$ho

mden43.ho<-mden43.hfstat.basicstats$Ho
mden43.hs<-mden43.hfstat.basicstats$Hs
mden43.ho<-mden43.hfstat.basicstats$ho

zczc125.ho.melt<-melt(zczc125.ho)
zczc125.ho<-zczc125.ho.melt[,2:3]
colnames(zczc125.ho)<-c("pop", "ho")

mden43.ho.melt<-melt(mden43.ho)
head(mden43.ho.melt)
mden43.ho<-mden43.ho.melt[,2:3]
colnames(mden43.ho)<-c("pop", "ho")
head(mden43.ho)

zczc125.ho<-zczc125.hfstat.basicstats$ho
mden43.ho<-mden43.hfstat.basicstats$ho

zczc125.ho.melt<-melt(zczc125.ho)
mden43.ho.melt<-melt(mden43.ho)

zczc125.ho.melt<-zczc125.ho.melt[,2:3]
colnames(zczc125.ho.melt)<-c("pop", "ho")

mden43.ho.melt<-mden43.ho.melt[,2:3]
colnames(mden43.ho.melt)<-c("pop", "ho")

kruskalTest(zczc125.ho$ho, zczc125.ho$pop, p.adjust.method = "bonferroni")
#Kruskal-Wallis test

#data:  zczc125.ho$ho and zczc125.ho$pop
#chi-squared = 6236.1, df = 2, p-value < 2.2e-16

kruskalTest(mden43.ho, p.adjust.method = "bonferroni")
#Kruskal-Wallis test

#data:  mden43.ho
#chi-squared = 43588, df = 1, p-value < 2.2e-16

#Posthoc tests
#Conover
#Dunn
#Nemenyi
#kwAllPairsDunnTest - Performs Dunn's non-parametric all-pairs comparison test for Kruskal-type ranked data. For all-pairs comparisons in an one-factorial layout with non-normally distributed residuals Dunn's non-parametric test can be performed. A total of m = k(k-1)/2 hypotheses can be tested. The null hypothesis H_{ij}: μ_i(x) = μ_j(x) is tested in the two-tailed test against the alternative A_{ij}: μ_i(x) \ne μ_j(x), ~~ i \ne j.The p-values are computed from the standard normal distribution using any of the p-adjustment methods as included in p.adjust. Originally, Dunn (1964) proposed Bonferroni's p-adjustment method.

kwAllPairsDunnTest(zczc125.ho$ho, zczc125.ho$pop, p.adjust.method = "bonferroni")
#               Atlantic Indopacific
# Indopacific   <2e-16   -          
# Mediterranean <2e-16   <2e-16 
kwAllPairsConoverTest(zczc125.ho$ho, zczc125.ho$pop, p.adjust.method = "bonferroni")
#               Atlantic Indopacific
# Indopacific   <2e-16   -          
# Mediterranean <2e-16   <2e-16
kwAllPairsNemenyiTest(zczc125.ho$ho, zczc125.ho$pop, p.adjust.method = "bonferroni")
#               Atlantic Indopacific
# Indopacific   <2e-16   -          
# Mediterranean <2e-16   <2e-16 

zczc125.ho.all<-zczc125.hfstat.basicstats$Ho
hist(zczc125.ho.all[,1])
hist(zczc125.ho.all[,2])
hist(zczc125.ho.all[,3])

head(mden43.ho)
kruskalTest(mden43.ho$ho, mden43.ho$pop, p.adjust.method = "bonferroni")
# Kruskal-Wallis test
# 
# data:  mden43.ho$ho and mden43.ho$pop
# chi-squared = 714.47, df = 1, p-value < 2.2e-16
kwAllPairsDunnTest(mden43.ho$ho, mden43.ho$pop, p.adjust.method = "bonferroni") #p-value < 2.2e-16
kwAllPairsConoverTest(mden43.ho$ho, mden43.ho$pop, p.adjust.method = "bonferroni") #p-value < 2.2e-16
kwAllPairsNemenyiTest(mden43.ho$ho, mden43.ho$pop, p.adjust.method = "bonferroni") #p-value < 2.2e-16

kwAllPairsConoverTest(zczc125.ho.melt$ho, zczc125.ho.melt$pop, p.adjust.method = "bonferroni")
#               Atlantic Indopacific
# Indopacific   2.2e-13  -          
# Mediterranean 1.1e-05  < 2e-16
kwAllPairsDunnTest(zczc125.ho.melt$ho, zczc125.ho.melt$pop, p.adjust.method = "bonferroni")
#               Atlantic Indopacific
# Indopacific   2.3e-13  -          
# Mediterranean 1.1e-05  < 2e-16 
kwAllPairsNemenyiTest(zczc125.ho.melt$ho, zczc125.ho.melt$pop, p.adjust.method = "bonferroni")
#               Atlantic Indopacific
# Indopacific   2.7e-13  -          
# Mediterranean 1.1e-05  3.0e-14 

kwAllPairsDunnTest(mden43.ho.melt$ho, mden43.ho.melt$pop, p.adjust.method = "bonferroni")  #0.0044
kwAllPairsConoverTest(mden43.ho.melt$ho, mden43.ho.melt$pop, p.adjust.method = "bonferroni") #0.0044
kwAllPairsNemenyiTest(mden43.ho.melt$ho, mden43.ho.melt$pop, p.adjust.method = "bonferroni") #0.0055


# Hs differences using Adegenet -------------------------------------------


#testing for differences in Hs from adegenet
library(adegenet)

zczc125.gi
zczc125.pop1.gi
zczc125.pop2.gi

mden43.gi
mden43.pop.gi

Hs.test(zczc125.gi[pop="Atlantic"], zczc125.gi[pop="Indopacific"])
# Monte-Carlo test
# Call: Hs.test(x = zczc125.gi[pop = "Atlantic"], y = zczc125.gi[pop = "Indopacific"])
# 
# Observation: 0.005007064 
# 
# Based on 999 replicates
# Simulated p-value: 0.001 
# Alternative hypothesis: two-sided 
# 
# Std.Obs.y  Expectation     Variance 
# 3.766721e+00 8.684779e-05 1.706245e-06 

Hs.test(zczc125.gi[pop="Atlantic"], zczc125.gi[pop="Mediterranean"])
# Monte-Carlo test
# Call: Hs.test(x = zczc125.gi[pop = "Atlantic"], y = zczc125.gi[pop = "Mediterranean"])
# 
# Observation: 0.03811084 
# 
# Based on 999 replicates
# Simulated p-value: 0.001 
# Alternative hypothesis: two-sided 
# 
# Std.Obs.y  Expectation     Variance 
# 1.312849e+01 4.696435e-04 8.220482e-06 

Hs.test(zczc125.gi[pop="Indopacific"], zczc125.gi[pop="Mediterranean"])
# Monte-Carlo test
# Call: Hs.test(x = zczc125.gi[pop = "Indopacific"], y = zczc125.gi[pop = "Mediterranean"])
# 
# Observation: 0.03310378 
# 
# Based on 999 replicates
# Simulated p-value: 0.001 
# Alternative hypothesis: two-sided 
# 
# Std.Obs.y  Expectation     Variance 
# 8.403159e+00 3.450145e-04 1.519742e-05 

Hs.test(mden43.gi[pop="Atlantic"], mden43.gi[pop="Indopacific"])
# Monte-Carlo test
# Call: Hs.test(x = mden43.gi[pop = "Atlantic"], y = mden43.gi[pop = "Indopacific"])
# 
# Observation: -0.02356708 
# 
# Based on 999 replicates
# Simulated p-value: 0.001 
# Alternative hypothesis: two-sided 
# 
# Std.Obs.y   Expectation      Variance 
# -3.715489e+00 -3.137009e-04  3.916881e-05 



# KW test for Ho, Hs and ho ----------------------------------------------
#KW tests for Ho, Hs and ho for populations:datafiles to work with
zczc125.pop1.basicstats
zczc125.pop2.basicstats
mden43.pop.basicstats

zczc125.pop1.ho<-zczc125.pop1.basicstats$Ho
zczc125.pop1.hs<-zczc125.pop1.basicstats$Hs
zczc125.pop1.fis<-zczc125.pop1.basicstats$Fis

zczc125.pop2.ho<-zczc125.pop2.basicstats$Ho
zczc125.pop2.hs<-zczc125.pop2.basicstats$Hs
zczc125.pop2.fis<-zczc125.pop2.basicstats$Fis

mden43.pop.ho<-mden43.pop.basicstats$Ho
mden43.pop.hs<-mden43.pop.basicstats$Hs
mden43.pop.fis<-mden43.pop.basicstats$Fis

mden43.pop.ho<-mden43.pop.ho[,1:6]
mden43.pop.hs<-mden43.pop.hs[,1:6]
mden43.pop.fis<-mden43.pop.fis[,1:6]

zczc125.pop1.ho.melt<-melt(zczc125.pop1.ho)
zczc125.pop1.hs.melt<-melt(zczc125.pop1.hs)
zczc125.pop1.fis.melt<-melt(zczc125.pop1.fis)

zczc125.pop2.ho.melt<-melt(zczc125.pop2.ho)
zczc125.pop2.hs.melt<-melt(zczc125.pop2.hs)
zczc125.pop2.fis.melt<-melt(zczc125.pop2.fis)

mden43.pop.ho.melt<-melt(mden43.pop.ho)
mden43.pop.hs.melt<-melt(mden43.pop.hs)
mden43.pop.fis.melt<-melt(mden43.pop.fis)

zczc125.pop1.ho<-zczc125.pop1.ho.melt[,2:3]
zczc125.pop1.hs<-zczc125.pop1.hs.melt[,2:3]
zczc125.pop1.fis<-zczc125.pop1.fis.melt[,2:3]

zczc125.pop2.ho<-zczc125.pop2.ho.melt[,2:3]
zczc125.pop2.hs<-zczc125.pop2.hs.melt[,2:3]
zczc125.pop2.fis<-zczc125.pop2.fis.melt[,2:3]

mden43.pop.ho<-mden43.pop.ho.melt[,2:3]
mden43.pop.hs<-mden43.pop.hs.melt[,2:3]
mden43.pop.fis<-mden43.pop.fis.melt[,2:3]

colnames(zczc125.pop1.ho)<-c("pop", "ho")
colnames(zczc125.pop1.hs)<-c("pop", "hs")
colnames(zczc125.pop1.fis)<-c("pop", "fis")

colnames(zczc125.pop2.ho)<-c("pop", "ho")
colnames(zczc125.pop2.hs)<-c("pop", "hs")
colnames(zczc125.pop2.fis)<-c("pop", "fis")

colnames(mden43.pop.ho)<-c("pop", "ho")
colnames(mden43.pop.hs)<-c("pop", "hs")
colnames(mden43.pop.fis)<-c("pop", "fis")

#run kruskal wallis tests on each dataset:
#zczc125 pop1
kwAllPairsDunnTest(zczc125.pop1.ho$ho, zczc125.pop1.ho$pop, p.adjust.method = "bonferroni")
kwAllPairsConoverTest(zczc125.pop1.ho$ho, zczc125.pop1.ho$pop, p.adjust.method = "bonferroni")
kwAllPairsNemenyiTest(zczc125.pop1.ho$ho, zczc125.pop1.ho$pop, p.adjust.method = "bonferroni")

kwAllPairsDunnTest(zczc125.pop1.hs$hs, zczc125.pop1.hs$pop, p.adjust.method = "bonferroni")
kwAllPairsConoverTest(zczc125.pop1.hs$hs, zczc125.pop1.hs$pop, p.adjust.method = "bonferroni")
kwAllPairsNemenyiTest(zczc125.pop1.hs$hs, zczc125.pop1.hs$pop, p.adjust.method = "bonferroni")

kwAllPairsDunnTest(zczc125.pop1.fis$fis, zczc125.pop1.fis$pop, p.adjust.method = "bonferroni")
kwAllPairsConoverTest(zczc125.pop1.fis$fis, zczc125.pop1.fis$pop, p.adjust.method = "bonferroni")
kwAllPairsNemenyiTest(zczc125.pop1.fis$fis, zczc125.pop1.fis$pop, p.adjust.method = "bonferroni")


#zczc125 pop2
kwAllPairsDunnTest(zczc125.pop2.ho$ho, zczc125.pop2.ho$pop, p.adjust.method = "bonferroni")
kwAllPairsConoverTest(zczc125.pop2.ho$ho, zczc125.pop2.ho$pop, p.adjust.method = "bonferroni")
kwAllPairsNemenyiTest(zczc125.pop2.ho$ho, zczc125.pop2.ho$pop, p.adjust.method = "bonferroni")

kwAllPairsDunnTest(zczc125.pop2.hs$hs, zczc125.pop2.hs$pop, p.adjust.method = "bonferroni")
kwAllPairsConoverTest(zczc125.pop2.hs$hs, zczc125.pop2.hs$pop, p.adjust.method = "bonferroni")
kwAllPairsNemenyiTest(zczc125.pop2.hs$hs, zczc125.pop2.hs$pop, p.adjust.method = "bonferroni")

kwAllPairsDunnTest(zczc125.pop2.fis$fis, zczc125.pop2.fis$pop, p.adjust.method = "bonferroni")
kwAllPairsConoverTest(zczc125.pop2.fis$fis, zczc125.pop2.fis$pop, p.adjust.method = "bonferroni")
kwAllPairsNemenyiTest(zczc125.pop2.fis$fis, zczc125.pop2.fis$pop, p.adjust.method = "bonferroni")

#mden43 pop
kwAllPairsDunnTest(mden43.pop.ho$ho, mden43.pop.ho$pop, p.adjust.method = "bonferroni")
kwAllPairsConoverTest(mden43.pop.ho$ho, mden43.pop.ho$pop, p.adjust.method = "bonferroni")
kwAllPairsNemenyiTest(mden43.pop.ho$ho, mden43.pop.ho$pop, p.adjust.method = "bonferroni")

kwAllPairsDunnTest(mden43.pop.hs$hs, mden43.pop.hs$pop, p.adjust.method = "bonferroni")
kwAllPairsConoverTest(mden43.pop.hs$hs, mden43.pop.hs$pop, p.adjust.method = "bonferroni")
kwAllPairsNemenyiTest(mden43.pop.hs$hs, mden43.pop.hs$pop, p.adjust.method = "bonferroni")

kwAllPairsDunnTest(mden43.pop.fis$fis, mden43.pop.fis$pop, p.adjust.method = "bonferroni")
kwAllPairsConoverTest(mden43.pop.fis$fis, mden43.pop.fis$pop, p.adjust.method = "bonferroni")
kwAllPairsNemenyiTest(mden43.pop.fis$fis, mden43.pop.fis$pop, p.adjust.method = "bonferroni")


# 95% CI for diversity indices means --------------------------------------
zczc125.pop1.ho
zczc125.pop1.hs
zczc125.pop1.fis

zczc125.pop2.ho
zczc125.pop2.hs
zczc125.pop2.fis

mden43.pop.ho
mden43.pop.hs
mden43.pop.fis

zczc125.pop1.ho.nona<-na.omit(zczc125.pop1.ho)

dim(zczc125.pop1.ho)
dim(zczc125.pop1.ho.nona)

# You can use the aggregate functions to calculate the
# mean and 95% CI.
# The formula for calculating a CI is: Z * standard error
# so in this case the 95% CI = 1.96 * (sd(x)/sqrt(length(x)))
# calculate the means of the populations
df = aggregate(list(ho = zczc125.pop1.ho.nona$ho),list(pop = zczc125.pop1.ho.nona$pop),mean)
# and the 95% CI (I've called them se in this case)
df1 = aggregate(list(se = zczc125.pop1.ho.nona$ho),list(pop = zczc125.pop1.ho.nona$pop),
                FUN = function(x){1.96*(sd(x)/sqrt(length(x)))})
# Add the CI's to the dataframe containing the means
df$se = df1$se
# and create your plot
ggplot(df,aes(x=pop, y=ho)) + geom_point() +
  geom_errorbar(aes(ymax = ho+se, ymin = ho-se,x=pop), width = .25)

zczc125.pop1.ho
zczc125.pop1.hs
zczc125.pop1.fis

zczc125.pop2.ho
zczc125.pop2.hs
zczc125.pop2.fis

mden43.pop.ho
mden43.pop.hs
mden43.pop.fis

# zczc125.pop1.ho.nona<-na.omit(zczc125.pop1.ho)
# 
# mymean = function(x, indices){return(mean(x[indices]))}
# 
# zczc125.pop1.ho.boot<-boot(zczc125.pop1.ho$ho, mymean, R=100)
# med.ho.boot.ci<-boot.ci(med.ho.boot, conf=0.95, type="basic")
# med.ho.boot.ci
# 
# aggregate(zczc125.pop1.ho, by=list(pop), FUN=mymean, na.rm=TRUE)
# 
# head(zczc125.pop1.ho$pop)

zczc125.pop1.basicstats
zczc125.pop2.basicstats
mden43.pop.basicstats

zczc125.ho<-zczc125.gl.basicstats$Ho
zczc125.pop1.ho<-zczc125.pop1.basicstats$Ho
zczc125.pop1.hs<-zczc125.pop1.basicstats$Hs
zczc125.pop1.fis<-zczc125.pop1.basicstats$Fis

zczc125.pop2.ho<-zczc125.pop2.basicstats$Ho
zczc125.pop2.hs<-zczc125.pop2.basicstats$Hs
zczc125.pop2.fis<-zczc125.pop2.basicstats$Fis

mden43.pop.ho<-mden43.pop.basicstats$Ho
mden43.pop.hs<-mden43.pop.basicstats$Hs
mden43.pop.fis<-mden43.pop.basicstats$Fis

#will go through each column and calculate the mean
apply(na.omit(zczc125.pop1.ho),2,mean)
apply(na.omit(zczc125.ho), 2, mean)

wilcox.test(zczc125.ho, conf.int=TRUE, conf.level=0.95)

apply(na.omit(zczc125.ho), 2, wilcox.test(zczc125.ho, conf.int=TRUE, conf.level=0.95))



# calculate confidence intervals for diversity indices data using rcompanion ----------------
library(rcompanion)

zczc125.pop1.ho<-na.omit(zczc125.pop1.ho)
zczc125.pop1.hs<-na.omit(zczc125.pop1.hs)
zczc125.pop1.fis<-na.omit(zczc125.pop1.fis)
zczc125.pop2.ho<-na.omit(zczc125.pop2.ho)
zczc125.pop2.hs<-na.omit(zczc125.pop2.hs)
zczc125.pop2.fis<-na.omit(zczc125.pop2.fis)
mden43.pop.ho<-na.omit(mden43.pop.ho)
mden43.pop.hs<-na.omit(mden43.pop.hs)
mden43.pop.fis<-na.omit(mden43.pop.fis)

zczc125.ho.melt<-melt(zczc125.ho)
zczc125.ho.melt<-zczc125.ho.melt[,2:3]
colnames(zczc125.ho.melt)<-c("pop", "ho")
zczc125.ho.melt

zczc125.hs.melt<-melt(zczc125.hs)
zczc125.hs.melt<-zczc125.hs.melt[,2:3]
colnames(zczc125.hs.melt)<-c("pop", "hs")
zczc125.hs.melt

zczc125.fis.melt<-melt(zczc125.fis)
zczc125.fis.melt<-zczc125.fis.melt[,2:3]
colnames(zczc125.fis.melt)<-c("pop", "fis")
zczc125.fis.melt

mden43.ho.melt<-mden43.ho

mden43.hs.melt<-melt(mden43.hs)
mden43.hs.melt<-mden43.hs.melt[,2:3]
colnames(mden43.hs.melt)<-c("pop", "hs")

mden43.fis.melt<-melt(mden43.fis)
mden43.fis.melt<-mden43.fis.melt[,2:3]
colnames(mden43.fis.melt)<-c("pop", "fis")

zczc125.ho.melt<-na.omit(zczc125.ho.melt)
zczc125.hs.melt<-na.omit(zczc125.hs.melt)
zczc125.fis.melt<-na.omit(zczc125.fis.melt)
mden43.ho.melt<-na.omit(mden43.ho.melt)
mden43.hs.melt<-na.omit(mden43.hs.melt)
mden43.fis.melt<-na.omit(mden43.fis.melt)

zczc125.pop1.ho
zczc125.pop1.hs
zczc125.pop1.fis
zczc125.pop2.ho
zczc125.pop2.hs
zczc125.pop2.fis
mden43.pop.ho
mden43.pop.hs
mden43.pop.fis

#zczc125 pop1 ho (plotting without population AE because there are only two samples[-5,])
Sum.zczc125.pop1.ho = groupwiseMean(ho ~ pop, data= zczc125.pop1.ho, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.pop1.ho
X = c(1:4,6-11)
Y = Sum.zczc125.pop1.ho$Percentile.upper + 0.2
png(filename="./zczc125_pop_ho.png", width=500, height=300)
ggplot(Sum.zczc125.pop1.ho[-5,], aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Observed Heterozygosity") + xlab(NULL)
dev.off()
# pop     n   Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1  Atl_AA 25056 0.1180       0.95     0.1160     0.1200           0.1160           0.1200
# 2  Atl_AB 25056 0.1100       0.95     0.1070     0.1120           0.1070           0.1120
# 3  Atl_AC 25056 0.1190       0.95     0.1170     0.1210           0.1170           0.1210
# 4  Atl_AD 25056 0.1190       0.95     0.1170     0.1210           0.1170           0.1210
# 5  Atl_AE 24790 0.1280       0.95     0.1230     0.1320           0.1230           0.1320
# 6  Indo_A 25056 0.1090       0.95     0.1060     0.1110           0.1060           0.1110
# 7  Indo_B 25056 0.1160       0.95     0.1140     0.1180           0.1140           0.1180
# 8  Indo_C 25056 0.1180       0.95     0.1160     0.1200           0.1160           0.1200
# 9   Med_A 25056 0.0886       0.95     0.0866     0.0907           0.0865           0.0905
# 10  Med_B 25056 0.0778       0.95     0.0757     0.0798           0.0758           0.0797
# 11  Med_C 25055 0.0860       0.95     0.0838     0.0881           0.0837           0.0882

#zczc125 pop1 hs
Sum.zczc125.pop1.hs = groupwiseMean(hs ~ pop, data= zczc125.pop1.hs, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.pop1.hs
X = 1:11
Y = Sum$Percentile.upper + 0.2
png(filename="./zczc125_pop_hs.png", width=500, height=300)
ggplot(Sum.zczc125.pop1.hs[-5,], aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Gene Diversity") + xlab(NULL)
dev.off()
# pop     n   Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1  Atl_AA 25056 0.1360       0.95     0.1340     0.1380           0.1340           0.1380
# 2  Atl_AB 25055 0.1390       0.95     0.1360     0.1410           0.1360           0.1410
# 3  Atl_AC 25056 0.1360       0.95     0.1340     0.1380           0.1340           0.1380
# 4  Atl_AD 25056 0.1360       0.95     0.1340     0.1380           0.1340           0.1380
# 5  Atl_AE 24621 0.0632       0.95     0.0611     0.0652           0.0611           0.0653
# 6  Indo_A 25052 0.1280       0.95     0.1260     0.1300           0.1260           0.1300
# 7  Indo_B 25056 0.1330       0.95     0.1310     0.1350           0.1310           0.1350
# 8  Indo_C 25056 0.1320       0.95     0.1300     0.1340           0.1300           0.1340
# 9   Med_A 25056 0.0969       0.95     0.0948     0.0990           0.0948           0.0991
# 10  Med_B 25051 0.0859       0.95     0.0838     0.0880           0.0837           0.0882
# 11  Med_C 25038 0.1080       0.95     0.1060     0.1100           0.1060           0.1100

#zczc125 pop1 fis
Sum.zczc125.pop1.fis = groupwiseMean(fis ~ pop, data= zczc125.pop1.fis, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.pop1.fis
X = 1:11
Y = Sum$Percentile.upper + 0.2
png(filename="./zczc125_pop_fis.png", width=500, height=300)
ggplot(Sum.zczc125.pop1.fis[-5,], aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Inbreeding Coefficient") + xlab(NULL)
dev.off()
# pop     n    Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1  Atl_AA 13978  0.0969       0.95     0.0916     0.1020           0.0922           0.1020
# 2  Atl_AB 10166  0.1420       0.95     0.1340     0.1500           0.1340           0.1490
# 3  Atl_AC 16426  0.0949       0.95     0.0905     0.0994           0.0906           0.0996
# 4  Atl_AD 17846  0.0969       0.95     0.0927     0.1010           0.0929           0.1010
# 5  Atl_AE  3110 -0.9840       0.95    -0.9880    -0.9790          -0.9880          -0.9790
# 6  Indo_A  9180  0.0984       0.95     0.0901     0.1070           0.0904           0.1070
# 7  Indo_B 14504  0.0961       0.95     0.0911     0.1010           0.0913           0.1020
# 8  Indo_C 16126  0.0899       0.95     0.0855     0.0943           0.0853           0.0942
# 9   Med_A  7889  0.0755       0.95     0.0689     0.0821           0.0686           0.0818
# 10  Med_B  6258  0.0760       0.95     0.0672     0.0847           0.0670           0.0846
# 11  Med_C  7347  0.1440       0.95     0.1340     0.1540           0.1340           0.1550

#zczc125 pop2 ho (plotting without population BD because there are only two samples)
Sum.zczc125.pop2.ho = groupwiseMean(ho ~ pop, data= zczc125.pop2.ho, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.pop2.ho
X = c(1:3,5-10)
Y = Sum.zczc125.pop2.ho$Percentile.upper + 0.2
png(filename="./zczc125_pop2_ho.png", width=500, height=300)
ggplot(Sum.zczc125.pop2.ho[-4,], aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Observed Heterozygosity") + xlab(NULL)
dev.off()
# pop     n   Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1  Atl_BA 25056 0.1160       0.95     0.1140     0.1170           0.1140           0.1180
# 2  Atl_BB 25056 0.1190       0.95     0.1170     0.1210           0.1170           0.1210
# 3  Atl_BC 25056 0.1190       0.95     0.1170     0.1210           0.1170           0.1210
# 4  Atl_BD 24790 0.1280       0.95     0.1230     0.1320           0.1230           0.1320
# 5  Indo_A 25056 0.1090       0.95     0.1060     0.1110           0.1060           0.1110
# 6  Indo_B 25056 0.1160       0.95     0.1140     0.1180           0.1140           0.1180
# 7  Indo_C 25056 0.1180       0.95     0.1160     0.1200           0.1160           0.1200
# 8   Med_A 25056 0.0886       0.95     0.0866     0.0907           0.0866           0.0907
# 9   Med_B 25056 0.0778       0.95     0.0757     0.0798           0.0757           0.0797
# 10  Med_C 25055 0.0860       0.95     0.0838     0.0881           0.0839           0.0882

#zczc125 pop2 hs (plotting without population BD because there are only two samples)
Sum.zczc125.pop2.hs = groupwiseMean(hs ~ pop, data= zczc125.pop2.hs, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.pop2.hs
X = c(1:3,5-10)
Y = Sum.zczc125.pop2.hs$Percentile.upper + 0.2
png(filename="./zczc125_pop2_hs.png", width=500, height=300)
ggplot(Sum.zczc125.pop2.hs[-4,], aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Gene Diversity") + xlab(NULL)
dev.off()
# pop     n   Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1  Atl_BA 25056 0.1370       0.95     0.1350     0.1390           0.1350           0.1390
# 2  Atl_BB 25056 0.1360       0.95     0.1340     0.1380           0.1340           0.1380
# 3  Atl_BC 25056 0.1360       0.95     0.1340     0.1380           0.1340           0.1380
# 4  Atl_BD 24621 0.0632       0.95     0.0611     0.0652           0.0610           0.0653
# 5  Indo_A 25052 0.1280       0.95     0.1260     0.1300           0.1260           0.1300
# 6  Indo_B 25056 0.1330       0.95     0.1310     0.1350           0.1310           0.1350
# 7  Indo_C 25056 0.1320       0.95     0.1300     0.1340           0.1300           0.1340
# 8   Med_A 25056 0.0969       0.95     0.0948     0.0990           0.0948           0.0991
# 9   Med_B 25051 0.0859       0.95     0.0838     0.0880           0.0837           0.0881
# 10  Med_C 25038 0.1080       0.95     0.1060     0.1100           0.1060           0.1110

#zczc125 pop2 fis(plotting without population BD because there are only two samples)
Sum.zczc125.pop2.fis = groupwiseMean(fis ~ pop, data= zczc125.pop2.fis, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.pop2.fis
X = c(1:3,5-10)
Y = Sum.zczc125.pop2.fis$Percentile.upper + 0.2
png(filename="./zczc125_pop2_fis.png", width=500, height=300)
ggplot(Sum.zczc125.pop2.fis[-4,], aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Inbreeding Coefficient") + xlab(NULL)
dev.off()
# pop     n    Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1  Atl_BA 16276  0.1250       0.95     0.1200     0.1300           0.1200           0.1290
# 2  Atl_BB 16426  0.0949       0.95     0.0905     0.0994           0.0903           0.0997
# 3  Atl_BC 17846  0.0969       0.95     0.0927     0.1010           0.0927           0.1010
# 4  Atl_BD  3110 -0.9840       0.95    -0.9880    -0.9790          -0.9880          -0.9790
# 5  Indo_A  9180  0.0984       0.95     0.0901     0.1070           0.0899           0.1070
# 6  Indo_B 14504  0.0961       0.95     0.0911     0.1010           0.0913           0.1020
# 7  Indo_C 16126  0.0899       0.95     0.0855     0.0943           0.0856           0.0943
# 8   Med_A  7889  0.0755       0.95     0.0689     0.0821           0.0688           0.0818
# 9   Med_B  6258  0.0760       0.95     0.0672     0.0847           0.0672           0.0850
# 10  Med_C  7347  0.1440       0.95     0.1340     0.1540           0.1340           0.1540

#mden43 pop ho(plotting without population BD because there are only two samples)
Sum.mden43.pop.ho = groupwiseMean(ho ~ pop, data= mden43.pop.ho, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.mden43.pop.ho
X = 6
Y = Sum.mden43.pop.ho$Percentile.upper + 0.2
png(filename="./mden43_pop_ho.png", width=500, height=300)
ggplot(Sum.mden43.pop.ho, aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Observed Heterozygosity") + xlab(NULL)
dev.off()
# pop     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1 Atl-Bahamas 13988 0.102       0.95      0.099      0.105           0.0991            0.105
# 2    Atl-East 13988 0.106       0.95      0.104      0.109           0.1040            0.109
# 3   Atl-Other 13988 0.106       0.95      0.103      0.110           0.1030            0.110
# 4 Indo-Africa 13988 0.125       0.95      0.122      0.128           0.1210            0.128
# 5 Indo-Hawaii 13987 0.128       0.95      0.125      0.131           0.1250            0.131
# 6  Indo-South 13987 0.127       0.95      0.124      0.131           0.1240            0.131

#mden43 pop hs(plotting without population BD because there are only two samples)
Sum.mden43.pop.hs = groupwiseMean(hs ~ pop, data= mden43.pop.hs, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.mden43.pop.hs
X = 6
Y = Sum.mden43.pop.hs$Percentile.upper + 0.2
png(filename="./mden43_pop_hs.png", width=500, height=300)
ggplot(Sum.mden43.pop.hs, aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Gene Diversity") + xlab(NULL)
dev.off()
# pop     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1 Atl-Bahamas 13988 0.110       0.95      0.107      0.113            0.107            0.113
# 2    Atl-East 13988 0.113       0.95      0.110      0.116            0.110            0.115
# 3   Atl-Other 13984 0.114       0.95      0.111      0.117            0.110            0.117
# 4 Indo-Africa 13988 0.142       0.95      0.139      0.145            0.139            0.145
# 5 Indo-Hawaii 13987 0.137       0.95      0.134      0.140            0.133            0.140
# 6  Indo-South 13956 0.138       0.95      0.135      0.142            0.135            0.142

#mden43 pop fis(plotting without population BD because there are only two samples)
Sum.mden43.pop.fis = groupwiseMean(fis ~ pop, data= mden43.pop.fis, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.mden43.pop.fis
X = 6
Y = Sum.mden43.pop.fis$Percentile.upper + 0.2
png(filename="./mden43_pop_fis.png", width=500, height=300)
ggplot(Sum.mden43.pop.fis, aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Inbreeding Coefficient") + xlab(NULL)
dev.off()
# pop    n   Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1 Atl-Bahamas 4758 0.0488       0.95    0.03820     0.0594          0.03760           0.0598
# 2    Atl-East 6278 0.0478       0.95    0.04110     0.0544          0.04080           0.0543
# 3   Atl-Other 4423 0.0327       0.95    0.02140     0.0440          0.02100           0.0449
# 4 Indo-Africa 6206 0.0764       0.95    0.06740     0.0855          0.06780           0.0852
# 5 Indo-Hawaii 6450 0.0372       0.95    0.02940     0.0450          0.02910           0.0449
# 6  Indo-South 4473 0.0208       0.95    0.00861     0.0329          0.00962           0.0336


####ocean basin confidence intervals
#zcav ho
Sum.zczc125.ho.melt = groupwiseMean(ho ~ pop, data= zczc125.ho.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.ho.melt
X = 3
Y = Sum.zczc125.ho.melt$Percentile.upper + 0.2
png(filename="./zczc125_ho.png", width=500, height=300)
ggplot(Sum.zczc125.ho.melt, aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Observed Heterozygosity") + xlab(NULL)
dev.off()
# pop     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1      Atlantic 25056 0.118       0.95     0.1170     0.1200            0.117           0.1200
# 2   Indopacific 25056 0.116       0.95     0.1140     0.1180            0.114           0.1180
# 3 Mediterranean 25056 0.085       0.95     0.0831     0.0868            0.083           0.0869

#zcav hs
Sum.zczc125.hs.melt = groupwiseMean(hs ~ pop, data= zczc125.hs.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.hs.melt
X = 3
Y = Sum.zczc125.hs.melt$Percentile.upper + 0.2
png(filename="./zczc125_hs.png", width=500, height=300)
ggplot(Sum.zczc125.hs.melt, aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Gene Diversity") + xlab(NULL)
dev.off()
# pop     n   Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1      Atlantic 25056 0.1370       0.95     0.1350      0.139           0.1350            0.139
# 2   Indopacific 25056 0.1330       0.95     0.1310      0.135           0.1310            0.135
# 3 Mediterranean 25056 0.0994       0.95     0.0974      0.101           0.0974            0.101

#zcav fis
Sum.zczc125.fis.melt = groupwiseMean(fis ~ pop, data= zczc125.fis.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.fis.melt
X = 3
Y = Sum.zczc125.fis.melt$Percentile.upper + 0.2
png(filename="./zczc125_fis.png", width=500, height=300)
ggplot(Sum.zczc125.fis.melt, aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Inbreeding Coefficient") + xlab(NULL)
dev.off()
# pop     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1      Atlantic 22444 0.126       0.95      0.123      0.130            0.123            0.130
# 2   Indopacific 19552 0.113       0.95      0.109      0.117            0.109            0.117
# 3 Mediterranean  9895 0.130       0.95      0.125      0.136            0.125            0.137

#mden43 ho
Sum.mden43.ho.melt = groupwiseMean(ho ~ pop, data= mden43.ho.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.mden43.ho.melt
X = 2
Y = Sum.mden43.ho.melt$Percentile.upper + 0.2
png(filename="./mden43_ho.png", width=500, height=300)
ggplot(Sum.mden43.ho.melt, aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Observed Heterozygosity") + xlab(NULL)
dev.off()
# pop     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1    Atlantic 13988 0.105       0.95      0.103      0.108            0.103            0.108
# 2 Indopacific 13988 0.127       0.95      0.125      0.129            0.125            0.129

#mden43 hs
Sum.mden43.hs.melt = groupwiseMean(hs ~ pop, data= mden43.hs.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.mden43.hs.melt
X = 2
Y = Sum.mden43.hs.melt$Percentile.upper + 0.2
png(filename="./mden43_hs.png", width=500, height=300)
ggplot(Sum.mden43.hs.melt, aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Gene Diversity") + xlab(NULL)
dev.off()
# pop     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1    Atlantic 13988 0.114       0.95      0.112      0.117            0.112            0.117
# 2 Indopacific 13988 0.141       0.95      0.138      0.143            0.138            0.143

#mden43 fis
Sum.mden43.fis.melt = groupwiseMean(fis ~ pop, data= mden43.fis.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.mden43.fis.melt
X = 2
Y = Sum.mden43.fis.melt$Percentile.upper + 0.2
png(filename="./mden43_fis.png", width=500, height=300)
ggplot(Sum.mden43.fis.melt, aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Inbreeding Coefficient") + xlab(NULL)
dev.off()
# pop     n   Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1    Atlantic  7723 0.0668       0.95     0.0613     0.0724           0.0616           0.0723
# 2 Indopacific 10703 0.0682       0.95     0.0631     0.0732           0.0631           0.0731

zczc.atl.wc<-read.table("~/Dropbox/Phd/Bioinformatics/bw_ddrad_global/Analysis Files/Demerelate/zczc/zczc_atl_fis.txt", header=TRUE)
head(zczc.atl.wc)
zczc.atl.wc[,3]<-"Atlantic"
colnames(zczc.atl.wc)<-c("locus", "fis", "pop")

zczc.indo.wc<-read.table("~/Dropbox/Phd/Bioinformatics/bw_ddrad_global/Analysis Files/Demerelate/zczc/zczc_indo_fis.txt", header=TRUE)
head(zczc.indo.wc)
zczc.indo.wc[,3]<-"Indopacific"
colnames(zczc.indo.wc)<-c("locus", "fis", "pop")

zczc.med.wc<-read.table("~/Dropbox/Phd/Bioinformatics/bw_ddrad_global/Analysis Files/Demerelate/zczc/zczc_med_fis.txt", header=TRUE)
head(zczc.med.wc)
zczc.med.wc[,3]<-"Mediterranean"
colnames(zczc.med.wc)<-c("locus", "fis", "pop")

zczc.all.wc<-rbind(zczc.atl.wc, zczc.indo.wc, zczc.med.wc)
head(zczc.all.wc)

Sum.zczc.all.wc = groupwiseMean(fis ~ pop, data= zczc.all.wc, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc.all.wc
X = c(1:4,6-11)
Y = Sum.zczc125.pop1.ho$Percentile.upper + 0.2
png(filename="./zczc125_pop_ho.png", width=500, height=300)
ggplot(Sum.zczc125.pop1.ho[-5,], aes(x = pop,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Observed Heterozygosity") + xlab(NULL)
dev.off()


# calculating confidence intervals for tajima's d data --------------------

#zczc125 tajima's d confidence intervals for ocean basins
zczc125.tajima.atl<-read.table("./zczc125.atl/zczc125.atl.100000.Tajima.D", header=TRUE)
zczc125.tajima.indo<-read.table("./zczc125.indo/zczc125.indo.100000.Tajima.D", header=TRUE)
zczc125.tajima.med<-read.table("./zczc125.med/zczc125.med.100000.Tajima.D", header=TRUE)

zczc125.tajima.atl[,5]<-"Atlantic"
zczc125.tajima.indo[,5]<-"Indopacific"
zczc125.tajima.med[,5]<-"Mediterranean"

zczc125.tajima<-rbind(zczc125.tajima.atl, zczc125.tajima.indo, zczc125.tajima.med)

Sum.zczc125.tajima = groupwiseMean(TajimaD ~ V5, data= zczc125.tajima, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.tajima
X = 3
Y = Sum.zczc125tajima$Percentile.upper + 0.2
png(filename="./zczc125_tajima.png", width=500, height=300)
ggplot(Sum.zczc125.tajima[-5,], aes(x = V5,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Tajima's D") + xlab(NULL)
dev.off()
# V5     n   Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1      Atlantic 20715 -0.223       0.95     -0.234     -0.211           -0.234           -0.212
# 2   Indopacific 18221 -0.221       0.95     -0.234     -0.209           -0.234           -0.208
# 3 Mediterranean  9512  0.253       0.95      0.233      0.274            0.233            0.274

#mden43 tajima's d confidence intervals for ocean basin
mden43.tajima.atl<-read.table("./mden43.atl/mden43.atl.100000.Tajima.D", header=TRUE)
mden43.tajima.indo<-read.table("./mden43.indo/mden43.indo.100000.Tajima.D", header=TRUE)

# mden43.tajima.atl.5<-read.table("./mden43.atl/mden43.atl.5000.Tajima.D", header=TRUE)
# mden43.tajima.indo.5<-read.table("./mden43.indo/mden43.indo.5000.Tajima.D", header=TRUE)

mden43.tajima.atl[,5]<-"Atlantic"
mden43.tajima.indo[,5]<-"Indopacific"
mden43.tajima<-na.omit(rbind(mden43.tajima.atl, mden43.tajima.indo))

# mden43.tajima.atl.5[,5]<-"Atlantic"
# mden43.tajima.indo.5[,5]<-"Indopacific"
# mden43.tajima.5<-na.omit(rbind(mden43.tajima.atl.5, mden43.tajima.indo.5))

Sum.mden43.tajima = groupwiseMean(TajimaD ~ V5, data= mden43.tajima, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.mden43.tajima
X = 2
Y = Sum.mden43.tajima$Percentile.upper + 0.2
png(filename="./mden43_tajima.png", width=500, height=300)
ggplot(Sum.mden43.tajima[-5,], aes(x = V5,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Tajima's D") + xlab(NULL)
dev.off()
# V5    n    Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1    Atlantic 6895 -0.0733       0.95    -0.0964    -0.0503          -0.0973          -0.0501
# 2 Indopacific 9153 -0.4500       0.95    -0.4680    -0.4320          -0.4680          -0.4310

# Sum.mden43.tajima.5 = groupwiseMean(TajimaD ~ V5, data= mden43.tajima.5, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
# Sum.mden43.tajima.5
# X = 2
# Y = Sum.mden43.tajima.5$Percentile.upper + 0.2
# png(filename="./mden43_tajima.5.png", width=500, height=300)
# ggplot(Sum.mden43.tajima.5, aes(x = V5,y = Mean)) +
#   geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
#   geom_point(shape = 15, size  = 4) +
#   theme_bw() +
#   theme(axis.title= element_text(face="bold")) +
#   ylab("Tajima's D") + xlab(NULL)
# dev.off()
# V5     n    Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1    Atlantic  7606 -0.0672       0.95    -0.0892    -0.0451          -0.0896           -0.046
# 2 Indopacific 10517 -0.4300       0.95    -0.4470    -0.4130          -0.4470           -0.413

#confidence intervals for zcav ocean basins
zczc125.pop1.atl.aa<-read.table("./zczc125.Atl_AA.100000.Tajima.D", header=TRUE)
zczc125.pop1.atl.ab<-read.table("./zczc125.Atl_AB.100000.Tajima.D", header=TRUE)
zczc125.pop1.atl.ac<-read.table("./zczc125.Atl_AC.100000.Tajima.D", header=TRUE)
zczc125.pop1.atl.ad<-read.table("./zczc125.Atl_AD.100000.Tajima.D", header=TRUE)
zczc125.pop1.atl.ae<-read.table("./zczc125.Atl_AE.100000.Tajima.D", header=TRUE)
zczc125.pop2.atl.ba<-read.table("./zczc125.Atl_BA.100000.Tajima.D", header=TRUE)
zczc125.pop2.atl.bb<-read.table("./zczc125.Atl_BB.100000.Tajima.D", header=TRUE)
zczc125.pop2.atl.bc<-read.table("./zczc125.Atl_BC.100000.Tajima.D", header=TRUE)
zczc125.pop2.atl.bd<-read.table("./zczc125.Atl_BD.100000.Tajima.D", header=TRUE)
zczc125.pop.indo.a<-read.table("./zczc125.Indo_A.100000.Tajima.D", header=TRUE)
zczc125.pop.indo.b<-read.table("./zczc125.Indo_B.100000.Tajima.D", header=TRUE)
zczc125.pop.indo.c<-read.table("./zczc125.Indo_C.100000.Tajima.D", header=TRUE)
zczc125.pop.med.a<-read.table("./zczc125.Med_A.100000.Tajima.D", header=TRUE)
zczc125.pop.med.b<-read.table("./zczc125.Med_B.100000.Tajima.D", header=TRUE)
zczc125.pop.med.c<-read.table("./zczc125.Med_C.100000.Tajima.D", header=TRUE)

zczc125.pop1.atl.aa[,5]<-"Atl_AA"
zczc125.pop1.atl.ab[,5]<-"Atl_AB"
zczc125.pop1.atl.ac[,5]<-"Atl_AC"
zczc125.pop1.atl.ad[,5]<-"Atl_AD"
zczc125.pop1.atl.ae[,5]<-"Atl_AE"
zczc125.pop2.atl.ba[,5]<-"Atl_BA"
zczc125.pop2.atl.bb[,5]<-"Atl_BB"
zczc125.pop2.atl.bc[,5]<-"Atl_BC"
zczc125.pop2.atl.bd[,5]<-"Atl_BD"
 zczc125.pop.indo.a[,5]<-"Indo_A"
 zczc125.pop.indo.b[,5]<-"Indo_B"
 zczc125.pop.indo.c[,5]<-"Indo_C"
  zczc125.pop.med.a[,5]<-"Med_A"
  zczc125.pop.med.b[,5]<-"Med_B"
  zczc125.pop.med.c[,5]<-"Med_C"

zczc125.pop1.tajima<-na.omit(rbind(zczc125.pop1.atl.aa,
                                   zczc125.pop1.atl.ab,
                                   zczc125.pop1.atl.ac,
                                   zczc125.pop1.atl.ad,
                                   zczc125.pop1.atl.ae,
                                   zczc125.pop.indo.a,
                                   zczc125.pop.indo.b,
                                   zczc125.pop.indo.c,
                                   zczc125.pop.med.a,
                                     zczc125.pop.med.b,
                                     zczc125.pop.med.c))
zczc125.pop2.tajima<-na.omit(rbind(zczc125.pop2.atl.ba,
                                   zczc125.pop2.atl.bb,
                                   zczc125.pop2.atl.bc,
                                   zczc125.pop2.atl.bd,
                                   zczc125.pop.indo.a,
                                     zczc125.pop.indo.b,
                                     zczc125.pop.indo.c,
                                     zczc125.pop.med.a,
                                     zczc125.pop.med.b,
                                     zczc125.pop.med.c))

Sum.zczc125.pop1.tajima = groupwiseMean(TajimaD ~ V5, data= zczc125.pop1.tajima, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.pop1.tajima
X = 11
Y = zczc125.pop1.tajima$Percentile.upper + 0.2
png(filename="./zczc125_pop1_tajima.png", width=500, height=300)
ggplot(Sum.zczc125.pop1.tajima[-5,], aes(x = V5,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower, ymax = Percentile.upper), width = .05, size  = .5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Tajima's D") + xlab(NULL)
dev.off()
# V5     n    Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1  Atl_AA 13305 -0.2050       0.95   -0.22100    -0.1890         -0.22100          -0.1890
# 2  Atl_AB  9776 -0.1670       0.95   -0.18500    -0.1480         -0.18600          -0.1470
# 3  Atl_AC 15446 -0.2380       0.95   -0.25200    -0.2240         -0.25100          -0.2240
# 4  Atl_AD 16742 -0.2800       0.95   -0.29400    -0.2670         -0.29300          -0.2660
# 5  Atl_AE  3151  1.6000       0.95    1.59000     1.6100          1.59000           1.6100
# 6  Indo_A  8901 -0.0885       0.95   -0.10800    -0.0686         -0.10800          -0.0686
# 7  Indo_B 13759 -0.2500       0.95   -0.26600    -0.2350         -0.26600          -0.2330
# 8  Indo_C 15191 -0.2090       0.95   -0.22400    -0.1950         -0.22400          -0.1950
# 9   Med_A  7642  0.4190       0.95    0.39700     0.4400          0.39800           0.4400
# 10  Med_B  6083  0.3700       0.95    0.34600     0.3940          0.34600           0.3950
# 11  Med_C  7121  0.0271       0.95    0.00407     0.0502          0.00389           0.0501

Sum.zczc125.pop2.tajima = groupwiseMean(TajimaD ~ V5, data= zczc125.pop2.tajima, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.pop2.tajima
X = 11
Y = zczc125.pop2.tajima$Percentile.upper + 0.2
png(filename="./zczc125_pop2_tajima.png", width=500, height=300)
ggplot(Sum.zczc125.pop2.tajima[-4,], aes(x = V5,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower, ymax = Percentile.upper), width = .05, size  = .5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Tajima's D") + xlab(NULL)
dev.off()
# V5     n    Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1  Atl_BA 15347 -0.2400       0.95   -0.25500    -0.2260         -0.25500          -0.2250
# 2  Atl_BB 15446 -0.2380       0.95   -0.25200    -0.2240         -0.25100          -0.2230
# 3  Atl_BC 16742 -0.2800       0.95   -0.29400    -0.2670         -0.29300          -0.2670
# 4  Atl_BD  3151  1.6000       0.95    1.59000     1.6100          1.59000           1.6100
# 5  Indo_A  8901 -0.0885       0.95   -0.10800    -0.0686         -0.11000          -0.0682
# 6  Indo_B 13759 -0.2500       0.95   -0.26600    -0.2350         -0.26600          -0.2340
# 7  Indo_C 15191 -0.2090       0.95   -0.22400    -0.1950         -0.22500          -0.1950
# 8   Med_A  7642  0.4190       0.95    0.39700     0.4400          0.39800           0.4410
# 9   Med_B  6083  0.3700       0.95    0.34600     0.3940          0.34700           0.3940
# 10  Med_C  7121  0.0271       0.95    0.00407     0.0502          0.00311           0.0492

#confidence intervals for mden ocean basins
mden43.atl.bah<-read.table("./mden43.atl-bahamas.100000.Tajima.D", header=TRUE)
mden43.atl.eas<-read.table("./mden43.atl-east.100000.Tajima.D", header=TRUE)
mden43.atl.oth<-read.table("./mden43.atl-other.100000.Tajima.D", header=TRUE)
mden43.indo.afr<-read.table("./mden43.indo-africa.100000.Tajima.D", header=TRUE)
mden43.indo.haw<-read.table("./mden43.indo-hawaii.100000.Tajima.D", header=TRUE)
mden43.indo.sou<-read.table("./mden43.indo-south.100000.Tajima.D", header=TRUE)

mden43.atl.bah[,5]<-"Atl_Bahamas"
mden43.atl.eas[,5]<-"Atl_East"
mden43.atl.oth[,5]<-"Atl_Other"
mden43.indo.afr[,5]<-"Indo_Africa"
mden43.indo.haw[,5]<-"Indo_Hawaii"
mden43.indo.sou[,5]<-"Indo_South"

mden43.pop.tajima<-na.omit(rbind(mden43.atl.bah, mden43.atl.eas, mden43.atl.oth, mden43.indo.afr, mden43.indo.haw, mden43.indo.sou))

Sum.mden43.pop.tajima = groupwiseMean(TajimaD ~ V5, data= mden43.pop.tajima, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.mden43.pop.tajima
X = 6
Y = Sum.mden43.pop.tajima$Percentile.upper + 0.2
png(filename="./mden43_pop_tajima.png", width=500, height=300)
ggplot(Sum.mden43.pop.tajima, aes(x = V5,y = Mean)) +
  geom_errorbar(aes(ymin = Percentile.lower,ymax = Percentile.upper), width = 0.05, size  = 0.5) +
  geom_point(shape = 15, size  = 4) +
  theme_bw() +
  theme(axis.title= element_text(face="bold")) +
  ylab("Tajima's D") + xlab(NULL)
dev.off()
# V5    n    Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
# 1 Atl_Bahamas 4405  0.0375       0.95    0.00858     0.0663          0.00696           0.0648
# 2    Atl_East 5702  0.0149       0.95   -0.01060     0.0404         -0.01180           0.0405
# 3   Atl_Other 4135 -0.0166       0.95   -0.04620     0.0131         -0.04600           0.0136
# 4 Indo_Africa 5625 -0.2850       0.95   -0.30900    -0.2600         -0.30900          -0.2600
# 5 Indo_Hawaii 5803 -0.2630       0.95   -0.28700    -0.2380         -0.28700          -0.2380
# 6  Indo_South 4902 -0.2170       0.95   -0.24400    -0.1900         -0.24500          -0.1910
# dapc --------------------------------------------------------------------
dapc1<-dapc(zczc125.gl)
scatter(dapc1)

dapc2<-dapc(zczc125.pop1.gl)
scatter(dapc2)

zczc125.pop1.nospain.gl<-gl.drop.pop(zczc125.pop1.gl, pop.list="Atl_AE")
zczc125.pop2.nospain.gl<-gl.drop.pop(zczc125.pop2.gl, pop.list="Atl_BD")
mden43.pop.nona.gl<-gl.drop.pop(mden43.pop.gl, pop.list="NA")
library(dartR)

dapc3<-dapc(zczc125.pop1.nospain.gl)
scatter(dapc3)

dapc4<-dapc(zczc125.pop2.nospain.gl)
scatter(dapc4)

dapc5<-dapc(mden43.gl)
scatter(dapc5)

dapc5<-dapc(mden43.pop.nona.gl)
scatter(dapc5)

#using optimisation
x<-mden43.gl
mat<-tab(x)
grp<-pop(x)
xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
xval
scatter(xval$DAPC)
assignplot(xval$DAPC) #shows both the membership probability (red rectangle) and the pre-assigned group membership (blue cross)

mden43.pop.gl
mat<-tab(mden43.pop.gl)
grp<-pop(mden43.pop.gl)
mden43.pop.xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
mden43.pop.xval
scatter(mden43.pop.xval$DAPC)
assignplot(mden43.pop.xval$DAPC) #shows both the membership probability (red rectangle) and the pre-assigned group membership (blue cross)

zczc125.gl
mat<-tab(zczc125.gl)
grp<-pop(zczc125.gl)
zczc125.xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
zczc125.xval
scatter(zczc125.xval$DAPC)
assignplot(zczc125.xval$DAPC) #shows both the membership probability (red rectangle) and the pre-assigned group membership (blue cross)

zczc125.pop1.gl
mat<-tab(zczc125.pop1.gl)
grp<-pop(zczc125.pop1.gl)
zczc125.pop1.xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
zczc125.pop1.xval
scatter(zczc125.pop1.xval$DAPC)
assignplot(zczc125.pop1.xval$DAPC) #shows both the membership probability (red rectangle) and the pre-assigned group membership (blue cross)

zczc125.pop2.gl
mat<-tab(zczc125.pop2.gl)
grp<-pop(zczc125.pop2.gl)
zczc125.pop2.xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
zczc125.pop2.xval
scatter(zczc125.pop2.xval$DAPC)
assignplot(zczc125.pop2.xval$DAPC) #shows both the membership probability (red rectangle) and the pre-assigned group membership (blue cross)

#2 spanish individuals messing up dapc. will remove and see how it looks
zczc125.pop1.gl<-gl.compliance.check(zczc125.pop1.gl)
zczc125.pop2.gl<-gl.compliance.check(zczc125.pop2.gl)

zcav.pop1.nosp.gl<-gl.drop.pop(zczc125.pop1.gl, "Atl_AE")
zcav.pop2.nosp.gl<-gl.drop.pop(zczc125.pop2.gl, "Atl_BD")

zcav.pop1.nosp.gl
mat<-tab(zcav.pop1.nosp.gl)
grp<-pop(zcav.pop1.nosp.gl)
zcav.pop1.nosp.xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
zcav.pop1.nosp.xval
scatter(zcav.pop1.nosp.xval$DAPC)
assignplot(zcav.pop1.nosp.xval$DAPC) #shows both the membership probability (red rectangle) and the pre-assigned group membership (blue cross)

zcav.pop2.nosp.gl
mat<-tab(zcav.pop2.nosp.gl)
grp<-pop(zcav.pop2.nosp.gl)
zcav.pop2.nosp.xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
zcav.pop2.nosp.xval
scatter(zcav.pop2.nosp.xval$DAPC)
assignplot(zcav.pop2.nosp.xval$DAPC) #shows both the membership probability (red rectangle) and the pre-assigned group membership (blue cross)

# Tess Results/Biogeography Maps ---------------------------------------------------------------

#data files
zczc125.coord #V2=long, V3=lat
tess3.zczc125
tess3.zczc125.atl
tess3.zczc125.indo
tess3.zczc125.med

mden43.coord
tess.mden43


qmatrix(tess.mden43, K=2)
for(i in 2:k) {
  Q.matrix<-qmatrix(tess.mden43, K = i)
  filename=paste0("Mden43.Q.matrix.k=",i,".csv")
  write.table(Q.matrix, file = filename)}

tess.mden43

# ## get the NOAA map (change lon and lat coordinates and resolution, see help)
# # use antimeridian = TRUE to center on antimeridian (Pacific Ocean)
# map.bathy <- marmap::getNOAA.bathy(lon1=-180, lon2=180, lat1= -90, lat2= 90, res = 10, keep=TRUE, antimeridian = FALSE)
# map.bathy2 <- marmap::getNOAA.bathy(lon1=-180, lon2=180, lat1= -90, lat2= 90, res = 10,keep=TRUE, antimeridian=TRUE)
# 
# # change sign (I think this inverts land/water color surface)
# map.bathy1 <-  - map.bathy
# plot(map.bathy1)
# 
# map.bathy2<--map.bathy2
# plot(map.bathy2, xlim=c(-50,300))
# 
# # convert bathy to raster (package = raster)
# asc.raster <- marmap::as.raster(map.bathy1)
# plot(asc.raster)
# 
# 
# #rewrite the modified raster in your working directory
# raster::writeRaster(asc.raster, "myraster.asc", overwrite=TRUE)
# raster::writeRaster(asc.raster2, "myraseter2.asc", overwrite=TRUE)
# 
# #check maps
# zczc125.coord #V2=long, V3=lat
# plot(zczc125.coord, pch = 19, col="blue", cex = .8,  
#      xlab = "Longitude (°E)", ylab = "Latitude (°N)", ylim=c(-90,90), main="Sample Distribution, Z. cavirostris (Z. cavirostris genome)")
# map(add = T)
# 
# mden43.coord #V2=long, V3=lat
# plot(mden43.coord, pch = 19, col="blue", cex = .8,  
#      xlab = "Longitude (°E)", ylab = "Latitude (°N)", ylim=c(-90,90), main="Sample Distribution, M. densirostris (M. bidens genome)")
# map(add = T)
# 
# q.matrix <- qmatrix(tess3.zczc125, K = 3)
# 
# plot(q.matrix, zczc125.coord, method = "map.max", cex = .5, raster.filename = "myraster.asc", interpol = FieldsKrigModel(10), main = "Ancestry coefficients, k=3", resolution = c(300, 300), col.palette = my.palette2, xlab = "Longitude", ylab = "Latitude", xlim=c(-180,180), ylim=c(-90,90))  
# 
# #convert the longitude data to 360°.
# coord_map <- as.data.frame(zczc125.coord)
# lon360 <- ifelse(coord_map$V2 < 0, 360 + coord_map$V2, coord_map$V2)
# coord_360 <- cbind(coord_map, lon360)
# coord_360$V2 <- NULL
# coord_360r <- coord_360[,c("lon360", "V3")]
# 
# plot(q.matrix, coord_360r, method = "map.max", cex = .5, raster.filename = "myraster.asc", interpol = FieldsKrigModel(10), main = "Ancestry coefficients, k=3", resolution = c(300, 300), col.palette = my.palette2, xlab = "Longitude", ylab = "Latitude") 
# 
# 
# library(ggplot2)
# library(rworldmap)
# 
# map.polygon <- getMap(resolution = "low")
# map.polygon<-map.bathy1
# ggtess3Q(q.matrix, zczc125.coord, raster.filename="myraster.asc", background=TRUE)
# 
# 
# plot(asc.raster)

sample.loc<-read.csv("./SampleLocationsPops.csv", head=TRUE)
head(sample.loc)
zc.loc<-subset(sample.loc, sample.loc$species=="zc")
zc.med.loc<-subset(zc.loc, zc.loc$ocean=="Mediterranean")
md.loc<-subset(sample.loc, sample.loc$species=="md")

#check locaitons on simple maps
#basic world map
wm<-map_data("world")
wm<-ggplot(wm, aes(long, lat, group=group)) +
  geom_polygon(fill="grey", colour="white") +
  coord_fixed(1.3) +
  theme_classic()
wm

med.map<-ggplot(wm, aes(long, lat, group=group)) +
  geom_polygon(fill="grey", colour="white") +
  coord_map(xlim = c(-5, 36),ylim = c(30, 45)) +
  theme_classic()
med.map

# png("zc.pop.map.png", width=6.75, height=4, units="in", res=300)
# dev.off()

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n=12, name="Paired")
paired.col<-brewer.pal(12, name="Paired")
paired.col
paired.col<-c("#A6CEE3" ,"#1F78B4" ,"#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#6A3D9A","#B15928")

png("zc.pop.map.png", width=6.75, height=4, units="in", res=300)
zc.pop.wm<-wm +
  geom_jitter(data=zc.loc, position=position_jitter(width = 1, height = 1),
              aes(x=long, y=lat, group=pop, colour=pop), size=1.5) +
  scale_color_manual(values=paired.col) +
  labs(colour=NULL, x="Longitude", y="Latitude") +
  theme(legend.position=c(1.02,0.58)) +
  theme(legend.background=element_blank()) +
  theme(legend.text=element_text(size=8)) +
  theme(plot.margin=unit(c(.5,1.5,.5,.5),"cm")) 
zc.pop.wm
dev.off()                   
                   
      
png("zc.med.map.png", width=6.75, height=4, units="in", res=300)
zc.med.pop.wm<-med.map +
  geom_jitter(data=zc.med.loc, position=position_jitter(width = 1, height = .1),
              aes(x=long, y=lat, group=pop, colour=pop), size=2.5) +
  scale_color_manual(values=paired.col[8:10])+
  labs(colour=NULL, x="Longitude", y="Latitude")
zc.med.pop.wm
dev.off()


png("md.pop.map.png", width=6.75, height=4, units="in", res=300)
md.pop.wm<-wm +
  geom_jitter(data=md.loc, position=position_jitter(width = 1, height = 1),
              aes(x=long, y=lat, group=pop, colour=pop), size=1.5) +
  scale_color_brewer(palette = "Paired") +
  labs(colour=NULL, x="Longitude", y="Latitude") +
  theme(legend.position=c(1.02,0.6)) +
  theme(legend.background=element_blank()) +
  theme(legend.text=element_text(size=8)) +
  theme(plot.margin=unit(c(.5,1.5,.5,.5),"cm"))
md.pop.wm
dev.off()

tess3.zczc125

plot(tess3.zczc125, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     main = "Cross-Entropy Score, Z. cavirostris",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score", type="l")
plot(tess3.zczc125.atl, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     main = "Cross-Entropy Score, Z. cavirostris, Atlantic",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score", type="l")
plot(tess3.zczc125.indo, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     main = "Cross-Entropy Score, Z. cavirostris, Indo-Pacific",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score", type="l")
plot(tess3.zczc125.med, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     main = "Cross-Entropy Score, Z. cavirostris, Mediterranean",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score", type="l")

plot(tess.mden43,crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     main = "Cross-Entropy Score, M. densirostris",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score3", type="l")
plot(tess3.mden.atl,crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     main = "Cross-Entropy Score, M. densirostris, Atlantic",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score", type="l")
plot(tess3.mden.indo, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     main = "Cross-Entropy Score, M. densirostris, Indo-Pacific",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score", type="l")


dim(mden43.geno)

dim(mden43.geno)
mden43.geno[1:3,1:3]

mden43.geno.atl<-mden43.geno[1:28,]
mden43.geno.indo<-mden43.geno[29:43,]

mden43.geno.atl[,1]
mden43.geno.atl[mden43.geno.atl=="-1"]<-NA

mden43.geno.indo[,1]
mden43.geno.indo[mden43.geno.indo=="-1"]<-NA

mden43.map
mden43.coord.atl<-subset(mden43.map, mden43.map$ocean=="Atlantic")
mden43.coord.atl<-mden43.coord.atl[,c(2,4,3)]
mden43.coord.atl<-mden43.coord.atl[,c(2,3)]

mden43.coord.indo<-subset(mden43.map, mden43.map$ocean=="Indopacific")
mden43.coord.indo<-mden43.coord.indo[,c(4,3)]

mden43.geno.atl
mden43.coord.atl<-as.matrix(mden43.coord.atl)

k<-10
tess3.mden.atl<-tess3(X=mden43.geno.atl, coord=mden43.coord.atl, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)
# rep = The number of times the algorithm will be repeated for each value of K ############(recommended = 10)
# max.iteration	= the maximum number of iterations of the optimization algorithm. ##########(recommend = 200 (default))
# keep = if "best", only the result with the lowest rmse score will be kept for each value of K. If "all", all results will be kept and returned for each value of K. The second option uses more space in memory.
# mask =  this numeric value is the proportion of masked data when computing the cross-validation criterion (default = 0).

#The plot function generates a plot for root mean-squared errors computed on a subset of loci used for cross-validation.
plot(tess3.mden.atl, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score",
     main="Cross validation scores, M. densirostris")
# specify crossentropy with error bars (appears to be very similar to "plot" above)
plot(tess3.mden.atl, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score", type="l")


mden43.geno.indo
mden43.coord.indo<-as.matrix(mden43.coord.indo)

k<-10
tess3.mden.indo<-tess3(X=mden43.geno.indo, coord=mden43.coord.indo, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)
# rep = The number of times the algorithm will be repeated for each value of K ############(recommended = 10)
# max.iteration	= the maximum number of iterations of the optimization algorithm. ##########(recommend = 200 (default))
# keep = if "best", only the result with the lowest rmse score will be kept for each value of K. If "all", all results will be kept and returned for each value of K. The second option uses more space in memory.
# mask =  this numeric value is the proportion of masked data when computing the cross-validation criterion (default = 0).

#The plot function generates a plot for root mean-squared errors computed on a subset of loci used for cross-validation.
plot(tess3.mden.indo, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score",
     main="Cross validation scores, M. densirostris")
# specify crossentropy with error bars (appears to be very similar to "plot" above)
plot(tess3.mden.indo, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score", type="l")

# AMOVA using Poppr -------------------------------------------------------------------


zczc125.gl$pop
zczc125.pop2.gl$pop

mden43.gl
mden43.pop.gl

zczc125.strata.gl<-zczc125.gl
mden43.strata.gl<-mden43.gl

#make a new data frame with strata info
zczc125.sample<-zczc125.gl$ind.names
zczc125.pop<-pop(zczc125.gl)
zczc125.subpop<-pop(zczc125.pop2.gl)

mden43.sample<-mden43.gl$ind.names
mden43.pop<-pop(mden43.gl)
mden43.subpop<-pop(mden43.pop.gl)

zczc125.sample<-as.data.frame(zczc125.sample)
zczc125.pop<-as.data.frame(zczc125.pop)
zczc125.subpop<-as.data.frame(zczc125.subpop)

mden43.sample<-as.data.frame(mden43.sample)
mden43.pop<-as.data.frame(mden43.pop)
mden43.subpop<-as.data.frame(mden43.subpop)

zczc125.strata<-cbind(zczc125.sample, zczc125.pop, zczc125.subpop)
zczc125.strata[,1]
zczc125.strata<-zczc125.strata[,2:3] #table with one column as ocean basin and the other as population
colnames(zczc125.strata)<-c("ocean", "pop")
head(zczc125.strata)

mden43.strata<-cbind(mden43.sample, mden43.pop, mden43.subpop)
mden43.strata
mden43.strata<-mden43.strata[,2:3]
colnames(mden43.strata)<-c("ocean", "pop")
head(mden43.strata)

#set population stratifiation in Adegenet
strata(zczc125.strata.gl)<-zczc125.strata
zczc125.strata.gl

strata(mden43.strata.gl)<-mden43.strata
mden43.strata.gl

#set hierarchy
strata(zczc125.strata.gl, ~ocean/pop)
strata(mden43.strata.gl, ~ocean/pop)

zczc125.amova<-poppr.amova(zczc125.strata.gl, hier=~ocean/pop, missing="asis")
zczc125.amova
zczc125.amova.sig<-randtest(zczc125.amova, nrepet = 999)
plot(zczc125.amova.sig)
zczc125.amova.sig

mden43.amova<-poppr.amova(mden43.strata.gl, hier=~ocean/pop, missing="asis")
mden43.amova
mden43.amova.sig<-randtest(mden43.amova, nrepet = 999)
plot(mden43.amova.sig)
mden43.amova.sig

zczc125.strata.gl
mden43.strata.gl

zczc125.strata.gi<-gl2gi(zczc125.strata.gl)
mden43.strata.gi<-gl2gi(mden43.strata.gl)

zczc125.strata.gt<-genind2gtypes(zczc125.strata.gi)
mden43.strata.gt<-genind2gtypes(mden43.strata.gi)

zczc125.strata.gt
mden43.strata.gt

rownames(zczc125.strata)<-indNames(zczc125.gl)
zczc125.strata<-cbind(zczc125.sample, zczc125.strata)
colnames(zczc125.strata)<-c("id", "ocean", "pop")
zczc125.strata

rownames(mden43.strata)<-indNames(mden43.gl)
mden43.strata<-cbind(mden43.sample, mden43.strata)
colnames(mden43.strata)<-c("id", "ocean", "pop")
mden43.strata

setSchemes(zczc125.strata.gt)<-zczc125.strata
zczc125.strata.gt@schemes
setSchemes(mden43.strata.gt)<-mden43.strata
mden43.strata.gt@schemes

arlequinWrite(zczc125.strata.gt, file="./zczc125")
arlequinWrite(mden43.strata.gt, file="./mden43")

zczc125.df<-as.data.frame(zczc125.gl)
dim(zczc125.df)
rownames(zczc125.df)

zczc125.df<-cbind(zczc125.strata, zczc125.df)
zczc125.df[1:5,1:5]

write.csv(zczc125.df, file="zczc125.df")

zczc125.pop2.gt

zczc125.pop2.gi
zczc125.pop.gt<-genind2gtypes(zczc125.pop2.gi)

arlequinWrite(zczc125.pop.gt, file="./zczc125.pop")

# Mitogenomes -------------------------------------------------------------

mito<-as.data.frame(read_xlsx("./mitogenome_map.xlsx", col_names = TRUE))
mito
mito$clade<-as.character(mito$clade)

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n=10, name="Paired")
paired.col<-brewer.pal(6, "Paired")
paired.col

png("zc.pop.map.png", width=6.75, height=4, units="in", res=300)

mito.ocean.wm<-wm +
  geom_jitter(data=mito, position=position_jitter(width = 1, height = 1),
              aes(x=long, y=lat, group=clade, colour=clade), size=1.5) +
  scale_color_brewer(palette = "Paired") +
  labs(colour=NULL, x="Longitude", y="Latitude") +
  theme(legend.position=c(1.02,0.58)) +
  theme(legend.background=element_blank()) +
  theme(legend.text=element_text(size=8)) +
  theme(plot.margin=unit(c(.5,1.5,.5,.5),"cm")) 
mito.ocean.wm
dev.off()

mito.pop.wm<-wm +
  geom_jitter(data=mito, position=position_jitter(width = 1, height = 1),
              aes(x=long, y=lat, group=subclade, colour=subclade), size=1.5) +
  scale_color_brewer(palette = "Paired") +
  labs(colour=NULL, x="Longitude", y="Latitude") +
  theme(legend.position=c(1.02,0.58)) +
  theme(legend.background=element_blank()) +
  theme(legend.text=element_text(size=8)) +
  theme(plot.margin=unit(c(.5,1.5,.5,.5),"cm")) 
mito.pop.wm
dev.off()


# adegenet trees with bootstraps; GOOD----------------------------------------------------------

zczc125.tree<-bionjs(dist(as.matrix(zczc125.gl)))
plot(zczc125.tree)
ggtree(zczc125.tree)

zczc.aboot<-aboot(zczc125.gi, tree=bionj, tip.label=NULL)
zczc.pop1.aboot<-aboot(zczc125.pop1.gi, tree=bionj)
mden.aboot<-aboot(mden43.gi, tree=bionj, tip.label=NULL)

#add ocean as tip label instead of individual
zczc.aboot$node.label
zczc.aboot$tip.label<-zczc125.gl$pop
zczc.aboot$tip.label<-as.character(zczc.aboot$tip.label)

zczc.pop1.aboot$tip.label<-zczc125.pop1.gl$pop
zczc.pop1.aboot$tip.label<-as.character(zczc.pop1.aboot$tip.label)

mden.aboot$node.label
mden.aboot$tip.label
mden.aboot$tip.label<-mden43.gl$pop
mden.aboot$tip.label<-as.character(mden.aboot$tip.label)

#plot tree
zczc.tibble
zczc.aboot.tibble2<-full_join(zczc.aboot.tib, d, by='label' )
zczc.aboot.tibble2
zczc.aboot.tibble3<-zczc.aboot.tibble2 %>% mutate(pop=recode(pop, "Indopacific"="Indo-Pacific"))
zczc.aboot.tibble3$pop
zczc.aboot.tree<-as.treedata(zczc.aboot.tibble3)
zczc.aboot.tree

zczc.pop1.aboot.tib<-as_tibble(zczc.pop1.aboot)
zczc.pop1.aboot.tibble2<-full_join(zczc.pop1.aboot.tib, d, by='label' )
zczc.pop1.aboot.tibble2
zczc.pop1.aboot.tibble2$pop
zczc.pop1.aboot.tree<-as.treedata(zczc.pop1.aboot.tibble2)
zczc.pop1.aboot.tree


mden.aboot.tibble<-as.tibble(mden.aboot)
mden.aboot.tibble
mden.aboot.tibble2<-full_join(mden.aboot.tibble, e, by='label')
mden.aboot.tibble2
mden.aboot.tibble3<-mden.aboot.tibble2 %>% mutate(pop=recode(pop, "Indopacific"="Indo-Pacific"))

mden.aboot.tree<-as.treedata(mden.aboot.tibble3)
mden.aboot.tree

#plot cladogram tree with tips coloured by ocean basin and add node labels (bootstrap values)
png("./zczc.boot.png", width=5, height=9, units="in", res=600)
ggtree(zczc.aboot.tree, branch.length = "none") +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=my_zissou2) +
  geom_label2(aes(label=label, subset=!is.na(as.numeric(label))), fill="white", size=1.5, label.padding = unit(0.1, "lines"), label.r = unit(0, "lines")) +
  theme(legend.position=c(.2,.9), legend.title = element_blank())
dev.off()

#plot phylogram tree with tips coloured by ocean basin and add node labels (bootstrap values)
png("./zczc.boot.phylo.png", width=5, height=9, units="in", res=1000)
ggtree(zczc.aboot.tree) +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=c("#7b848f", "#c85200", "#5fa2ce")) +
  geom_label2(aes(label=label, subset=as.numeric(label)>50), fill='white', size=2, label.padding = unit(0.1, "lines"), label.r = unit(0, "lines"), nudge_x=-.001, show.legend=NULL) +
  theme_tree2() +
  theme(legend.position=c(.9,.9), legend.title = element_blank())
dev.off()

+
  theme_tree2()+
  theme(legend.position=c(.9,.9), legend.title = element_blank())


png("./zczc.tree.august.png", width=5, height=10, units="in", res=1000)
plot(zczc.aboot.tree@phylo, cex=.5)
dev.off()
zczc.aboot.tree@info

png("./mden.boot.png", width=4, height=6, units="in", res=600)
m<-ggtree(mden.aboot.tree, branch.length = "none") +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=my_zissou2) +
  geom_label2(aes(label=label, subset=!is.na(as.numeric(label))), fill="white", size=2.5, label.padding = unit(0.1, "lines"), label.r = unit(0, "lines")) +
  theme(legend.position=c(.2,.9), legend.title = element_blank())
gridExtra::grid.arrange(rotate(m, 48) %>% rotate(49) %>% rotate(63))
dev.off()

png("./mden.boot.phylo.png", width=5, height=9, units="in", res=1000)
m<-ggtree(mden.aboot.tree) +
  geom_tippoint(aes(col=factor(pop)), size=3) +
  scale_colour_manual(values=c("#7b848f", "#c85200", "#5fa2ce")) +
  geom_label2(aes(label=label, subset=as.numeric(label)>50), fill='white', size=3, label.padding = unit(0.1, "lines"), label.r = unit(0, "lines"), nudge_x=-.001, show.legend=NULL) +
  theme_tree2() +
  theme(legend.position=c(.9,.5), legend.title = element_blank())
m
gridExtra::grid.arrange(rotate(m, 48) %>% rotate(49) %>% rotate(63))
dev.off()

#fix colours so that Red=Atlantic, Yellow=Indo-Pacific, Blue=Mediterranean

library(wesanderson)
wes_palette("Zissou1")
wes_palette("Darjeeling1")

my_zissou<-wes_palette("Zissou1")[c(2,3,5,1,4)]
my_zissou2<-wes_palette("Zissou1")[c(5,3,1)]
wes_palette("Zissou1")[c(2,3,5,1,4)]

my_zissou

table(zczc125.gl@n.loc)

zczc125.hfstat$X1005_37

#fix colours to be colorblind friendly
tableau_color_pal(palette="Color Blind")

# assigner for fis --------------------------------------------------------
library(assigner)
library(radiator)

#turn genlight file into a tidy dataframe
zczc125.tidy<-radiator::genomic_converter(zczc125.gl, output="tidy")
zczc125.tidy<-radiator::tidy_genlight(zczc125.gl, gds = FALSE)
#calculate w and c fis 

# STRATAG for dA ----------------------------------------------------------
library(strataG)
data(dloop.g)
dloop.g
zcav.arp<-arlequinRead("./NEW_282bp_6Pops.arp")
zcav.g<-arp2gtypes(zcav.arp, avoid.dups = TRUE)
zcav.g

zcav.da<-nucleotideDivergence(zcav.g)
zcav.da

zcav.nei.da<-neiDa(zcav.g)
zcav.nei.da

zcav.g.fst<-popStructTest(zcav.g, nrep=1000, stats="fst", type="both")
zcav.g.fst

zczc125.gt<-genind2gtypes(zczc125.gi)
zczc.d<-popStructTest(zczc125.gt, nrep=1000, stats="d", type="pairwise")
#Population structure results:
#                                        D            D.p.val
#Atlantic (55) v. Indopacific (36)      0.0003194935 0.000999001
#Atlantic (55) v. Mediterranean (34)    0.0040695368 0.000999001
#Indopacific (36) v. Mediterranean (34) 0.0044979618 0.000999001

zczc125.pop1.gt<-genind2gtypes(zczc125.pop1.gi)
zczc125.pop2.gt<-genind2gtypes(zczc125.pop2.gi)
zczcpop1.d<-popStructTest(zczc125.pop1.gt, nrep=1000, stats="d", type="pairwise")
#Population structure results:
#  D     D.p.val
#Atl_AA (11) v. Atl_AB (5)  0.0006353102 0.033966034
#Atl_AA (11) v. Atl_AC (17) 0.0004472821 0.000999001
#Atl_AA (11) v. Atl_AD (20) 0.0003751885 0.000999001
#Atl_AA (11) v. Atl_AE (2)  0.0045541756 0.011988012
#Atl_AA (11) v. Indo_A (5)  0.0011048012 0.000999001
#Atl_AA (11) v. Indo_B (12) 0.0006254962 0.000999001
#Atl_AA (11) v. Indo_C (19) 0.0006973782 0.000999001
#Atl_AA (11) v. Med_A (19)  0.0049972329 0.000999001
#Atl_AA (11) v. Med_B (10)  0.0055000550 0.000999001
#Atl_AA (11) v. Med_C (5)   0.0036761390 0.001998002
#Atl_AB (5) v. Atl_AC (17)  0.0006671647 0.012987013
#Atl_AB (5) v. Atl_AD (20)  0.0005944834 0.003996004
#Atl_AB (5) v. Atl_AE (2)   0.0045887923 0.043956044
#Atl_AB (5) v. Indo_A (5)   0.0011760515 0.007992008
#Atl_AB (5) v. Indo_B (12)  0.0007648280 0.000999001
#Atl_AB (5) v. Indo_C (19)  0.0008508739 0.000999001
#Atl_AB (5) v. Med_A (19)   0.0029369413 0.000999001
#Atl_AB (5) v. Med_B (10)   0.0029433423 0.005994006
#Atl_AB (5) v. Med_C (5)    0.0018913952 0.139860140
#Atl_AC (17) v. Atl_AD (20) 0.0002638527 0.000999001
#Atl_AC (17) v. Atl_AE (2)  0.0044960721 0.005994006
#Atl_AC (17) v. Indo_A (5)  0.0009705584 0.000999001
#Atl_AC (17) v. Indo_B (12) 0.0005004311 0.000999001
#Atl_AC (17) v. Indo_C (19) 0.0005443222 0.000999001
#Atl_AC (17) v. Med_A (19)  0.0046612693 0.000999001
#Atl_AC (17) v. Med_B (10)  0.0050977498 0.000999001
#Atl_AC (17) v. Med_C (5)   0.0033936088 0.001998002
#Atl_AD (20) v. Atl_AE (2)  0.0042624120 0.004995005
#Atl_AD (20) v. Indo_A (5)  0.0008863707 0.000999001
#Atl_AD (20) v. Indo_B (12) 0.0003726065 0.000999001
#Atl_AD (20) v. Indo_C (19) 0.0004592980 0.000999001
#Atl_AD (20) v. Med_A (19)  0.0045357318 0.000999001
#Atl_AD (20) v. Med_B (10)  0.0049622268 0.000999001
#Atl_AD (20) v. Med_C (5)   0.0032970732 0.000999001
#Atl_AE (2) v. Indo_A (5)   0.0053311068 0.044955045
#Atl_AE (2) v. Indo_B (12)  0.0046072515 0.009990010
#Atl_AE (2) v. Indo_C (19)  0.0050478872 0.003996004
#Atl_AE (2) v. Med_A (19)   0.0102282213 0.005994006
#Atl_AE (2) v. Med_B (10)   0.0108374479 0.010989011
#Atl_AE (2) v. Med_C (5)    0.0087235751 0.062937063
#Indo_A (5) v. Indo_B (12)  0.0007089331 0.000999001
#Indo_A (5) v. Indo_C (19)  0.0006463077 0.000999001
#Indo_A (5) v. Med_A (19)   0.0051637689 0.000999001
#Indo_A (5) v. Med_B (10)   0.0055861751 0.002997003
#Indo_A (5) v. Med_C (5)    0.0040281371 0.007992008
#Indo_B (12) v. Indo_C (19) 0.0003212235 0.000999001
#Indo_B (12) v. Med_A (19)  0.0046225813 0.000999001
#Indo_B (12) v. Med_B (10)  0.0050534787 0.000999001
#Indo_B (12) v. Med_C (5)   0.0034624673 0.000999001
#Indo_C (19) v. Med_A (19)  0.0050349933 0.000999001
#Indo_C (19) v. Med_B (10)  0.0054882091 0.000999001
#Indo_C (19) v. Med_C (5)   0.0038374449 0.000999001
#Med_A (19) v. Med_B (10)   0.0010465984 0.000999001
#Med_A (19) v. Med_C (5)    0.0009911861 0.000999001
#Med_B (10) v. Med_C (5)    0.0002736392 0.262737263
zczcpop2.d<-popStructTest(zczc125.pop2.gt, nrep=1000, stats="d", type="pairwise")
#Population structure results:
#  D     D.p.val
#Atl_BA (16) v. Atl_BB (17) 0.0003749781 0.000999001
#Atl_BA (16) v. Atl_BC (20) 0.0002866787 0.000999001
#Atl_BA (16) v. Atl_BD (2)  0.0044585916 0.004995005
#Atl_BA (16) v. Indo_A (5)  0.0010324528 0.000999001
#Atl_BA (16) v. Indo_B (12) 0.0005395685 0.000999001
#Atl_BA (16) v. Indo_C (19) 0.0005847121 0.000999001
#Atl_BA (16) v. Med_A (19)  0.0043066847 0.000999001
#Atl_BA (16) v. Med_B (10)  0.0046269876 0.000999001
#Atl_BA (16) v. Med_C (5)   0.0030101154 0.001998002
#Atl_BB (17) v. Atl_BC (20) 0.0002638527 0.000999001
#Atl_BB (17) v. Atl_BD (2)  0.0044960721 0.005994006
#Atl_BB (17) v. Indo_A (5)  0.0009705584 0.000999001
#Atl_BB (17) v. Indo_B (12) 0.0005004311 0.000999001
#Atl_BB (17) v. Indo_C (19) 0.0005443222 0.000999001
#Atl_BB (17) v. Med_A (19)  0.0046612693 0.000999001
#Atl_BB (17) v. Med_B (10)  0.0050977498 0.000999001
#Atl_BB (17) v. Med_C (5)   0.0033936088 0.000999001
#Atl_BC (20) v. Atl_BD (2)  0.0042624120 0.005994006
#Atl_BC (20) v. Indo_A (5)  0.0008863707 0.000999001
#Atl_BC (20) v. Indo_B (12) 0.0003726065 0.000999001
#Atl_BC (20) v. Indo_C (19) 0.0004592980 0.000999001
#Atl_BC (20) v. Med_A (19)  0.0045357318 0.000999001
#Atl_BC (20) v. Med_B (10)  0.0049622268 0.000999001
#Atl_BC (20) v. Med_C (5)   0.0032970732 0.000999001
#Atl_BD (2) v. Indo_A (5)   0.0053311068 0.041958042
#Atl_BD (2) v. Indo_B (12)  0.0046072515 0.014985015
#Atl_BD (2) v. Indo_C (19)  0.0050478872 0.004995005
#Atl_BD (2) v. Med_A (19)   0.0102282213 0.004995005
#Atl_BD (2) v. Med_B (10)   0.0108374479 0.015984016
#Atl_BD (2) v. Med_C (5)    0.0087235751 0.059940060
#Indo_A (5) v. Indo_B (12)  0.0007089331 0.000999001
#Indo_A (5) v. Indo_C (19)  0.0006463077 0.000999001
#Indo_A (5) v. Med_A (19)   0.0051637689 0.000999001
#Indo_A (5) v. Med_B (10)   0.0055861751 0.000999001
#Indo_A (5) v. Med_C (5)    0.0040281371 0.010989011
#Indo_B (12) v. Indo_C (19) 0.0003212235 0.000999001
#Indo_B (12) v. Med_A (19)  0.0046225813 0.000999001
#Indo_B (12) v. Med_B (10)  0.0050534787 0.000999001
#Indo_B (12) v. Med_C (5)   0.0034624673 0.000999001
#Indo_C (19) v. Med_A (19)  0.0050349933 0.000999001
#Indo_C (19) v. Med_B (10)  0.0054882091 0.000999001
#Indo_C (19) v. Med_C (5)   0.0038374449 0.000999001
#Med_A (19) v. Med_B (10)   0.0010465984 0.000999001
#Med_A (19) v. Med_C (5)    0.0009911861 0.000999001
#Med_B (10) v. Med_C (5)    0.0002736392 0.266733267

mden.gt<-genind2gtypes(mden43.gi)
mden.pop.gt<-genind2gtypes(mden43.pop.gi)
mden.d<-popStructTest(mden.gt, nrep=1000, stats="d")
#Population structure results:
#  estimate       p.val
#D 0.002539642 0.000999001
mden.pop.d<-popStructTest(mden.pop.gt, nrep=1000, stats="d", type="pairwise")
#Population structure results:
#  D     D.p.val
#Atl-Bahamas (7) v. Atl-East (16)   0.0007304617 0.000999001
#Atl-Bahamas (7) v. Atl-Other (5)   0.0005096556 0.228771229
#Atl-Bahamas (7) v. Indo-Africa (5) 0.0021496983 0.004995005
#Atl-Bahamas (7) v. Indo-Hawaii (6) 0.0036290272 0.000999001
#Atl-Bahamas (7) v. Indo-South (3)  0.0033229895 0.010989011
#Atl-East (16) v. Atl-Other (5)     0.0004595851 0.005994006
#Atl-East (16) v. Indo-Africa (5)   0.0020853483 0.000999001
#Atl-East (16) v. Indo-Hawaii (6)   0.0035367964 0.000999001
#Atl-East (16) v. Indo-South (3)    0.0032471023 0.002997003
#Atl-Other (5) v. Indo-Africa (5)   0.0016911061 0.044955045
#Atl-Other (5) v. Indo-Hawaii (6)   0.0029224989 0.002997003
#Atl-Other (5) v. Indo-South (3)    0.0026777653 0.019980020
#Indo-Africa (5) v. Indo-Hawaii (6) 0.0008923798 0.002997003


#Indo-Africa (5) v. Indo-South (3)  0.0008788438 0.165834166
#Indo-Hawaii (6) v. Indo-South (3)  0.0008375751 0.114885115

# Plotting FST ------------------------------------------------------------
#Plotting confidence intervals and Fst point estimates for within ocean-basin pairs that were significant
fst.all<-read.csv("./fst.csv")
fst.all<-fst.all[,1:10]
fst.all

fst.zcav<-subset(fst.all, fst.all$species=="zcav")
fst.zcav
fst.mden<-subset(fst.all, fst.all$species=="mden")
fst.mden


zcav.fst.plot<-ggplot(fst.zcav, aes(PopComb, Fst, colour=Ocean1)) + 
  scale_color_manual(values=c("#8c96a0", "#be6620", "#67abd1")) +
  facet_grid(cols=vars(Ocean1), scales="free", space="free") +
  geom_linerange(aes(ymin=LowerBound, ymax=UpperBound)) +
  geom_point(size=1) +
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="none") 
zcav.fst.plot
ggsave("zcav.fst.png", plot=zcav.fst.plot)

mden.fst.plot<-ggplot(fst.mden, aes(PopComb, Fst, colour=Ocean1)) + 
  scale_color_manual(values=c("#8c96a0", "#be6620")) +
  facet_grid(cols=vars(Ocean1), scales="free", space="free") +
  geom_linerange(aes(ymin=LowerBound, ymax=UpperBound)) +
  geom_point() +
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="none") 
mden.fst.plot
ggsave("mden.fst.png", plot=mden.fst.plot)

# IQtree plotting ---------------------------------------------------------


