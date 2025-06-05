#Genotyping western leatherback gtseq data using microhaplot

#Based on Mnov.microhaplot.analysis.R by Karen Martien
devtools::install_github("ngthomas/microhaplot", build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))

library(vcfR)
library(microhaplot)
library(tidyverse)
runShinyHaplot(app.path)
microhaplot::mvShinyHaplot("~/Documents/GitHub/Shiny")
app.path <- "~/Documents/GitHub/Shiny/microhaplot/" #path to where the shiny app is located
runShinyHaplot(app.path)

#Prepare sample IDs
library(swfscMisc)
inds<-read.csv("data-raw/mplot_labels/w-leatherback-inds.txt", sep = "\t", header=FALSE)
head(inds)
dim(inds)
ind.name<-inds$V2
ind.name
ind.name.norep<-ind.name %<>%
  gsub("b", "", .) %>%
  gsub("c", "", .) 
ind.name.norep.pad<-zero.pad(as.numeric(ind.name.norep))
ind.name.norep.pad.z<-paste("z0", ind.name.norep.pad, sep="")
ind.name.norep.pad.z
new.inds<-cbind(inds, ind.name.norep.pad.z)
write.table(new.inds, "data-raw/mplot_labels/w-leatherback-inds-zpad.txt", sep="\t", col.names =FALSE, row.names = FALSE)
#Manually go into file and add the "b" and "c" for duplicates. 
#Manually reformat into required format for microhaploty: filename, ID, pop (NA here)

#Preparing paths for running microhaplot
run.label <- "wpac.sams"

#Microhaplot: Customize the following paths
sam.path <- "data-raw/sam.files/all/" #path to where all of the sam files are located
label.path <- "data-raw/mplot_labels/w-leatherback-inds-zpad.txt" #path to where the sample id file is located
vcf.path <- "data-raw/vcf/Dcor_DcPanel_205_maf.targetSNPs_012224.recode.vcf" #path to where the vcf file with snp loci is located
out.path <- "results-R/microhaplot/" #path to where you want the results to be located
app.path <- "~/Documents/GitHub/Shiny/microhaplot/" #path to where the shiny app is located

#read in vcf file
vcf <- read.vcfR(vcf.path)
vcf

#Run Microhaplot
haplo.read.tbl <- prepHaplotFiles(run.label = run.label,
                                  sam.path = sam.path,
                                  out.path = out.path,
                                  label.path = label.path,
                                  vcf.path = vcf.path,
                                  app.path = app.path,
                                  n.jobs = 4) # use all the cores!
runShinyHaplot(app.path) 

#Call and evaluate mplot Genotypes
#Based on call.and.evaluate.mplot.genos.R by Karen Martien
setwd("/Users/aonoufriou/Documents/GitHub/Turtle_GTseq_SNPs")

#devtools::install_github("kmartien/vcftoolsR")
library(vcfR)
library(vcftoolsR) #Karen's package
library(tidyverse)
library(dplyr)
source("R/functions/mplot2tgt.R") #MUST UPDATE WITH PATH TO YOUR OWN FILES
source("R/functions/Compare.replicates.R")

project <- "wpac.sams"

AB.min.het <- 3/7
AB.max.homo <- 2/8
min.AR.het <- 3/10
max.AR.homo <- 2/10
min.read.depth <- 20
num.locs <- 185 #from CHROM in vcf file
min.genos.per.ind <- 111 #num.locs * 0.6

tgt <- mplot2tgt(project = project, AB.min.het = AB.min.het, AB.max.homo = AB.max.homo,
                 min.read.depth = min.read.depth)
saveRDS(tgt, file = paste0('results-R/tgt.', project, '.rds'))
#tgt.after.merge<-readRDS("results-R/tgt.wpac.sams.rds")
#tgt<-tgt.after.merge

# compare replicates (only run if replicates present)
LABIDs <- unique(tgt$Indiv) %>% substr(start = 1, stop = 8)
replicates <- LABIDs[duplicated(LABIDs)]
mismatches.to.check <-do.call('rbind',lapply(replicates, function(r){
  rep.tgt <- tgt[grep(substr(r, start = 1, stop = 8), tgt$Indiv),]
  #  rep.tgt <- filter(tgt, Indiv %in% c(r,paste0(r,"b")))
  mismatches <- compare.replicates(rep.tgt)
}))
if(nrow(mismatches.to.check > 0)) {
  print("Some replicates have mismatched genotypes")
  print(paste0("Mismatches saved to results-R/", project, ".genotype.mismatches.rda"))
  save(mismatches.to.check, file = paste0("results-R/", project, ".genotype.mismatches.rda"))
}
mismatches.to.check
nrow(mismatches.to.check) #438 mismatched genotypes
write.csv(mismatches.to.check, "results-raw/mismatches.to.check.csv")

# compare replicates-ONLY NEEDED IF YOU DON"T DO THE ZPAD STEP
#LABIDs <- unique(tgt$Indiv)
#length(unique(LABIDs)) #227
#LABIDs
#LABIDs %<>%
#  gsub("b", "", .) %>%
#  gsub("c", "_", .) 
#summary(duplicated(LABIDs)) #51 TRUE
#replicates <- LABIDs[duplicated(LABIDs)]
#replicates #list of individuals that are duplicated (or triplicated)
#length(replicates) #51
#
##mismatches.to.check <-do.call('rbind',lapply(replicates, function(r){
##  rep.tgt <- tgt[grep(substr(r, start = 1, stop = 8), tgt$Indiv),]
##  #  rep.tgt <- filter(tgt, Indiv %in% c(r,paste0(r,"b")))
##  mismatches <- compare.replicates(rep.tgt)
##}))
#
##Subset the tgt file to contain only the duplicated individuals
#  a<-subset(tgt, tgt$Indiv %in% replicates) #tgt subset based on the "replicates" list from above
#  b<-tgt[grep('*b', tgt$Indiv),] #tgt subset based on any samples that have a "b" in Indiv column
#  c<-tgt[grep('*c', tgt$Indiv),] #tgt subset based on any samples that have a "c" in Indiv column
#
#head(a)
#head(b)  
#head(c)
#
#length(unique(a$Indiv)) #51
#length(unique(b$Indiv)) #51
#length(unique(c$Indiv)) #1
#
#rep.tgt<-rbind(a, b, c)
#length(unique(rep.tgt$Indiv)) #103  
#unique(sort(rep.tgt$Indiv)) #Checking to make sure that the replicated samples are all here
#
###Remove the "b" and "c" from replicated sample IDs
##rep.tgt$Indiv %<>%
##  gsub("b", "", .) %>%
##  gsub("c", "_", .) 
##unique(sort(rep.tgt$Indiv))
#
#for (i in replicates){
#  replicate.tgt <- rep.tgt[grepl("i*", rep.tgt$Indiv),]
#  replicatedID.comp <- compare.replicates(replicate.tgt)
#  mismatches <- rbind(replicatedID.comp)
#  }
#mismatches
#length(unique(mismatches$Indiv)) #103
#
#
#if(nrow(mismatches > 0)) {
#  print("Some replicates have mismatched genotypes")
#  print(paste0("Mismatches saved to results-R/", project, ".genotype.mismatches.rda"))
#  save(mismatches, file = paste0("results-R/", project, ".genotype.mismatches.rda"))
#}

####I have 476 rows of mismatched genotypes. 
#There still may be some problem loci and individuals in the dataset that need to be removed before comparing replicated genotypes. 
#I will move forward with identifying and removing problematic loci/individuals before comparing replicated genotypes

# Identify samples with Ns or Xs in their haplotypes or more than 2 haplotypes
questionable.hap <- sapply(1:nrow(tgt), function(i){
  ifelse(length(grep("N", tgt$gt[i],)) > 0 || length(grep("X", tgt$gt[i],)) > 0 
         || tgt$num.haps[i] > 2, TRUE, FALSE) 
})
genos.to.check <- filter(tgt, questionable.hap == TRUE)
write.csv(genos.to.check, "results-raw/genos.to.check.csv")

table(genos.to.check$locus)
#locus016 locus019 locus069 
#3       44        1


############################################################################
###### NEED TO DECIDE HOW TO DEAL WITH QUESTIONABLE HAPLOTYPES AT THIS POINT
############################################################################

tgt <- filter(tgt, questionable.hap == FALSE)
if(nrow(genos.to.check > 0)) {
  print("Some samples with Ns or Xs in their haplotypes or more than 2 haplotypes")
  print(paste0("Questionable genotypes saved to results-R/", project, ".genos.to.check.rda"))
  save(genos.to.check, file = paste0("results-R/", project, ".genos.to.check.rda"))
}

# summarize individual data
missing.data.ind <- data.frame(table(tgt$Indiv[!is.na(tgt$gt)])) %>%
  mutate(missing = num.locs-Freq)
names(missing.data.ind) <- c("labID", "genos", "missing")
length(which(missing.data.ind$genos >= min.genos.per.ind))
inds.2.keep <- filter(missing.data.ind, genos >= min.genos.per.ind) |> 
  pull(labID)
tgt <- filter(tgt, Indiv %in% inds.2.keep)
num.inds <- length(unique(tgt$Indiv))
write.csv(missing.data.ind, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.num.genos.per.ind.csv"))
dim(missing.data.ind)

# summarize locus data
tgt_long <- tgt |> 
  select(locus, Indiv, gt, depth.1, depth.2) |> 
  separate_wider_delim(
    cols = gt, 
    delim = '/', 
    names = c('haplo.1', 'haplo.2'), 
    cols_remove = FALSE
  ) |> 
  mutate(depth.2 = ifelse(haplo.2 == haplo.1, NA, depth.2))
#loc.sum <- tgt_long %>%
#  pivot_longer(cols = c(haplo.1, haplo.2), names_to = 'hap') |> 
#  mutate(tmp = strsplit(as.character(value), "")) %>%
#  unnest(tmp) %>%
#  group_by(locus, Indiv, hap) %>%
#  mutate(name = 1:n()) %>%
#  pivot_wider(id_cols = c(locus, Indiv, gt, hap), values_from = tmp, names_prefix = 'snp') |> 
#  ungroup() |> 
#  filter(!is.na(gt)) |> 
#  group_by(locus) |> 
#  summarise(
#    inds.genoed = n() / 2,
#    num.unique.genos = length(unique(gt)),
#    num.alleles.pos1 = length(unique(snp1)),
#    num.alleles.pos2 = length(unique(snp2)),
#    num.alleles.pos3 = length(unique(snp3)),
#    num.alleles.pos4 = length(unique(snp4)),
#    num.alleles.pos5 = length(unique(snp5))
#  )
loc.sum <- tgt_long %>%
  pivot_longer(cols = c(haplo.1, haplo.2), names_to = 'hap') |> 
  mutate(tmp = strsplit(as.character(value), "")) %>%
  unnest(tmp) %>%
  group_by(locus, Indiv, hap) %>%
  mutate(name = 1:n()) %>%
  pivot_wider(id_cols = c(locus, Indiv, gt, hap), values_from = tmp, names_prefix = 'snp') |> 
  ungroup() |> 
  filter(!is.na(gt)) |> 
  group_by(locus) |> 
  summarise(
    inds.genoed = n() / 2,
    num.unique.genos = length(unique(gt)),
    num.alleles.pos1 = length(unique(snp1)),
    num.alleles.pos2 = length(unique(snp2)),
    num.alleles.pos3 = length(unique(snp3))
  )
write.csv(loc.sum, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.locus.summary.csv"))

tgt <- filter(tgt,
              locus %in% 
                (filter(loc.sum,
                        inds.genoed >= num.inds * 0.5) |> 
                   pull(locus)))
geno.table <- tgt.2.geno.table(tgt) 

save(geno.table, tgt, loc.sum, file = paste0("results-R/", project, ".", min.read.depth, "readsMin.geno.eval.rda"))


# Karen's notes
# I went through all of the loci and eliminated ones that were questionable - multiple
# SNPs, but only 2 alleles, loci where one SNP had a MAC < 10, etc. I made notes
# in the .locus.summary.csv document generated ~20 lines above and changed the name
# to 'results_raw/RunMS58.locus.notes.csv. I then printed out all locus names
# from the vcf file used to generate the tgt object above, and deleted SNPs that
# I had deemed untrustworthy (see results_raw/RunMS58.locus.notes.csv for reason)
# The next lines remove all SNPs not deemed trustworthy re-do the microhaplot calling

##########################
#   Filter
##########################

###Pre-filter tgt
#load tgt table from previous step if the tgt has been modified in any way
tgt.prefilter<-tgt
tgt.after.merge.before.filter<-tgt
#####First Filter
#Individuals to remove: z0012526, z0216912, z0216919, z0216911, z0216922, z0216917
###These all had >90% missing data

#Loci to remove: locus111, locus155, locus201, locus271
###These were all genotyped in <10% of individuals
###Also remove locus019 (questionable haplotypes across 67 individuals)

#Delete bad individuals and loci and make a new tgt
#tgt<-tgt.prefilter
tgt<-tgt.after.merge.before.filter
tgt

tgt %>% filter(!(locus %in% c("locus019", "locus111", "locus155", "locus201", "locus271")) & !(Indiv %in% c("z0012526", "z0216912", "z0216919", "z0216911", "z0216922", "z0216917")))
sort(unique(tgt$locus))
sort(unique(tgt$Indiv))

saveRDS(tgt, file = paste0('results-R/tgt_filter3.', project, '.rds'))

#compare replicates after filter
LABIDs <- unique(tgt$Indiv) %>% substr(start = 1, stop = 8)
replicates <- LABIDs[duplicated(LABIDs)]
mismatches.to.check <-do.call('rbind',lapply(replicates, function(r){
  rep.tgt <- tgt[grep(substr(r, start = 1, stop = 8), tgt$Indiv),]
  #  rep.tgt <- filter(tgt, Indiv %in% c(r,paste0(r,"b")))
  mismatches <- compare.replicates(rep.tgt)
}))
if(nrow(mismatches.to.check > 0)) {
  print("Some replicates have mismatched genotypes")
  print(paste0("Mismatches saved to results-R/", project, ".genotype.mismatches.filter2.rda"))
  save(mismatches.to.check, file = paste0("results-R/", project, ".genotype.mismatches.filter2.rda"))
}
mismatches.to.check
nrow(mismatches.to.check) #438 mismatched genotypes
write.csv(mismatches.to.check, "results-raw/mismatches.to.check.filter2.csv")


#caluclate questionable haplotypes after filter 1
questionable.hap <- sapply(1:nrow(tgt), function(i){
  ifelse(length(grep("N", tgt$gt[i],)) > 0 || length(grep("X", tgt$gt[i],)) > 0 
         || tgt$num.haps[i] > 2, TRUE, FALSE) 
})
genos.to.check <- filter(tgt, questionable.hap == TRUE)
write.csv(genos.to.check, "results-raw/genos.to.check.filter3.csv")

table(genos.to.check$locus) #EMPTY

#calcuate locus summaries after filter 
tgt_long <- tgt |> 
  select(locus, Indiv, gt, depth.1, depth.2) |> 
  separate_wider_delim(
    cols = gt, 
    delim = '/', 
    names = c('haplo.1', 'haplo.2'), 
    cols_remove = FALSE
  ) |> 
  mutate(depth.2 = ifelse(haplo.2 == haplo.1, NA, depth.2))

loc.sum <- tgt_long %>%
  pivot_longer(cols = c(haplo.1, haplo.2), names_to = 'hap') |> 
  mutate(tmp = strsplit(as.character(value), "")) %>%
  unnest(tmp) %>%
  group_by(locus, Indiv, hap) %>%
  mutate(name = 1:n()) %>%
  pivot_wider(id_cols = c(locus, Indiv, gt, hap), values_from = tmp, names_prefix = 'snp') |> 
  ungroup() |> 
  filter(!is.na(gt)) |> 
  group_by(locus) |> 
  summarise(
    inds.genoed = n() / 2,
    num.unique.genos = length(unique(gt)),
    num.alleles.pos1 = length(unique(snp1)),
    num.alleles.pos2 = length(unique(snp2)),
    num.alleles.pos3 = length(unique(snp3))
  )
#write.csv(loc.sum, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.locus.summary.filter1.csv"))
#NOT writing this file yet. Want to re-calculate after the individual summary section since some individuals may be dropped
#HOWEVER, need to keep this here before individual summaries since it also has a filtering stip for na genotypes

#calculate individual summaries after filter 1(calculating after locus summaries since that step can remove NA loci)
num.locs<-length(unique(tgt$locus))
missing.data.ind <- data.frame(table(tgt$Indiv[!is.na(tgt$gt)])) %>%
  mutate(missing = num.locs-Freq)
names(missing.data.ind) <- c("labID", "genos", "missing")
length(which(missing.data.ind$genos >= min.genos.per.ind))
inds.2.keep <- filter(missing.data.ind, genos >= min.genos.per.ind) |> 
  pull(labID)
tgt <- filter(tgt, Indiv %in% inds.2.keep)
num.inds <- length(unique(tgt$Indiv))
write.csv(missing.data.ind, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.num.genos.per.ind.filter3.csv"))
missing.data.ind

#recalculate locus summaries with updated sample list after filter 1
tgt_long <- tgt |> 
  select(locus, Indiv, gt, depth.1, depth.2) |> 
  separate_wider_delim(
    cols = gt, 
    delim = '/', 
    names = c('haplo.1', 'haplo.2'), 
    cols_remove = FALSE
  ) |> 
  mutate(depth.2 = ifelse(haplo.2 == haplo.1, NA, depth.2))

loc.sum <- tgt_long %>%
  pivot_longer(cols = c(haplo.1, haplo.2), names_to = 'hap') |> 
  mutate(tmp = strsplit(as.character(value), "")) %>%
  unnest(tmp) %>%
  group_by(locus, Indiv, hap) %>%
  mutate(name = 1:n()) %>%
  pivot_wider(id_cols = c(locus, Indiv, gt, hap), values_from = tmp, names_prefix = 'snp') |> 
  ungroup() |> 
  filter(!is.na(gt)) |> 
  group_by(locus) |> 
  summarise(
    inds.genoed = n() / 2,
    num.unique.genos = length(unique(gt)),
    num.alleles.pos1 = length(unique(snp1)),
    num.alleles.pos2 = length(unique(snp2)),
    num.alleles.pos3 = length(unique(snp3))
  )
write.csv(loc.sum, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.locus.summary.filter3.csv"))

#Checking what happened to the mismatched genotypes from replicates
#geno.table <- tgt.2.geno.table(tgt) 
#geno.table %>% filter(Indiv=="z0012519") %>% pull(locus165) #Genotype is NA
#geno.table %>% filter(Indiv=="z0006799") %>% pull(locus314) #Genotype is NA
#geno.table %>% filter(Indiv=="z0006802") %>% pull(locus314) #Genotype is NA
#geno.table %>% filter(Indiv=="z0134098") %>% pull(locus314) #Genotype is NA

save(geno.table, tgt, loc.sum, file = paste0("results-R/", project, ".", min.read.depth, "readsMin.geno.eval.rda"))

#####Histograms of data#####
hist(loc.sum$inds.genoed)
hist(missing.data.ind$genos)

