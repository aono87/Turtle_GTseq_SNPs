#Genotyping western leatherback gtseq data using microhaplot

#Based on Mnov.microhaplot.analysis.R by Karen Martien
devtools::install_github("ngthomas/microhaplot", build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))

library(vcfR)
library(microhaplot)
library(tidyverse)
runShinyHaplot(app.path)
microhaplot::mvShinyHaplot("~/Documents/GitHub/Shiny")
app.path <- "~/Documents/GitHub/Shiny/microhaplot/" #path to where the shiny app is located

#Preparing paths for running microhaplot
run.label <- "dcor.wpac.test"

#Microhaplot: Customize the following paths
sam.path <- "~/Documents/GitHub/Turtle_GTseq_SNPs/data-raw/sam.files/all/" #path to where all of the sam files are located
label.path <- "~/Documents/GitHub/Turtle_GTseq_SNPs/data-raw/mplot_labels/w-leatherback-inds-zpad.txt" #path to where the sample id file is located
vcf.path <- "~/Documents/GitHub/Turtle_GTseq_SNPs/data-raw/vcf/Dcor_DcPanel_205_maf.targetSNPs_012224.recode.vcf" #path to where the vcf file with snp loci is located
out.path <- "results-R/microhaplot/" #path to where you want the results to be located
app.path <- "~/Documents/GitHub/Shiny/microhaplot/" #path to where the shiny app is located

#read in vcf file
vcf <- read.vcfR(vcf.path)
vcf #185 loci

gi<-vcfR2genind(vcf)
gi #231 loci
loci.list<-gi@all.names

gi@all.names

write.csv(loci.list, file="data/loci.list.csv")

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

#devtools::install_github("kmartien/vcftoolsR")
library(vcfR)
library(vcftoolsR) #Karen's package
library(tidyverse)
library(dplyr)
source("R/functions/mplot2tgt.R") #MUST UPDATE WITH PATH TO YOUR OWN FILES
source("R/functions/Compare.replicates.R")

project <- "dcor.wpac.test"

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


# Identify samples with Ns or Xs in their haplotypes or more than 2 haplotypes
questionable.hap <- sapply(1:nrow(tgt), function(i){
  ifelse(length(grep("N", tgt$gt[i],)) > 0 || length(grep("X", tgt$gt[i],)) > 0 
         || tgt$num.haps[i] > 2, TRUE, FALSE) 
})
genos.to.check <- filter(tgt, questionable.hap == TRUE)
table(genos.to.check$locus)
#locus016 locus019 locus069 
#3       44        1

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
tgt <- filter(tgt, Indiv %in% inds.2.keep) #filters out individuals that have fewer genotypes than "min.genos.per.ind"
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

##########################
#   Filter
##########################

#####First Filter
#Individuals to remove: z0012526, z0216912, z0216919, z0216911, z0216922, z0216917
###These all had >90% missing data

#Loci to remove: locus111, locus155, locus201, locus271
###These were all genotyped in <10% of individuals
###Also remove locus019 (questionable haplotypes across 67 individuals)

#Delete bad individuals and loci and make a new tgt
tgt<-tgt %>% filter(!(locus %in% c("locus019", "locus111", "locus155", "locus201", "locus271")) & !(Indiv %in% c("z0012526", "z0216912", "z0216919", "z0216911", "z0216922", "z0216917")))
sort(unique(tgt$locus))
sort(unique(tgt$Indiv))

saveRDS(tgt, file = paste0('results-R/tgt.', project, '.rds'))


#caluclate questionable haplotypes after filter 1
questionable.hap <- sapply(1:nrow(tgt), function(i){
  ifelse(length(grep("N", tgt$gt[i],)) > 0 || length(grep("X", tgt$gt[i],)) > 0 
         || tgt$num.haps[i] > 2, TRUE, FALSE) 
})
genos.to.check <- filter(tgt, questionable.hap == TRUE)
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
write.csv(missing.data.ind, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.num.genos.per.ind.filter.csv"))

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
write.csv(loc.sum, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.locus.summary.filter.csv"))

geno.table<-tgt.2.geno.table(tgt)

#Checking what happened to the mismatched genotypes from replicates
save(geno.table, tgt, loc.sum, file = paste0("results-R/", project, ".", min.read.depth, "readsMin.geno.eval.rda"))

#####Histograms of data#####
hist(loc.sum$inds.genoed)
hist(missing.data.ind$genos)

