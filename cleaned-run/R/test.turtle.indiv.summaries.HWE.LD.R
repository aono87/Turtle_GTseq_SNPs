library(tidyverse)
library(strataG)
library(dplyr)
#library(swfscMisc)

project <- 'dcor.wpac.test'
min.reads <- 20

load(paste0('data/gtypes_', project, '_minReads.', min.reads, '.rda'))
g #if going straight from other scripts
strata2keep <- c("CentralCA", "INDONESIA", "MALAYSIA", "PAPUA_NEW_GUINEA", "SOLOMON_ISLANDS")


### REMOVE INDIVIDUALS THAT ARE OUTLILERS IN TERMS OF HIGH HOMOZYGOSITY ####
ind.summary <- summarizeInds(g)
hist(ind.summary$pct.loci.homozygous)
high.homo.samps <- filter(ind.summary, pct.loci.homozygous > 0.7) |> 
  pull(id)
high.homo.samps #none in my dataset
#g <- g[-which(getIndNames(g) %in% high.homo.samps),,]

### HARDY-WEINBERG EQUILIBRIUM ###############################
##Karen originally calculated a simple Bonferroni correction: 0.05/num populations
##I also calculated the Sequential bonferron correction using p.adjust()
pop.g <- stratify(g, "sampling.location")
pop.g

hwe.list <- lapply(c("CentralCA", "INDONESIA", "MALAYSIA", "PAPUA_NEW_GUINEA", "SOLOMON_ISLANDS"), function(s){
  x <- hweTest(pop.g[,,s]) %>% data.frame() %>% rownames_to_column()
  names(x) <- c("locus", s)
  return(x)
})

hwe.res <- hwe.list |> reduce(full_join, by = 'locus') 

#Karen's original code to plot the p-value and Bonferroni correction (p-value/number of strata)
hwe.res.p.bon <- hwe.res |> 
  rowwise() %>%
  mutate(
    num.sig.p = sum(
      c_across('CentralCA':'SOLOMON_ISLANDS') < 0.05),
    num.sig.bon = sum(
        c_across('CentralCA':'SOLOMON_ISLANDS') < (0.05/5) #have updated to 5 to reflect the number of Dcor populations
      ))
head(hwe.res.p.bon)

#Calculate adjusted p-value using sequential bonferroni correction ("holm")
hwe.res.adj<-hwe.res %>% column_to_rownames("locus")
hwe.res.adj <- apply(hwe.res.adj, 1, function(x) p.adjust(x, method = "holm"))
hwe.res.adj<-t(hwe.res.adj)
hwe.res.adj<-as.data.frame(hwe.res.adj)
hwe.res.adj<-hwe.res.adj %>% rownames_to_column("locus")
head(hwe.res.adj)

#Add column and count number of strata with adjusted p-value < 0.05
hwe.res.seq <-hwe.res.adj |> 
  rowwise() %>%
  mutate(
    num.sig.seq = sum(
      c_across('CentralCA':'SOLOMON_ISLANDS') < 0.05
    ))

head(hwe.res.p.bon)
head(hwe.res.seq)

#Combine the different p-values into table
HWE.sig<- full_join(hwe.res.p.bon, hwe.res.seq, by="locus") %>%
  select(locus, CentralCA.x, INDONESIA.x, MALAYSIA.x, PAPUA_NEW_GUINEA.x, SOLOMON_ISLANDS.x, num.sig.p, num.sig.bon, num.sig.seq)
head(HWE.sig)

write.csv(HWE.sig, file = paste0("results-raw/",project, ".hwe.results.csv"))

### SUMMARIZE REMAINING LOCI ################################################

loc.sum <- summarizeLoci(g)
loc.sums.by.strat <- summarizeLoci(pop.g, by.strata = TRUE) |> 
  filter(stratum %in% c("CentralCA", "INDONESIA", "MALAYSIA", "PAPUA_NEW_GUINEA", "SOLOMON_ISLANDS")) |> 
  mutate(obs_minus_exp = exptd.het - obsvd.het)

loc.sum <- left_join(loc.sum, 
g@data |> select(c(locus, allele)) |> 
  distinct() |> 
  filter(!is.na(allele)) |> 
  mutate(num.snps = nchar(allele)) |> 
  select(c(locus, num.snps)) |> 
  distinct()
)
write.csv(loc.sum, file = 'results-raw/final.loc.sum.csv')
save(loc.sums.by.strat, loc.sum, file = 'data/final.loc.sum.rda')
