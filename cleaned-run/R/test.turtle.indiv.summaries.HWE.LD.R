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

#Look at homozygosity within populations
strats.to.analyze<-strataSplit(g)
summary.list <- lapply(1:length(strats.to.analyze), function(s){ ##Takes a while to run
  print('next')
  x <- summarizeInds(strats.to.analyze[[s]])# %>% data.frame() %>% rownames_to_column()
  #names(x) <- c("locus", names(strats.to.analyze)[s])
  return(x)
})
names(summary.list) <- names(strats.to.analyze)

table(ind.summary$stratum)
#CentralCA        INDONESIA         MALAYSIA PAPUA_NEW_GUINEA  SOLOMON_ISLANDS 
#7               63                7               12               54 
summary(ind.summary$stratum)

par(mfrow=c(3,2))
hist(summary.list$CentralCA$pct.loci.homozygous, main="Central CA (n=7)", xlab="Percent homozygous loci per individual", xlim=c(0.35,0.70), breaks=8)
hist(summary.list$INDONESIA$pct.loci.homozygous, main="Indonesia (n=63)", xlab="Percent homozygous loci per individual", xlim=c(0.35,0.70), breaks=8)
hist(summary.list$MALAYSIA$pct.loci.homozygous, main="Malaysia (n=7)", xlab="Percent homozygous loci per individual", xlim=c(0.35,0.70), breaks=8)
hist(summary.list$PAPUA_NEW_GUINEA$pct.loci.homozygous, main="Papua New Guinea (n=12)", xlab="Percent homozygous loci per individual", xlim=c(0.35,0.70), breaks=8)
hist(summary.list$SOLOMON_ISLANDS$pct.loci.homozygous, main="Solomon Islands (n=54)", xlab="Percent homozygous loci per individual", xlim=c(0.35,0.70), breaks=8)
hist(ind.summary$pct.loci.homozygous, main="All populations combined (n=143)", xlab="Percent homozygous loci per individual", xlim=c(0.35,0.70), breaks=8)
dev.off()

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

############LINKAGE DISEQUILIBRIUM#####################
#Create a list of gtypes files for each stratum
strats.to.analyze<-strataSplit(g, remove.sequences = TRUE)
length(unique(strats.to.analyze$CentralCA@data$locus)) #163 unique loci
length(unique(strats.to.analyze$INDONESIA@data$locus)) #163 unique loci
length(unique(strats.to.analyze$MALAYSIA@data$locus)) #163 unique loci
length(unique(strats.to.analyze$PAPUA_NEW_GUINEA@data$locus)) #163 unique loci
length(unique(strats.to.analyze$SOLOMON_ISLANDS@data$locus)) #163 unique loci

#calculate LD for all individuals together
ld.overall <- LDgenepop(g) ##Takes a while to run
length(unique(ld.overall$Locus.1))
length(unique(ld.overall$Locus.2))

#caluclate LD for each stratum individually
ld.list <- lapply(1:length(strats.to.analyze), function(s){ ##Takes a while to run
  print('next')
  x <- LDgenepop(strats.to.analyze[[s]])# %>% data.frame() %>% rownames_to_column()
  #names(x) <- c("locus", names(strats.to.analyze)[s])
  return(x)
})
names(ld.list) <- names(strats.to.analyze)

ld.sig.res <- lapply(ld.list, function(s){
  filter(s, p.value < 0.05) |> 
    select(c(Locus.1, Locus.2, p.value))
}) 

for (i in 1:length(ld.sig.res)){
  names(ld.sig.res[[i]])[3] <- paste0('p.val.',names(ld.list[i]))
}
 
ld.sig.res <- ld.sig.res |> reduce(full_join) |> 
  rowwise() %>%
  mutate(
    num.sig = sum(
      c_across('p.val.CentralCA':'p.val.SOLOMON_ISLANDS') < 0.05, na.rm = TRUE
    ),
    num.sig.after.correction = sum(
      c_across('p.val.CentralCA':'p.val.SOLOMON_ISLANDS') < (0.05/5), na.rm = TRUE
    )
  ) %>%
  ungroup()

#add overall values
ld.list2<-ld.list
ld.list2$all <- ld.overall
ld.list<-ld.list2

#Finding some NA's in the locus columns. 
unique(ld.sig.res$Locus.1)
unique(ld.sig.res$Locus.2)

length(unique(ld.list$CentralCA$Locus.1))
length(unique(ld.list$CentralCA$Locus.2))

length(unique(ld.list$INDONESIA$Locus.1))
length(unique(ld.list$INDONESIA$Locus.2))

length(unique(ld.list$MALAYSIA$Locus.1))
length(unique(ld.list$MALAYSIA$Locus.2))

length(unique(ld.list$PAPUA_NEW_GUINEA$Locus.1))
length(unique(ld.list$PAPUA_NEW_GUINEA$Locus.2))

unique(ld.list$SOLOMON_ISLANDS$Locus.1)
unique(ld.list$SOLOMON_ISLANDS$Locus.2)

length(unique(ld.list$all$Locus.1))
length(unique(ld.list$all$Locus.2))

#Will run the LD command on each stratum-specific gtypes file separately
g.california<-strats.to.analyze$CentralCA
g.indonesia<-strats.to.analyze$INDONESIA
g.malaysia<-strats.to.analyze$MALAYSIA
g.png<-strats.to.analyze$PAPUA_NEW_GUINEA
g.solomon<-strats.to.analyze$SOLOMON_ISLANDS

sum.california<-summarizeLoci(g.california)
sum.indonesia<-summarizeLoci(g.indonesia)
sum.malaysia<-summarizeLoci(g.malaysia)
sum.png<-summarizeLoci(g.png)
sum.solomon<-summarizeLoci(g.solomon)

california.loci.to.remove<- sum.california %>%
  filter(prop.genotyped < 0.5) %>%
  pull(locus)
california.loci.to.remove #locus110
g.california <- g.california[,-which(getLociNames(g.california) %in% california.loci.to.remove),]

indonesia.loci.to.remove<- sum.indonesia %>%
  filter(prop.genotyped < 0.5) %>%
  pull(locus)
indonesia.loci.to.remove #none
g.indonesia <- g.indonesia[,-which(getLociNames(g.indonesia) %in% indonesia.loci.to.remove),]

malaysia.loci.to.remove<- sum.malaysia %>%
  filter(prop.genotyped < 0.5) %>%
  pull(locus)
malaysia.loci.to.remove #"Dc00544"  "Dc05864"  "Dc06503"  "Dc10878"  "Dc22883"  "Dc27955"  "Dc31007"  "Dc36767"  "Dc41234" "locus027" "locus076" "locus093" "locus128" "locus231" "locus255" "locus349"
g.malaysia <- g.malaysia[,-which(getLociNames(g.malaysia) %in% malaysia.loci.to.remove),]

png.loci.to.remove<- sum.png %>%
  filter(prop.genotyped < 0.5) %>%
  pull(locus)
png.loci.to.remove #"Dc00544"  "Dc10003"  "locus255"
g.png <- g.png[,-which(getLociNames(g.png) %in% png.loci.to.remove),]

solomon.loci.to.remove<- sum.solomon %>%
  filter(prop.genotyped < 0.5) %>%
  pull(locus)
solomon.loci.to.remove #none
g.solomon <- g.solomon[,-which(getLociNames(g.solomon) %in% solomon.loci.to.remove),]

#cleaned so that no locus was genotyped in <50% of individuals
ld.california<-LDgenepop(g.california) 
ld.indonesia<-LDgenepop(g.indonesia)
ld.malaysia<-LDgenepop(g.malaysia)
ld.png<-LDgenepop(g.png)
ld.solomon<-LDgenepop(g.solomon)
