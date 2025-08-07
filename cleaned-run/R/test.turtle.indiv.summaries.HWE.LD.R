library(tidyverse)
library(strataG)
library(dplyr)
#library(swfscMisc)

project <- 'dcor.wpac.test'
min.reads <- 20

load(paste0('data/gtypes_', project, '_minReads.', min.reads, '.rda'))
g #if going straight from other scripts


### REMOVE INDIVIDUALS THAT ARE OUTLILERS IN TERMS OF HIGH HOMOZYGOSITY ####
ind.summary <- summarizeInds(g)
hist(ind.summary$pct.loci.homozygous)
high.homo.samps <- filter(ind.summary, pct.loci.homozygous > 0.7) |> 
  pull(id)
high.homo.samps #none in my dataset
#g <- g[-which(getIndNames(g) %in% high.homo.samps),,]

#Look at homozygosity within populations
strata2keep <- c("Bird's Head-Summer", "Bird's Head-Winter", "HaevoSI-Summer", "HaevoSI-Winter", "IsabelSI-South", "Malaysia", "PNG")
all.strats.g<-strataSplit(g)
strats.to.analyze <- all.strats.g[names(all.strats.g) %in% strata2keep]
strats.to.analyze

summary.list <- lapply(1:length(strats.to.analyze), function(s){ ##Takes a while to run
  print('next')
  x <- summarizeInds(strats.to.analyze[[s]])# %>% data.frame() %>% rownames_to_column()
  #names(x) <- c("locus", names(strats.to.analyze)[s])
  return(x)
})
names(summary.list) <- names(strats.to.analyze)
summary.list

#par(mfrow=c(3,2))
#hist(summary.list$CentralCA$pct.loci.homozygous, main="Central CA (n=7)", xlab="Percent homozygous loci per individual", xlim=c(0.35,0.70), breaks=8)
#hist(summary.list$INDONESIA$pct.loci.homozygous, main="Indonesia (n=63)", xlab="Percent homozygous loci per individual", xlim=c(0.35,0.70), breaks=8)
#hist(summary.list$MALAYSIA$pct.loci.homozygous, main="Malaysia (n=7)", xlab="Percent homozygous loci per individual", xlim=c(0.35,0.70), breaks=8)
#hist(summary.list$PAPUA_NEW_GUINEA$pct.loci.homozygous, main="Papua New Guinea (n=12)", xlab="Percent homozygous loci per individual", xlim=c(0.35,0.70), breaks=8)
#hist(summary.list$SOLOMON_ISLANDS$pct.loci.homozygous, main="Solomon Islands (n=54)", xlab="Percent homozygous loci per individual", xlim=c(0.35,0.70), breaks=8)
#hist(ind.summary$pct.loci.homozygous, main="All populations combined (n=143)", xlab="Percent homozygous loci per individual", xlim=c(0.35,0.70), breaks=8)
#dev.off()

### HARDY-WEINBERG EQUILIBRIUM ###############################
##Karen originally calculated a simple Bonferroni correction: 0.05/num populations
##I also calculated the Sequential bonferron correction using p.adjust()
pop.g <- stratify(g, "Stratum_ABO")
pop.g

#hwe.list <- lapply(c("CentralCA", "INDONESIA", "MALAYSIA", "PAPUA_NEW_GUINEA", "SOLOMON_ISLANDS"), function(s){
#  x <- hweTest(pop.g[,,s]) %>% data.frame() %>% rownames_to_column()
#  names(x) <- c("locus", s)
#  return(x)
#})

hwe.list <- lapply(strata2keep, function(s){
  x <- hweTest(pop.g[,,s]) %>% data.frame() %>% rownames_to_column()
  names(x) <- c("locus", s)
  return(x)
})

hwe.res <- hwe.list |> reduce(full_join, by = 'locus') 

#Karen's original code to plot the p-value and Bonferroni correction (p-value/number of strata)
#hwe.res.p.bon <- hwe.res |> 
#  rowwise() %>%
#  mutate(
#    num.sig.p = sum(
#      c_across('CentralCA':'SOLOMON_ISLANDS') < 0.05),
#    num.sig.bon = sum(
#        c_across('CentralCA':'SOLOMON_ISLANDS') < (0.05/5) #have updated to 5 to reflect the number of Dcor populations
#      ))
#head(hwe.res.p.bon)

#Unadjusted p-values
hwe.res.unadj <-hwe.res |>
  rowwise() |>
  mutate(
    num.sig.p=sum(c_across(!locus)<0.05),
    num.sig.after.correction=sum(c_across(!locus)<0.05/7))
head(hwe.res.unadj)

#Calculate adjusted p-value using sequential bonferroni correction ("holm")
hwe.res.bon<-hwe.res %>% column_to_rownames("locus")
hwe.res.bon <- apply(hwe.res.bon, 1, function(x) p.adjust(x, method = "holm"))
hwe.res.bon<-t(hwe.res.bon)
hwe.res.bon<-as.data.frame(hwe.res.bon)
hwe.res.bon<-hwe.res.bon %>% rownames_to_column("locus")
head(hwe.res.bon)

#Add column and count number of strata with adjusted p-value < 0.05
hwe.res.adj <- hwe.res.bon |> 
  rowwise() |> 
  mutate(
    num.sig.p = sum(c_across(!locus) < 0.05, na.rm = TRUE)) |>
      arrange(desc(num.sig.p))

head(hwe.res.unadj) 
head(hwe.res.adj)

#Write CSV files for unadjusted and adjusted p-values
write.csv(hwe.res.unadj, file = paste0("results-raw/",project, ".hwe.results.unadjusted.csv"))
write.csv(hwe.res.adj, file = paste0("results-raw/",project, ".hwe.results.adjusted.csv"))

save(hwe.res.unadj, hwe.res.adj, file = 'data/hwe.sig.res.rda')

### SUMMARIZE REMAINING LOCI ################################################
loc.sum <- summarizeLoci(g)
loc.sums.by.strat <- summarizeLoci(pop.g, by.strata = TRUE) |> 
#  filter(stratum %in% strata2keep) |> 
  mutate(obs_minus_exp = exptd.het - obsvd.het)

loc.sum <- left_join(loc.sum, g@data |> select(c(locus, allele)) |> 
  distinct() |> 
  filter(!is.na(allele)) |> 
  mutate(num.snps = nchar(allele)) |> 
  select(c(locus, num.snps)) |> 
  distinct()
)

loc.sums.by.strat <- left_join(loc.sums.by.strat, pop.g@data |> select(c(locus, allele)) |> 
                       distinct() |> 
                       filter(!is.na(allele)) |> 
                       mutate(num.snps = nchar(allele)) |> 
                       select(c(locus, num.snps)) |> 
                       distinct()
)

write.csv(loc.sum, file = 'results-raw/final.loc.sum.csv')
write.csv(loc.sums.by.strat, file='results-raw/final.loc.sums.by.strat.csv')
save(loc.sums.by.strat, loc.sum, file = 'data/final.loc.sum.rda')

############LINKAGE DISEQUILIBRIUM#####################
#Create a list of gtypes files for each stratum
strats.to.analyze<-strataSplit(g, remove.sequences = TRUE)
strata2keep
length(unique(strats.to.analyze$`Bird's Head-Summer`@data$locus)) #163 unique loci
length(unique(strats.to.analyze$`Bird's Head-Winter`@data$locus)) #163 unique loci
length(unique(strats.to.analyze$`HaevoSI-Summer`@data$locus)) #163 unique loci
length(unique(strats.to.analyze$`HaevoSI-Winter`@data$locus)) #163 unique loci
length(unique(strats.to.analyze$`IsabelSI-South`@data$locus)) #163 unique loci
length(unique(strats.to.analyze$Malaysia@data$locus)) #163 unique loci
length(unique(strats.to.analyze$PNG@data$locus)) #163 unique loci

#calculate LD for all individuals together
ld.overall <- LDgenepop(g) ##Takes a while to run
length(unique(ld.overall$Locus.1))
length(unique(ld.overall$Locus.2))

#caluclate LD for each stratum individually
#Will run the LD command on each stratum-specific gtypes file separately
g.bhs  <-strats.to.analyze$`Bird's Head-Summer`
g.bhw  <-strats.to.analyze$`Bird's Head-Winter`
g.hues <-strats.to.analyze$`HaevoSI-Summer`
g.huew <-strats.to.analyze$`HaevoSI-Winter`
g.sisa <-strats.to.analyze$`IsabelSI-South`
g.mal  <-strats.to.analyze$Malaysia
g.png  <-strats.to.analyze$PNG

sum.bhs <-summarizeLoci(g.bhs )
sum.bhw <-summarizeLoci(g.bhw )
sum.hues<-summarizeLoci(g.hues)
sum.huew<-summarizeLoci(g.huew)
sum.sisa<-summarizeLoci(g.sisa)
sum.mal <-summarizeLoci(g.mal )
sum.png <-summarizeLoci(g.png )

bhs.loci.to.remove<- sum.bhs %>%
  filter(prop.genotyped < 0.5) %>%
  pull(locus)
bhs.loci.to.remove #none
#g.bhs <- g.bhs[,-which(getLociNames(g.bhs) %in% bhs.loci.to.remove),]

bhw.loci.to.remove<- sum.bhw %>%
  filter(prop.genotyped < 0.5) %>%
  pull(locus)
bhw.loci.to.remove #none
#g.bhw <- g.bhw[,-which(getLociNames(g.bhw) %in% bhw.loci.to.remove),]

hues.loci.to.remove<- sum.hues %>%
  filter(prop.genotyped < 0.5) %>%
  pull(locus)
hues.loci.to.remove #"locus083" "locus110"
g.hues <- g.hues[,-which(getLociNames(g.hues) %in% hues.loci.to.remove),]

huew.loci.to.remove<- sum.huew %>%
  filter(prop.genotyped < 0.5) %>%
  pull(locus)
huew.loci.to.remove # [1] "locus002" "locus061" "locus077" "locus083" "locus110" "locus126" "locus130" "locus134" "locus165" "locus176" "locus192" "locus210" "locus215" "locus349" "locus360"
g.huew <- g.huew[,-which(getLociNames(g.huew) %in% huew.loci.to.remove),]

sisa.loci.to.remove<- sum.sisa %>%
  filter(prop.genotyped < 0.5) %>%
  pull(locus)
sisa.loci.to.remove # locus255
g.sisa <- g.sisa[,-which(getLociNames(g.sisa) %in% sisa.loci.to.remove),]

mal.loci.to.remove<- sum.mal %>%
  filter(prop.genotyped < 0.5) %>%
  pull(locus)
mal.loci.to.remove #  [1] "Dc00544"  "Dc05864"  "Dc06503"  "Dc10878"  "Dc22883"  "Dc27955"  "Dc31007"  "Dc36767"  "Dc41234"  "locus027" "locus076" "locus093" "locus128" "locus231" "locus255" "locus349"
g.mal <- g.mal[,-which(getLociNames(g.mal) %in% mal.loci.to.remove),]

png.loci.to.remove<- sum.png %>%
  filter(prop.genotyped < 0.5) %>%
  pull(locus)
png.loci.to.remove #  "Dc00544"  "Dc10003"  "locus255"
g.png <- g.png[,-which(getLociNames(g.png) %in% png.loci.to.remove),]

#cleaned so that no locus was genotyped in <50% of individuals
ld.california<-LDgenepop(g.california) 
ld.indonesia<-LDgenepop(g.indonesia)
ld.malaysia<-LDgenepop(g.malaysia)
ld.png<-LDgenepop(g.png)
ld.solomon<-LDgenepop(g.solomon)

ld.bhs <-LDgenepop(g.bhs )
ld.bhw <-LDgenepop(g.bhw )
ld.hues<-LDgenepop(g.hues)
ld.huew<-LDgenepop(g.huew)
ld.sisa<-LDgenepop(g.sisa)
ld.mal <-LDgenepop(g.mal )
ld.png <-LDgenepop(g.png )


ld.strata.list<-list(ld.bhs, ld.bhw, ld.hues, ld.huew, ld.sisa, ld.mal, ld.png)
names(ld.strata.list) <- strata2keep
names(ld.strata.list)

ld.sig.res <- lapply(ld.strata.list, function(s){
  #filter(s, p.value < 0.05) |> 
  s %>% select(c(Locus.1, Locus.2, p.value))
  }) 
  for (i in 1:length(ld.sig.res)){
  names(ld.sig.res[[i]])[3] <- paste0('p.val.',names(ld.strata.list[i]))
  }

ld.sig.res <- ld.sig.res %>%
  map(~ .x %>% drop_na(Locus.1, Locus.2)) %>%
  reduce(full_join)
#There are some loci in individual populations that dont have enough data for LD equations and you end up with NAs in the locus columns
#this removes any row that has an NA in either Locus 1 or Locus 2
head(ld.sig.res)

ld.sig.res <- ld.sig.res |>
  rowwise() %>%
  mutate(
    num.sig.p = sum(c_across(!starts_with("Locus")) < 0.05, na.rm = TRUE),
    num.sig.after.correction = sum(c_across(!starts_with(c("Locus", "num"))) < (0.05/7), na.rm = TRUE)
  ) %>%
  ungroup()
head(ld.sig.res)

#Calculate adjusted p-value using sequential bonferroni correction ("holm")
ld.res.bon<-ld.sig.res 

#%>%
#mutate(locus.pair=paste(Locus.1, Locus.2, sep="_") ) %>%
#  relocate(locus.pair) %>%
#  column_to_rownames(var="locus.pair") %>%
#  select(-c(Locus.1, Locus.2, num.sig.p, num.sig.after.correction)) 

head(ld.res.bon)

ld.res.bon<-ld.res.bon %>% rename_with(.cols=!starts_with(c("Locus", "num")), ~paste0(.x, '.bon.corrected'))

# 1. Calculate the adjusted p-values for the specified columns
# This creates a matrix with just the new values
bon.pvals <- t(apply(ld.res.bon[,-c(1,2, 10, 11)], 1, function(x) {
  p.adjust(x, method = "holm")
}))
head(bon.pvals)

# 2. Merge the dropped columns with the caluclated bonferonni values
# This keeps all other columns in place
ld.res.bon<-cbind(ld.res.bon[,c(1, 2)], bon.pvals)
head(ld.res.bon)

#Add column and count number of strata with adjusted p-value < 0.05
ld.res.adj <- ld.res.bon %>% 
  rowwise() %>%
  mutate(num.sig.p.after.bon.correction = sum(c_across(!starts_with(c("Locus", "num"))) < 0.05, na.rm = TRUE)) %>%
  ungroup()
head(ld.res.adj)


head(ld.sig.res) #LD p-values, and count of p-values that are less than 0.05 AND 0.05/7 (no. of strata)
head(ld.res.adj) #bonferroni corrected p-values and count of p-values that are less than 0.05


write.csv(ld.sig.res, file = 'results-raw/final.uncorrected.ld.sum.csv')
write.csv(ld.res.adj, file = 'results-raw/final.bonferroni-corrected.ld.sum.csv')

save(ld.sig.res, ld.res.adj, ld.strata.list, file = 'data/ld.sig.res.rda')


###Plotting Missingness vs Homozygosity
ind.summary
head(ind.summary)
plot(ind.summary$pct.loci.missing.genotypes, ind.summary$pct.loci.homozygous)
