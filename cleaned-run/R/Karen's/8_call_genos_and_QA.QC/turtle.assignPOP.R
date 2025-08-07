library(assignPOP)
library(strataG)
library(genepop)
library(tidyverse)
library(ggplot2)
source("R/functions/genepopWrite.KKM.R")
load("data/gtypes_final.sams.no.dupes_minReads.20.rda")
g

strat.scheme <- "sampling.location"

g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)

g.n7 <- g.stratified[,,which(getNumInd(g.stratified, by.strata = TRUE)$num.ind > 6)] #changed to 7 as this is my smallest group size
pws.struct <- pairwiseSummary(pairwiseTest(g.n7, nrep = 10))

strata2keep <- c("CentralCA", "INDONESIA", "MALAYSIA", "PAPUA_NEW_GUINEA", "SOLOMON_ISLANDS")

#g.westcoast <- g.stratified[,,c("CentAm-CA.OR.WA", "MnMx-CA.OR.WA", "HI-WA.SBC")]

#g.pop.infile <- genepopWrite.KKM(g.westcoast, path = "data-raw", filename = "genepop.westcoast.txt")
g.pop.infile <- genepopWrite(g.stratified)

genepop.pop <- read.Genepop(g.pop.infile$fname, pop.names = strata2keep, haploid = FALSE)
genepop.pop$SampleID <- gsub(" ", "_", genepop.pop$SampleID)

genepop.pop.Rd <- reduce.allele(genepop.pop, p = 0.95)

mdl <- "randomForest"
res.dir <- paste0("results-raw/assignPOP/MC.",mdl,"/")
# MCMC Assignment
assign.MC( genepop.pop.Rd, dir=res.dir, train.inds=c(0.7, 0.8, 0.9),
           train.loci=c(0.5,1), loci.sample="fst", iterations=50,
           model=mdl )

accuMC <- accuracy.MC(dir = res.dir)

#plot and save population-level plots
plot_list<-lapply(strata2keep, function(pop_name){
accuracy.plot( accuMC, pop=pop_name) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=1/length(strata2keep),yend=1/length(strata2keep),colour="red",size=1) +
  #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci for", pop_name)+
  #Add a plot title
  theme(plot.title = element_text(size=20, face="bold"))
})
for (i in seq_along(strata2keep)) {
  ggsave(
    filename = paste0("results-raw/assignPOP/MC-cross-validation-plot-", i, ".png"),
    plot = plot_list[[i]],
    width = 8,
    height = 6
  )
}

#Single mc plot with all pops and overall
accuracy.plot( accuMC, pop=c("all", strata2keep)) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=1/length(strata2keep),yend=1/length(strata2keep),colour="red",size=1) +
  #Add a red horizontal line at y = 0.2 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loc")+
  #Add a plot title
  theme(plot.title = element_text(size=20, face="bold"))

# K-fold assignment
assign.kfold( genepop.pop.Rd, k.fold=c(3,4,5), train.loci=c(0.8, 0.9, 1),
              loci.sample="fst", dir=paste0("results-raw/assignPOP/kfold.lda/"), model="lda" )

#Single k-fold plot with all pops and overall
accuMC.kfold.switch <- accuracy.kfold(dir = "results-raw/assignPOP/kfold.lda/")
#names(accuMC.kfold.switch)[1:2] <- c("train.loci", "KF")
#accuMC.kfold.switch <- accuMC.kfold.switch %>% relocate(KF, .before = train.loci)
#strata2keep
accuracy.plot(accuMC.kfold.switch, pop=c("all", strata2keep)) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=1/length(strata2keep),yend=1/length(strata2keep),colour="red",size=1) +
  #Add a red horizontal line at y = 0.2 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci")+
  #Add a plot title
  theme(plot.title = element_text(size=20, face="bold"))

