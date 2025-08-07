library(adegenet)
library(strataG)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(viridis)
library(RColorBrewer)
source('R/functions/DAPC.fit.and.predict.R')

# Load data
load("data/gtypes_final.sams.no.dupes_minReads.20.rda")
g #can use this if continuing along
project <- 'final_sams.no.dupes'
minReads <- 20

#############################################
# stratified by herds
strat.scheme <- 'sampling.location'
g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)

# choose which strata I want to include
#g.stratified <- g.stratified[,,which(getStrataNames(g.stratified) %in% c('CentAm.CA.OR.WA', 'HI.SEAK.NBC', 'MnMx.CA.OR.WA'))]
g.stratified <- g.stratified[,,
                             filter(getNumInd(g.stratified, by.strata = TRUE), num.ind >= 7) |> #changed to 7 as this is the size of my smallest pops 
                               pull(stratum)]

genind.strat <- gtypes2genind(g.stratified)

# cross-validation to choose number of PCs and execute DAPC with the chosen number
mat <- tab(genind.strat, NA.method = 'mean')
grp <- pop(genind.strat)
xval.pop <- xvalDapc(mat, grp, n.pca.max = 200, 
                      training.set = 0.9, result = 'groupMean', center = TRUE,
                      scale = FALSE, n.pca = NULL, n.rep = 30, xval.plot = TRUE)
res <- xval.pop$DAPC
res
scatter.dapc(res)

#predicted stratum based on DAPC
dapc.loov <- predictAllIDsDAPC(genind.strat, n.da = 4, n.pca = 40)
dapc.loov
save(dapc.loov, file = 'results-R/dapc_loov_res.rda')
write.csv(dapc.loov, file='results-raw/dapc_loov_res.csv')

# print results
jpeg(filename = paste0('results-raw/DAPC_scatter_', strat.scheme, '.jpg'))
scatter.dapc(res, scree.da = FALSE)
dev.off()
jpeg(filename = paste0('results-raw/DAPC_assign_', strat.scheme, '.jpg'))
table.value(table(res$assign, grp), col.lab = levels(grp))
dev.off()
mean(as.character(res$assign) == grp) #0.8450704
