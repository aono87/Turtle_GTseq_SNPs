###Global Analysis - Example Script for all steps

##Updated Z. cavirostris samples with all duplicates removed (removed individuals using vcftools

##Import new vcf into R for use in adegenet
library(adegenet)
library(vcfR)


# Importing data  ----------------------------------------------

#convert vcf file to vcfR file
zczc125.vcfR<-read.vcfR(("~/Dropbox/Phd/Bioinformatics/bw_ddrad/SNP_files/no_populations/one_snp_per_locus/zczc_snp_files/zczc125.snps.vcf"), verbose = TRUE)
mden43.vcfR<-read.vcfR("~/Dropbox/Phd/Bioinformatics/bw_ddrad_global/SNP_files/no_populations/one_snp_per_locus/mden_snp_files/mden43.snps.vcf", verbose=TRUE)


#convert vcfR files to genlight files
zczc125.gl<-vcfR2genlight(zczc125.vcfR)
#this just shows that it worked and counts missing data, 
zczc125.gl

mden43.gl<-vcfR2genlight(mden43.vcfR)
mden43.gl

pop(mden43.gl)<-as.factor(c("Atlantic",
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
                            "Indopacific"))

pop(mden43.gl)

##Adding population data to zczc127.gl
pop(zczc125.gl)<-as.factor(c("Mediterranean", 
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

popNames(zczc125.gl)
#[1] "Atlantic"      "Indopacific"   "Mediterranean"
nPop(zczc125.gl)
#[1] 3

zczc125.pop.gl<-zczc125.gl
pop(zczc125.pop.gl)<-as.factor(c("Med-Italy", 
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

popNames(zczc125.pop.gl)
# [1] "Atl-Bahamas"  "Atl-East"     "Atl-Madeira" 
# [4] "Atl-Other"    "Atl-Spain"    "Indo-Central"
# [7] "Indo-East"    "Indo-Mexico"  "Indo-South"  
#[10] "Indo-Spac"    "Med-East"     "Med-France"  
#[13] "Med-Italy"   
nPop(zczc125.pop.gl)
#[1] 13

mden43.pop.gl$pop<-as.factor(c("Atl_East",
                               "Atl_Oth",
                               "Atl_Bah",
                               "Atl_Bah",
                               "Atl_Bah",
                               "Atl_Bah",
                               "Atl_Bah",
                               "Atl_East",
                               "Atl_East",
                               "Atl_East",
                               "Atl_Oth",
                               "Atl_East",
                               "Atl_East",
                               "Atl_East",
                               "Atl_Oth",
                               "Atl_East",
                               "Atl_Oth",
                               "Atl_East",
                               "Atl_Bah",
                               "Atl_East",
                               "Atl_East",
                               "Atl_East",
                               "Atl_East",
                               "Atl_East",
                               "Atl_East",
                               "Atl_Bah",
                               "Atl_East",
                               "Atl_Oth",
                               "Indo_Haw",
                               "Indo_Haw",
                               "Indo_Haw",
                               "Indo_Haw",
                               "Indo_Haw",
                               "Indo_Haw",
                               "Indo_Afr",
                               "Indo_Afr",
                               "Indo_Afr",
                               "Indo_Afr",
                               "Indo_Afr",
                               "Indo_Sou",
                               "Indo_Sou",
                               "Indo_Sou",
                               "Indo_Sou"))

# Tess3R ------------------------------------------------------------------


##########taken from Phil Morin's github

###Genotypes Data File
#should be diploid, codominant, format:
#1st column=lab id
#then two columns/locus

###Stratification Data File
#columns: lab id, latitude, longitude, strata names
devtools::install_github("bcm-uga/TESS3_encho_sen")
library(tess3r)
install.packages(c("fields", "RColorBrewer", "mapplots"))
library(fields)
library(RColorBrewer)  
library(mapplots)
BiocManager::install("LEA")
library(LEA)

 #import 012 matrix as genotype file
zczc125.geno<-read.table("~/Dropbox/Phd/Bioinformatics/bw_ddrad/SNP_files/no_populations/one_snp_per_locus/zczc_snp_files/zczc125.snps.012")
zczc125.geno[1:10, 1:5] #first column is generic sample number, need to replace with sample name

zczc125.names<-read.table("~/Dropbox/Phd/Bioinformatics/bw_ddrad/SNP_files/no_populations/one_snp_per_locus/zczc_snp_files/zczc125.snps.012.indv")
zczc125.names #just a column with sample names

zczc125.geno.names<-cbind(zczc125.names, zczc125.geno) #binds column of sample names to genotypes file.
zczc125.geno.names[1:10,1:10]

zczc125.geno.names<-zczc125.geno.names[,-2] #removes column with generic sample id numbers
zczc125.geno.names[1:10,1:10]

zczc125.geno<-zczc125.geno.names #rename file back to zczc125.geno for ease later
dim(zczc125.geno) # 125 25057
rownames(zczc125.geno)<-zczc125.geno[,1] #changing file so that the sample names are the row names, and not a separate column
zczc125.geno[1:10,1:10]
zczc125.geno<-zczc125.geno[,-1]
zczc125.geno[1:10,1:10]
##missing data must be coded as NA. currently missing data stored as a -1
zczc125.geno[zczc125.geno=="-1"]<-NA
zczc125.geno[1:10,1:10]

#Prepare lat/long file
zczc125.coord<-read.table("~/Dropbox/Phd/Bioinformatics/bw_ddrad/SNP_files/no_populations/one_snp_per_locus/zczc_snp_files/zczc125.coord.csv", sep=",")
head(zczc125.coord)
rownames(zczc125.coord)<-zczc125.coord[,1]
head(zczc125.coord)
zczc125.coord<-zczc125.coord[,-1]
head(zczc125.coord)
dim(zczc125.coord)
#make sure coord file is a matrix
is.matrix(zczc125.coord)
zczc125.coord<-as.matrix(zczc125.coord)
#test its right by mapping coordinates on a basic world map
plot(zczc125.coord)
map(add=T, interior=F)

##two files needed:
zczc125.geno
zczc125.coord


# tess-all samples --------------------------------------------------------
#####Running Tess and testing parameters-->changed repatitions from 20 to 100
k<-10
tess.zczc.test<-tess3(X=zczc125.geno, coord=zczc125.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=100, max.iteration=200, keep="best", mask=0)
plot(tess.zczc.test, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
# specify crossentropy with error bars (appears to be very similar to "plot" above)
plot(tess.zczc.test, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score")
for(i in 2:k) {
  Q.matrix <- qmatrix(tess.zczc.test, K = i)
  barplot(Q.matrix, sort.by.Q = FALSE, 
          border = NA, space = 0,
          col.palette = my.palette2,
          xlab = "Individuals", ylab = "Ancestry coefficients") -> bp
  axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .3)
  write.table(Q.matrix, file=paste0("zczc125.qmatrix.k",i,""))
}
for(i in 2:k) {
  Q.matrix <- qmatrix(tess.zczc.test, K = i)
  barplot(Q.matrix, sort.by.Q = TRUE, 
          border = NA, space = 0,
          col.palette = my.palette2,
          xlab = "Individuals", ylab = "Ancestry coefficients") -> bp
  axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .3)
  write.table(Q.matrix, file=paste0("zczc125.qmatrix.test.k",i,""))
}

#####Running Tess with same parameters as before
k<-10
tess3.zczc125<-tess3(X=zczc125.geno, coord=zczc125.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)
# rep = The number of times the algorithm will be repeated for each value of K ############(recommended = 10)
# max.iteration	= the maximum number of iterations of the optimization algorithm. ##########(recommend = 200 (default))
# keep = if "best", only the result with the lowest rmse score will be kept for each value of K. If "all", all results will be kept and returned for each value of K. The second option uses more space in memory.
# mask =  this numeric value is the proportion of masked data when computing the cross-validation criterion (default = 0).

#The plot function generates a plot for root mean-squared errors computed on a subset of loci used for cross-validation.
# specify crossentropy with error bars 
plot(tess3.zczc125, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score")

# Plot cross-validation errors for all values of number of ancestral populations. 
# From Olivier Francçois: For the tess parameters, rep is the number of repetitions. Tess3r is a local optimizer, so increasing rep increases the chance of getting a better local optimum for ancestry coefficients, just like STRUCTURE or other programs. Usually 10 reps are ok, but the program is fast and you could explore more repetitions.  For max.iteration and mask, just use the default values. 

##############visualizing admixture proportions as barplot for various ks

# Generate Colour Palette########
library(wesanderson)
my.col<-wes_palette("Zissou1", 10, type=c("continuous"))
my.palette2<-CreatePalette(my.col)

####Generate barplots and write q.matrix to files (working directory: "~/Dropbox/Phd/Bioinformatics/bw_ddrad/Results/R-summaries/zczc125.tess")
##with samples in order
for(i in 2:k) {
  Q.matrix <- qmatrix(tess3.zczc125, K = i)
  barplot(Q.matrix, sort.by.Q = FALSE, 
          border = NA, space = 0,
          col.palette = my.palette2,
          xlab = "Individuals", ylab = "Ancestry coefficients") -> bp
  axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .3)
  write.table(Q.matrix, file=paste0("zczc125.qmatrix.k",i,""))
}

##with samples in order
for(i in 2:k) {
  Q.matrix.i <- qmatrix(tess3.zczc125, K = i)
  barplot(Q.matrix, sort.by.Q = TRUE, 
          border = NA, space = 0,
          col.palette = my.palette2,
          xlab = "Individuals", ylab = "Ancestry coefficients") -> bp
  axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .3)
}

# Tess for ocean basins##########
#Atlantc
#Indo-Pacific
#Mediterranean
#add ocean data and combine into one table to later subset by ocean basin
zczc125.all<-cbind(zczc125.coord, zczc125.geno)
zczc125.all[1:10,1:10]
zczc125.ocean<-read.table("~/Dropbox/Phd/Bioinformatics/bw_ddrad/SNP_files/no_populations/one_snp_per_locus/zczc_snp_files/zczc125.ocean.txt")
zczc125.ocean
zczc125.ocean.all<-cbind(zczc125.ocean, zczc125.all)
zczc125.ocean.all[1:10,1:10]

#Subdivide dataset by ocean basin
zczc125.ocean.atl<-subset(zczc125.ocean.all, zczc125.ocean.all$V1=="Atlantic")
zczc125.ocean.indo<-subset(zczc125.ocean.all, zczc125.ocean.all$V1=="Indopacific")
zczc125.ocean.med<-subset(zczc125.ocean.all, zczc125.ocean.all$V1=="Mediterranean")

#split back into separate coord/geno files
zczc125.ocean.atl[1:10,1:10]
zczc125.ocean.indo[1:10,1:10]
zczc125.ocean.med[1:10,1:10]

zczc125.atl.geno<-zczc125.ocean.atl[,-(1:3)]
zczc125.atl.coord<-zczc125.ocean.atl[,2:3]
zczc125.atl.coord<-as.matrix(zczc125.atl.coord)

zczc125.indo.geno<-zczc125.ocean.indo[,-(1:3)]
zczc125.indo.coord<-zczc125.ocean.indo[,2:3]
zczc125.indo.coord<-as.matrix(zczc125.indo.coord)

zczc125.med.geno<-zczc125.ocean.med[,-(1:3)]
zczc125.med.coord<-zczc125.ocean.med[,2:3]
zczc125.med.coord<-as.matrix(zczc125.med.coord)


# Atlantic Tess -----------------------------------------------------------


tess3.zczc125.atl<-tess3(X=zczc125.atl.geno, coord=zczc125.atl.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)
tess3.zczc125.atl.100<-tess3(X=zczc125.atl.geno, coord=zczc125.atl.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=100, max.iteration=200, keep="best", mask=0, verbose=F)

par(mfrow=c(3,1))

plot(tess3.zczc125.atl, crossvalid=FALSE, crossentropy=TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score", 
     main = "Z. cavirostris, Atlantic")



for(i in 2:k) {
  Q.matrix <- qmatrix(tess3.zczc125.atl.100, K = i)
  barplot(Q.matrix, sort.by.Q = TRUE, 
          border = NA, space = 0,
          col.palette = my.palette2,
          xlab = "Individuals", ylab = "Ancestry coefficients") -> bp
  axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .3)
  write.table(Q.matrix, file=paste0("zczc125.atl.100.qmatrix.k",i,""))
}


# Indo-Pacific Tess -------------------------------------------------------
tess3.zczc125.indo<-tess3(X=zczc125.indo.geno, coord=zczc125.indo.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)
tess3.zczc125.indo.100<-tess3(X=zczc125.indo.geno, coord=zczc125.indo.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=100, max.iteration=200, keep="best", mask=0, verbose=F)

plot(tess3.zczc125.indo, crossvalid=FALSE, crossentropy=TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score", 
     main = "Z. cavirostris, Indo-Pacific")

for(i in 2:k) {
  Q.matrix <- qmatrix(tess3.zczc125.indo.100, K = i)
  barplot(Q.matrix, sort.by.Q = FALSE, 
          border = NA, space = 0,
          col.palette = my.palette2,
          xlab = "Individuals", ylab = "Ancestry coefficients") -> bp
  axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .3)
  write.table(Q.matrix, file=paste0("zczc125.indo.100.qmatrix.k",i,""))
}


# Mediterranean Tess ------------------------------------------------------


tess3.zczc125.med<-tess3(X=zczc125.med.geno, coord=zczc125.med.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)
tess3.zczc125.med.100<-tess3(X=zczc125.med.geno, coord=zczc125.med.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=100, max.iteration=200, keep="best", mask=0, verbose=F)


plot(tess3.zczc125.med, crossvalid=FALSE, crossentropy=TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score", 
     main = "Z. cavirostris, Mediterranean")



for(i in 2:k) {
  Q.matrix <- qmatrix(tess3.zczc125.med, K = i)
  barplot(Q.matrix, sort.by.Q = FALSE, 
          border = NA, space = 0,
          col.palette = my.palette2,
          xlab = "Individuals", ylab = "Ancestry coefficients") -> bp
  axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .3)
  write.table(Q.matrix, file=paste0("zczc125.med.qmatrix.k",i,""))
}


# Mediterranean-East Tess -------------------------------------------------

####Split italian samples from med and re-run tess
library(maps)
plot(zczc125.atl.coord, pch=19, cex=0.75, col="red",
     xlab="Longitude",
     ylab="Latitude")
map(add=T, interior=F)

plot(zczc125.indo.coord, pch=19, cex=0.75, col="blue",
     xlab="Longitude",
     ylab="Latitude")
map(add=T, interior=F)

plot(zczc125.med.coord, pch=19, cex=0.75, col="green",
     xlab="Longitude",
     ylab="Latitude")
map(add=T, interior=F)

zczc125.med.map<-read_excel("~/Dropbox/Phd/Bioinformatics/bw_ddrad/Results/R-summaries/zczc125.tess/med.maps.xlsx", sheet="med")
head(zczc125.med.map)

#map of mediterranean samples coloured by the cluster assigned with tess and k=3
plot(zczc125.med.map$long, zczc125.med.map$lat, pch=19, cex=0.75, col=zczc125.med.map$pop)
map(add=T, interior=F)


zczc125.ocean.med[1:10,1:10]
zczc125.med.pop<-read.csv("~/Dropbox/Phd/Bioinformatics/bw_ddrad/Results/R-summaries/zczc125.tess/zczc125.med.pop.txt")
zczc125.med.pop
zczc125.pop.med<-cbind(zczc125.med.pop, zczc125.ocean.med)
zczc125.pop.med[1:10,1:10]

zczc125.pop.east.med<-subset(zczc125.pop.med, zczc125.pop.med$pop!=1)
zczc125.pop.east.med[1:10,1:10]

zczc125.med.east.geno<-zczc125.pop.east.med[,-(1:4)]
zczc125.med.east.geno[1:10,1:10]

zczc125.med.east.coord<-zczc125.pop.east.med[,3:4]
zczc125.med.east.coord[1:10,1:2]

zczc125.med.east.coord<-as.matrix(zczc125.med.east.coord)
plot(zczc125.med.east.coord)
map(add=T, interior=F)

tess3.zczc125.med.east<-tess3(X=zczc125.med.east.geno, coord=zczc125.med.east.coord, K=1:6, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)


plot(tess3.zczc125.med.east, crossvalid=FALSE, crossentropy=TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score", 
     main = "Z. cavirostris, Mediterranean-East")


for(i in 2:6) {
  Q.matrix <- qmatrix(tess3.zczc125.med.east, K = i)
  barplot(Q.matrix, sort.by.Q = FALSE, 
          border = NA, space = 0,
          col.palette = my.palette2,
          xlab = "Individuals", ylab = "Ancestry coefficients") -> bp
  axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .3)
  write.table(Q.matrix, file=paste0("zczc125.med.east.qmatrix.k",i,""))
}

zczc125.med.map<-read_excel("~/Dropbox/Phd/Bioinformatics/bw_ddrad/Results/R-summaries/zczc125.tess/med.maps.xlsx", sheet="k=3")
head(zczc125.med.map)
zczc125.med.map

#map of mediterranean samples coloured by the cluster assigned with tess and k=3
plot(zczc125.med.map$long, zczc125.med.map$lat, pch=19, cex=0.75, col=zczc125.med.map$pop)
map(add=T, interior=F)

# VCFTOOLS Codes #####################################
##Need to export data from vcftools as a 012 file for tess3r (this is the kind of data matrix we need
##vcftools code (in linux, not R)
##change directory to where the snp file is
## vcftools --vcf zczc125.snps.vcf --out zczc125.snps --012 --recode-INFO-all

##heterozygosity and F in vcftools
#  vcftools --vcf zczc125.snps.vcf --out zczc125.snps --het --recode-INFO-all

#Explore vcftools results: depth, heterozygosities and F (inbreeding coefficient)
sum<-read.table("~/Dropbox/Phd/Bioinformatics/bw_ddrad/Results/R-summaries/r-working-summaries.csv", sep=",", header=TRUE)
sum.zc<-subset(sum, sum$species=="Zcav")

sum.zc$ocean.basin<-as.factor(sum.zc$ocean.basin)
sum.zc$pop1<-as.character(sum.zc$pop1)
sum.zc$pop1<-as.factor(sum.zc$pop1)
sum.zc$pop2<-as.character(sum.zc$pop2)
sum.zc$pop2<-as.factor(sum.zc$pop2)

boxplot(sum.zc$f ~ sum.zc$pop1)
boxplot(sum.zc$f ~ sum.zc$pop2)
boxplot(sum.zc$f ~ sum.zc$ocean.basin)


summary(sum.zc$f ~ sum.zc$pop1)
summary(sum.zc$f ~ sum.zc$pop2)
summary(sum.zc$f ~ sum.zc$ocean.basin)

library(psych)
describeBy(sum.zc$f, sum.zc$ocean.basin, mat=TRUE, digits=3)
describeBy(sum.zc$f, sum.zc$pop1, mat=TRUE, digits=3)
describeBy(sum.zc$f, sum.zc$pop2, mat=TRUE, digits=3)

obs.het.ind<-(1-sum.zc$obs.homo.snps/sum.zc$no.snps.ind)
sum.zc<-cbind(sum.zc, obs.het.ind)

describeBy(sum.zc$obs.het.ind, sum.zc$ocean.basin, mat=TRUE, digits=3)
describeBy(sum.zc$obs.het.ind, sum.zc$pop1, mat=TRUE, digits=3)
describeBy(sum.zc$obs.het.ind, sum.zc$pop2, mat=TRUE, digits=3)

obs.het.tot<-(1-sum.zc$obs.homo.snps/sum.zc$no.snps.total)
sum.zc<-cbind(sum.zc, obs.het.tot)

describeBy(sum.zc$obs.het.tot, sum.zc$ocean.basin, mat=TRUE, digits=3)
describeBy(sum.zc$obs.het.tot, sum.zc$pop1, mat=TRUE, digits=3)
describeBy(sum.zc$obs.het.tot, sum.zc$pop2, mat=TRUE, digits=3)

describeBy(sum.zc$no.snps.ind, sum.zc$ocean.basin, mat=TRUE, digits=3)
describeBy(sum.zc$no.snps.ind, sum.zc$pop1, mat=TRUE, digits=3)
describeBy(sum.zc$no.snps.ind, sum.zc$pop2, mat=TRUE, digits=3)

describeBy(sum.zc$mean.site.depth, sum.zc$ocean.basin, mat=TRUE, digits=3)
describeBy(sum.zc$mean.site.depth, sum.zc$pop1, mat=TRUE, digits=3)
describeBy(sum.zc$mean.site.depth, sum.zc$pop2, mat=TRUE, digits=3)

describeBy(sum.zc$prop.missing.snps, sum.zc$ocean.basin, mat=TRUE, digits=3)
describeBy(sum.zc$prop.missing.snps, sum.zc$pop1, mat=TRUE, digits=3)
describeBy(sum.zc$prop.missing.snps, sum.zc$pop2, mat=TRUE, digits=3)

# BASIC Trees -------------------------------------------------------------------

###NJ tree
library(adegenet)
library(ape)
library(tidyverse)
BiocManager::install("ggtree")
library(ggtree)
library(ggplot2)

glPlot(zczc125.gl)

zczc125.tree<-nj(dist(as.matrix(zczc125.gl)))
plot(zczc125.tree, typ="fan", cex=0.7, show.tip.label=FALSE, use.edge.length=TRUE, no.margin=TRUE, )


ggplot(zczc125.tree) +
  geom_tree() +
  theme_tree()

ggtree(zczc125.tree, layout="fan")+
  geom_point(aes(color=pop(zczc125.gl)))



zczc125.tree.tib<-as_tibble(zczc125.tree)
zczc125.tree.tib

d<-tibble(ocean.basin=pop(zczc125.gl))
zczc.tib.tree<- add_column(zczc125.tree.tib, d)

plot.phylo(x=zczc125.tree, typ="fan", show.tip=FALSE)
zczc125.pops<-zczc125.gl$pop
tiplabels(text=indNames(zczc125.gl), bg=(col=zczc125.pops), cex=0.4, col="white")

zczc125.pop.tree<-nj(dist(as.matrix(zczc125.pop.gl)))

plot.phylo(x=zczc125.pop.tree, typ="fan", show.tip=FALSE)
zczc125.pop.pops<-zczc125.pop.gl$pop
tiplabels(text=indNames(zczc125.pop.gl), bg=(col=zczc125.pop.pops), cex=0.4, col="white")

# Tess plots --------------------------------------------------------------
librar(readxl)
zczc125.q.data <- read_excel("Results/global_paper/structure-plots/q-matrix-data.xlsx", sheet = "zczc125")
library(wesanderson)


zc.facets.ocean<-c("Atlantic", "Indo-Pacific", "Mediterranean")
names(zc.facets.ocean)<-c("Atlantic", "Indopacific", "Mediterranean")
md.facets.ocean<-c("Atlantic", "Indo-Pacific")
names(md.facets.ocean)<-c("Atlantic", "Indopacific")


dpi=600    #pixels per square inch
png("zczc125.qplot.png", width=6*dpi, height=3*dpi, res=dpi)
ggplot(zczc125.q.data, aes( x = reorder(filename, -sort), y = q.val.ocean, fill=factor(q.num.ocean)) ) + 
  geom_bar( stat = "identity", position = "stack", width=1 ) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Ancestry Coefficient") +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none") +
  facet_grid(~ocean, scales="free", space="free_x", labeller=labeller(ocean=zc.facets.ocean)) +
  scale_fill_manual(values=my_zissou)
dev.off()

# Preparing populations-------------------------------------------------------------------------

##Need to make new .pop.gl files for each possible format

#3 possible atlantic configurations
#Atl1-5 pops
#Atl2-6 pops
#Atl3-5 pops

#3 possible indopacific configurations
#Indo1-4 pops
#Indo2-3 pops
#Indo3-4 pops

#3 possible mediterranean configurations
#Med1-3 pops, remove 1 individual
#Med2-4 pops, remove 1 individual
#Med3-3 pops


# StrataG -----------------------------------------------------------------
#Need files to be in genind format
library(dartR)
library(strataG)

#convert ocean file to genind
zczc125.gi<-gl2gi(zczc125.gl)
zczc125.pop1.gi<-gl2gi(zczc125.pop1.gl)
zczc125.pop2.gi<-gl2gi(zczc125.pop2.gl)

# zczc125.opt.a.gi<-gl2gi(zczc125.pop.opt.a.gl)
# zczc125.opt.b.gi<-gl2gi(zczc125.pop.opt.b.gl)
# zczc125.opt.c.gi<-gl2gi(zczc125.pop.opt.c.gl)

#convert from genind to genotypes file

zczc125.gt<-genind2gtypes(zczc125.gi)
zczc125.pop1.gt<-genind2gtypes(zczc125.pop1.gi)
zczc125.pop2.gt<-genind2gtypes(zczc125.pop2.gi)

# zczc125.opt.a.gt<-genind2gtypes(zczc125.opt.a.gi)
# zczc125.opt.b.gt<-genind2gtypes(zczc125.opt.b.gi)
# zczc125.opt.c.gt<-genind2gtypes(zczc125.opt.c.gi)

getNumInd(zczc125.gt)
getNumLoci(zczc125.gt)
getNumStrata(zczc125.gt)
getStrataNames(zczc125.gt)
getStrata(zczc125.gt)

#Adding strata
new.strata<-c("Med_A",
              "Med_A",
              "Med_A",
              "Med_A",
              "Med_A",
              "Med_A",
              "Med_A",
              "Med_C",
              "Med_A",
              "Med_C",
              "Med_A",
              "Med_C",
              "Med_A",
              "Med_A",
              "Med_B",
              "Indo_A",
              "Indo_A",
              "Atl_BC",
              "Atl_BB",
              "Med_B",
              "Atl_BC",
              "Med_B",
              "Atl_BC",
              "Indo_C",
              "Med_B",
              "Med_A",
              "Indo_C",
              "Med_B",
              "Med_C",
              "Med_C",
              "Med_B",
              "Med_B",
              "Med_B",
              "Atl_BB",
              "Atl_BB",
              "Atl_BB",
              "Atl_BB",
              "Atl_BB",
              "Atl_BC",
              "Atl_BB",
              "Atl_BB",
              "Atl_BC",
              "Atl_BC",
              "Atl_BC",
              "Atl_BC",
              "Atl_BC",
              "Atl_BB",
              "Atl_BC",
              "Atl_BA",
              "Atl_BC",
              "Atl_BB",
              "Atl_BD",
              "Atl_BC",
              "Atl_BA",
              "Atl_BD",
              "Atl_BB",
              "Atl_BB",
              "Atl_BC",
              "Atl_BC",
              "Atl_BC",
              "Atl_BB",
              "Atl_BA",
              "Atl_BA",
              "Atl_BB",
              "Indo_C",
              "Indo_C",
              "Atl_BA",
              "Med_B",
              "Indo_B",
              "Atl_BA",
              "Indo_C",
              "Atl_BA",
              "Atl_BC",
              "Atl_BB",
              "Med_A",
              "Med_A",
              "Indo_B",
              "Med_A",
              "Med_B",
              "Med_A",
              "Atl_BB",
              "Atl_BA",
              "Atl_BA",
              "Atl_BA",
              "Atl_BC",
              "Atl_BC",
              "Atl_BC",
              "Atl_BB",
              "Atl_BA",
              "Atl_BA",
              "Atl_BA",
              "Atl_BA",
              "Atl_BC",
              "Indo_C",
              "Indo_C",
              "Med_A",
              "Med_A",
              "Indo_C",
              "Indo_C",
              "Atl_BA",
              "Atl_BA",
              "Indo_C",
              "Indo_C",
              "Indo_A",
              "Indo_C",
              "Indo_C",
              "Indo_C",
              "Indo_C",
              "Indo_C",
              "Indo_A",
              "Indo_C",
              "Indo_C",
              "Indo_A",
              "Indo_B",
              "Indo_B",
              "Indo_B",
              "Indo_B",
              "Indo_B",
              "Indo_B",
              "Indo_B",
              "Med_A",
              "Indo_C",
              "Indo_B",
              "Indo_B",
              "Indo_B")
names(new.strata)<-getIndNames(zczc125.gt)
zczc125.fine.gt<-zczc125.gt
setStrata(zczc125.fine.gt)<-new.strata
getStrata(zczc125.fine.gt)
getSchemes(zczc125.fine.gt)




zczc125.Struct<-popStructTest(zczc125.gt, stats=c(statFst), nrep=100)
# N
# Atlantic      55
# Indopacific   36
# Mediterranean 34
# 
# Population structure results:
#   estimate      p.val
# Fst 0.1163269 0.00990099
# 
# Population structure results:
#   pair.label        Fst  Fst.p.val
# 1      Atlantic (55) v. Indopacific (36) 0.01783218 0.00990099
# 2    Atlantic (55) v. Mediterranean (34) 0.16940166 0.00990099
# 3 Indopacific (36) v. Mediterranean (34) 0.19097749 0.00990099

zczc125.pop1.Struct<-popStructTest(zczc125.pop1.gt, stats=c(statFst), nrep=100)
# N
# Atl_AA 11
# Atl_AB  5
# Atl_AC 17
# Atl_AD 20
# Atl_AE  2
# Indo_A  5
# Indo_B 12
# Indo_C 19
# Med_A  19
# Med_B  10
# Med_C   5
# 
# Population structure results:
#   estimate p.val
# Fst 0.109461     1
# 
# Population structure results:
#   pair.label         Fst  Fst.p.val
# 1   Atl_AA (11) v. Atl_AB (5) 0.006053763 0.04123711
# 2  Atl_AA (11) v. Atl_AC (17) 0.012350897 0.00990099
# 3  Atl_AA (11) v. Atl_AD (20) 0.012525243 0.00990099
# 4   Atl_AA (11) v. Atl_AE (2) 0.156176759 0.75000000
# 5   Atl_AA (11) v. Indo_A (5) 0.035350704 0.01020408
# 6  Atl_AA (11) v. Indo_B (12) 0.021663719 0.00990099
# 7  Atl_AA (11) v. Indo_C (19) 0.034876468 0.00990099
# 8   Atl_AA (11) v. Med_A (19) 0.215430571 0.00990099
# 9   Atl_AA (11) v. Med_B (10) 0.233128742 0.00990099
# 10   Atl_AA (11) v. Med_C (5) 0.142783897 0.01010101
# 11  Atl_AB (5) v. Atl_AC (17) 0.012355496 0.02197802
# 12  Atl_AB (5) v. Atl_AD (20) 0.010047364 0.01098901
# 13   Atl_AB (5) v. Atl_AE (2) 0.170494461 0.55555556
# 14   Atl_AB (5) v. Indo_A (5) 0.027942568 0.03260870
# 15  Atl_AB (5) v. Indo_B (12) 0.017939332 0.01052632
# 16  Atl_AB (5) v. Indo_C (19) 0.033259359 0.01111111
# 17   Atl_AB (5) v. Med_A (19) 0.161725820 0.01639344
# 18   Atl_AB (5) v. Med_B (10) 0.169309355 0.01562500
# 19    Atl_AB (5) v. Med_C (5) 0.070547931 0.07692308
# 20 Atl_AC (17) v. Atl_AD (20) 0.005441310 0.00990099
# 21  Atl_AC (17) v. Atl_AE (2) 0.149693209 1.00000000
# 22  Atl_AC (17) v. Indo_A (5) 0.031951483 0.01030928
# 23 Atl_AC (17) v. Indo_B (12) 0.016423868 0.00990099
# 24 Atl_AC (17) v. Indo_C (19) 0.029577293 0.00990099
# 25  Atl_AC (17) v. Med_A (19) 0.200165372 0.00990099
# 26  Atl_AC (17) v. Med_B (10) 0.215272885 0.00990099
# 27   Atl_AC (17) v. Med_C (5) 0.134161394 0.01020408
# 28  Atl_AD (20) v. Atl_AE (2) 0.142523371 1.00000000
# 29  Atl_AD (20) v. Indo_A (5) 0.028864413 0.01111111
# 30 Atl_AD (20) v. Indo_B (12) 0.011147368 0.00990099
# 31 Atl_AD (20) v. Indo_C (19) 0.024311146 0.00990099
# 32  Atl_AD (20) v. Med_A (19) 0.194030427 0.00990099
# 33  Atl_AD (20) v. Med_B (10) 0.208544224 0.00990099
# 34   Atl_AD (20) v. Med_C (5) 0.131194896 0.01052632
# 35   Atl_AE (2) v. Indo_A (5) 0.217182458 0.46666667
# 36  Atl_AE (2) v. Indo_B (12) 0.164388856 1.00000000
# 37  Atl_AE (2) v. Indo_C (19) 0.171880123 1.00000000
# 38   Atl_AE (2) v. Med_A (19) 0.366194796 1.00000000
# 39   Atl_AE (2) v. Med_B (10) 0.420776161 1.00000000
# 40    Atl_AE (2) v. Med_C (5) 0.327572194 0.55555556
# 41  Indo_A (5) v. Indo_B (12) 0.016592605 0.01010101
# 42  Indo_A (5) v. Indo_C (19) 0.017759121 0.01041667
# 43   Indo_A (5) v. Med_A (19) 0.243172741 0.01470588
# 44   Indo_A (5) v. Med_B (10) 0.269851202 0.01315789
# 45    Indo_A (5) v. Med_C (5) 0.166374501 0.02970297
# 46 Indo_B (12) v. Indo_C (19) 0.010024243 0.00990099
# 47  Indo_B (12) v. Med_A (19) 0.208691507 0.00990099
# 48  Indo_B (12) v. Med_B (10) 0.225771151 0.00990099
# 49   Indo_B (12) v. Med_C (5) 0.141220666 0.01020408
# 50  Indo_C (19) v. Med_A (19) 0.211902730 0.00990099
# 51  Indo_C (19) v. Med_B (10) 0.227383922 0.00990099
# 52   Indo_C (19) v. Med_C (5) 0.151797471 0.01030928
# 53   Med_A (19) v. Med_B (10) 0.088840500 0.00990099
# 54    Med_A (19) v. Med_C (5) 0.072676384 0.01923077
# 55    Med_B (10) v. Med_C (5) 0.007056218 0.75609756

zczc125.pop2.Struct<-popStructTest(zczc125.pop2.gt, stats=c(statFst), nrep=100)
# N
# Atl_BA 16
# Atl_BB 17
# Atl_BC 20
# Atl_BD  2
# Indo_A  5
# Indo_B 12
# Indo_C 19
# Med_A  19
# Med_B  10
# Med_C   5
# 
# Population structure results:
#   estimate p.val
# Fst 0.1091136     1
# Population structure results:
#   pair.label         Fst  Fst.p.val
# 1  Atl_BA (16) v. Atl_BB (17) 0.010328401 0.00990099
# 2  Atl_BA (16) v. Atl_BC (20) 0.009738826 0.00990099
# 3   Atl_BA (16) v. Atl_BD (2) 0.144353124 0.50000000
# 4   Atl_BA (16) v. Indo_A (5) 0.030173927 0.01010101
# 5  Atl_BA (16) v. Indo_B (12) 0.018116823 0.00990099
# 6  Atl_BA (16) v. Indo_C (19) 0.031506636 0.00990099
# 7   Atl_BA (16) v. Med_A (19) 0.185598626 0.00990099
# 8   Atl_BA (16) v. Med_B (10) 0.196417339 0.00990099
# 9    Atl_BA (16) v. Med_C (5) 0.114863834 0.01000000
# 10 Atl_BB (17) v. Atl_BC (20) 0.005441310 0.00990099
# 11  Atl_BB (17) v. Atl_BD (2) 0.149693209 1.00000000
# 12  Atl_BB (17) v. Indo_A (5) 0.031951483 0.01000000
# 13 Atl_BB (17) v. Indo_B (12) 0.016423868 0.00990099
# 14 Atl_BB (17) v. Indo_C (19) 0.029577293 0.00990099
# 15  Atl_BB (17) v. Med_A (19) 0.200165372 0.00990099
# 16  Atl_BB (17) v. Med_B (10) 0.215272885 0.00990099
# 17   Atl_BB (17) v. Med_C (5) 0.134161394 0.01063830
# 18  Atl_BC (20) v. Atl_BD (2) 0.142523371 1.00000000
# 19  Atl_BC (20) v. Indo_A (5) 0.028864413 0.01086957
# 20 Atl_BC (20) v. Indo_B (12) 0.011147368 0.00990099
# 21 Atl_BC (20) v. Indo_C (19) 0.024311146 0.00990099
# 22  Atl_BC (20) v. Med_A (19) 0.194030427 0.00990099
# 23  Atl_BC (20) v. Med_B (10) 0.208544224 0.00990099
# 24   Atl_BC (20) v. Med_C (5) 0.131194896 0.01149425
# 25   Atl_BD (2) v. Indo_A (5) 0.217182458 0.36363636
# 26  Atl_BD (2) v. Indo_B (12) 0.164388856 1.00000000
# 27  Atl_BD (2) v. Indo_C (19) 0.171880123 1.00000000
# 28   Atl_BD (2) v. Med_A (19) 0.366194796 1.00000000
# 29   Atl_BD (2) v. Med_B (10) 0.420776161 1.00000000
# 30    Atl_BD (2) v. Med_C (5) 0.327572194 0.56250000
# 31  Indo_A (5) v. Indo_B (12) 0.016592605 0.01052632
# 32  Indo_A (5) v. Indo_C (19) 0.017759121 0.01176471
# 33   Indo_A (5) v. Med_A (19) 0.243172741 0.01449275
# 34   Indo_A (5) v. Med_B (10) 0.269851202 0.01190476
# 35    Indo_A (5) v. Med_C (5) 0.166374501 0.02040816
# 36 Indo_B (12) v. Indo_C (19) 0.010024243 0.00990099
# 37  Indo_B (12) v. Med_A (19) 0.208691507 0.00990099
# 38  Indo_B (12) v. Med_B (10) 0.225771151 0.00990099
# 39   Indo_B (12) v. Med_C (5) 0.141220666 0.01020408
# 40  Indo_C (19) v. Med_A (19) 0.211902730 0.00990099
# 41  Indo_C (19) v. Med_B (10) 0.227383922 0.00990099
# 42   Indo_C (19) v. Med_C (5) 0.151797471 0.01041667
# 43   Med_A (19) v. Med_B (10) 0.088840500 0.01000000
# 44    Med_A (19) v. Med_C (5) 0.072676384 0.02439024
# 45    Med_B (10) v. Med_C (5) 0.007056218 0.48717949

##Trying Fis calcs with Stratag
zczc125.stratag.fis<-statFis(zczc125.gt)
zczc125.stratag.fis$result

#Per locus Heterozygosity
exptdHet(zczc125.gt)
obsvdHet(zczc125.gt)

lapply(seppop(zczc125.gi), function(e) mean(summary(e)$Hobs, na.rm = TRUE))

# Convert to Genepop files ------------------------------------------------


genomic_converter(zczc125.gl, output="genepop")

zczc125.genepop<-readGenepop(infile=("~/Dropbox/Phd/Bioinformatics/bw_ddrad/266_radiator_genomic_converter_20190714@2329/radiator_data_20190714@2329_genepop.gen"))

zczc125.genepop$pop_names<-c("Atlantic", "Indopacific", "Mediterranean")
zczc125.genepop$pop_list

fastDivPart(infile = "~/Dropbox/Phd/Bioinformatics/bw_ddrad/266_radiator_genomic_converter_20190714@2329/radiator_data_20190714@2329_genepop.gen", outfile = "./zczc125.gp.results")

zczc125.hier.fst<-pairwise.fst(zczc125.gi,res.type=c("dist", "matrix"))
zczc125.hier.fst

zczc125.diff<-diffCalc(infile=zczc125.genepop, fst=TRUE, pairwise=TRUE, boots=100)


# OLD StrataG results ---------------------------------------------------------


zczc125.Struct<-popStructTest(zczc125.gt, stats=c(statFst), nrep=1000)
# << gtypes created on 2019-07-04 19:21:26 >>>
#   2019-07-04 19:49:23 : Overall test : 1000 permutations
# 
# N
# Atlantic      55
# Indopacific   36
# Mediterranean 34
# 
# Population structure results:
#   estimate       p.val
# Fst 0.1163269 0.000999001
# 
# 
# <<< gtypes created on 2019-07-04 19:21:26 >>>
#   2019-07-04 20:06:46 : Pairwise tests : 1000 permutations
# 2019-07-04 20:06:46 : Atlantic v. Indopacific 
# 2019-07-04 20:20:16 : Atlantic v. Mediterranean 
# 2019-07-04 20:33:08 : Indopacific v. Mediterranean 
# 
# Population structure results:
#   pair.label        Fst   Fst.p.val
# 1      Atlantic (55) v. Indopacific (36) 0.01783218 0.000999001
# 2    Atlantic (55) v. Mediterranean (34) 0.16940166 0.000999001
# 3 Indopacific (36) v. Mediterranean (34) 0.19097749 0.000999001

zczc125.Struct

# $overall
# $overall$strata.freq
# 
# Atlantic   Indopacific Mediterranean 
# 55            36            34 
# 
# $overall$result
# estimate       p.val
# Fst 0.1163269 0.000999001
# 
# $overall$null.dist
# NULL
# 
# 
# $pairwise
# $pairwise$result
# pair.label    strata.1      strata.2 n.1 n.2        Fst   Fst.p.val
# 1      Atlantic (55) v. Indopacific (36)    Atlantic   Indopacific  55  36 0.01783218 0.000999001
# 2    Atlantic (55) v. Mediterranean (34)    Atlantic Mediterranean  55  34 0.16940166 0.000999001
# 3 Indopacific (36) v. Mediterranean (34) Indopacific Mediterranean  36  34 0.19097749 0.000999001
# 
# $pairwise$pair.mat
# $pairwise$pair.mat$Fst
# Atlantic Indopacific Mediterranean
# Atlantic              NA 0.000999001   0.000999001
# Indopacific   0.01783218          NA   0.000999001
# Mediterranean 0.16940166 0.190977491            NA
# 
# 
# $pairwise$null.dist
# NULL


zczc125.Struct.optA<-popStructTest(zczc125.opt.a.gt, stats=c(statFst), nrep=1000)

# <<< gtypes created on 2019-07-04 19:28:45 >>>
#   2019-07-04 20:51:12 : Overall test : 1000 permutations
# 
# N
# Atl_A  11
# Atl_B   5
# Atl_C  17
# Atl_D  20
# Atl_E   2
# Indo_A  5
# Indo_B  9
# Indo_C  3
# Indo_D 19
# Med_A  19
# Med_B   9
# Med_C   5
# Med_NA  1
# 
# Population structure results:
#   estimate p.val
# Fst 0.1021232     1
# 
# 
# <<< gtypes created on 2019-07-04 19:28:45 >>>
#   2019-07-04 21:05:32 : Pairwise tests : 1000 permutations
# 2019-07-04 21:05:32 : Atl_A v. Atl_B 
# 2019-07-04 21:08:40 : Atl_A v. Atl_C 
# 2019-07-04 21:13:24 : Atl_A v. Atl_D 
# 2019-07-04 21:18:37 : Atl_A v. Atl_E 
# 2019-07-04 21:21:28 : Atl_A v. Indo_A 
# 2019-07-04 21:24:35 : Atl_A v. Indo_B 
# 2019-07-04 21:28:15 : Atl_A v. Indo_C 
# 2019-07-04 21:31:15 : Atl_A v. Indo_D 
# 2019-07-04 21:36:19 : Atl_A v. Med_A 
# 2019-07-04 21:40:39 : Atl_A v. Med_B 
# 2019-07-04 21:44:05 : Atl_A v. Med_C 
# 2019-07-04 21:47:18 : Atl_A v. Med_NA 
# 2019-07-04 21:49:31 : Atl_B v. Atl_C 
# 2019-07-04 21:53:24 : Atl_B v. Atl_D 
# 2019-07-04 21:57:48 : Atl_B v. Atl_E 
# 2019-07-04 21:59:49 : Atl_B v. Indo_A 
# 2019-07-04 22:02:11 : Atl_B v. Indo_B 
# 2019-07-04 22:04:52 : Atl_B v. Indo_C 
# 2019-07-04 22:06:55 : Atl_B v. Indo_D 
# 2019-07-04 22:11:03 : Atl_B v. Med_A 
# 2019-07-04 22:14:12 : Atl_B v. Med_B 
# 2019-07-04 22:16:31 : Atl_B v. Med_C 
# 2019-07-04 22:18:53 : Atl_B v. Med_NA 
# 2019-07-04 22:20:27 : Atl_C v. Atl_D 
# 2019-07-04 22:26:26 : Atl_C v. Atl_E 
# 2019-07-04 22:30:02 : Atl_C v. Indo_A 
# 2019-07-04 22:33:58 : Atl_C v. Indo_B 
# 2019-07-04 22:38:29 : Atl_C v. Indo_C 
# 2019-07-04 22:42:23 : Atl_C v. Indo_D 
# 2019-07-04 22:48:16 : Atl_C v. Med_A 
# 2019-07-04 22:53:35 : Atl_C v. Med_B 
# 2019-07-04 22:57:52 : Atl_C v. Med_C 
# 2019-07-04 23:02:01 : Atl_C v. Med_NA 
# 2019-07-04 23:04:55 : Atl_D v. Atl_E 
# 2019-07-04 23:09:02 : Atl_D v. Indo_A 
# 2019-07-04 23:13:28 : Atl_D v. Indo_B 
# 2019-07-04 23:18:26 : Atl_D v. Indo_C 
# 2019-07-04 23:22:46 : Atl_D v. Indo_D 
# 2019-07-04 23:29:10 : Atl_D v. Med_A 
# 2019-07-04 23:35:01 : Atl_D v. Med_B 
# 2019-07-04 23:39:47 : Atl_D v. Med_C 
# 2019-07-04 23:44:22 : Atl_D v. Med_NA 
# 2019-07-04 23:47:38 : Atl_E v. Indo_A 
# 2019-07-04 23:49:33 : Atl_E v. Indo_B 
# 2019-07-04 23:52:05 : Atl_E v. Indo_C 
# 2019-07-04 23:53:49 : Atl_E v. Indo_D 
# 2019-07-04 23:57:36 : Atl_E v. Med_A 
# 2019-07-05 00:00:22 : Atl_E v. Med_B 
# 2019-07-05 00:02:21 : Atl_E v. Med_C 
# 2019-07-05 00:04:08 : Atl_E v. Med_NA 
# 2019-07-05 00:05:13 : Indo_A v. Indo_B 
# 2019-07-05 00:07:47 : Indo_A v. Indo_C 
# 2019-07-05 00:09:57 : Indo_A v. Indo_D 
# 2019-07-05 00:13:52 : Indo_A v. Med_A 
# 2019-07-05 00:17:06 : Indo_A v. Med_B 
# 2019-07-05 00:19:27 : Indo_A v. Med_C 
# 2019-07-05 00:21:48 : Indo_A v. Med_NA 
# 2019-07-05 00:23:20 : Indo_B v. Indo_C 
# 2019-07-05 00:25:58 : Indo_B v. Indo_D 
# 2019-07-05 00:30:29 : Indo_B v. Med_A 
# 2019-07-05 00:34:32 : Indo_B v. Med_B 
# 2019-07-05 00:37:38 : Indo_B v. Med_C 
# 2019-07-05 00:40:36 : Indo_B v. Med_NA 
# 2019-07-05 00:42:30 : Indo_C v. Indo_D 
# 2019-07-05 00:46:12 : Indo_C v. Med_A 
# 2019-07-05 00:49:13 : Indo_C v. Med_B 
# 2019-07-05 00:51:33 : Indo_C v. Med_C 
# 2019-07-05 00:53:38 : Indo_C v. Med_NA 
# 2019-07-05 00:54:59 : Indo_D v. Med_A 
# 2019-07-05 01:00:30 : Indo_D v. Med_B 
# 2019-07-05 01:04:54 : Indo_D v. Med_C 
# 2019-07-05 01:09:12 : Indo_D v. Med_NA 
# 2019-07-05 01:12:12 : Med_A v. Med_B 
# 2019-07-05 01:14:53 : Med_A v. Med_C 
# 2019-07-05 01:17:47 : Med_A v. Med_NA 
# 2019-07-05 01:19:38 : Med_B v. Med_C 
# 2019-07-05 01:21:46 : Med_B v. Med_NA 
# 2019-07-05 01:23:09 : Med_C v. Med_NA 
# 
# Population structure results:
#   pair.label           Fst   Fst.p.val
# 1    Atl_A (11) v. Atl_B (5)  0.0060537626 0.033023736
# 2   Atl_A (11) v. Atl_C (17)  0.0123508969 0.000999001
# 3   Atl_A (11) v. Atl_D (20)  0.0125252432 0.000999001
# 4    Atl_A (11) v. Atl_E (2)  0.1561767587 0.617647059
# 5   Atl_A (11) v. Indo_A (5)  0.0353507037 0.001031992
# 6   Atl_A (11) v. Indo_B (9)  0.0207338315 0.001000000
# 7   Atl_A (11) v. Indo_C (3)  0.0246246807 0.029411765
# 8  Atl_A (11) v. Indo_D (19)  0.0348764680 0.000999001
# 9   Atl_A (11) v. Med_A (19)  0.2154305711 0.000999001
# 10   Atl_A (11) v. Med_B (9)  0.2339863115 0.000999001
# 11   Atl_A (11) v. Med_C (5)  0.1427838972 0.001019368
# 12  Atl_A (11) v. Med_NA (1) -0.0057109451 1.000000000
# 13   Atl_B (5) v. Atl_C (17)  0.0123554964 0.011603376
# 14   Atl_B (5) v. Atl_D (20)  0.0100473637 0.001104972
# 15    Atl_B (5) v. Atl_E (2)  0.1704944606 0.395973154
# 16   Atl_B (5) v. Indo_A (5)  0.0279425680 0.012698413
# 17   Atl_B (5) v. Indo_B (9)  0.0169356351 0.001022495
# 18   Atl_B (5) v. Indo_C (3)  0.0118612630 0.163841808
# 19  Atl_B (5) v. Indo_D (19)  0.0332593592 0.001106195
# 20   Atl_B (5) v. Med_A (19)  0.1617258197 0.001666667
# 21    Atl_B (5) v. Med_B (9)  0.1704527941 0.004126547
# 22    Atl_B (5) v. Med_C (5)  0.0705479308 0.108671789
# 23   Atl_B (5) v. Med_NA (1) -0.1162041167 1.000000000
# 24  Atl_C (17) v. Atl_D (20)  0.0054413103 0.000999001
# 25   Atl_C (17) v. Atl_E (2)  0.1496932089 1.000000000
# 26  Atl_C (17) v. Indo_A (5)  0.0319514832 0.001034126
# 27  Atl_C (17) v. Indo_B (9)  0.0155802166 0.000999001
# 28  Atl_C (17) v. Indo_C (3)  0.0202361282 0.015748031
# 29 Atl_C (17) v. Indo_D (19)  0.0295772929 0.000999001
# 30  Atl_C (17) v. Med_A (19)  0.2001653718 0.000999001
# 31   Atl_C (17) v. Med_B (9)  0.2166605028 0.000999001
# 32   Atl_C (17) v. Med_C (5)  0.1341613936 0.001025641
# 33  Atl_C (17) v. Med_NA (1) -0.0114294505 1.000000000
# 34   Atl_D (20) v. Atl_E (2)  0.1425233712 1.000000000
# 35  Atl_D (20) v. Indo_A (5)  0.0288644132 0.001067236
# 36  Atl_D (20) v. Indo_B (9)  0.0099772541 0.000999001
# 37  Atl_D (20) v. Indo_C (3)  0.0160492715 0.030303030
# 38 Atl_D (20) v. Indo_D (19)  0.0243111460 0.000999001
# 39  Atl_D (20) v. Med_A (19)  0.1940304270 0.000999001
# 40   Atl_D (20) v. Med_B (9)  0.2100720643 0.000999001
# 41   Atl_D (20) v. Med_C (5)  0.1311948961 0.001047120
# 42  Atl_D (20) v. Med_NA (1) -0.0148096016 1.000000000
# 43   Atl_E (2) v. Indo_A (5)  0.2171824580 0.364238411
# 44   Atl_E (2) v. Indo_B (9)  0.1713939032 1.000000000
# 45   Atl_E (2) v. Indo_C (3)  0.2382027484 0.473988439
# 46  Atl_E (2) v. Indo_D (19)  0.1718801226 0.538461538
# 47   Atl_E (2) v. Med_A (19)  0.3661947959 1.000000000
# 48    Atl_E (2) v. Med_B (9)  0.4226191699 1.000000000
# 49    Atl_E (2) v. Med_C (5)  0.3275721940 0.556603774
# 50   Atl_E (2) v. Med_NA (1)  0.5976764450 0.318681319
# 51  Indo_A (5) v. Indo_B (9)  0.0197375417 0.001024590
# 52  Indo_A (5) v. Indo_C (3)  0.0127285537 0.120192308
# 53 Indo_A (5) v. Indo_D (19)  0.0177591205 0.001085776
# 54  Indo_A (5) v. Med_A (19)  0.2431727410 0.001483680
# 55   Indo_A (5) v. Med_B (9)  0.2707443513 0.001141553
# 56   Indo_A (5) v. Med_C (5)  0.1663745010 0.012257406
# 57  Indo_A (5) v. Med_NA (1)  0.0153360476 1.000000000
# 58  Indo_B (9) v. Indo_C (3)  0.0042176745 0.103448276
# 59 Indo_B (9) v. Indo_D (19)  0.0118140808 0.000999001
# 60  Indo_B (9) v. Med_A (19)  0.2136973345 0.001000000
# 61   Indo_B (9) v. Med_B (9)  0.2339840081 0.000999001
# 62   Indo_B (9) v. Med_C (5)  0.1429066926 0.002081165
# 63  Indo_B (9) v. Med_NA (1) -0.0038264395 1.000000000
# 64 Indo_C (3) v. Indo_D (19)  0.0103358117 0.056603774
# 65  Indo_C (3) v. Med_A (19)  0.2472006399 0.031250000
# 66   Indo_C (3) v. Med_B (9)  0.2843033232 0.026490066
# 67   Indo_C (3) v. Med_C (5)  0.1654352704 0.034951456
# 68  Indo_C (3) v. Med_NA (1) -0.0005747696 1.000000000
# 69 Indo_D (19) v. Med_A (19)  0.2119027305 0.000999001
# 70  Indo_D (19) v. Med_B (9)  0.2290537634 0.001000000
# 71  Indo_D (19) v. Med_C (5)  0.1517974708 0.001064963
# 72 Indo_D (19) v. Med_NA (1)  0.0125756976 1.000000000
# 73   Med_A (19) v. Med_B (9)  0.0899794412 0.001029866
# 74   Med_A (19) v. Med_C (5)  0.0726763837 0.002341920
# 75  Med_A (19) v. Med_NA (1)  0.0046267414 1.000000000
# 76    Med_B (9) v. Med_C (5)  0.0077735819 0.508108108
# 77   Med_B (9) v. Med_NA (1) -0.0183255653 1.000000000
# 78   Med_C (5) v. Med_NA (1) -0.1528296141 1.000000000
# 
# There were 43 warnings (use warnings() to see them)

zczc125.Struct.optA

# $overall
# $overall$strata.freq
# 
# Atl_A  Atl_B  Atl_C  Atl_D  Atl_E Indo_A Indo_B Indo_C Indo_D  Med_A  Med_B  Med_C Med_NA 
# 11      5     17     20      2      5      9      3     19     19      9      5      1 
# 
# $overall$result
# estimate p.val
# Fst 0.1021232     1
# 
# $overall$null.dist
# NULL
# 
# 
# $pairwise
# $pairwise$result
# pair.label strata.1 strata.2 n.1 n.2           Fst   Fst.p.val
# 1    Atl_A (11) v. Atl_B (5)    Atl_A    Atl_B  11   5  0.0060537626 0.033023736
# 2   Atl_A (11) v. Atl_C (17)    Atl_A    Atl_C  11  17  0.0123508969 0.000999001
# 3   Atl_A (11) v. Atl_D (20)    Atl_A    Atl_D  11  20  0.0125252432 0.000999001
# 4    Atl_A (11) v. Atl_E (2)    Atl_A    Atl_E  11   2  0.1561767587 0.617647059
# 5   Atl_A (11) v. Indo_A (5)    Atl_A   Indo_A  11   5  0.0353507037 0.001031992
# 6   Atl_A (11) v. Indo_B (9)    Atl_A   Indo_B  11   9  0.0207338315 0.001000000
# 7   Atl_A (11) v. Indo_C (3)    Atl_A   Indo_C  11   3  0.0246246807 0.029411765
# 8  Atl_A (11) v. Indo_D (19)    Atl_A   Indo_D  11  19  0.0348764680 0.000999001
# 9   Atl_A (11) v. Med_A (19)    Atl_A    Med_A  11  19  0.2154305711 0.000999001
# 10   Atl_A (11) v. Med_B (9)    Atl_A    Med_B  11   9  0.2339863115 0.000999001
# 11   Atl_A (11) v. Med_C (5)    Atl_A    Med_C  11   5  0.1427838972 0.001019368
# 12  Atl_A (11) v. Med_NA (1)    Atl_A   Med_NA  11   1 -0.0057109451 1.000000000
# 13   Atl_B (5) v. Atl_C (17)    Atl_B    Atl_C   5  17  0.0123554964 0.011603376
# 14   Atl_B (5) v. Atl_D (20)    Atl_B    Atl_D   5  20  0.0100473637 0.001104972
# 15    Atl_B (5) v. Atl_E (2)    Atl_B    Atl_E   5   2  0.1704944606 0.395973154
# 16   Atl_B (5) v. Indo_A (5)    Atl_B   Indo_A   5   5  0.0279425680 0.012698413
# 17   Atl_B (5) v. Indo_B (9)    Atl_B   Indo_B   5   9  0.0169356351 0.001022495
# 18   Atl_B (5) v. Indo_C (3)    Atl_B   Indo_C   5   3  0.0118612630 0.163841808
# 19  Atl_B (5) v. Indo_D (19)    Atl_B   Indo_D   5  19  0.0332593592 0.001106195
# 20   Atl_B (5) v. Med_A (19)    Atl_B    Med_A   5  19  0.1617258197 0.001666667
# 21    Atl_B (5) v. Med_B (9)    Atl_B    Med_B   5   9  0.1704527941 0.004126547
# 22    Atl_B (5) v. Med_C (5)    Atl_B    Med_C   5   5  0.0705479308 0.108671789
# 23   Atl_B (5) v. Med_NA (1)    Atl_B   Med_NA   5   1 -0.1162041167 1.000000000
# 24  Atl_C (17) v. Atl_D (20)    Atl_C    Atl_D  17  20  0.0054413103 0.000999001
# 25   Atl_C (17) v. Atl_E (2)    Atl_C    Atl_E  17   2  0.1496932089 1.000000000
# 26  Atl_C (17) v. Indo_A (5)    Atl_C   Indo_A  17   5  0.0319514832 0.001034126
# 27  Atl_C (17) v. Indo_B (9)    Atl_C   Indo_B  17   9  0.0155802166 0.000999001
# 28  Atl_C (17) v. Indo_C (3)    Atl_C   Indo_C  17   3  0.0202361282 0.015748031
# 29 Atl_C (17) v. Indo_D (19)    Atl_C   Indo_D  17  19  0.0295772929 0.000999001
# 30  Atl_C (17) v. Med_A (19)    Atl_C    Med_A  17  19  0.2001653718 0.000999001
# 31   Atl_C (17) v. Med_B (9)    Atl_C    Med_B  17   9  0.2166605028 0.000999001
# 32   Atl_C (17) v. Med_C (5)    Atl_C    Med_C  17   5  0.1341613936 0.001025641
# 33  Atl_C (17) v. Med_NA (1)    Atl_C   Med_NA  17   1 -0.0114294505 1.000000000
# 34   Atl_D (20) v. Atl_E (2)    Atl_D    Atl_E  20   2  0.1425233712 1.000000000
# 35  Atl_D (20) v. Indo_A (5)    Atl_D   Indo_A  20   5  0.0288644132 0.001067236
# 36  Atl_D (20) v. Indo_B (9)    Atl_D   Indo_B  20   9  0.0099772541 0.000999001
# 37  Atl_D (20) v. Indo_C (3)    Atl_D   Indo_C  20   3  0.0160492715 0.030303030
# 38 Atl_D (20) v. Indo_D (19)    Atl_D   Indo_D  20  19  0.0243111460 0.000999001
# 39  Atl_D (20) v. Med_A (19)    Atl_D    Med_A  20  19  0.1940304270 0.000999001
# 40   Atl_D (20) v. Med_B (9)    Atl_D    Med_B  20   9  0.2100720643 0.000999001
# 41   Atl_D (20) v. Med_C (5)    Atl_D    Med_C  20   5  0.1311948961 0.001047120
# 42  Atl_D (20) v. Med_NA (1)    Atl_D   Med_NA  20   1 -0.0148096016 1.000000000
# 43   Atl_E (2) v. Indo_A (5)    Atl_E   Indo_A   2   5  0.2171824580 0.364238411
# 44   Atl_E (2) v. Indo_B (9)    Atl_E   Indo_B   2   9  0.1713939032 1.000000000
# 45   Atl_E (2) v. Indo_C (3)    Atl_E   Indo_C   2   3  0.2382027484 0.473988439
# 46  Atl_E (2) v. Indo_D (19)    Atl_E   Indo_D   2  19  0.1718801226 0.538461538
# 47   Atl_E (2) v. Med_A (19)    Atl_E    Med_A   2  19  0.3661947959 1.000000000
# 48    Atl_E (2) v. Med_B (9)    Atl_E    Med_B   2   9  0.4226191699 1.000000000
# 49    Atl_E (2) v. Med_C (5)    Atl_E    Med_C   2   5  0.3275721940 0.556603774
# 50   Atl_E (2) v. Med_NA (1)    Atl_E   Med_NA   2   1  0.5976764450 0.318681319
# 51  Indo_A (5) v. Indo_B (9)   Indo_A   Indo_B   5   9  0.0197375417 0.001024590
# 52  Indo_A (5) v. Indo_C (3)   Indo_A   Indo_C   5   3  0.0127285537 0.120192308
# 53 Indo_A (5) v. Indo_D (19)   Indo_A   Indo_D   5  19  0.0177591205 0.001085776
# 54  Indo_A (5) v. Med_A (19)   Indo_A    Med_A   5  19  0.2431727410 0.001483680
# 55   Indo_A (5) v. Med_B (9)   Indo_A    Med_B   5   9  0.2707443513 0.001141553
# 56   Indo_A (5) v. Med_C (5)   Indo_A    Med_C   5   5  0.1663745010 0.012257406
# 57  Indo_A (5) v. Med_NA (1)   Indo_A   Med_NA   5   1  0.0153360476 1.000000000
# 58  Indo_B (9) v. Indo_C (3)   Indo_B   Indo_C   9   3  0.0042176745 0.103448276
# 59 Indo_B (9) v. Indo_D (19)   Indo_B   Indo_D   9  19  0.0118140808 0.000999001
# 60  Indo_B (9) v. Med_A (19)   Indo_B    Med_A   9  19  0.2136973345 0.001000000
# 61   Indo_B (9) v. Med_B (9)   Indo_B    Med_B   9   9  0.2339840081 0.000999001
# 62   Indo_B (9) v. Med_C (5)   Indo_B    Med_C   9   5  0.1429066926 0.002081165
# 63  Indo_B (9) v. Med_NA (1)   Indo_B   Med_NA   9   1 -0.0038264395 1.000000000
# 64 Indo_C (3) v. Indo_D (19)   Indo_C   Indo_D   3  19  0.0103358117 0.056603774
# 65  Indo_C (3) v. Med_A (19)   Indo_C    Med_A   3  19  0.2472006399 0.031250000
# 66   Indo_C (3) v. Med_B (9)   Indo_C    Med_B   3   9  0.2843033232 0.026490066
# 67   Indo_C (3) v. Med_C (5)   Indo_C    Med_C   3   5  0.1654352704 0.034951456
# 68  Indo_C (3) v. Med_NA (1)   Indo_C   Med_NA   3   1 -0.0005747696 1.000000000
# 69 Indo_D (19) v. Med_A (19)   Indo_D    Med_A  19  19  0.2119027305 0.000999001
# 70  Indo_D (19) v. Med_B (9)   Indo_D    Med_B  19   9  0.2290537634 0.001000000
# 71  Indo_D (19) v. Med_C (5)   Indo_D    Med_C  19   5  0.1517974708 0.001064963
# 72 Indo_D (19) v. Med_NA (1)   Indo_D   Med_NA  19   1  0.0125756976 1.000000000
# 73   Med_A (19) v. Med_B (9)    Med_A    Med_B  19   9  0.0899794412 0.001029866
# 74   Med_A (19) v. Med_C (5)    Med_A    Med_C  19   5  0.0726763837 0.002341920
# 75  Med_A (19) v. Med_NA (1)    Med_A   Med_NA  19   1  0.0046267414 1.000000000
# 76    Med_B (9) v. Med_C (5)    Med_B    Med_C   9   5  0.0077735819 0.508108108
# 77   Med_B (9) v. Med_NA (1)    Med_B   Med_NA   9   1 -0.0183255653 1.000000000
# 78   Med_C (5) v. Med_NA (1)    Med_C   Med_NA   5   1 -0.1528296141 1.000000000
# 
# $pairwise$pair.mat
# $pairwise$pair.mat$Fst
# Atl_A       Atl_B        Atl_C        Atl_D     Atl_E      Indo_A       Indo_B
# Atl_A            NA  0.03302374  0.000999001  0.000999001 0.6176471 0.001031992  0.001000000
# Atl_B   0.006053763          NA  0.011603376  0.001104972 0.3959732 0.012698413  0.001022495
# Atl_C   0.012350897  0.01235550           NA  0.000999001 1.0000000 0.001034126  0.000999001
# Atl_D   0.012525243  0.01004736  0.005441310           NA 1.0000000 0.001067236  0.000999001
# Atl_E   0.156176759  0.17049446  0.149693209  0.142523371        NA 0.364238411  1.000000000
# Indo_A  0.035350704  0.02794257  0.031951483  0.028864413 0.2171825          NA  0.001024590
# Indo_B  0.020733831  0.01693564  0.015580217  0.009977254 0.1713939 0.019737542           NA
# Indo_C  0.024624681  0.01186126  0.020236128  0.016049272 0.2382027 0.012728554  0.004217674
# Indo_D  0.034876468  0.03325936  0.029577293  0.024311146 0.1718801 0.017759121  0.011814081
# Med_A   0.215430571  0.16172582  0.200165372  0.194030427 0.3661948 0.243172741  0.213697335
# Med_B   0.233986311  0.17045279  0.216660503  0.210072064 0.4226192 0.270744351  0.233984008
# Med_C   0.142783897  0.07054793  0.134161394  0.131194896 0.3275722 0.166374501  0.142906693
# Med_NA -0.005710945 -0.11620412 -0.011429450 -0.014809602 0.5976764 0.015336048 -0.003826439
# Indo_C      Indo_D       Med_A        Med_B        Med_C    Med_NA
# Atl_A   0.0294117647 0.000999001 0.000999001  0.000999001  0.001019368 1.0000000
# Atl_B   0.1638418079 0.001106195 0.001666667  0.004126547  0.108671789 1.0000000
# Atl_C   0.0157480315 0.000999001 0.000999001  0.000999001  0.001025641 1.0000000
# Atl_D   0.0303030303 0.000999001 0.000999001  0.000999001  0.001047120 1.0000000
# Atl_E   0.4739884393 0.538461538 1.000000000  1.000000000  0.556603774 0.3186813
# Indo_A  0.1201923077 0.001085776 0.001483680  0.001141553  0.012257406 1.0000000
# Indo_B  0.1034482759 0.000999001 0.001000000  0.000999001  0.002081165 1.0000000
# Indo_C            NA 0.056603774 0.031250000  0.026490066  0.034951456 1.0000000
# Indo_D  0.0103358117          NA 0.000999001  0.001000000  0.001064963 1.0000000
# Med_A   0.2472006399 0.211902730          NA  0.001029866  0.002341920 1.0000000
# Med_B   0.2843033232 0.229053763 0.089979441           NA  0.508108108 1.0000000
# Med_C   0.1654352704 0.151797471 0.072676384  0.007773582           NA 1.0000000
# Med_NA -0.0005747696 0.012575698 0.004626741 -0.018325565 -0.152829614        NA
# 
# 
# $pairwise$null.dist
# NULL


zczc125.Struct.optB<-popStructTest(zczc125.opt.b.gt, stats=c(statFst), nrep=1000)
# 
# <<< gtypes created on 2019-07-04 19:44:57 >>>
#   2019-07-05 01:28:10 : Overall test : 1000 permutations
# 
# N
# Atl_A  11
# Atl_B   5
# Atl_C  17
# Atl_D   2
# Atl_E  18
# Atl_F   2
# Indo_A  5
# Indo_B 12
# Indo_C 19
# Med_A  19
# Med_B   4
# Med_C   5
# Med_D   5
# Med_NA  1
# 
# Population structure results:
#   estimate p.val
# Fst 0.1002582     1
# 
# 
# <<< gtypes created on 2019-07-04 19:44:57 >>>
#   2019-07-05 01:42:06 : Pairwise tests : 1000 permutations
# 2019-07-05 01:42:06 : Atl_A v. Atl_B 
# 2019-07-05 01:44:53 : Atl_A v. Atl_C 
# 2019-07-05 01:49:27 : Atl_A v. Atl_D 
# 2019-07-05 01:52:20 : Atl_A v. Atl_E 
# 2019-07-05 01:57:13 : Atl_A v. Atl_F 
# 2019-07-05 01:59:52 : Atl_A v. Indo_A 
# 2019-07-05 02:02:50 : Atl_A v. Indo_B 
# 2019-07-05 02:06:40 : Atl_A v. Indo_C 
# 2019-07-05 02:11:43 : Atl_A v. Med_A 
# 2019-07-05 02:16:03 : Atl_A v. Med_B 
# 2019-07-05 02:19:03 : Atl_A v. Med_C 
# 2019-07-05 02:22:14 : Atl_A v. Med_D 
# 2019-07-05 02:25:19 : Atl_A v. Med_NA 
# 2019-07-05 02:27:33 : Atl_B v. Atl_C 
# 2019-07-05 02:31:30 : Atl_B v. Atl_D 
# 2019-07-05 02:33:40 : Atl_B v. Atl_E 
# 2019-07-05 02:37:45 : Atl_B v. Atl_F 
# 2019-07-05 02:39:48 : Atl_B v. Indo_A 
# 2019-07-05 02:42:02 : Atl_B v. Indo_B 
# 2019-07-05 02:45:22 : Atl_B v. Indo_C 
# 2019-07-05 02:49:34 : Atl_B v. Med_A 
# 2019-07-05 02:52:48 : Atl_B v. Med_B 
# 2019-07-05 02:55:00 : Atl_B v. Med_C 
# 2019-07-05 02:57:19 : Atl_B v. Med_D 
# 2019-07-05 02:59:32 : Atl_B v. Med_NA 
# 2019-07-05 03:01:03 : Atl_C v. Atl_D 
# 2019-07-05 03:04:46 : Atl_C v. Atl_E 
# 2019-07-05 03:10:25 : Atl_C v. Atl_F 
# 2019-07-05 03:14:02 : Atl_C v. Indo_A 
# 2019-07-05 03:17:56 : Atl_C v. Indo_B 
# 2019-07-05 03:22:50 : Atl_C v. Indo_C 
# 2019-07-05 03:28:41 : Atl_C v. Med_A 
# 2019-07-05 03:33:58 : Atl_C v. Med_B 
# 2019-07-05 03:37:54 : Atl_C v. Med_C 
# 2019-07-05 03:42:03 : Atl_C v. Med_D 
# 2019-07-05 03:46:07 : Atl_C v. Med_NA 
# 2019-07-05 03:49:00 : Atl_D v. Atl_E 
# 2019-07-05 03:52:52 : Atl_D v. Atl_F 
# 2019-07-05 03:54:29 : Atl_D v. Indo_A 
# 2019-07-05 03:56:24 : Atl_D v. Indo_B 
# 2019-07-05 03:59:23 : Atl_D v. Indo_C 
# 2019-07-05 04:03:17 : Atl_D v. Med_A 
# 2019-07-05 04:06:11 : Atl_D v. Med_B 
# 2019-07-05 04:08:02 : Atl_D v. Med_C 
# 2019-07-05 04:09:59 : Atl_D v. Med_D 
# 2019-07-05 04:11:54 : Atl_D v. Med_NA 
# 2019-07-05 04:13:09 : Atl_E v. Atl_F 
# 2019-07-05 04:16:59 : Atl_E v. Indo_A 
# 2019-07-05 04:21:08 : Atl_E v. Indo_B 
# 2019-07-05 04:26:14 : Atl_E v. Indo_C 
# 2019-07-05 04:32:19 : Atl_E v. Med_A 
# 2019-07-05 04:37:49 : Atl_E v. Med_B 
# 2019-07-05 04:41:59 : Atl_E v. Med_C 
# 2019-07-05 04:46:19 : Atl_E v. Med_D 
# 2019-07-05 04:50:36 : Atl_E v. Med_NA 
# 2019-07-05 04:53:35 : Atl_F v. Indo_A 
# 2019-07-05 04:55:34 : Atl_F v. Indo_B 
# 2019-07-05 04:58:25 : Atl_F v. Indo_C 
# 2019-07-05 05:02:14 : Atl_F v. Med_A 
# 2019-07-05 05:04:57 : Atl_F v. Med_B 
# 2019-07-05 05:06:39 : Atl_F v. Med_C 
# 2019-07-05 05:08:29 : Atl_F v. Med_D 
# 2019-07-05 05:10:14 : Atl_F v. Med_NA 
# 2019-07-05 05:11:25 : Indo_A v. Indo_B 
# 2019-07-05 05:14:35 : Indo_A v. Indo_C 
# 2019-07-05 05:18:28 : Indo_A v. Med_A 
# 2019-07-05 05:21:27 : Indo_A v. Med_B 
# 2019-07-05 05:23:28 : Indo_A v. Med_C 
# 2019-07-05 05:25:37 : Indo_A v. Med_D 
# 2019-07-05 05:27:44 : Indo_A v. Med_NA 
# 2019-07-05 05:29:19 : Indo_B v. Indo_C 
# 2019-07-05 05:34:13 : Indo_B v. Med_A 
# 2019-07-05 05:38:41 : Indo_B v. Med_B 
# 2019-07-05 05:41:51 : Indo_B v. Med_C 
# 2019-07-05 05:45:16 : Indo_B v. Med_D 
# 2019-07-05 05:48:38 : Indo_B v. Med_NA 
# 2019-07-05 05:50:58 : Indo_C v. Med_A 
# 2019-07-05 05:56:27 : Indo_C v. Med_B 
# 2019-07-05 06:00:33 : Indo_C v. Med_C 
# 2019-07-05 06:04:51 : Indo_C v. Med_D 
# 2019-07-05 06:09:06 : Indo_C v. Med_NA 
# 2019-07-05 06:12:05 : Med_A v. Med_B 
# 2019-07-05 06:14:43 : Med_A v. Med_C 
# 2019-07-05 06:17:38 : Med_A v. Med_D 
# 2019-07-05 06:20:22 : Med_A v. Med_NA 
# 2019-07-05 06:22:16 : Med_B v. Med_C 
# 2019-07-05 06:24:10 : Med_B v. Med_D 
# 2019-07-05 06:25:49 : Med_B v. Med_NA 
# 2019-07-05 06:27:00 : Med_C v. Med_D 
# 2019-07-05 06:28:56 : Med_C v. Med_NA 
# 2019-07-05 06:30:17 : Med_D v. Med_NA 
# 
# Population structure results:
#   pair.label          Fst   Fst.p.val
# 1     Atl_A (11) v. Atl_B (5)  0.006053763 0.037190083
# 2    Atl_A (11) v. Atl_C (17)  0.012350897 0.000999001
# 3     Atl_A (11) v. Atl_D (2)  0.020192907 1.000000000
# 4    Atl_A (11) v. Atl_E (18)  0.012921358 0.000999001
# 5     Atl_A (11) v. Atl_F (2)  0.156176759 0.571428571
# 6    Atl_A (11) v. Indo_A (5)  0.035350704 0.002068252
# 7   Atl_A (11) v. Indo_B (12)  0.021663719 0.000999001
# 8   Atl_A (11) v. Indo_C (19)  0.034876468 0.000999001
# 9    Atl_A (11) v. Med_A (19)  0.215430571 0.000999001
# 10    Atl_A (11) v. Med_B (4)  0.208103261 0.001199041
# 11    Atl_A (11) v. Med_C (5)  0.142783897 0.001022495
# 12    Atl_A (11) v. Med_D (5)  0.217791173 0.001037344
# 13   Atl_A (11) v. Med_NA (1) -0.005710945 1.000000000
# 14    Atl_B (5) v. Atl_C (17)  0.012355496 0.013669821
# 15     Atl_B (5) v. Atl_D (2)  0.006868253 1.000000000
# 16    Atl_B (5) v. Atl_E (18)  0.010195723 0.001109878
# 17     Atl_B (5) v. Atl_F (2)  0.170494461 0.330985915
# 18    Atl_B (5) v. Indo_A (5)  0.027942568 0.006309148
# 19   Atl_B (5) v. Indo_B (12)  0.017939332 0.002127660
# 20   Atl_B (5) v. Indo_C (19)  0.033259359 0.001118568
# 21    Atl_B (5) v. Med_A (19)  0.161725820 0.001605136
# 22     Atl_B (5) v. Med_B (4)  0.137318679 0.014120668
# 23     Atl_B (5) v. Med_C (5)  0.070547931 0.082039911
# 24     Atl_B (5) v. Med_D (5)  0.150951343 0.009735744
# 25    Atl_B (5) v. Med_NA (1) -0.116204117 1.000000000
# 26    Atl_C (17) v. Atl_D (2)  0.013369699 1.000000000
# 27   Atl_C (17) v. Atl_E (18)  0.005897420 0.000999001
# 28    Atl_C (17) v. Atl_F (2)  0.149693209 1.000000000
# 29   Atl_C (17) v. Indo_A (5)  0.031951483 0.001043841
# 30  Atl_C (17) v. Indo_B (12)  0.016423868 0.000999001
# 31  Atl_C (17) v. Indo_C (19)  0.029577293 0.000999001
# 32   Atl_C (17) v. Med_A (19)  0.200165372 0.000999001
# 33    Atl_C (17) v. Med_B (4)  0.194587394 0.001305483
# 34    Atl_C (17) v. Med_C (5)  0.134161394 0.001033058
# 35    Atl_C (17) v. Med_D (5)  0.203209740 0.001047120
# 36   Atl_C (17) v. Med_NA (1) -0.011429450 1.000000000
# 37    Atl_D (2) v. Atl_E (18)  0.013569443 1.000000000
# 38     Atl_D (2) v. Atl_F (2)  0.289846678 0.318681319
# 39    Atl_D (2) v. Indo_A (5)  0.044737851 1.000000000
# 40   Atl_D (2) v. Indo_B (12)  0.025761880 1.000000000
# 41   Atl_D (2) v. Indo_C (19)  0.044268028 1.000000000
# 42    Atl_D (2) v. Med_A (19)  0.264287596 1.000000000
# 43     Atl_D (2) v. Med_B (4)  0.280939973 1.000000000
# 44     Atl_D (2) v. Med_C (5)  0.177060026 1.000000000
# 45     Atl_D (2) v. Med_D (5)  0.298208477 1.000000000
# 46    Atl_D (2) v. Med_NA (1) -0.007809595 1.000000000
# 47    Atl_E (18) v. Atl_F (2)  0.145305011 1.000000000
# 48   Atl_E (18) v. Indo_A (5)  0.029355485 0.001070664
# 49  Atl_E (18) v. Indo_B (12)  0.011233763 0.000999001
# 50  Atl_E (18) v. Indo_C (19)  0.024186678 0.000999001
# 51   Atl_E (18) v. Med_A (19)  0.196244556 0.000999001
# 52    Atl_E (18) v. Med_B (4)  0.191286751 0.001367989
# 53    Atl_E (18) v. Med_C (5)  0.132120938 0.001053741
# 54    Atl_E (18) v. Med_D (5)  0.199654119 0.001088139
# 55   Atl_E (18) v. Med_NA (1) -0.013107383 1.000000000
# 56    Atl_F (2) v. Indo_A (5)  0.217182458 0.319148936
# 57   Atl_F (2) v. Indo_B (12)  0.164388856 1.000000000
# 58   Atl_F (2) v. Indo_C (19)  0.171880123 0.666666667
# 59    Atl_F (2) v. Med_A (19)  0.366194796 1.000000000
# 60     Atl_F (2) v. Med_B (4)  0.440930633 0.510948905
# 61     Atl_F (2) v. Med_C (5)  0.327572194 0.404040404
# 62     Atl_F (2) v. Med_D (5)  0.443923657 0.542553191
# 63    Atl_F (2) v. Med_NA (1)  0.597676445 0.327672328
# 64  Indo_A (5) v. Indo_B (12)  0.016592605 0.001035197
# 65  Indo_A (5) v. Indo_C (19)  0.017759121 0.001107420
# 66   Indo_A (5) v. Med_A (19)  0.243172741 0.001574803
# 67    Indo_A (5) v. Med_B (4)  0.242991747 0.005202914
# 68    Indo_A (5) v. Med_C (5)  0.166374501 0.010193680
# 69    Indo_A (5) v. Med_D (5)  0.255365302 0.007486631
# 70   Indo_A (5) v. Med_NA (1)  0.015336048 1.000000000
# 71 Indo_B (12) v. Indo_C (19)  0.010024243 0.000999001
# 72  Indo_B (12) v. Med_A (19)  0.208691507 0.000999001
# 73   Indo_B (12) v. Med_B (4)  0.203731198 0.001186240
# 74   Indo_B (12) v. Med_C (5)  0.141220666 0.001035197
# 75   Indo_B (12) v. Med_D (5)  0.212067202 0.001042753
# 76  Indo_B (12) v. Med_NA (1) -0.004458035 1.000000000
# 77  Indo_C (19) v. Med_A (19)  0.211902730 0.000999001
# 78   Indo_C (19) v. Med_B (4)  0.210306006 0.001375516
# 79   Indo_C (19) v. Med_C (5)  0.151797471 0.001069519
# 80   Indo_C (19) v. Med_D (5)  0.215992250 0.001064963
# 81  Indo_C (19) v. Med_NA (1)  0.012575698 1.000000000
# 82    Med_A (19) v. Med_B (4)  0.081045093 0.012500000
# 83    Med_A (19) v. Med_C (5)  0.072676384 0.002392344
# 84    Med_A (19) v. Med_D (5)  0.092801587 0.005291005
# 85   Med_A (19) v. Med_NA (1)  0.004626741 1.000000000
# 86     Med_B (4) v. Med_C (5) -0.007708296 0.932489451
# 87     Med_B (4) v. Med_D (5)  0.009370310 0.071428571
# 88    Med_B (4) v. Med_NA (1) -0.028441764 1.000000000
# 89     Med_C (5) v. Med_D (5)  0.001557007 0.299492386
# 90    Med_C (5) v. Med_NA (1) -0.152829614 1.000000000
# 91    Med_D (5) v. Med_NA (1) -0.002227887 1.000000000

zczc125.Struct.optB

# $overall
# $overall$strata.freq
# 
# Atl_A  Atl_B  Atl_C  Atl_D  Atl_E  Atl_F Indo_A Indo_B Indo_C  Med_A  Med_B  Med_C  Med_D Med_NA 
# 11      5     17      2     18      2      5     12     19     19      4      5      5      1 
# 
# $overall$result
# estimate p.val
# Fst 0.1002582     1
# 
# $overall$null.dist
# NULL
# 
# 
# $pairwise
# $pairwise$result
# pair.label strata.1 strata.2 n.1 n.2          Fst   Fst.p.val
# 1     Atl_A (11) v. Atl_B (5)    Atl_A    Atl_B  11   5  0.006053763 0.037190083
# 2    Atl_A (11) v. Atl_C (17)    Atl_A    Atl_C  11  17  0.012350897 0.000999001
# 3     Atl_A (11) v. Atl_D (2)    Atl_A    Atl_D  11   2  0.020192907 1.000000000
# 4    Atl_A (11) v. Atl_E (18)    Atl_A    Atl_E  11  18  0.012921358 0.000999001
# 5     Atl_A (11) v. Atl_F (2)    Atl_A    Atl_F  11   2  0.156176759 0.571428571
# 6    Atl_A (11) v. Indo_A (5)    Atl_A   Indo_A  11   5  0.035350704 0.002068252
# 7   Atl_A (11) v. Indo_B (12)    Atl_A   Indo_B  11  12  0.021663719 0.000999001
# 8   Atl_A (11) v. Indo_C (19)    Atl_A   Indo_C  11  19  0.034876468 0.000999001
# 9    Atl_A (11) v. Med_A (19)    Atl_A    Med_A  11  19  0.215430571 0.000999001
# 10    Atl_A (11) v. Med_B (4)    Atl_A    Med_B  11   4  0.208103261 0.001199041
# 11    Atl_A (11) v. Med_C (5)    Atl_A    Med_C  11   5  0.142783897 0.001022495
# 12    Atl_A (11) v. Med_D (5)    Atl_A    Med_D  11   5  0.217791173 0.001037344
# 13   Atl_A (11) v. Med_NA (1)    Atl_A   Med_NA  11   1 -0.005710945 1.000000000
# 14    Atl_B (5) v. Atl_C (17)    Atl_B    Atl_C   5  17  0.012355496 0.013669821
# 15     Atl_B (5) v. Atl_D (2)    Atl_B    Atl_D   5   2  0.006868253 1.000000000
# 16    Atl_B (5) v. Atl_E (18)    Atl_B    Atl_E   5  18  0.010195723 0.001109878
# 17     Atl_B (5) v. Atl_F (2)    Atl_B    Atl_F   5   2  0.170494461 0.330985915
# 18    Atl_B (5) v. Indo_A (5)    Atl_B   Indo_A   5   5  0.027942568 0.006309148
# 19   Atl_B (5) v. Indo_B (12)    Atl_B   Indo_B   5  12  0.017939332 0.002127660
# 20   Atl_B (5) v. Indo_C (19)    Atl_B   Indo_C   5  19  0.033259359 0.001118568
# 21    Atl_B (5) v. Med_A (19)    Atl_B    Med_A   5  19  0.161725820 0.001605136
# 22     Atl_B (5) v. Med_B (4)    Atl_B    Med_B   5   4  0.137318679 0.014120668
# 23     Atl_B (5) v. Med_C (5)    Atl_B    Med_C   5   5  0.070547931 0.082039911
# 24     Atl_B (5) v. Med_D (5)    Atl_B    Med_D   5   5  0.150951343 0.009735744
# 25    Atl_B (5) v. Med_NA (1)    Atl_B   Med_NA   5   1 -0.116204117 1.000000000
# 26    Atl_C (17) v. Atl_D (2)    Atl_C    Atl_D  17   2  0.013369699 1.000000000
# 27   Atl_C (17) v. Atl_E (18)    Atl_C    Atl_E  17  18  0.005897420 0.000999001
# 28    Atl_C (17) v. Atl_F (2)    Atl_C    Atl_F  17   2  0.149693209 1.000000000
# 29   Atl_C (17) v. Indo_A (5)    Atl_C   Indo_A  17   5  0.031951483 0.001043841
# 30  Atl_C (17) v. Indo_B (12)    Atl_C   Indo_B  17  12  0.016423868 0.000999001
# 31  Atl_C (17) v. Indo_C (19)    Atl_C   Indo_C  17  19  0.029577293 0.000999001
# 32   Atl_C (17) v. Med_A (19)    Atl_C    Med_A  17  19  0.200165372 0.000999001
# 33    Atl_C (17) v. Med_B (4)    Atl_C    Med_B  17   4  0.194587394 0.001305483
# 34    Atl_C (17) v. Med_C (5)    Atl_C    Med_C  17   5  0.134161394 0.001033058
# 35    Atl_C (17) v. Med_D (5)    Atl_C    Med_D  17   5  0.203209740 0.001047120
# 36   Atl_C (17) v. Med_NA (1)    Atl_C   Med_NA  17   1 -0.011429450 1.000000000
# 37    Atl_D (2) v. Atl_E (18)    Atl_D    Atl_E   2  18  0.013569443 1.000000000
# 38     Atl_D (2) v. Atl_F (2)    Atl_D    Atl_F   2   2  0.289846678 0.318681319
# 39    Atl_D (2) v. Indo_A (5)    Atl_D   Indo_A   2   5  0.044737851 1.000000000
# 40   Atl_D (2) v. Indo_B (12)    Atl_D   Indo_B   2  12  0.025761880 1.000000000
# 41   Atl_D (2) v. Indo_C (19)    Atl_D   Indo_C   2  19  0.044268028 1.000000000
# 42    Atl_D (2) v. Med_A (19)    Atl_D    Med_A   2  19  0.264287596 1.000000000
# 43     Atl_D (2) v. Med_B (4)    Atl_D    Med_B   2   4  0.280939973 1.000000000
# 44     Atl_D (2) v. Med_C (5)    Atl_D    Med_C   2   5  0.177060026 1.000000000
# 45     Atl_D (2) v. Med_D (5)    Atl_D    Med_D   2   5  0.298208477 1.000000000
# 46    Atl_D (2) v. Med_NA (1)    Atl_D   Med_NA   2   1 -0.007809595 1.000000000
# 47    Atl_E (18) v. Atl_F (2)    Atl_E    Atl_F  18   2  0.145305011 1.000000000
# 48   Atl_E (18) v. Indo_A (5)    Atl_E   Indo_A  18   5  0.029355485 0.001070664
# 49  Atl_E (18) v. Indo_B (12)    Atl_E   Indo_B  18  12  0.011233763 0.000999001
# 50  Atl_E (18) v. Indo_C (19)    Atl_E   Indo_C  18  19  0.024186678 0.000999001
# 51   Atl_E (18) v. Med_A (19)    Atl_E    Med_A  18  19  0.196244556 0.000999001
# 52    Atl_E (18) v. Med_B (4)    Atl_E    Med_B  18   4  0.191286751 0.001367989
# 53    Atl_E (18) v. Med_C (5)    Atl_E    Med_C  18   5  0.132120938 0.001053741
# 54    Atl_E (18) v. Med_D (5)    Atl_E    Med_D  18   5  0.199654119 0.001088139
# 55   Atl_E (18) v. Med_NA (1)    Atl_E   Med_NA  18   1 -0.013107383 1.000000000
# 56    Atl_F (2) v. Indo_A (5)    Atl_F   Indo_A   2   5  0.217182458 0.319148936
# 57   Atl_F (2) v. Indo_B (12)    Atl_F   Indo_B   2  12  0.164388856 1.000000000
# 58   Atl_F (2) v. Indo_C (19)    Atl_F   Indo_C   2  19  0.171880123 0.666666667
# 59    Atl_F (2) v. Med_A (19)    Atl_F    Med_A   2  19  0.366194796 1.000000000
# 60     Atl_F (2) v. Med_B (4)    Atl_F    Med_B   2   4  0.440930633 0.510948905
# 61     Atl_F (2) v. Med_C (5)    Atl_F    Med_C   2   5  0.327572194 0.404040404
# 62     Atl_F (2) v. Med_D (5)    Atl_F    Med_D   2   5  0.443923657 0.542553191
# 63    Atl_F (2) v. Med_NA (1)    Atl_F   Med_NA   2   1  0.597676445 0.327672328
# 64  Indo_A (5) v. Indo_B (12)   Indo_A   Indo_B   5  12  0.016592605 0.001035197
# 65  Indo_A (5) v. Indo_C (19)   Indo_A   Indo_C   5  19  0.017759121 0.001107420
# 66   Indo_A (5) v. Med_A (19)   Indo_A    Med_A   5  19  0.243172741 0.001574803
# 67    Indo_A (5) v. Med_B (4)   Indo_A    Med_B   5   4  0.242991747 0.005202914
# 68    Indo_A (5) v. Med_C (5)   Indo_A    Med_C   5   5  0.166374501 0.010193680
# 69    Indo_A (5) v. Med_D (5)   Indo_A    Med_D   5   5  0.255365302 0.007486631
# 70   Indo_A (5) v. Med_NA (1)   Indo_A   Med_NA   5   1  0.015336048 1.000000000
# 71 Indo_B (12) v. Indo_C (19)   Indo_B   Indo_C  12  19  0.010024243 0.000999001
# 72  Indo_B (12) v. Med_A (19)   Indo_B    Med_A  12  19  0.208691507 0.000999001
# 73   Indo_B (12) v. Med_B (4)   Indo_B    Med_B  12   4  0.203731198 0.001186240
# 74   Indo_B (12) v. Med_C (5)   Indo_B    Med_C  12   5  0.141220666 0.001035197
# 75   Indo_B (12) v. Med_D (5)   Indo_B    Med_D  12   5  0.212067202 0.001042753
# 76  Indo_B (12) v. Med_NA (1)   Indo_B   Med_NA  12   1 -0.004458035 1.000000000
# 77  Indo_C (19) v. Med_A (19)   Indo_C    Med_A  19  19  0.211902730 0.000999001
# 78   Indo_C (19) v. Med_B (4)   Indo_C    Med_B  19   4  0.210306006 0.001375516
# 79   Indo_C (19) v. Med_C (5)   Indo_C    Med_C  19   5  0.151797471 0.001069519
# 80   Indo_C (19) v. Med_D (5)   Indo_C    Med_D  19   5  0.215992250 0.001064963
# 81  Indo_C (19) v. Med_NA (1)   Indo_C   Med_NA  19   1  0.012575698 1.000000000
# 82    Med_A (19) v. Med_B (4)    Med_A    Med_B  19   4  0.081045093 0.012500000
# 83    Med_A (19) v. Med_C (5)    Med_A    Med_C  19   5  0.072676384 0.002392344
# 84    Med_A (19) v. Med_D (5)    Med_A    Med_D  19   5  0.092801587 0.005291005
# 85   Med_A (19) v. Med_NA (1)    Med_A   Med_NA  19   1  0.004626741 1.000000000
# 86     Med_B (4) v. Med_C (5)    Med_B    Med_C   4   5 -0.007708296 0.932489451
# 87     Med_B (4) v. Med_D (5)    Med_B    Med_D   4   5  0.009370310 0.071428571
# 88    Med_B (4) v. Med_NA (1)    Med_B   Med_NA   4   1 -0.028441764 1.000000000
# 89     Med_C (5) v. Med_D (5)    Med_C    Med_D   5   5  0.001557007 0.299492386
# 90    Med_C (5) v. Med_NA (1)    Med_C   Med_NA   5   1 -0.152829614 1.000000000
# 91    Med_D (5) v. Med_NA (1)    Med_D   Med_NA   5   1 -0.002227887 1.000000000
# 
# $pairwise$pair.mat
# $pairwise$pair.mat$Fst
# Atl_A        Atl_B        Atl_C        Atl_D        Atl_E     Atl_F      Indo_A
# Atl_A            NA  0.037190083  0.000999001  1.000000000  0.000999001 0.5714286 0.002068252
# Atl_B   0.006053763           NA  0.013669821  1.000000000  0.001109878 0.3309859 0.006309148
# Atl_C   0.012350897  0.012355496           NA  1.000000000  0.000999001 1.0000000 0.001043841
# Atl_D   0.020192907  0.006868253  0.013369699           NA  1.000000000 0.3186813 1.000000000
# Atl_E   0.012921358  0.010195723  0.005897420  0.013569443           NA 1.0000000 0.001070664
# Atl_F   0.156176759  0.170494461  0.149693209  0.289846678  0.145305011        NA 0.319148936
# Indo_A  0.035350704  0.027942568  0.031951483  0.044737851  0.029355485 0.2171825          NA
# Indo_B  0.021663719  0.017939332  0.016423868  0.025761880  0.011233763 0.1643889 0.016592605
# Indo_C  0.034876468  0.033259359  0.029577293  0.044268028  0.024186678 0.1718801 0.017759121
# Med_A   0.215430571  0.161725820  0.200165372  0.264287596  0.196244556 0.3661948 0.243172741
# Med_B   0.208103261  0.137318679  0.194587394  0.280939973  0.191286751 0.4409306 0.242991747
# Med_C   0.142783897  0.070547931  0.134161394  0.177060026  0.132120938 0.3275722 0.166374501
# Med_D   0.217791173  0.150951343  0.203209740  0.298208477  0.199654119 0.4439237 0.255365302
# Med_NA -0.005710945 -0.116204117 -0.011429450 -0.007809595 -0.013107383 0.5976764 0.015336048
# Indo_B      Indo_C       Med_A        Med_B        Med_C        Med_D    Med_NA
# Atl_A   0.000999001 0.000999001 0.000999001  0.001199041  0.001022495  0.001037344 1.0000000
# Atl_B   0.002127660 0.001118568 0.001605136  0.014120668  0.082039911  0.009735744 1.0000000
# Atl_C   0.000999001 0.000999001 0.000999001  0.001305483  0.001033058  0.001047120 1.0000000
# Atl_D   1.000000000 1.000000000 1.000000000  1.000000000  1.000000000  1.000000000 1.0000000
# Atl_E   0.000999001 0.000999001 0.000999001  0.001367989  0.001053741  0.001088139 1.0000000
# Atl_F   1.000000000 0.666666667 1.000000000  0.510948905  0.404040404  0.542553191 0.3276723
# Indo_A  0.001035197 0.001107420 0.001574803  0.005202914  0.010193680  0.007486631 1.0000000
# Indo_B           NA 0.000999001 0.000999001  0.001186240  0.001035197  0.001042753 1.0000000
# Indo_C  0.010024243          NA 0.000999001  0.001375516  0.001069519  0.001064963 1.0000000
# Med_A   0.208691507 0.211902730          NA  0.012500000  0.002392344  0.005291005 1.0000000
# Med_B   0.203731198 0.210306006 0.081045093           NA  0.932489451  0.071428571 1.0000000
# Med_C   0.141220666 0.151797471 0.072676384 -0.007708296           NA  0.299492386 1.0000000
# Med_D   0.212067202 0.215992250 0.092801587  0.009370310  0.001557007           NA 1.0000000
# Med_NA -0.004458035 0.012575698 0.004626741 -0.028441764 -0.152829614 -0.002227887        NA
# 
# 
# $pairwise$null.dist
# NULL


zczc125.Struct.optC<-popStructTest(zczc125.opt.c.gt, stats=c(statFst))

# <<< gtypes created on 2019-07-04 19:35:47 >>>
#   2019-07-05 10:20:54 : Overall test : 1000 permutations
# 
# N
# Atl_A  16
# Atl_B  17
# Atl_C   2
# Atl_D  18
# Atl_E   2
# Indo_A  5
# Indo_B 10
# Indo_C  2
# Indo_D 19
# Med_A  19
# Med_B  10
# Med_C   5
# 
# Population structure results:
#   estimate p.val
# Fst 0.1100663     1
# 
# 
# <<< gtypes created on 2019-07-04 19:35:47 >>>
#   2019-07-05 10:38:42 : Pairwise tests : 1000 permutations
# 2019-07-05 10:38:42 : Atl_A v. Atl_B 
# 2019-07-05 10:44:08 : Atl_A v. Atl_C 
# 2019-07-05 10:47:46 : Atl_A v. Atl_D 
# 2019-07-05 10:53:27 : Atl_A v. Atl_E 
# 2019-07-05 10:56:59 : Atl_A v. Indo_A 
# 2019-07-05 11:00:52 : Atl_A v. Indo_B 
# 2019-07-05 11:05:24 : Atl_A v. Indo_C 
# 2019-07-05 11:08:58 : Atl_A v. Indo_D 
# 2019-07-05 11:14:44 : Atl_A v. Med_A 
# 2019-07-05 11:19:46 : Atl_A v. Med_B 
# 2019-07-05 11:23:58 : Atl_A v. Med_C 
# 2019-07-05 11:27:54 : Atl_B v. Atl_C 
# 2019-07-05 11:31:21 : Atl_B v. Atl_D 
# 2019-07-05 11:36:58 : Atl_B v. Atl_E 
# 2019-07-05 11:40:37 : Atl_B v. Indo_A 
# 2019-07-05 11:44:32 : Atl_B v. Indo_B 
# 2019-07-05 11:49:12 : Atl_B v. Indo_C 
# 2019-07-05 11:52:53 : Atl_B v. Indo_D 
# 2019-07-05 11:58:45 : Atl_B v. Med_A 
# 2019-07-05 12:04:02 : Atl_B v. Med_B 
# 2019-07-05 12:08:30 : Atl_B v. Med_C 
# 2019-07-05 12:12:38 : Atl_C v. Atl_D 
# 2019-07-05 12:16:38 : Atl_C v. Atl_E 
# 2019-07-05 12:18:16 : Atl_C v. Indo_A 
# 2019-07-05 12:20:21 : Atl_C v. Indo_B 
# 2019-07-05 12:23:05 : Atl_C v. Indo_C 
# 2019-07-05 12:24:47 : Atl_C v. Indo_D 
# 2019-07-05 12:28:41 : Atl_C v. Med_A 
# 2019-07-05 12:31:39 : Atl_C v. Med_B 
# 2019-07-05 12:33:48 : Atl_C v. Med_C 
# 2019-07-05 12:35:44 : Atl_D v. Atl_E 
# 2019-07-05 12:39:35 : Atl_D v. Indo_A 
# 2019-07-05 12:43:41 : Atl_D v. Indo_B 
# 2019-07-05 12:48:33 : Atl_D v. Indo_C 
# 2019-07-05 12:52:27 : Atl_D v. Indo_D 
# 2019-07-05 12:58:34 : Atl_D v. Med_A 
# 2019-07-05 13:04:10 : Atl_D v. Med_B 
# 2019-07-05 13:08:49 : Atl_D v. Med_C 
# 2019-07-05 13:13:10 : Atl_E v. Indo_A 
# 2019-07-05 13:15:07 : Atl_E v. Indo_B 
# 2019-07-05 13:17:46 : Atl_E v. Indo_C 
# 2019-07-05 13:19:22 : Atl_E v. Indo_D 
# 2019-07-05 13:23:11 : Atl_E v. Med_A 
# 2019-07-05 13:25:54 : Atl_E v. Med_B 
# 2019-07-05 13:27:59 : Atl_E v. Med_C 
# 2019-07-05 13:29:51 : Indo_A v. Indo_B 
# 2019-07-05 13:32:43 : Indo_A v. Indo_C 
# 2019-07-05 13:34:45 : Indo_A v. Indo_D 
# 2019-07-05 13:38:41 : Indo_A v. Med_A 
# 2019-07-05 13:41:57 : Indo_A v. Med_B 
# 2019-07-05 13:44:26 : Indo_A v. Med_C 
# 2019-07-05 13:46:47 : Indo_B v. Indo_C 
# 2019-07-05 13:49:30 : Indo_B v. Indo_D 
# 2019-07-05 13:54:14 : Indo_B v. Med_A 
# 2019-07-05 13:58:23 : Indo_B v. Med_B 
# 2019-07-05 14:01:41 : Indo_B v. Med_C 
# 2019-07-05 14:04:43 : Indo_C v. Indo_D 
# 2019-07-05 14:08:34 : Indo_C v. Med_A 
# 2019-07-05 14:11:31 : Indo_C v. Med_B 
# 2019-07-05 14:13:42 : Indo_C v. Med_C 
# 2019-07-05 14:15:37 : Indo_D v. Med_A 
# 2019-07-05 14:21:06 : Indo_D v. Med_B 
# 2019-07-05 14:25:34 : Indo_D v. Med_C 
# 2019-07-05 14:29:51 : Med_A v. Med_B 
# 2019-07-05 14:32:37 : Med_A v. Med_C 
# 2019-07-05 14:35:34 : Med_B v. Med_C 
# 
# Population structure results:
#   pair.label         Fst   Fst.p.val
# 1    Atl_A (16) v. Atl_B (17) 0.010328401 0.000999001
# 2     Atl_A (16) v. Atl_C (2) 0.014813275 1.000000000
# 3    Atl_A (16) v. Atl_D (18) 0.009890580 0.000999001
# 4     Atl_A (16) v. Atl_E (2) 0.144353124 0.263157895
# 5    Atl_A (16) v. Indo_A (5) 0.030173927 0.001025641
# 6   Atl_A (16) v. Indo_B (10) 0.017587706 0.000999001
# 7    Atl_A (16) v. Indo_C (2) 0.021662207 1.000000000
# 8   Atl_A (16) v. Indo_D (19) 0.031506636 0.000999001
# 9    Atl_A (16) v. Med_A (19) 0.185598626 0.000999001
# 10   Atl_A (16) v. Med_B (10) 0.196417339 0.000999001
# 11    Atl_A (16) v. Med_C (5) 0.114863834 0.002051282
# 12    Atl_B (17) v. Atl_C (2) 0.013369699 1.000000000
# 13   Atl_B (17) v. Atl_D (18) 0.005897420 0.000999001
# 14    Atl_B (17) v. Atl_E (2) 0.149693209 1.000000000
# 15   Atl_B (17) v. Indo_A (5) 0.031951483 0.001044932
# 16  Atl_B (17) v. Indo_B (10) 0.016225445 0.000999001
# 17   Atl_B (17) v. Indo_C (2) 0.024324807 1.000000000
# 18  Atl_B (17) v. Indo_D (19) 0.029577293 0.000999001
# 19   Atl_B (17) v. Med_A (19) 0.200165372 0.000999001
# 20   Atl_B (17) v. Med_B (10) 0.215272885 0.000999001
# 21    Atl_B (17) v. Med_C (5) 0.134161394 0.001035197
# 22    Atl_C (2) v. Atl_D (18) 0.013569443 1.000000000
# 23     Atl_C (2) v. Atl_E (2) 0.289846678 0.314685315
# 24    Atl_C (2) v. Indo_A (5) 0.044737851 1.000000000
# 25   Atl_C (2) v. Indo_B (10) 0.026775738 1.000000000
# 26    Atl_C (2) v. Indo_C (2) 0.026338453 0.358641359
# 27   Atl_C (2) v. Indo_D (19) 0.044268028 1.000000000
# 28    Atl_C (2) v. Med_A (19) 0.264287596 1.000000000
# 29    Atl_C (2) v. Med_B (10) 0.307958399 1.000000000
# 30     Atl_C (2) v. Med_C (5) 0.177060026 1.000000000
# 31    Atl_D (18) v. Atl_E (2) 0.145305011 1.000000000
# 32   Atl_D (18) v. Indo_A (5) 0.029355485 0.001062699
# 33  Atl_D (18) v. Indo_B (10) 0.010732590 0.000999001
# 34   Atl_D (18) v. Indo_C (2) 0.022133534 1.000000000
# 35  Atl_D (18) v. Indo_D (19) 0.024186678 0.000999001
# 36   Atl_D (18) v. Med_A (19) 0.196244556 0.000999001
# 37   Atl_D (18) v. Med_B (10) 0.211015074 0.000999001
# 38    Atl_D (18) v. Med_C (5) 0.132120938 0.001059322
# 39    Atl_E (2) v. Indo_A (5) 0.217182458 0.335616438
# 40   Atl_E (2) v. Indo_B (10) 0.168673200 1.000000000
# 41    Atl_E (2) v. Indo_C (2) 0.313620522 0.316683317
# 42   Atl_E (2) v. Indo_D (19) 0.171880123 0.500000000
# 43    Atl_E (2) v. Med_A (19) 0.366194796 1.000000000
# 44    Atl_E (2) v. Med_B (10) 0.420776161 1.000000000
# 45     Atl_E (2) v. Med_C (5) 0.327572194 0.476744186
# 46  Indo_A (5) v. Indo_B (10) 0.019300819 0.001036269
# 47   Indo_A (5) v. Indo_C (2) 0.017122789 1.000000000
# 48  Indo_A (5) v. Indo_D (19) 0.017759121 0.001092896
# 49   Indo_A (5) v. Med_A (19) 0.243172741 0.001562500
# 50   Indo_A (5) v. Med_B (10) 0.269851202 0.001269036
# 51    Indo_A (5) v. Med_C (5) 0.166374501 0.005096840
# 52  Indo_B (10) v. Indo_C (2) 0.011594961 1.000000000
# 53 Indo_B (10) v. Indo_D (19) 0.011070054 0.000999001
# 54  Indo_B (10) v. Med_A (19) 0.212378983 0.000999001
# 55  Indo_B (10) v. Med_B (10) 0.230566955 0.000999001
# 56   Indo_B (10) v. Med_C (5) 0.143121121 0.001039501
# 57  Indo_C (2) v. Indo_D (19) 0.018465180 1.000000000
# 58   Indo_C (2) v. Med_A (19) 0.264142233 1.000000000
# 59   Indo_C (2) v. Med_B (10) 0.309148794 1.000000000
# 60    Indo_C (2) v. Med_C (5) 0.177172862 1.000000000
# 61  Indo_D (19) v. Med_A (19) 0.211902730 0.000999001
# 62  Indo_D (19) v. Med_B (10) 0.227383922 0.000999001
# 63   Indo_D (19) v. Med_C (5) 0.151797471 0.001068376
# 64   Med_A (19) v. Med_B (10) 0.088840500 0.001006036
# 65    Med_A (19) v. Med_C (5) 0.072676384 0.002439024
# 66    Med_B (10) v. Med_C (5) 0.007056218 0.632530120
# 
# There were 39 warnings (use warnings() to see them)

zczc125.Struct.optC

# $overall
# $overall$strata.freq
# 
# Atl_A  Atl_B  Atl_C  Atl_D  Atl_E Indo_A Indo_B Indo_C Indo_D  Med_A  Med_B  Med_C 
# 16     17      2     18      2      5     10      2     19     19     10      5 
# 
# $overall$result
# estimate p.val
# Fst 0.1100663     1
# 
# $overall$null.dist
# NULL
# 
# 
# $pairwise
# $pairwise$result
# pair.label strata.1 strata.2 n.1 n.2         Fst   Fst.p.val
# 1    Atl_A (16) v. Atl_B (17)    Atl_A    Atl_B  16  17 0.010328401 0.000999001
# 2     Atl_A (16) v. Atl_C (2)    Atl_A    Atl_C  16   2 0.014813275 1.000000000
# 3    Atl_A (16) v. Atl_D (18)    Atl_A    Atl_D  16  18 0.009890580 0.000999001
# 4     Atl_A (16) v. Atl_E (2)    Atl_A    Atl_E  16   2 0.144353124 0.263157895
# 5    Atl_A (16) v. Indo_A (5)    Atl_A   Indo_A  16   5 0.030173927 0.001025641
# 6   Atl_A (16) v. Indo_B (10)    Atl_A   Indo_B  16  10 0.017587706 0.000999001
# 7    Atl_A (16) v. Indo_C (2)    Atl_A   Indo_C  16   2 0.021662207 1.000000000
# 8   Atl_A (16) v. Indo_D (19)    Atl_A   Indo_D  16  19 0.031506636 0.000999001
# 9    Atl_A (16) v. Med_A (19)    Atl_A    Med_A  16  19 0.185598626 0.000999001
# 10   Atl_A (16) v. Med_B (10)    Atl_A    Med_B  16  10 0.196417339 0.000999001
# 11    Atl_A (16) v. Med_C (5)    Atl_A    Med_C  16   5 0.114863834 0.002051282
# 12    Atl_B (17) v. Atl_C (2)    Atl_B    Atl_C  17   2 0.013369699 1.000000000
# 13   Atl_B (17) v. Atl_D (18)    Atl_B    Atl_D  17  18 0.005897420 0.000999001
# 14    Atl_B (17) v. Atl_E (2)    Atl_B    Atl_E  17   2 0.149693209 1.000000000
# 15   Atl_B (17) v. Indo_A (5)    Atl_B   Indo_A  17   5 0.031951483 0.001044932
# 16  Atl_B (17) v. Indo_B (10)    Atl_B   Indo_B  17  10 0.016225445 0.000999001
# 17   Atl_B (17) v. Indo_C (2)    Atl_B   Indo_C  17   2 0.024324807 1.000000000
# 18  Atl_B (17) v. Indo_D (19)    Atl_B   Indo_D  17  19 0.029577293 0.000999001
# 19   Atl_B (17) v. Med_A (19)    Atl_B    Med_A  17  19 0.200165372 0.000999001
# 20   Atl_B (17) v. Med_B (10)    Atl_B    Med_B  17  10 0.215272885 0.000999001
# 21    Atl_B (17) v. Med_C (5)    Atl_B    Med_C  17   5 0.134161394 0.001035197
# 22    Atl_C (2) v. Atl_D (18)    Atl_C    Atl_D   2  18 0.013569443 1.000000000
# 23     Atl_C (2) v. Atl_E (2)    Atl_C    Atl_E   2   2 0.289846678 0.314685315
# 24    Atl_C (2) v. Indo_A (5)    Atl_C   Indo_A   2   5 0.044737851 1.000000000
# 25   Atl_C (2) v. Indo_B (10)    Atl_C   Indo_B   2  10 0.026775738 1.000000000
# 26    Atl_C (2) v. Indo_C (2)    Atl_C   Indo_C   2   2 0.026338453 0.358641359
# 27   Atl_C (2) v. Indo_D (19)    Atl_C   Indo_D   2  19 0.044268028 1.000000000
# 28    Atl_C (2) v. Med_A (19)    Atl_C    Med_A   2  19 0.264287596 1.000000000
# 29    Atl_C (2) v. Med_B (10)    Atl_C    Med_B   2  10 0.307958399 1.000000000
# 30     Atl_C (2) v. Med_C (5)    Atl_C    Med_C   2   5 0.177060026 1.000000000
# 31    Atl_D (18) v. Atl_E (2)    Atl_D    Atl_E  18   2 0.145305011 1.000000000
# 32   Atl_D (18) v. Indo_A (5)    Atl_D   Indo_A  18   5 0.029355485 0.001062699
# 33  Atl_D (18) v. Indo_B (10)    Atl_D   Indo_B  18  10 0.010732590 0.000999001
# 34   Atl_D (18) v. Indo_C (2)    Atl_D   Indo_C  18   2 0.022133534 1.000000000
# 35  Atl_D (18) v. Indo_D (19)    Atl_D   Indo_D  18  19 0.024186678 0.000999001
# 36   Atl_D (18) v. Med_A (19)    Atl_D    Med_A  18  19 0.196244556 0.000999001
# 37   Atl_D (18) v. Med_B (10)    Atl_D    Med_B  18  10 0.211015074 0.000999001
# 38    Atl_D (18) v. Med_C (5)    Atl_D    Med_C  18   5 0.132120938 0.001059322
# 39    Atl_E (2) v. Indo_A (5)    Atl_E   Indo_A   2   5 0.217182458 0.335616438
# 40   Atl_E (2) v. Indo_B (10)    Atl_E   Indo_B   2  10 0.168673200 1.000000000
# 41    Atl_E (2) v. Indo_C (2)    Atl_E   Indo_C   2   2 0.313620522 0.316683317
# 42   Atl_E (2) v. Indo_D (19)    Atl_E   Indo_D   2  19 0.171880123 0.500000000
# 43    Atl_E (2) v. Med_A (19)    Atl_E    Med_A   2  19 0.366194796 1.000000000
# 44    Atl_E (2) v. Med_B (10)    Atl_E    Med_B   2  10 0.420776161 1.000000000
# 45     Atl_E (2) v. Med_C (5)    Atl_E    Med_C   2   5 0.327572194 0.476744186
# 46  Indo_A (5) v. Indo_B (10)   Indo_A   Indo_B   5  10 0.019300819 0.001036269
# 47   Indo_A (5) v. Indo_C (2)   Indo_A   Indo_C   5   2 0.017122789 1.000000000
# 48  Indo_A (5) v. Indo_D (19)   Indo_A   Indo_D   5  19 0.017759121 0.001092896
# 49   Indo_A (5) v. Med_A (19)   Indo_A    Med_A   5  19 0.243172741 0.001562500
# 50   Indo_A (5) v. Med_B (10)   Indo_A    Med_B   5  10 0.269851202 0.001269036
# 51    Indo_A (5) v. Med_C (5)   Indo_A    Med_C   5   5 0.166374501 0.005096840
# 52  Indo_B (10) v. Indo_C (2)   Indo_B   Indo_C  10   2 0.011594961 1.000000000
# 53 Indo_B (10) v. Indo_D (19)   Indo_B   Indo_D  10  19 0.011070054 0.000999001
# 54  Indo_B (10) v. Med_A (19)   Indo_B    Med_A  10  19 0.212378983 0.000999001
# 55  Indo_B (10) v. Med_B (10)   Indo_B    Med_B  10  10 0.230566955 0.000999001
# 56   Indo_B (10) v. Med_C (5)   Indo_B    Med_C  10   5 0.143121121 0.001039501
# 57  Indo_C (2) v. Indo_D (19)   Indo_C   Indo_D   2  19 0.018465180 1.000000000
# 58   Indo_C (2) v. Med_A (19)   Indo_C    Med_A   2  19 0.264142233 1.000000000
# 59   Indo_C (2) v. Med_B (10)   Indo_C    Med_B   2  10 0.309148794 1.000000000
# 60    Indo_C (2) v. Med_C (5)   Indo_C    Med_C   2   5 0.177172862 1.000000000
# 61  Indo_D (19) v. Med_A (19)   Indo_D    Med_A  19  19 0.211902730 0.000999001
# 62  Indo_D (19) v. Med_B (10)   Indo_D    Med_B  19  10 0.227383922 0.000999001
# 63   Indo_D (19) v. Med_C (5)   Indo_D    Med_C  19   5 0.151797471 0.001068376
# 64   Med_A (19) v. Med_B (10)    Med_A    Med_B  19  10 0.088840500 0.001006036
# 65    Med_A (19) v. Med_C (5)    Med_A    Med_C  19   5 0.072676384 0.002439024
# 66    Med_B (10) v. Med_C (5)    Med_B    Med_C  10   5 0.007056218 0.632530120
# 
# $pairwise$pair.mat
# $pairwise$pair.mat$Fst
# Atl_A       Atl_B      Atl_C       Atl_D     Atl_E      Indo_A      Indo_B     Indo_C
# Atl_A          NA 0.000999001 1.00000000 0.000999001 0.2631579 0.001025641 0.000999001 1.00000000
# Atl_B  0.01032840          NA 1.00000000 0.000999001 1.0000000 0.001044932 0.000999001 1.00000000
# Atl_C  0.01481327 0.013369699         NA 1.000000000 0.3146853 1.000000000 1.000000000 0.35864136
# Atl_D  0.00989058 0.005897420 0.01356944          NA 1.0000000 0.001062699 0.000999001 1.00000000
# Atl_E  0.14435312 0.149693209 0.28984668 0.145305011        NA 0.335616438 1.000000000 0.31668332
# Indo_A 0.03017393 0.031951483 0.04473785 0.029355485 0.2171825          NA 0.001036269 1.00000000
# Indo_B 0.01758771 0.016225445 0.02677574 0.010732590 0.1686732 0.019300819          NA 1.00000000
# Indo_C 0.02166221 0.024324807 0.02633845 0.022133534 0.3136205 0.017122789 0.011594961         NA
# Indo_D 0.03150664 0.029577293 0.04426803 0.024186678 0.1718801 0.017759121 0.011070054 0.01846518
# Med_A  0.18559863 0.200165372 0.26428760 0.196244556 0.3661948 0.243172741 0.212378983 0.26414223
# Med_B  0.19641734 0.215272885 0.30795840 0.211015074 0.4207762 0.269851202 0.230566955 0.30914879
# Med_C  0.11486383 0.134161394 0.17706003 0.132120938 0.3275722 0.166374501 0.143121121 0.17717286
# Indo_D       Med_A       Med_B       Med_C
# Atl_A  0.000999001 0.000999001 0.000999001 0.002051282
# Atl_B  0.000999001 0.000999001 0.000999001 0.001035197
# Atl_C  1.000000000 1.000000000 1.000000000 1.000000000
# Atl_D  0.000999001 0.000999001 0.000999001 0.001059322
# Atl_E  0.500000000 1.000000000 1.000000000 0.476744186
# Indo_A 0.001092896 0.001562500 0.001269036 0.005096840
# Indo_B 0.000999001 0.000999001 0.000999001 0.001039501
# Indo_C 1.000000000 1.000000000 1.000000000 1.000000000
# Indo_D          NA 0.000999001 0.000999001 0.001068376
# Med_A  0.211902730          NA 0.001006036 0.002439024
# Med_B  0.227383922 0.088840500          NA 0.632530120
# Med_C  0.151797471 0.072676384 0.007056218          NA
# 
# 
# $pairwise$null.dist
# NULL



# DartR: Fst and SUm Stats ------------------------------------------------
library(dartR)
# This script calculates pairwise fst values based on the implementation in the StAMPP package (?stamppFst). It allows to run bootstrap to estimate probability of fst values to be different from zero. For detailed information please check the help pages (?stamppFst). 
#From STAMPP: This function calculates pairwise Fst values along with confidence intervals and p-values between populations according to the method proposed by Wright(1949) and updated by Weir and Cockerham (1984).

#genlight files to use
zczc125.gl
zczc125.pop1.gl
zczc125.pop2.gl

zczc125.fst.ci<-gl.fst.pop(zczc125.gl, nboots=100, percent=95)
zczc125.fst.ci
zczc125.fst.ci$Fsts
#               Mediterranean Indopacific Atlantic
# Mediterranean            NA          NA       NA
# Indopacific       0.1909775          NA       NA
# Atlantic          0.1694017  0.01783218       NA
zczc125.fst.ci$Pvalues
#               Mediterranean Indopacific Atlantic
# Mediterranean            NA          NA       NA
# Indopacific               0          NA       NA
# Atlantic                  0           0       NA
zczc125.fst.boots<-zczc125.fst.ci$Bootstraps[,c(1:2,103:106)]
zczc125.fst.boots
# Population1 Population2 Lower bound CI limit Upper bound CI limit p-value        Fst
# 1 Mediterranean Indopacific           0.18673098           0.19562435       0 0.19097749
# 2 Mediterranean    Atlantic           0.16436397           0.17384722       0 0.16940166
# 3   Indopacific    Atlantic           0.01697198           0.01887862       0 0.01783218

zczc125.pop1.fst.ci<-gl.fst.pop(zczc125.pop1.gl, nboots=100, percent=95)
zczc125.pop1.fst.ci
zczc125.pop1.fst.ci$Fsts
#            Med_A       Med_C     Med_B     Indo_A     Atl_AD     Atl_AC     Indo_C      Atl_AA    Atl_AE      Atl_AB  Indo_B
#Med_A          NA          NA        NA         NA         NA         NA         NA          NA        NA          NA     NA
#Med_C  0.07267638          NA        NA         NA         NA         NA         NA          NA        NA          NA     NA
#Med_B  0.08884050 0.007056218        NA         NA         NA         NA         NA          NA        NA          NA     NA
#Indo_A 0.24317274 0.166374501 0.2698512         NA         NA         NA         NA          NA        NA          NA     NA
#Atl_AD 0.19403043 0.131194896 0.2085442 0.02886441         NA         NA         NA          NA        NA          NA     NA
#Atl_AC 0.20016537 0.134161394 0.2152729 0.03195148 0.00544131         NA         NA          NA        NA          NA     NA
#Indo_C 0.21190273 0.151797471 0.2273839 0.01775912 0.02431115 0.02957729         NA          NA        NA          NA     NA
#Atl_AA 0.21543057 0.142783897 0.2331287 0.03535070 0.01252524 0.01235090 0.03487647          NA        NA          NA     NA
#Atl_AE 0.36619480 0.327572194 0.4207762 0.21718246 0.14252337 0.14969321 0.17188012 0.156176759        NA          NA     NA
#Atl_AB 0.16172582 0.070547931 0.1693094 0.02794257 0.01004736 0.01235550 0.03325936 0.006053763 0.1704945          NA     NA
#Indo_B 0.20869151 0.141220666 0.2257712 0.01659260 0.01114737 0.01642387 0.01002424 0.021663719 0.1643889  0.01793933     NA

zczc125.pop1.fst.ci$Pvalues
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
zczc125.fst.pop1.boots<-zczc125.pop1.fst.ci$Bootstraps[,c(1:2,103:106)]
zczc125.fst.pop1.boots
# Population1 Population2 Lower bound CI limit Upper bound CI limit p-value         Fst
# 1        Med_A       Med_C          0.068673874          0.076971006       0 0.072676384
# 2        Med_A       Med_B          0.085733124          0.093149802       0 0.088840500
# 3        Med_A      Indo_A          0.237272807          0.248864500       0 0.243172741
# 4        Med_A      Atl_AD          0.189081480          0.199912982       0 0.194030427
# 5        Med_A      Atl_AC          0.195057253          0.205369969       0 0.200165372
# 6        Med_A      Indo_C          0.206734542          0.217605745       0 0.211902730
# 7        Med_A      Atl_AA          0.208603387          0.222243778       0 0.215430571
# 8        Med_A      Atl_AE          0.357187443          0.374031307       0 0.366194796
# 9        Med_A      Atl_AB          0.156203818          0.166242503       0 0.161725820
# 10       Med_A      Indo_B          0.202598656          0.213748346       0 0.208691507
# 11       Med_C       Med_B          0.003778618          0.010022805       0 0.007056218
# 12       Med_C      Indo_A          0.159647400          0.172020479       0 0.166374501
# 13       Med_C      Atl_AD          0.125422719          0.138247701       0 0.131194896
# 14       Med_C      Atl_AC          0.129148303          0.140382994       0 0.134161394
# 15       Med_C      Indo_C          0.146138462          0.158968115       0 0.151797471
# 16       Med_C      Atl_AA          0.136716788          0.148757836       0 0.142783897
# 17       Med_C      Atl_AE          0.317070618          0.335548004       0 0.327572194
# 18       Med_C      Atl_AB          0.065665073          0.075612662       0 0.070547931
# 19       Med_C      Indo_B          0.135639640          0.148627936       0 0.141220666
# 20       Med_B      Indo_A          0.261837454          0.278956704       0 0.269851202
# 21       Med_B      Atl_AD          0.201161461          0.216097480       0 0.208544224
# 22       Med_B      Atl_AC          0.207468510          0.222527247       0 0.215272885
# 23       Med_B      Indo_C          0.219783283          0.236256995       0 0.227383922
# 24       Med_B      Atl_AA          0.223740940          0.240380359       0 0.233128742
# 25       Med_B      Atl_AE          0.411658511          0.428181964       0 0.420776161
# 26       Med_B      Atl_AB          0.162696802          0.175384268       0 0.169309355
# 27       Med_B      Indo_B          0.217528380          0.233445079       0 0.225771151
# 28      Indo_A      Atl_AD          0.025573177          0.031101659       0 0.028864413
# 29      Indo_A      Atl_AC          0.029044549          0.034372308       0 0.031951483
# 30      Indo_A      Indo_C          0.015129641          0.020513543       0 0.017759121
# 31      Indo_A      Atl_AA          0.031614129          0.038085391       0 0.035350704
# 32      Indo_A      Atl_AE          0.208588390          0.227028545       0 0.217182458
# 33      Indo_A      Atl_AB          0.024318325          0.031327625       0 0.027942568
# 34      Indo_A      Indo_B          0.013605556          0.019691454       0 0.016592605
# 35      Atl_AD      Atl_AC          0.004600690          0.006423556       0 0.005441310
# 36      Atl_AD      Indo_C          0.022778856          0.025766163       0 0.024311146
# 37      Atl_AD      Atl_AA          0.011042390          0.013790777       0 0.012525243
# 38      Atl_AD      Atl_AE          0.134324421          0.150102504       0 0.142523371
# 39      Atl_AD      Atl_AB          0.008213877          0.011964281       0 0.010047364
# 40      Atl_AD      Indo_B          0.009934122          0.012447656       0 0.011147368
# 41      Atl_AC      Indo_C          0.028183954          0.031579208       0 0.029577293
# 42      Atl_AC      Atl_AA          0.011054409          0.013419995       0 0.012350897
# 43      Atl_AC      Atl_AE          0.140720356          0.157534188       0 0.149693209
# 44      Atl_AC      Atl_AB          0.009978662          0.014920124       0 0.012355496
# 45      Atl_AC      Indo_B          0.014969722          0.017724272       0 0.016423868
# 46      Indo_C      Atl_AA          0.032978228          0.036702986       0 0.034876468
# 47      Indo_C      Atl_AE          0.162951094          0.178682820       0 0.171880123
# 48      Indo_C      Atl_AB          0.030411273          0.036109622       0 0.033259359
# 49      Indo_C      Indo_B          0.008890458          0.011453113       0 0.010024243
# 50      Atl_AA      Atl_AE          0.146595540          0.165169269       0 0.156176759
# 51      Atl_AA      Atl_AB          0.003593094          0.008403983       0 0.006053763
# 52      Atl_AA      Indo_B          0.019778201          0.023696853       0 0.021663719
# 53      Atl_AE      Atl_AB          0.159422681          0.178341479       0 0.170494461
# 54      Atl_AE      Indo_B          0.156383216          0.171449796       0 0.164388856
# 55      Atl_AB      Indo_B          0.014956186          0.019990766       0 0.017939332

zczc125.pop2.fst.ci<-gl.fst.pop(zczc125.pop2.gl, nboots=100, percent=95)
zczc125.pop2.fst.ci
zczc125.pop2.fst.ci$Fsts
#             Med_A       Med_C     Med_B     Indo_A      Atl_BC     Atl_BB     Indo_C     Atl_BA    Atl_BD Indo_B
# Med_A          NA          NA        NA         NA          NA         NA         NA         NA        NA     NA
# Med_C  0.07267638          NA        NA         NA          NA         NA         NA         NA        NA     NA
# Med_B  0.08884050 0.007056218        NA         NA          NA         NA         NA         NA        NA     NA
# Indo_A 0.24317274 0.166374501 0.2698512         NA          NA         NA         NA         NA        NA     NA
# Atl_BC 0.19403043 0.131194896 0.2085442 0.02886441          NA         NA         NA         NA        NA     NA
# Atl_BB 0.20016537 0.134161394 0.2152729 0.03195148 0.005441310         NA         NA         NA        NA     NA
# Indo_C 0.21190273 0.151797471 0.2273839 0.01775912 0.024311146 0.02957729         NA         NA        NA     NA
# Atl_BA 0.18559863 0.114863834 0.1964173 0.03017393 0.009738826 0.01032840 0.03150664         NA        NA     NA
# Atl_BD 0.36619480 0.327572194 0.4207762 0.21718246 0.142523371 0.14969321 0.17188012 0.14435312        NA     NA
# Indo_B 0.20869151 0.141220666 0.2257712 0.01659260 0.011147368 0.01642387 0.01002424 0.01811682 0.1643889     NA
zczc125.pop2.fst.ci$Pvalues
#         Med_A Med_C Med_B Indo_A Atl_BC Atl_BB Indo_C Atl_BA Atl_BD Indo_B
# Med_A     NA    NA    NA     NA     NA     NA     NA     NA     NA     NA
# Med_C      0    NA    NA     NA     NA     NA     NA     NA     NA     NA
# Med_B      0     0    NA     NA     NA     NA     NA     NA     NA     NA
# Indo_A     0     0     0     NA     NA     NA     NA     NA     NA     NA
# Atl_BC     0     0     0      0     NA     NA     NA     NA     NA     NA
# Atl_BB     0     0     0      0      0     NA     NA     NA     NA     NA
# Indo_C     0     0     0      0      0      0     NA     NA     NA     NA
# Atl_BA     0     0     0      0      0      0      0     NA     NA     NA
# Atl_BD     0     0     0      0      0      0      0      0     NA     NA
# Indo_B     0     0     0      0      0      0      0      0      0     NA
zczc125.fst.pop2.boots<-zczc125.pop2.fst.ci$Bootstraps[,c(1:2,103:106)]
zczc125.fst.pop2.boots
# Population1 Population2 Lower bound CI limit Upper bound CI limit p-value         Fst
# 1        Med_A       Med_C          0.068289171          0.077193946       0 0.072676384
# 2        Med_A       Med_B          0.084774906          0.093285087       0 0.088840500
# 3        Med_A      Indo_A          0.237310665          0.248680731       0 0.243172741
# 4        Med_A      Atl_BC          0.187773147          0.198600168       0 0.194030427
# 5        Med_A      Atl_BB          0.193813468          0.204962362       0 0.200165372
# 6        Med_A      Indo_C          0.205935549          0.216620673       0 0.211902730
# 7        Med_A      Atl_BA          0.179820920          0.190477216       0 0.185598626
# 8        Med_A      Atl_BD          0.356563746          0.374776309       0 0.366194796
# 9        Med_A      Indo_B          0.203533132          0.213472892       0 0.208691507
# 10       Med_C       Med_B          0.004041340          0.010071019       0 0.007056218
# 11       Med_C      Indo_A          0.158146873          0.171222070       0 0.166374501
# 12       Med_C      Atl_BC          0.125624537          0.135881402       0 0.131194896
# 13       Med_C      Atl_BB          0.128530301          0.139102126       0 0.134161394
# 14       Med_C      Indo_C          0.145826902          0.156375467       0 0.151797471
# 15       Med_C      Atl_BA          0.109421290          0.119754072       0 0.114863834
# 16       Med_C      Atl_BD          0.315894288          0.335725608       0 0.327572194
# 17       Med_C      Indo_B          0.136119963          0.146517004       0 0.141220666
# 18       Med_B      Indo_A          0.263039361          0.275726794       0 0.269851202
# 19       Med_B      Atl_BC          0.202209444          0.213040097       0 0.208544224
# 20       Med_B      Atl_BB          0.209052032          0.221736200       0 0.215272885
# 21       Med_B      Indo_C          0.221291406          0.232699658       0 0.227383922
# 22       Med_B      Atl_BA          0.190692502          0.201800435       0 0.196417339
# 23       Med_B      Atl_BD          0.410663623          0.430111264       0 0.420776161
# 24       Med_B      Indo_B          0.219895646          0.231584502       0 0.225771151
# 25      Indo_A      Atl_BC          0.026530090          0.031420477       0 0.028864413
# 26      Indo_A      Atl_BB          0.029048636          0.034683843       0 0.031951483
# 27      Indo_A      Indo_C          0.015816401          0.019700008       0 0.017759121
# 28      Indo_A      Atl_BA          0.027866463          0.033009294       0 0.030173927
# 29      Indo_A      Atl_BD          0.208651032          0.224499544       0 0.217182458
# 30      Indo_A      Indo_B          0.014201735          0.019205849       0 0.016592605
# 31      Atl_BC      Atl_BB          0.004474818          0.006159607       0 0.005441310
# 32      Atl_BC      Indo_C          0.022644580          0.025841099       0 0.024311146
# 33      Atl_BC      Atl_BA          0.008278836          0.010834848       0 0.009738826
# 34      Atl_BC      Atl_BD          0.133949640          0.148183276       0 0.142523371
# 35      Atl_BC      Indo_B          0.009428582          0.012522177       0 0.011147368
# 36      Atl_BB      Indo_C          0.027746675          0.030886059       0 0.029577293
# 37      Atl_BB      Atl_BA          0.009028759          0.011329402       0 0.010328401
# 38      Atl_BB      Atl_BD          0.141396549          0.155643718       0 0.149693209
# 39      Atl_BB      Indo_B          0.015150212          0.017850369       0 0.016423868
# 40      Indo_C      Atl_BA          0.029538354          0.033129029       0 0.031506636
# 41      Indo_C      Atl_BD          0.164493464          0.177222322       0 0.171880123
# 42      Indo_C      Indo_B          0.008554546          0.011194612       0 0.010024243
# 43      Atl_BA      Atl_BD          0.136042238          0.149864756       0 0.144353124
# 44      Atl_BA      Indo_B          0.016596512          0.019531219       0 0.018116823
# 45      Atl_BD      Indo_B          0.155980755          0.170724327       0 0.164388856


##These give same results as Hierfstat basic.stats
zczc125.gl.basicstats<-gl.basic.stats(zczc125.gl)
zczc125.pop1.basicstats<-gl.basic.stats(zczc125.pop1.gl)
zczc125.pop2.basicstats<-gl.basic.stats(zczc125.pop2.gl)

summary(zczc125.gl.basicstats$Ho)
summary(zczc125.gl.basicstats$Hs)
summary(zczc125.gl.basicstats$Fis)

summary(zczc125.pop1.basicstats$Ho)
summary(zczc125.pop1.basicstats$Hs)
summary(zczc125.pop1.basicstats$Fis)

summary(zczc125.pop2.basicstats$Ho)
summary(zczc125.pop2.basicstats$Hs)
summary(zczc125.pop2.basicstats$Fis)


# STAMPP ------------------------------------------------------------------

library(StAMPP)

#calculate fst with confidence intervals
zczc125.stampp<-stamppFst(zczc125.gl, nboots=100, percent=95)
zczc125.stampp$Fsts
zczc125.stampp$Pvalues
zczc125.stampp$Bootstraps

#calcualte Neis Genetic Distannce for AMOVA
zczc125.stampp.nei.pop<-stamppNeisD(zczc125.gl, pop=TRUE)
zczc125.stampp.nei.ind<-stamppNeisD(zczc125.gl, pop=FALSE)

#amova
zczc125.stampp.amova<-stamppAmova(zczc125.stampp.nei.ind, zczc125.gl, perm=100)
zczc125.stampp.amova$tab
zczc125.stampp.amova$varcoef
zczc125.stampp.amova$varcomp
zczc125.stampp.amova$call


# DiveRsity ---------------------------------------------------------------
devtools::install_github('wrengels/HWxtest', subdir='pkg')
library(HWxtest)
install.packages("diveRsity", dependencies = TRUE)
devtools::install_github("kkeenan02/diveRsity")
library(diveRsity)
#dependencies and optional support packagees

install.packages("plotrix", dependencies = TRUE)
install.packages("shiny", dependencies = TRUE)
install.packages("xlsx", dependencies = TRUE)
install.packages("sendplot", dependencies = TRUE)
install.packages("doParallel", dependencies = TRUE)
install.packages("parallel", dependencies = TRUE)
install.packages("foreach", dependencies = TRUE)
install.packages("iterators", dependencies = TRUE)
devtools::install_github("rstudio/shiny-incubator")
install.packages("doSNOW", dependencies=TRUE)


library(plotrix)
library(shiny)
library(xlsx)
library(sendplot)
library(doParallel)
library(parallel)
library(foreach)
library(iterators)
library(shinyIncubator)
library(doSNOW)


runApp(getwd())

zczc125.fstWC<-fstWC(infile = "./zczc125.genepop.txt")

zczc125.basicstats<-basicStats(infile = "./zczc125_genepop.txt", outfile="./zczc125.basicstats.txt", fis_ci = TRUE, fis_boots = 100)

zczc125.divBasic<-divBasic(infile = "./zczc125_genepop.txt")

#Confidence Intervals for Fst values based on Ocean basins
zczc125.diffCalc<-diffCalc(infile="./zczc125.genepop.txt", outfile="./zczc125.fststats.noloc.ind.txt", fst = TRUE, pairwise=TRUE, para=TRUE, bs_locus = TRUE, boots=100, bs_pairwise=TRUE)
tail(zczc125.diffCalc$std_stats)
names(zczc125.diffCalc)

#Confidence Intervals for Fst values based on populations 1
diffCalc(infile="./zczc125.pop1.genepop.txt", outfile="./zczc125.fststats.pop1.txt", fst = TRUE, pairwise=TRUE, bs_locus=TRUE, para=TRUE, boots=100, bs_pairwise=TRUE)

#Confidence Intervals for Fst values based on populations 2
diffCalc(infile="./zczc125.pop2.genepop.txt", outfile="./zczc125.fststats.pop2.txt", fst = TRUE, pairwise=TRUE, bs_locus=TRUE, para=TRUE, boots=100, bs_pairwise=TRUE)


###This is Emma's code that worked

dif.diversity<-diffCalc(infile = "/Users/ecar026/Dropbox/Bioinformatics/Cluster_UOA/Output_high_quality/withrep/SRWNoRep_populations.snps.genepop", outfile = "fst.calc.pairwise", fst = TRUE, pairwise = TRUE, 
                        bs_locus = TRUE, para = TRUE, boots = 999, bs_pairwise = TRUE)


# Hierfstat/Pegas ---------------------------------------------------------------
install_github("jgx65/hierfstat")
library("hierfstat")

#compute Fit, Fst and Fis for each locus, using formulae in Weir and Cockerham (1984) for each allele and then averaged within each locus over the different alleles as suggested by these authors. Will generate a table, one column with locus, then three with each of the F-statistics 
zczc125.fsttab<-Fst(as.loci(zczc125.gi), pop=zczc125.gi$pop)
zczc125.fsttab[1:3,1:3]
head(zczc125.fsttab)

#convert genind file to hierfstat formatted file, first column is population name then the remaining columns are loci
zczc125.hfstat<-genind2hierfstat(zczc125.gi)
dim(zczc125.hfstat)
zczc125.hfstat[1:3, 1:3]
zczc125.pop1<-zczc125.pop1.gi$pop
zczc125.pop2<-zczc125.pop2.gi$pop
zczc125.hfstat1<-cbind(zczc125.pop1,zczc125.hfstat)
zczc125.hfstat2<-cbind(zczc125.pop2, zczc125.hfstat)

zczc125.atl1.hfstat<-subset(zczc125.hfstat1, zczc125.hfstat1$pop=="Atlantic")
zczc125.atl2.hfstat<-subset(zczc125.hfstat2, zczc125.hfstat2$pop=="Atlantic")
zczc125.atl1.hfstat[1:5,1:5]
zczc125.atl2.hfstat[1:5,1:5]
zczc125.atl1.hfstat<-zczc125.atl1.hfstat[,-2]
zczc125.atl2.hfstat<-zczc125.atl2.hfstat[,-2]


zczc125.atl1.hfstat<-zczc125.atl1.hfstat[order(zczc125.atl1.hfstat$zczc125.pop1),]
zczc125.atl1.hfstat[1:3,1:3]
colnames(zczc125.atl1.hfstat)[1]<-"pop"
dim(zczc125.atl1.hfstat)


zczc125.indo1.hfstat<-subset(zczc125.hfstat1, zczc125.hfstat1$pop=="Indopacific")
zczc125.indo2.hfstat<-subset(zczc125.hfstat1, zczc125.hfstat2$pop=="Indopacific")
zczc125.indo1.hfstat<-zczc125.indo1.hfstat[,-2]
zczc125.indo2.hfstat<-zczc125.indo2.hfstat[,-2]


zczc125.indo1.hfstat<-subset(zczc125.hfstat1, zczc125.hfstat1$pop=="Indopacific")
zczc125.indo2.hfstat<-subset(zczc125.hfstat1, zczc125.hfstat2$pop=="Indopacific")
zczc125.indo1.hfstat<-zczc125.indo1.hfstat[,-2]
zczc125.indo2.hfstat<-zczc125.indo2.hfstat[,-2]


#calculate weir and cockerham f statistics
zczc125.wc<-wc(zczc125.hfstat)
zczc125.wc$FIS

zczc125.atl1.levels<-zczc125.atl1.hfstat$pop
zczc125.atl1.loci<-zczc125.atl1.hfstat[,-1]
dim(zczc125.atl1.loci)

varcomp.glob(zczc125.atl1.levels, zczc125.atl1.loci)

zczc125.indo.wc<-wc(zczc125.indo.hfstat)


##For some reason, just this code isnt working! Atlantic and Indo-pacific worked just fine!
zczc125.med.wc<-wc(zczc125.med.hfstat[,-1])
nrows(zczc125.med.wc)


sum(is.na(zczc125.atl.hfstat)) #33462
sum(is.na(zczc125.indo.hfstat)) #21895
sum(is.na(zczc125.med.hfstat)) #39389


#Calculate individual counts, allelic frequencies, observed heterozygosities and genetic diversities per locus and per population. and mean observed heterozygosities, mean gene diversities and within population Hs. Based on Nei 1987
zczc125.hfstat.basicstats<-basic.stats(zczc125.hfstat, diploid=TRUE)
zczc125.hfstat.basicstats$overall
print(zczc125.hfstat.basicstats)

#If not including populations, need to remove that column with [,-1]
zczc125.atl.basicstats<-basic.stats(zczc125.atl.hfstat[,-1], diploid=TRUE)
zczc125.indo.basicstats<-basic.stats(zczc125.indo.hfstat[,-1], diploid=TRUE)
zczc125.med.basicstats<-basic.stats(zczc125.med.hfstat[,-1], diploid=TRUE)

zczc125.atl.basicstats$overall
# Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
# 0.1227  0.1137  0.1273  0.0135  0.1407  0.0273  0.1063  0.1943 -0.0786  0.0309 
zczc125.indo.basicstats$overall
# Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
# 0.1163  0.1339  0.1338 -0.0001  0.1337 -0.0002 -0.0007 -0.0015  0.1318 -0.0002 
zczc125.med.basicstats$overall
# Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
# 0.0850 0.0994 0.0994 0.0000    NaN    NaN 0.0000    NaN 0.1451    NaN 

zczc125.hfstat.ho<-zczc125.hfstat.basicstats$Ho
zczc125.hfstat.hs<-zczc125.hfstat.basicstats$Hs
zczc125.hfstat.Fis<-zczc125.hfstat.basicstats$Fis

zczc125.hfstat.ho<-as.data.frame(zczc125.hfstat.ho)
zczc125.hfstat.hs<-as.data.frame(zczc125.hfstat.hs)
zczc125.hfstat.Fis<-as.data.frame(zczc125.hfstat.Fis)

boxplot(zczc125.hfstat.ho)
boxplot(zczc125.hfstat.hs)
boxplot(zczc125.hfstat.Fis)

summary(zczc125.hfstat.ho)
summary(zczc125.hfstat.hs)
summary(zczc125.hfstat.Fis)

zczc125.hfstat.atl<-subset(zczc125.hfstat, zczc125.hfstat$pop=="Atlantic")
dim(zczc125.hfstat.atl)
zczc125.hfstat.atl[1:5,1:5]

zczc125.hfstat.indo<-subset(zczc125.hfstat, zczc125.hfstat$pop=="Indopacific")

zczc125.hfstat.med<-subset(zczc125.hfstat, zczc125.hfstat$pop=="Mediterranean")

zczc125.hfstat.basicstats.atl<-basic.stats(zczc125.hfstat.atl, diploid=TRUE)
zczc125.hfstat.basicstats.indo<-basic.stats(zczc125.hfstat.indo, diploid=TRUE)
zczc125.hfstat.basicstats.med<-basic.stats(zczc125.hfstat.med, diploid=TRUE)


zczc125.gt[1:5,1:5]

wc(zczc125.atl.hfstat)

# bootstrapping hierfstat results -----------------------------------------


install.packages("boot")
library(boot)
zczc125.ho.atlantic.boot<-boot(zczc125.hfstat.ho$Atlantic, function(x,i) mean(x[i]), R=10000)
# ORDINARY NONPARAMETRIC BOOTSTRAP
# Call:
#   boot(data = zczc125.hfstat.ho$Atlantic, statistic = function(x,i) mean(x[i]), R = 100)
# Bootstrap Statistics :
#   original       bias     std. error
# t1* 0.118471 0.0001347708 0.0007286473
zczc125.ho.atlantic.boot.mean<-mean(zczc125.ho.atlantic.boot$t[,1])
#0.1184771
hist(zczc125.ho.atlantic.boot$t[,1])
zczc125.ho.atlantic.boot.ci<-boot.ci(zczc125.ho.atlantic.boot, conf=0.95, type = "norm")
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 10000 bootstrap replicates
# 
# CALL : 
#   boot.ci(boot.out = zczc125.ho.boot, conf = 0.95, type = "norm")
# 
# Intervals : 
#   Level      Normal        
# 95%   ( 0.1168,  0.1201 )  
# Calculations and Intervals on Original Scale
                                                                        
zczc125.ho.indopacific.boot<-boot(zczc125.hfstat.ho$Indopacific, function(x,i) mean(x[i]), R=10000)
# ORDINARY NONPARAMETRIC BOOTSTRAP
# 
# 
# Call:
#   boot(data = zczc125.hfstat.ho$Indopacific, statistic = function(x, 
#                                                                   i) mean(x[i]), R = 10000)
# 
# 
# Bootstrap Statistics :
#   original       bias     std. error
# t1* 0.1159951 4.792618e-06 0.0008672516
zczc125.ho.indopacific.boot.mean<-mean(zczc125.ho.indopacific.boot$t[,1])
#0.1159999
hist(zczc125.ho.indopacific.boot$t[,1])
zczc125.ho.indopacific.boot.ci<-boot.ci(zczc125.ho.indopacific.boot, conf=0.95, type = "norm")
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 10000 bootstrap replicates
# 
# CALL : 
#   boot.ci(boot.out = zczc125.ho.indopacific.boot, conf = 0.95, 
#           type = "norm")
# 
# Intervals : 
#   Level      Normal        
# 95%   ( 0.1143,  0.1177 )  
# Calculations and Intervals on Original Scale

zczc125.ho.mediterranean.boot<-boot(zczc125.hfstat.ho$Mediterranean, function(x,i) mean(x[i]), R=10000)
# ORDINARY NONPARAMETRIC BOOTSTRAP
# 
# 
# Call:
#   boot(data = zczc125.hfstat.ho$Mediterranean, statistic = function(x, 
#                                                                     i) mean(x[i]), R = 10000)
# 
# 
# Bootstrap Statistics :
#   original       bias     std. error
# t1* 0.08499457 7.572997e-06 0.0009300494
zczc125.ho.mediterranean.boot.mean<-mean(zczc125.ho.mediterranean.boot$t[,1])
#0.08500214
hist(zczc125.ho.mediterranean.boot$t[,1])
zczc125.ho.mediterranean.boot.ci<-boot.ci(zczc125.ho.mediterranean.boot, conf=0.95, type = "norm")
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 10000 bootstrap replicates
# 
# CALL : 
#   boot.ci(boot.out = zczc125.ho.mediterranean.boot, conf = 0.95, 
#           type = "norm")
# 
# Intervals : 
#   Level      Normal        
# 95%   ( 0.0832,  0.0868 )  
# Calculations and Intervals on Original Scale

zczc125.hs.atlantic.boot<-boot(zczc125.hfstat.hs$Atlantic, function(x,i) mean(x[i]), R=10000)
# ORDINARY NONPARAMETRIC BOOTSTRAP
# 
# 
# Call:
#   boot(data = zczc125.hfstat.hs$Atlantic, statistic = function(x, 
#                                                                i) mean(x[i]), R = 10000)
# 
# 
# Bootstrap Statistics :
#   original       bias     std. error
# t1* 0.1371626 -3.36902e-06 0.0009066062
zczc125.hs.atlantic.boot.mean<-mean(zczc125.hs.atlantic.boot$t[,1])
#0.1371592
hist(zczc125.hs.atlantic.boot$t[,1])
zczc125.hs.atlantic.boot.ci<-boot.ci(zczc125.hs.atlantic.boot, conf=0.95, type = "norm")
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 10000 bootstrap replicates
# 
# CALL : 
#   boot.ci(boot.out = zczc125.hs.atlantic.boot, conf = 0.95, type = "norm")
# 
# Intervals : 
#   Level      Normal        
# 95%   ( 0.1354,  0.1389 )  
# Calculations and Intervals on Original Scale

zczc125.hs.indopacific.boot<-boot(zczc125.hfstat.hs$Indopacific, function(x,i) mean(x[i]), R=10000)
# ORDINARY NONPARAMETRIC BOOTSTRAP
# 
# 
# Call:
#   boot(data = zczc125.hfstat.hs$Indopacific, statistic = function(x, 
#                                                                   i) mean(x[i]), R = 10000)
# 
# 
# Bootstrap Statistics :
#   original       bias     std. error
# t1* 0.1328643 2.730111e-06 0.0009477275
zczc125.hs.indopacific.boot.mean<-mean(zczc125.hs.indopacific.boot$t[,1])
#0.132867
hist(zczc125.hs.indopacific.boot$t[,1])
zczc125.hs.indopacific.boot.ci<-boot.ci(zczc125.hs.indopacific.boot, conf=0.95, type = "norm")
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 10000 bootstrap replicates
# 
# CALL : 
#   boot.ci(boot.out = zczc125.hs.indopacific.boot, conf = 0.95, 
#           type = "norm")
# 
# Intervals : 
#   Level      Normal        
# 95%   ( 0.1310,  0.1347 )  
# Calculations and Intervals on Original Scale

zczc125.hs.mediterranean.boot<-boot(zczc125.hfstat.hs$Mediterranean, function(x,i) mean(x[i]), R=10000)
# RDINARY NONPARAMETRIC BOOTSTRAP
# 
# 
# Call:
#   boot(data = zczc125.hfstat.hs$Mediterranean, statistic = function(x, 
#                                                                     i) mean(x[i]), R = 10000)
# 
# 
# Bootstrap Statistics :
#   original        bias    std. error
# t1* 0.09941738 -6.238131e-06 0.001023246
zczc125.hs.mediterranean.boot.mean<-mean(zczc125.hs.mediterranean.boot$t[,1])
#0.09941114
hist(zczc125.hs.mediterranean.boot$t[,1])
zczc125.hs.mediterranean.boot.ci<-boot.ci(zczc125.hs.mediterranean.boot, conf=0.95, type = "norm")
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 10000 bootstrap replicates
# 
# CALL : 
#   boot.ci(boot.out = zczc125.hs.mediterranean.boot, conf = 0.95, 
#           type = "norm")
# 
# Intervals : 
#   Level      Normal        
# 95%   ( 0.0974,  0.1014 )  
# Calculations and Intervals on Original Scale

#Need to remove NAs from Fis data before bootstrapping
length(zczc125.hfstat.Fis$Atlantic)
zczc125.hfstat.Fis.Atlantic<-na.omit(zczc125.hfstat.Fis$Atlantic)
length(zczc125.hfstat.Fis.Atlantic)

zczc125.Fis.atlantic.boot<-boot(zczc125.hfstat.Fis.Atlantic, function(x,i) mean(x[i]), R=10000)
# ORDINARY NONPARAMETRIC BOOTSTRAP
# 
# 
# Call:
#   boot(data = zczc125.hfstat.Fis.Atlantic, statistic = function(x, 
#                                                                 i) mean(x[i]), R = 10000)
# 
# 
# Bootstrap Statistics :
#   original       bias    std. error
# t1* 0.126398 2.648219e-05 0.001751465
zczc125.Fis.atlantic.boot.mean<-mean(zczc125.Fis.atlantic.boot$t[,1])
#0.1264245
hist(zczc125.Fis.atlantic.boot$t[,1])
zczc125.Fis.atlantic.boot.ci<-boot.ci(zczc125.Fis.atlantic.boot, conf=0.95, type = "norm")
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 10000 bootstrap replicates
# 
# CALL : 
#   boot.ci(boot.out = zczc125.Fis.atlantic.boot, conf = 0.95, type = "norm")
# 
# Intervals : 
#   Level      Normal        
# 95%   ( 0.1229,  0.1298 )  
# Calculations and Intervals on Original Scale

zczc125.hfstat.Fis.Indopacific<-na.omit(zczc125.hfstat.Fis$Indopacific)
zczc125.Fis.indopacific.boot<-boot(zczc125.hfstat.Fis.Indopacific, function(x,i) mean(x[i]), R=10000)
# ORDINARY NONPARAMETRIC BOOTSTRAP
# 
# 
# Call:
#   boot(data = zczc125.hfstat.Fis.Indopacific, statistic = function(x, 
#                                                                    i) mean(x[i]), R = 10000)
# 
# 
# Bootstrap Statistics :
#   original       bias    std. error
# t1* 0.1128735 1.821679e-05 0.001921171
zczc125.Fis.indopacific.boot.mean<-mean(zczc125.Fis.indopacific.boot$t[,1])
#0.1128917
hist(zczc125.Fis.indopacific.boot$t[,1])
zczc125.Fis.indopacific.boot.ci<-boot.ci(zczc125.Fis.indopacific.boot, conf=0.95, type = "norm")
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 10000 bootstrap replicates
# 
# CALL : 
#   boot.ci(boot.out = zczc125.Fis.indopacific.boot, conf = 0.95, 
#           type = "norm")
# 
# Intervals : 
#   Level      Normal        
# 95%   ( 0.1091,  0.1166 )  
# Calculations and Intervals on Original Scale


zczc125.hfstat.Fis.Mediterranean<-na.omit(zczc125.hfstat.Fis$Mediterranean)
zczc125.Fis.mediterranean.boot<-boot(zczc125.hfstat.Fis.Mediterranean, function(x,i) mean(x[i]), R=10000)
# RDINARY NONPARAMETRIC BOOTSTRAP
# 
# 
# Call:
#   boot(data = zczc125.hfstat.hs$Mediterranean, statistic = function(x, 
#                                                                     i) mean(x[i]), R = 10000)
# 
# 
# Bootstrap Statistics :
#   original        bias    std. error
# t1* 0.09941738 -6.238131e-06 0.001023246
zczc125.hs.mediterranean.boot.mean<-mean(zczc125.hs.mediterranean.boot$t[,1])
#0.09941114
hist(zczc125.hs.mediterranean.boot$t[,1])
zczc125.hs.mediterranean.boot.ci<-boot.ci(zczc125.hs.mediterranean.boot, conf=0.95, type = "norm")
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 10000 bootstrap replicates
# 
# CALL : 
#   boot.ci(boot.out = zczc125.hs.mediterranean.boot, conf = 0.95, 
#           type = "norm")
# 
# Intervals : 
#   Level      Normal        
# 95%   ( 0.0974,  0.1014 )  
# Calculations and Intervals on Original Scale

#bootstrapping over loci of populations Fis. 
zczc125.hfstat.bootfis<-boot.ppfis(dat=zczc125.hfstat, nboot=100, quant=c(0.025,0.975),diploid=TRUE)
zczc125.hfstat.bootfis$fis.ci



#PCA
zczc125.hfstat.pca<-indpca(zczc125.hfstat, ind.labels=indNames(zczc125.gl))
zczc125.hfstat.pca

plot(zczc125.hfstat.pca, col=zczc125.gl$pop)

#Bootstrap confidence intervals (over loci) for gene diversity
zczc125.hfstat.booths<-boot.vc(levels=zczc125.hfstat[,1], loci=zczc125.hfstat[,-1], nboot=100)
zczc125.hfstat.booths$ci


dim(zczc125.hfstat)
zczc125.hfstat[,1]


# making opt gl files -----------------------------------------------------


zczc125.pop.opt.a.gl<-zczc125.gl
pop(zczc125.pop.opt.a.gl)<-as.factor(c("Med_A",
                                      "Med_A",
                                      "Med_A",
                                      "Med_A",
                                      "Med_A",
                                      "Med_A",
                                      "Med_A",
                                      "Med_C",
                                      "Med_A",
                                      "Med_C",
                                      "Med_A",
                                      "Med_C",
                                      "Med_A",
                                      "Med_A",
                                      "Med_B",
                                      "Indo_A",
                                      "Indo_A",
                                      "Atl_D",
                                      "Atl_C",
                                      "Med_B",
                                      "Atl_D",
                                      "Med_B",
                                      "Atl_D",
                                      "Indo_D",
                                      "Med_B",
                                      "Med_A",
                                      "Indo_D",
                                      "Med_B",
                                      "Med_C",
                                      "Med_C",
                                      "Med_NA",
                                      "Med_B",
                                      "Med_B",
                                      "Atl_C",
                                      "Atl_C",
                                      "Atl_C",
                                      "Atl_C",
                                      "Atl_C",
                                      "Atl_D",
                                      "Atl_C",
                                      "Atl_C",
                                      "Atl_D",
                                      "Atl_D",
                                      "Atl_D",
                                      "Atl_D",
                                      "Atl_D",
                                      "Atl_C",
                                      "Atl_D",
                                      "Atl_A",
                                      "Atl_D",
                                      "Atl_C",
                                      "Atl_E",
                                      "Atl_D",
                                      "Atl_A",
                                      "Atl_E",
                                      "Atl_C",
                                      "Atl_C",
                                      "Atl_D",
                                      "Atl_D",
                                      "Atl_D",
                                      "Atl_C",
                                      "Atl_B",
                                      "Atl_A",
                                      "Atl_C",
                                      "Indo_D",
                                      "Indo_D",
                                      "Atl_A",
                                      "Med_B",
                                      "Indo_C",
                                      "Atl_A",
                                      "Indo_D",
                                      "Atl_A",
                                      "Atl_D",
                                      "Atl_C",
                                      "Med_A",
                                      "Med_A",
                                      "Indo_B",
                                      "Med_A",
                                      "Med_B",
                                      "Med_A",
                                      "Atl_C",
                                      "Atl_B",
                                      "Atl_B",
                                      "Atl_A",
                                      "Atl_D",
                                      "Atl_D",
                                      "Atl_D",
                                      "Atl_C",
                                      "Atl_A",
                                      "Atl_A",
                                      "Atl_A",
                                      "Atl_A",
                                      "Atl_D",
                                      "Indo_D",
                                      "Indo_D",
                                      "Med_A",
                                      "Med_A",
                                      "Indo_D",
                                      "Indo_D",
                                      "Atl_B",
                                      "Atl_B",
                                      "Indo_D",
                                      "Indo_D",
                                      "Indo_A",
                                      "Indo_D",
                                      "Indo_D",
                                      "Indo_D",
                                      "Indo_D",
                                      "Indo_D",
                                      "Indo_A",
                                      "Indo_D",
                                      "Indo_D",
                                      "Indo_A",
                                      "Indo_B",
                                      "Indo_B",
                                      "Indo_C",
                                      "Indo_C",
                                      "Indo_B",
                                      "Indo_B",
                                      "Indo_B",
                                      "Med_A",
                                      "Indo_D",
                                      "Indo_B",
                                      "Indo_B",
                                      "Indo_B"))
popNames(zczc125.pop.opt.a.gl)

zczc125.pop.opt.b.gl<-zczc125.gl
pop(zczc125.pop.opt.b.gl)<-as.factor(c("Med_A",
                                       "Med_A",
                                       "Med_A",
                                       "Med_A",
                                       "Med_A",
                                       "Med_A",
                                       "Med_A",
                                       "Med_C",
                                       "Med_A",
                                       "Med_C",
                                       "Med_A",
                                       "Med_C",
                                       "Med_A",
                                       "Med_A",
                                       "Med_B",
                                       "Indo_A",
                                       "Indo_A",
                                       "Atl_E",
                                       "Atl_C",
                                       "Med_D",
                                       "Atl_E",
                                       "Med_B",
                                       "Atl_D",
                                       "Indo_C",
                                       "Med_B",
                                       "Med_A",
                                       "Indo_C",
                                       "Med_D",
                                       "Med_C",
                                       "Med_C",
                                       "Med_NA",
                                       "Med_D",
                                       "Med_D",
                                       "Atl_C",
                                       "Atl_C",
                                       "Atl_C",
                                       "Atl_C",
                                       "Atl_C",
                                       "Atl_E",
                                       "Atl_C",
                                       "Atl_C",
                                       "Atl_E",
                                       "Atl_E",
                                       "Atl_E",
                                       "Atl_E",
                                       "Atl_D",
                                       "Atl_C",
                                       "Atl_E",
                                       "Atl_A",
                                       "Atl_E",
                                       "Atl_C",
                                       "Atl_F",
                                       "Atl_E",
                                       "Atl_A",
                                       "Atl_F",
                                       "Atl_C",
                                       "Atl_C",
                                       "Atl_E",
                                       "Atl_E",
                                       "Atl_E",
                                       "Atl_C",
                                       "Atl_B",
                                       "Atl_A",
                                       "Atl_C",
                                       "Indo_C",
                                       "Indo_C",
                                       "Atl_A",
                                       "Med_D",
                                       "Indo_B",
                                       "Atl_A",
                                       "Indo_C",
                                       "Atl_A",
                                       "Atl_E",
                                       "Atl_C",
                                       "Med_A",
                                       "Med_A",
                                       "Indo_B",
                                       "Med_A",
                                       "Med_B",
                                       "Med_A",
                                       "Atl_C",
                                       "Atl_B",
                                       "Atl_B",
                                       "Atl_A",
                                       "Atl_E",
                                       "Atl_E",
                                       "Atl_E",
                                       "Atl_C",
                                       "Atl_A",
                                       "Atl_A",
                                       "Atl_A",
                                       "Atl_A",
                                       "Atl_E",
                                       "Indo_C",
                                       "Indo_C",
                                       "Med_A",
                                       "Med_A",
                                       "Indo_C",
                                       "Indo_C",
                                       "Atl_B",
                                       "Atl_B",
                                       "Indo_C",
                                       "Indo_C",
                                       "Indo_A",
                                       "Indo_C",
                                       "Indo_C",
                                       "Indo_C",
                                       "Indo_C",
                                       "Indo_C",
                                       "Indo_A",
                                       "Indo_C",
                                       "Indo_C",
                                       "Indo_A",
                                       "Indo_B",
                                       "Indo_B",
                                       "Indo_B",
                                       "Indo_B",
                                       "Indo_B",
                                       "Indo_B",
                                       "Indo_B",
                                       "Med_A",
                                       "Indo_C",
                                       "Indo_B",
                                       "Indo_B",
                                       "Indo_B"))
popNames(zczc125.pop.opt.b.gl)

zczc125.pop.opt.c.gl<-zczc125.gl
pop(zczc125.pop.opt.c.gl)<-as.factor(c("Med_A",
                                       "Med_A",
                                       "Med_A",
                                       "Med_A",
                                       "Med_A",
                                       "Med_A",
                                       "Med_A",
                                       "Med_C",
                                       "Med_A",
                                       "Med_C",
                                       "Med_A",
                                       "Med_C",
                                       "Med_A",
                                       "Med_A",
                                       "Med_B",
                                       "Indo_A",
                                       "Indo_A",
                                       "Atl_D",
                                       "Atl_B",
                                       "Med_B",
                                       "Atl_D",
                                       "Med_B",
                                       "Atl_C",
                                       "Indo_D",
                                       "Med_B",
                                       "Med_A",
                                       "Indo_D",
                                       "Med_B",
                                       "Med_C",
                                       "Med_C",
                                       "Med_B",
                                       "Med_B",
                                       "Med_B",
                                       "Atl_B",
                                       "Atl_B",
                                       "Atl_B",
                                       "Atl_B",
                                       "Atl_B",
                                       "Atl_D",
                                       "Atl_B",
                                       "Atl_B",
                                       "Atl_D",
                                       "Atl_D",
                                       "Atl_D",
                                       "Atl_D",
                                       "Atl_C",
                                       "Atl_B",
                                       "Atl_D",
                                       "Atl_A",
                                       "Atl_D",
                                       "Atl_B",
                                       "Atl_E",
                                       "Atl_D",
                                       "Atl_A",
                                       "Atl_E",
                                       "Atl_B",
                                       "Atl_B",
                                       "Atl_D",
                                       "Atl_D",
                                       "Atl_D",
                                       "Atl_B",
                                       "Atl_A",
                                       "Atl_A",
                                       "Atl_B",
                                       "Indo_D",
                                       "Indo_D",
                                       "Atl_A",
                                       "Med_B",
                                       "Indo_C",
                                       "Atl_A",
                                       "Indo_D",
                                       "Atl_A",
                                       "Atl_D",
                                       "Atl_B",
                                       "Med_A",
                                       "Med_A",
                                       "Indo_B",
                                       "Med_A",
                                       "Med_B",
                                       "Med_A",
                                       "Atl_B",
                                       "Atl_A",
                                       "Atl_A",
                                       "Atl_A",
                                       "Atl_D",
                                       "Atl_D",
                                       "Atl_D",
                                       "Atl_B",
                                       "Atl_A",
                                       "Atl_A",
                                       "Atl_A",
                                       "Atl_A",
                                       "Atl_D",
                                       "Indo_D",
                                       "Indo_D",
                                       "Med_A",
                                       "Med_A",
                                       "Indo_D",
                                       "Indo_D",
                                       "Atl_A",
                                       "Atl_A",
                                       "Indo_D",
                                       "Indo_D",
                                       "Indo_A",
                                       "Indo_D",
                                       "Indo_D",
                                       "Indo_D",
                                       "Indo_D",
                                       "Indo_D",
                                       "Indo_A",
                                       "Indo_D",
                                       "Indo_D",
                                       "Indo_A",
                                       "Indo_B",
                                       "Indo_B",
                                       "Indo_C",
                                       "Indo_B",
                                       "Indo_B",
                                       "Indo_B",
                                       "Indo_B",
                                       "Med_A",
                                       "Indo_D",
                                       "Indo_B",
                                       "Indo_B",
                                       "Indo_B"))
popNames(zczc125.pop.opt.c.gl)

zczc125.pop1.gl<-zczc125.gl
pop(zczc125.pop1.gl)<-as.factor(c("Med_A",
                                 "Med_A",
                                 "Med_A",
                                 "Med_A",
                                 "Med_A",
                                 "Med_A",
                                 "Med_A",
                                 "Med_C",
                                 "Med_A",
                                 "Med_C",
                                 "Med_A",
                                 "Med_C",
                                 "Med_A",
                                 "Med_A",
                                 "Med_B",
                                 "Indo_A",
                                 "Indo_A",
                                 "Atl_AD",
                                 "Atl_AC",
                                 "Med_B",
                                 "Atl_AD",
                                 "Med_B",
                                 "Atl_AD",
                                 "Indo_C",
                                 "Med_B",
                                 "Med_A",
                                 "Indo_C",
                                 "Med_B",
                                 "Med_C",
                                 "Med_C",
                                 "Med_B",
                                 "Med_B",
                                 "Med_B",
                                 "Atl_AC",
                                 "Atl_AC",
                                 "Atl_AC",
                                 "Atl_AC",
                                 "Atl_AC",
                                 "Atl_AD",
                                 "Atl_AC",
                                 "Atl_AC",
                                 "Atl_AD",
                                 "Atl_AD",
                                 "Atl_AD",
                                 "Atl_AD",
                                 "Atl_AD",
                                 "Atl_AC",
                                 "Atl_AD",
                                 "Atl_AA",
                                 "Atl_AD",
                                 "Atl_AC",
                                 "Atl_AE",
                                 "Atl_AD",
                                 "Atl_AA",
                                 "Atl_AE",
                                 "Atl_AC",
                                 "Atl_AC",
                                 "Atl_AD",
                                 "Atl_AD",
                                 "Atl_AD",
                                 "Atl_AC",
                                 "Atl_AB",
                                 "Atl_AA",
                                 "Atl_AC",
                                 "Indo_C",
                                 "Indo_C",
                                 "Atl_AA",
                                 "Med_B",
                                 "Indo_B",
                                 "Atl_AA",
                                 "Indo_C",
                                 "Atl_AA",
                                 "Atl_AD",
                                 "Atl_AC",
                                 "Med_A",
                                 "Med_A",
                                 "Indo_B",
                                 "Med_A",
                                 "Med_B",
                                 "Med_A",
                                 "Atl_AC",
                                 "Atl_AB",
                                 "Atl_AB",
                                 "Atl_AA",
                                 "Atl_AD",
                                 "Atl_AD",
                                 "Atl_AD",
                                 "Atl_AC",
                                 "Atl_AA",
                                 "Atl_AA",
                                 "Atl_AA",
                                 "Atl_AA",
                                 "Atl_AD",
                                 "Indo_C",
                                 "Indo_C",
                                 "Med_A",
                                 "Med_A",
                                 "Indo_C",
                                 "Indo_C",
                                 "Atl_AB",
                                 "Atl_AB",
                                 "Indo_C",
                                 "Indo_C",
                                 "Indo_A",
                                 "Indo_C",
                                 "Indo_C",
                                 "Indo_C",
                                 "Indo_C",
                                 "Indo_C",
                                 "Indo_A",
                                 "Indo_C",
                                 "Indo_C",
                                 "Indo_A",
                                 "Indo_B",
                                 "Indo_B",
                                 "Indo_B",
                                 "Indo_B",
                                 "Indo_B",
                                 "Indo_B",
                                 "Indo_B",
                                 "Med_A",
                                 "Indo_C",
                                 "Indo_B",
                                 "Indo_B",
                                 "Indo_B"))
popNames(zczc125.pop1.gl)

zczc125.pop2.gl<-zczc125.gl
pop(zczc125.pop2.gl)<-as.factor(c("Med_A",
                                  "Med_A",
                                  "Med_A",
                                  "Med_A",
                                  "Med_A",
                                  "Med_A",
                                  "Med_A",
                                  "Med_C",
                                  "Med_A",
                                  "Med_C",
                                  "Med_A",
                                  "Med_C",
                                  "Med_A",
                                  "Med_A",
                                  "Med_B",
                                  "Indo_A",
                                  "Indo_A",
                                  "Atl_BC",
                                  "Atl_BB",
                                  "Med_B",
                                  "Atl_BC",
                                  "Med_B",
                                  "Atl_BC",
                                  "Indo_C",
                                  "Med_B",
                                  "Med_A",
                                  "Indo_C",
                                  "Med_B",
                                  "Med_C",
                                  "Med_C",
                                  "Med_B",
                                  "Med_B",
                                  "Med_B",
                                  "Atl_BB",
                                  "Atl_BB",
                                  "Atl_BB",
                                  "Atl_BB",
                                  "Atl_BB",
                                  "Atl_BC",
                                  "Atl_BB",
                                  "Atl_BB",
                                  "Atl_BC",
                                  "Atl_BC",
                                  "Atl_BC",
                                  "Atl_BC",
                                  "Atl_BC",
                                  "Atl_BB",
                                  "Atl_BC",
                                  "Atl_BA",
                                  "Atl_BC",
                                  "Atl_BB",
                                  "Atl_BD",
                                  "Atl_BC",
                                  "Atl_BA",
                                  "Atl_BD",
                                  "Atl_BB",
                                  "Atl_BB",
                                  "Atl_BC",
                                  "Atl_BC",
                                  "Atl_BC",
                                  "Atl_BB",
                                  "Atl_BA",
                                  "Atl_BA",
                                  "Atl_BB",
                                  "Indo_C",
                                  "Indo_C",
                                  "Atl_BA",
                                  "Med_B",
                                  "Indo_B",
                                  "Atl_BA",
                                  "Indo_C",
                                  "Atl_BA",
                                  "Atl_BC",
                                  "Atl_BB",
                                  "Med_A",
                                  "Med_A",
                                  "Indo_B",
                                  "Med_A",
                                  "Med_B",
                                  "Med_A",
                                  "Atl_BB",
                                  "Atl_BA",
                                  "Atl_BA",
                                  "Atl_BA",
                                  "Atl_BC",
                                  "Atl_BC",
                                  "Atl_BC",
                                  "Atl_BB",
                                  "Atl_BA",
                                  "Atl_BA",
                                  "Atl_BA",
                                  "Atl_BA",
                                  "Atl_BC",
                                  "Indo_C",
                                  "Indo_C",
                                  "Med_A",
                                  "Med_A",
                                  "Indo_C",
                                  "Indo_C",
                                  "Atl_BA",
                                  "Atl_BA",
                                  "Indo_C",
                                  "Indo_C",
                                  "Indo_A",
                                  "Indo_C",
                                  "Indo_C",
                                  "Indo_C",
                                  "Indo_C",
                                  "Indo_C",
                                  "Indo_A",
                                  "Indo_C",
                                  "Indo_C",
                                  "Indo_A",
                                  "Indo_B",
                                  "Indo_B",
                                  "Indo_B",
                                  "Indo_B",
                                  "Indo_B",
                                  "Indo_B",
                                  "Indo_B",
                                  "Med_A",
                                  "Indo_C",
                                  "Indo_B",
                                  "Indo_B",
                                  "Indo_B"))
popNames(zczc125.pop2.gl)
# Phylogenetic trees ------------------------------------------------------

#VERY basic neighbor joining tree
zczc125.tree<-nj(dist(as.matrix(zczc125.gl)))
plot(zczc125.tree)

#modification where tree is a fan and the tips have the individual name and are coloured by ocean basin. OK but hard to read all the indiviual names. 
plot.phylo(x=zczc125.tree, typ="fan", show.tip=FALSE)
zczc125.pops<-zczc125.gl$pop
tiplabels(text=indNames(zczc125.gl), bg=(col=zczc125.pops), cex=0.4, col="white")

#better and includes bootstrap scores, but still touch to read the sample names
library(poppr)
library(RColorBrewer)
zczc125.tree.boot<-aboot(zczc125.gl, tree="upgma", distance = bitwise.dist, sample=100, cutoff=50)
cols <- brewer.pal(n = nPop(zczc125.gl), name = "Dark2")
plot.phylo(zczc125.tree.boot, cex = 0.5, font = 2, adj = 0, tip.color =  cols[pop(zczc125.gl)])
nodelabels(zczc125.tree.boot$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 3, xpd = TRUE)
legend('bottomleft', legend = c("ATL", "INDO-PAC", "MED"), bty="n", fill = cols, border = FALSE, cex =1)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different) - Z. cavirostris")


plot(zczc125.tree, typ="fan", show.tip=FALSE)
tiplabels(pch=20, col=my_zissou, cex=1.5)
legend('bottomleft', legend = c("ATL", "INDO-PAC", "MED"), bty="n", fill = cols, border = FALSE, cex =1)

##Using ggtree and tibble

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
  geom_tippoint(aes(col=factor(pop)), size=3) +
  scale_colour_manual(values=my_zissou[c(3,2,1)]) +
  theme(legend.position=c(0.05,.2), legend.title = element_blank())



ggtree(zczc125.tibtree) +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=my_zissou) +
  theme(legend.position=c(.9,.2), legend.title = element_blank())

ggtree(zczc125.tibtree, layout="slanted") +
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=my_zissou) +
  theme(legend.position=c(.9,.2), legend.title = element_blank())

ggtree(zczc125.tibtree, layout="daylight") + #######not good
  geom_tippoint(aes(col=factor(pop)), size=1.5) +
  scale_colour_manual(values=my_zissou) +
  theme(legend.position=c(0.2,0), legend.title = element_blank())

#make df with sample name and ocean basin to pass into tibble
zczc125.tree$tip.label
zczc125.ocean
zczc125.names
zczc125.ind<-cbind(zczc125.names, zczc125.ocean)
colnames(zczc125.ind)<-c("ind", "ocean")  zczc125.ind



# maps --------------------------------------------------------------------
#basic world map
wm<-map_data("world")
wm<-ggplot(wm, aes(long, lat, group=group)) +
  geom_polygon(fill="grey", colour="white") +
  coord_fixed(1.3) +
  theme_classic()
wm

#make a mediterranean map
install.packages("oceanmap")
library(oceanmap)
library(rworldmap)
library(mapproj)

plotmap("med4", main =" Mediterranean Sea ")

zczc125.map<-read.csv("~/Dropbox/Phd/Bioinformatics/bw_ddrad/SNP_files/no_populations/one_snp_per_locus/zczc_snp_files/zczc125.coord.csv", sep=",", header=TRUE)
head(zczc125.map)


mm<-wm + coord_map(xlim = c(-10,40),ylim = c(30,47))
mm



is.factor(zczc125.map$pop1)
is.numeric(zczc125.map$lat)
is.numeric(zczc125.map$long)

zczc125.atl.map<-subset(zczc125.map, zczc125.map$ocean=="Atlantic")
zczc125.indo.map<-subset(zczc125.map, zczc125.map$ocean=="Indopacific")
zczc125.med.map<-subset(zczc125.map, zczc125.map$ocean=="Mediterranean")

png("zcav.allddRAD.map.png", width = 6, height = 4, units = 'in', res = 300)
zcavwm<-wm + 
  ##this is the information needed to add the data to the plot. 
  geom_jitter(data=zczc125.map, position=position_jitter(width=1, height=1),
              aes(x=long, y=lat, group=ocean, col=ocean), size=1.5) +
  labs(colour=NULL, x="Longitude", y="Latitude") +
  theme(legend.title = element_blank())
##this will show the new map
zcavwm
dev.off()

zcav.atlA.wm<-wm + 
  ##this is the information needed to add the data to the plot. 
  geom_jitter(data=zczc125.atl.map, position=position_jitter(width=1, height=1),
              aes(x=long, y=lat, group=pop1, col=pop1), size=2) +
  labs(colour=NULL, x="Longitude", y="Latitude") +
  theme(legend.title = element_blank())
##this will show the new map
zcav.atlA.wm

zcav.atlB.wm<-wm + 
  ##this is the information needed to add the data to the plot. 
  geom_jitter(data=zczc125.atl.map, position=position_jitter(width=1, height=1),
              aes(x=long, y=lat, group=pop2, col=pop2), size=2) +
  labs(colour=NULL, x="Longitude", y="Latitude") +
  theme(legend.title = element_blank())
##this will show the new map
zcav.atlB.wm

zcav.indo.wm<-wm + 
  ##this is the information needed to add the data to the plot. 
  geom_jitter(data=zczc125.indo.map, position=position_jitter(width=1, height=1),
              aes(x=long, y=lat, group=pop1, col=pop1), size=2) +
  labs(colour=NULL, x="Longitude", y="Latitude") +
  theme(legend.title = element_blank())
##this will show the new map
zcav.indo.wm

zcav.med.wm<-mm + 
  ##this is the information needed to add the data to the plot. 
  geom_jitter(data=zczc125.med.map, position=position_jitter(width=.1, height=.1),
              aes(x=long, y=lat, group=pop1, col=pop1), size=2) +
  labs(colour=NULL, x="Longitude", y="Latitude") +
  theme(legend.title = element_blank())
##this will show the new map
zcav.med.wm

# Genepop -----------------------------------------------------------------

install.packages("genepop")
library(genepop)
#used pgdspider to convert from vcf to genepop format

#Calculates basic stats for each locus. Allele and genotype frequencies per locus. Not very informative because there are SO many loci and alleles 
basic_info("./zczc125_genepop.txt", outputFile="./zczc125.genepop.basicstats.txt")


#Genotypic Differentiation. Exact conditional contingency-table tests for genotypic differentiation (due to genic=false command). A single test for all populations or distinct tests for all pairs of populations (pairs=TRUE). 
test_diff("./zczc125.genepop.txt", genic=FALSE, pairs=TRUE, outputFile = "./zczc125.genepop.diff.txt")

#
zczc125.gp.Fis<-genedivFis("./zczc125_genepop.txt")
zczc125.pop1.gp.Fis<-genedivFis("./zczc125.pop1.genepop.txt")
zczc125.pop2.gp.Fis<-genedivFis("./zczc125.pop2.genepop.txt")

###Update genepop files with population as sample name:
devtools::install_github("rystanley/genepopedit") 
library(genepopedit) # load the library

#
genepop_detective("./zczc125.genepop.txt", variable = "Pops")
#[1] "MED"  "INDO" "ATL" 
genepop_detective("./zczc125.genepop.txt", variable = "PopNum")
# Population Freq
# 1        ATL   55
# 2       INDO   36
# 3        MED   34

#Population pair               Chi2      df   P-Value
#Zcav20182_L7- & Zcav20183_L12 >103328.0229542888Highly sign.
#Zcav20182_L7- & Zcav20182_L9- >115734.1542546710Highly sign.
#Zcav20183_L12 & Zcav20182_L9- >67159.6120748816Highly sign.


#Fst: evaluates Fst or related measures based on allele sizes, for all populations or for all pairs of populations
Fst("./zczc125.genepop.txt", pairs=TRUE, outputFile = "./zczc125.genepop.fst.txt", dataType = "Diploid")

#Fis and Gene Diversity
genedivFis("./zczc125.genepop.txt", outputFile = "./zczc125.genepop.fis.txt", dataType = "Diploid")

#Estimation of Nm by private allele method
Nm_private("./zczc125.genepop.txt", outputFile = "./zczc125.genepop.nmprivate.txt", dataType = "Diploid")

# Demerelate ---------------------------------------------------------------
zczc125.deme<-zczc125.geno
zczc125.deme[1:5, 1:5]
zczc125.deme<-cbind(zczc125.ind, zczc125.deme)
zczc125.deme[1:5, 1:5]
colnames(zczc125.deme)[1:2]<-c("Individual", "population")
zczc125.deme[1:5, 1:5]

zczc125.demerelate<-Demerelate(inputdata=zczc125.deme, Fis=TRUE, iteration=100)

zczc125.deme<-gl2demerelate(zczc125.gl)
zczc125.demerelate<-Demerelate(inputdata=zczc125.deme, Fis=TRUE, iteration=100)



# Adegenet ----------------------------------------------------------------

zczc125.gi
zczc125.gi.sum<-summary(zczc125.gi)
boxplot(zczc125.gi.sum$Hobs)

#make a genpop file in adegenet
zczc125.gp<-genind2genpop(zczc125.gi)
zczc125.pop1.gp<-genind2genpop(zczc125.pop1.gi)
zczc125.pop2.gp<-genind2genpop(zczc125.pop2.gi)

zczc125.gl.sum<-summary(zczc125.gl)
zczc125.gp.sum<-summary(zczc125.gp)


#Calculate Hs for each ocean basin using genpop file
Hs(zczc125.gp)
# Atlantic   Indopacific Mediterranean 
# 0.13569395    0.13068689    0.09758311 
Hs(zczc125.pop1.gp)
# Atl_AA     Atl_AB     Atl_AC     Atl_AD     Atl_AE     Indo_A     Indo_B     Indo_C      Med_A      Med_B      Med_C 
# 0.12909228 0.12075041 0.13128402 0.13173343 0.07395933 0.11250089 0.12674538 0.12796312 0.09385036 0.08047348 0.09431465 
Hs(zczc125.pop2.gp)
# Atl_BA     Atl_BB     Atl_BC     Atl_BD     Indo_A     Indo_B     Indo_C      Med_A      Med_B      Med_C 
# 0.13202948 0.13128402 0.13173343 0.07395933 0.11250089 0.12674538 0.12796312 0.09385036 0.08047348 0.09431465 

#Observed heterozygosity
lapply(seppop(zczc125.gi), function(e) mean(summary(e)$Hobs, na.rm = TRUE))
# $Atlantic
# [1] 0.1184717
# 
# $Indopacific
# [1] 0.115993
# 
# $Mediterranean
# [1] 0.08499716

lapply(seppop(zczc125.pop1.gi), function(e) mean(summary(e)$Hobs, na.rm = TRUE))
# $Atl_AA
# [1] 0.1183948
# 
# $Atl_AB
# [1] 0.1095007
# 
# $Atl_AC
# [1] 0.119318
# 
# $Atl_AD
# [1] 0.1192011
# 
# $Atl_AE
# [1] 0.1275313
# 
# $Indo_A
# [1] 0.1085522
# 
# $Indo_B
# [1] 0.1162764
# 
# $Indo_C
# [1] 0.1178388
# 
# $Med_A
# [1] 0.08862283
# 
# $Med_B
# [1] 0.07775722
# 
# $Med_C
# [1] 0.0859682

lapply(seppop(zczc125.pop2.gi), function(e) mean(summary(e)$Hobs, na.rm = TRUE))
# $Atl_BA
# [1] 0.1156869
# 
# $Atl_BB
# [1] 0.119318
# 
# $Atl_BC
# [1] 0.1192011
# 
# $Atl_BD
# [1] 0.1275313
# 
# $Indo_A
# [1] 0.1085522
# 
# $Indo_B
# [1] 0.1162764
# 
# $Indo_C
# [1] 0.1178388
# 
# $Med_A
# [1] 0.08862283
# 
# $Med_B
# [1] 0.07775722
# 
# $Med_C
# [1] 0.0859682


# poppr -------------------------------------------------------------------

#calculate diversity statistics on a matrix containing counts of multilocus genotypes per population. Gives same gene diversity value as hierfstat!
zczc125.poppr<-poppr(zczc125.gi)
# Pop   N MLG eMLG       SE    H   G lambda E.5   Hexp    Ia   rbarD       File
# 1      Atlantic  55  55   34 1.12e-06 4.01  55  0.982   1 0.1370  47.0 0.00240 zczc125.gi
# 2   Indopacific  36  36   34 3.94e-08 3.58  36  0.972   1 0.1326  41.5 0.00240 zczc125.gi
# 3 Mediterranean  34  34   34 0.00e+00 3.53  34  0.971   1 0.0992 198.5 0.02251 zczc125.gi
# 4         Total 125 125   34 0.00e+00 4.83 125  0.992   1 0.1362 136.0 0.00645 zczc125.gi

zczc125.pop1.poppr<-poppr(zczc125.pop1.gi)
# Pop   N MLG eMLG       SE     H   G lambda E.5   Hexp     Ia    rbarD            File
# 1  Atl_AA  11  11   10 0.00e+00 2.398  11  0.909   1 0.1354   2.58 0.000199 zczc125.pop1.gi
# 2  Atl_AB   5   5    5 0.00e+00 1.609   5  0.800   1 0.1351  88.62 0.009492 zczc125.pop1.gi
# 3  Atl_AC  17  17   10 2.31e-07 2.833  17  0.941   1 0.1354  27.74 0.001849 zczc125.pop1.gi
# 4  Atl_AD  20  20   10 0.00e+00 2.996  20  0.950   1 0.1353  70.93 0.004395 zczc125.pop1.gi
# 5  Atl_AE   2   2    2 0.00e+00 0.693   2  0.500   1 0.0864     NA       NA zczc125.pop1.gi
# 6  Indo_A   5   5    5 0.00e+00 1.609   5  0.800   1 0.1257  40.60 0.004797 zczc125.pop1.gi
# 7  Indo_B  12  12   10 0.00e+00 2.485  12  0.917   1 0.1324   2.76 0.000206 zczc125.pop1.gi
# 8  Indo_C  19  19   10 2.51e-07 2.944  19  0.947   1 0.1316  62.04 0.004219 zczc125.pop1.gi
# 9   Med_A  19  19   10 2.51e-07 2.944  19  0.947   1 0.0966  37.25 0.005043 zczc125.pop1.gi
# 10  Med_B  10  10   10 0.00e+00 2.303  10  0.900   1 0.0853 158.01 0.026901 zczc125.pop1.gi
# 11  Med_C   5   5    5 0.00e+00 1.609   5  0.800   1 0.1054 707.75 0.104795 zczc125.pop1.gi
# 12  Total 125 125   10 9.88e-06 4.828 125  0.992   1 0.1362 135.97 0.006453 zczc125.pop1.gi

zczc125.pop2.poppr<-poppr(zczc125.pop2.gi)
# Pop   N MLG eMLG       SE     H   G lambda E.5   Hexp     Ia    rbarD            File
# 1  Atl_BA  16  16   10 0.00e+00 2.773  16  0.938   1 0.1364  35.50 0.002389 zczc125.pop2.gi
# 2  Atl_BB  17  17   10 2.31e-07 2.833  17  0.941   1 0.1354  27.74 0.001849 zczc125.pop2.gi
# 3  Atl_BC  20  20   10 0.00e+00 2.996  20  0.950   1 0.1353  70.93 0.004395 zczc125.pop2.gi
# 4  Atl_BD   2   2    2 0.00e+00 0.693   2  0.500   1 0.0864     NA       NA zczc125.pop2.gi
# 5  Indo_A   5   5    5 0.00e+00 1.609   5  0.800   1 0.1257  40.60 0.004797 zczc125.pop2.gi
# 6  Indo_B  12  12   10 0.00e+00 2.485  12  0.917   1 0.1324   2.76 0.000206 zczc125.pop2.gi
# 7  Indo_C  19  19   10 2.51e-07 2.944  19  0.947   1 0.1316  62.04 0.004219 zczc125.pop2.gi
# 8   Med_A  19  19   10 2.51e-07 2.944  19  0.947   1 0.0966  37.25 0.005043 zczc125.pop2.gi
# 9   Med_B  10  10   10 0.00e+00 2.303  10  0.900   1 0.0853 158.01 0.026901 zczc125.pop2.gi
# 10  Med_C   5   5    5 0.00e+00 1.609   5  0.800   1 0.1054 707.75 0.104795 zczc125.pop2.gi
# 11  Total 125 125   10 9.88e-06 4.828 125  0.992   1 0.1362 135.97 0.006453 zczc125.pop2.gi

#make poppr file format from genlight file
#zczc125.mlg.table<-mlg.table(zczc125.gi)
#diversity functions in poppr give same results as the poppr command
# zczc125.poppr.diversity<-diversity_stats(zczc125.mlg.table)
# zczc125.poppr.diversity.ci<-diversity_ci(zczc125.mlg.table)

##Private loci - NOT WORKING
# private_alleles(zczc125.gi, locus ~ zczc125.gi$pop, level="population")




# Tajimas D ---------------------------------------------------------------

#Ocean Basins
zczc123.tajima.atl<-read.table("./ocean_basins/zczc123/zczc123.atl.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.atl$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.59765 -0.80569 -0.60031 -0.27938  0.05382  2.56639 

png("./ocean_basins/zczc123.tajima.atl.png", width = 6, height = 4, units = 'in', res = 300)
zczc123.tajima.atl.plot<-ggplot(zczc123.tajima.atl, aes(x=zczc123.tajima.atl$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc123.tajima.atl$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc123.tajima.atl$TajimaD, na.rm=TRUE)) +1, y=2500, 
           label=(paste("Z. cavirostris, Atlantic, mean Tajima's D=", round(mean(zczc123.tajima.atl$TajimaD, na.rm=TRUE), 4)))) +
  annotate("text", x=(mean(zczc123.tajima.atl$TajimaD, na.rm=TRUE)) +1, y=2300, 
           label=(paste("n=55"))) +
  xlab("Tajima's D")
zczc123.tajima.atl.plot
dev.off()

zczc123.tajima.indo<-read.table("./ocean_basins/zczc123/zczc123.pac.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.indo$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.65323 -0.90334 -0.59550 -0.30579  0.09248  2.39375 
png("./ocean_basins/zczc123.tajima.indo.png", width = 6, height = 4, units = 'in', res = 300)
zczc123.tajima.indo.plot<-ggplot(zczc123.tajima.indo, aes(x=zczc123.tajima.indo$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc123.tajima.indo$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc123.tajima.indo$TajimaD, na.rm=TRUE) +1.25), y=4000, 
           label=(paste("Z. cavirostris, Indo-Pacific, mean Tajima's D=", round(mean(zczc123.tajima.indo$TajimaD, na.rm=TRUE), 4)))) +
  annotate("text", x=(mean(zczc123.tajima.indo$TajimaD, na.rm=TRUE) +1.25), y=3800, 
           label=(paste("n=36"))) +
  xlab("Tajima's D")
zczc123.tajima.indo.plot
dev.off()

zczc123.tajima.med<-read.table("./ocean_basins/zczc123/zczc123.med.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.med$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.4656 -0.4056  0.4502  0.4430  1.3694  2.3661 
png("./ocean_basins/zczc123.tajima.med.png", width = 6, height = 4, units = 'in', res = 300)
zczc123.tajima.med.plot<-ggplot(zczc123.tajima.med, aes(x=zczc123.tajima.med$TajimaD)) +
  geom_histogram(binwidth=.25, alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept=mean(zczc123.tajima.med$TajimaD, na.rm=TRUE), linetype="dashed", size=1) + 
  annotate("text", x=(mean(zczc123.tajima.med$TajimaD, na.rm=TRUE)), y=1800, 
           label=(paste("Z. cavirostris, Mediterranean, mean Tajima's D=", round(mean(zczc123.tajima.med$TajimaD, na.rm=TRUE), 4)))) +
  annotate("text", x=(mean(zczc123.tajima.med$TajimaD, na.rm=TRUE)), y=1500, 
           label=(paste("n=34"))) +
  xlab("Tajima's D")
zczc123.tajima.med.plot
dev.off()

#Populations
zczc123.tajima.Atl_CanIs<-read.table("./populations/zcav/zczc123.pops/zczc123.Atl_CanIs.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.Atl_CanIs$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.85641 -1.11603 -0.52823 -0.32975  0.08512  2.17739 

zczc123.tajima.Atl_Carib<-read.table("./populations/zcav/zczc123.pops/zczc123.Atl_Carib.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.Atl_Carib$TajimaD)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.7318 -1.1470 -0.6515 -0.3039  0.2157  2.1490 

zczc123.tajima.Atl_France<-read.table("./populations/zcav/zczc123.pops/zczc123.Atl_France.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.Atl_France$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.5622 -1.1117 -0.1952 -0.2496  0.8048  1.8443 

zczc123.tajima.Atl_NE<-read.table("./populations/zcav/zczc123.pops/zczc123.Atl_NE.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.Atl_NE$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.73324 -1.15142 -0.74099 -0.37413  0.08083  2.16958 

zczc123.tajima.Atl_Spain<-read.table("./populations/zcav/zczc123.pops/zczc123.Atl_Spain.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.Atl_Spain$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.6124  1.6330  1.6330  1.6117  1.6330  2.0119

# zczc125.tajima.Atl_BA<-read.table("./populations/zcav/zczc123.pop/zczc125.Atl_BA.100000.Tajima.D", header=TRUE)
# summary(zczc125.tajima.Atl_BA$TajimaD)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # -1.7295 -1.1177 -0.4483 -0.2404  0.4071  2.1592 
# 
# zczc125.tajima.Atl_BB<-read.table("./populations/zcav/zczc123.pop/zczc125.Atl_BB.100000.Tajima.D", header=TRUE)
# summary(zczc125.tajima.Atl_BB$TajimaD)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # -1.7267 -1.1160 -0.4827 -0.2379  0.3363  2.1774 
# 
# zczc125.tajima.Atl_BC<-read.table("./populations/zcav/zczc123.pop/zczc125.Atl_BC.100000.Tajima.D", header=TRUE)
# summary(zczc125.tajima.Atl_BC$TajimaD)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # -1.7004 -1.1086 -0.5635 -0.2803  0.3024  2.2234 
# 
# zczc125.tajima.Atl_BD<-read.table("./populations/zcav/zczc125.pop2/zczc125.Atl_BD.100000.Tajima.D", header=TRUE)
#summary(zczc125.tajima.Atl_BD$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.6124  1.6330  1.6330  1.6002  1.6330  1.8931 

zczc123.tajima.Indo_Cent<-read.table("./populations/zcav/zczc123.pops/zczc123.Indo_Cent.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.Indo_Cent$TajimaD)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.56222 -1.11173  0.01499 -0.17730  0.81980  1.84427 

zczc123.tajima.Indo_Mix<-read.table("./populations/zcav/zczc123.pops/zczc123.Indo_Mix.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.Indo_Mix$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.2331 -0.9330 -0.9330 -0.1522  0.8506  1.7532 

zczc123.tajima.Indo_NE<-read.table("./populations/zcav/zczc123.pops/zczc123.Indo_NE.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.Indo_NE$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.7200 -1.0921 -0.5397 -0.3021  0.2138  2.3930  

zczc123.tajima.Indo_South<-read.table("./populations/zcav/zczc123.pops/zczc123.Indo_South.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.Indo_South$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.7130 -1.1647 -0.5290 -0.3365  0.3147  1.9760 

zczc123.tajima.Med_West<-read.table("./populations/zcav/zczc123.pops/zczc123.Med_West.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.Med_West$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.4911 -0.4766  0.4303  0.4072  1.2879  2.4023 

zczc123.tajima.Med_East<-read.table("./populations/zcav/zczc123.pops/zczc123.Med_East.100000.Tajima.D", header=TRUE)
summary(zczc123.tajima.Med_East$TajimaD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.5106 -0.3634  0.5724  0.3855  1.3184  2.1294 

#zczc125.tajima.Med_C<-read.table("./populations/zcav/zczc123.pop/zczc125.Med_C.100000.Tajima.D", header=TRUE)
#summary(zczc125.tajima.Med_C$TajimaD)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.41614 -1.11173  0.01499  0.02714  0.81980  1.84427 

zczc123.tajima.atl
zczc123.tajima.indo
zczc123.tajima.med

zczc123.tajima.Atl_CanIs
zczc123.tajima.Atl_Carib
zczc123.tajima.Atl_France
zczc123.tajima.Atl_NE
zczc123.tajima.Atl_Spain
zczc123.tajima.Indo_Cent
zczc123.tajima.Indo_Mix
zczc123.tajima.Indo_NE
zczc123.tajima.Indo_South
zczc123.tajima.Med_East
zczc123.tajima.Med_West

quantileCI(zczc123.tajima.atl$TajimaD, level=0.95, method="normal", digits=3)#-0.6   -0.591
quantileCI(zczc123.tajima.indo$TajimaD, digits=3)# -0.596   -0.596
quantileCI(zczc123.tajima.med$TajimaD, method="normal", digits=3)#
quantileCI(zczc123.tajima.Atl_CanIs, method="normal", digits=3)#
quantileCI(zczc123.tajima.Atl_Carib, method="normal", digits=3)#
quantileCI(zczc123.tajima.Atl_France, method="normal", digits=3)#
quantileCI(zczc123.tajima.Atl_NE, method="normal", digits=3)#
quantileCI(zczc123.tajima.Atl_Spain, method="normal", digits=3)#
quantileCI(zczc123.tajima.Indo_Cent, method="normal", digits=3)#
quantileCI(zczc123.tajima.Indo_Mix, method="normal", digits=3)#
quantileCI(zczc123.tajima.Indo_NE, method="normal", digits=3)#
quantileCI(zczc123.tajima.Indo_South, method="normal", digits=3)#
quantileCI(zczc123.tajima.Med_East, method="normal", digits=3)#
quantileCI(zczc123.tajima.Med_West, method="normal", digits=3)#

groupwiseMean(TajimaD~1, data=zczc123.tajima.atl, percentile = TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.indo, percentile=TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.med, percentile=TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.Atl_CanIs, percentile=TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.Atl_Carib, percentile=TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.Atl_France, percentile=TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.Atl_NE, percentile=TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.Atl_Spain, percentile=TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.Indo_Cent, percentile=TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.Indo_Mix, percentile=TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.Indo_NE, percentile=TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.Indo_South, percentile=TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.Med_East, percentile=TRUE)
groupwiseMean(TajimaD~1, data=zczc123.tajima.Med_West, percentile=TRUE)

if(!require(Rmisc)){install.packages("Rmisc")}
if(!require(DescTools)){install.packages("DescTools")}
if(!require(plyr)){install.packages("plyr")}
if(!require(boot)){install.packages("boot")}
if(!require(rcompanion)){install.packages("rcompanion")}






# zczc123 tess3r ----------------------------------------------------------
zczc123.gl

#"/Users/aubrieonoufriou/Dropbox/Phd/Bioinformatics/bw_ddrad_global/Results/Tess3r"

#all 123 zcavs
zczc123.geno<-as.matrix(zczc123.gl)
zczc125.coord
#remove the rows with these names:
#Zcav20181_L3-30, row 10
#Zcav20181_L4-40, row 62
zczc123.coord<-zczc125.coord[c(1:9,11:61, 63:125),]

zczc123.geno
zczc123.coord

k<-10
tess3.zczc123<-tess3(X=zczc123.geno, coord=zczc123.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)

#The plot function generates a plot for root mean-squared errors computed on a subset of loci used for cross-validation.
# specify crossentropy with error bars 
plot(tess3.zczc123, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score")

for(i in 2:k) {
  ppi<-300
  png(paste0("zczc123.barplots.sorted.k",i,".png"), width=10*ppi,height=6*ppi, res=ppi)
  Q.matrix<-qmatrix(tess3.zczc123, K = i)
  rownames(Q.matrix)<-rownames(zczc123.geno)
  
  par(mar=c(10,5,5,0))
  barplot(Q.matrix, sort.by.Q = TRUE, 
          border = NA, space = 0,
          col.palette = my.palette2,
          las=3,
          cex.axis=0.8,
          ylab = "Ancestry coefficients", main=paste0("Admixture Proportions, Z. cavirostris, k=",i,""))
  
  filename=paste0("zczc123.Q.matrix.k=",i,".csv")
  write.table(Q.matrix, file = filename)
  dev.off()
}

#atlantic only
zczc123.atl.gi
zczc123.atl.gl<-gi2gl(zczc123.atl.gi)

zczc123.atl.gl
zczc123.atl.geno<-as.matrix(zczc123.atl.gl)
zczc125.atl.coord

cbind(rownames(zczc123.atl.geno), rownames(zczc125.atl.coord))
#remove row 33 from coord file
zczc123.atl.coord<-zczc125.atl.coord[c(1:32,34:55),]

cbind(rownames(zczc123.atl.geno), rownames(zczc123.atl.coord))

zczc123.atl.coord
zczc123.atl.geno

k<-10
tess3.zczc123.atl<-tess3(X=zczc123.atl.geno, coord=zczc123.atl.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)

#The plot function generates a plot for root mean-squared errors computed on a subset of loci used for cross-validation.
# specify crossentropy with error bars 
plot(tess3.zczc123.atl, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score")

for(i in 2:k) {
  ppi<-300
  png(paste0("zczc123.atl.barplots.sorted.k",i,".png"), width=10*ppi,height=6*ppi, res=ppi)
  Q.matrix<-qmatrix(tess3.zczc123.atl, K = i)
  rownames(Q.matrix)<-rownames(zczc123.atl.geno)
  
  par(mar=c(10,5,5,0))
  barplot(Q.matrix, sort.by.Q = TRUE, 
          border = NA, space = 0,
          col.palette = my.palette2,
          las=3,
          cex.axis=0.8,
          ylab = "Ancestry coefficients", main=paste0("Admixture Proportions, Z. cavirostris, Atlantic, k=",i,""))
  
  filename=paste0("zczc123.atl.Q.matrix.k=",i,".csv")
  write.table(Q.matrix, file = filename)
  dev.off()
}

#Med only

zczc123.med.gi
zczc123.med.gl<-gi2gl(zczc123.med.gi)

zczc123.med.gl
zczc123.med.geno<-as.matrix(zczc123.med.gl)
zczc125.med.coord

cbind(rownames(zczc123.med.geno), rownames(zczc125.med.coord))
#remove row 10 from coord file
zczc123.med.coord<-zczc125.med.coord[c(1:9,11:34),]

cbind(rownames(zczc123.med.geno), rownames(zczc123.med.coord))

zczc123.med.coord
zczc123.med.geno

k<-10
tess3.zczc123.med<-tess3(X=zczc123.med.geno, coord=zczc123.med.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)

#The plot function generates a plot for root mean-squared errors computed on a subset of loci used for cross-validation.
# specify crossentropy with error bars 
plot(tess3.zczc123.med, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score")

for(i in 2:k) {
  ppi<-300
  png(paste0("zczc123.med.barplots.sorted.k",i,".png"), width=10*ppi,height=6*ppi, res=ppi)
  Q.matrix<-qmatrix(tess3.zczc123.med, K = i)
  rownames(Q.matrix)<-rownames(zczc123.med.geno)
  
  par(mar=c(10,5,5,0))
  barplot(Q.matrix, sort.by.Q = TRUE, 
          border = NA, space = 0,
          col.palette = my.palette2,
          las=3,
          cex.axis=0.8,
          ylab = "Ancestry coefficients", main=paste0("Admixture Proportions, Z. cavirostris, medantic, k=",i,""))
  
  filename=paste0("zczc123.med.Q.matrix.k=",i,".csv")
  write.table(Q.matrix, file = filename)
  dev.off()
}

#pac only

zczc123.pac.gi
zczc123.pac.gl<-gi2gl(zczc123.pac.gi)

zczc123.pac.gl
zczc123.pac.geno<-as.matrix(zczc123.pac.gl)
zczc123.pac.coord<-zczc125.indo.coord

zczc123.pac.coord
zczc123.pac.geno

k<-10
tess3.zczc123.pac<-tess3(X=zczc123.pac.geno, coord=zczc123.pac.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)

#The plot function generates a plot for root mean-squared errors computed on a subset of loci used for cross-validation.
# specify crossentropy with error bars 
plot(tess3.zczc123.pac, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score")

for(i in 2:k) {
  ppi<-300
  png(paste0("zczc123.pac.barplots.sorted.k",i,".png"), width=10*ppi,height=6*ppi, res=ppi)
  Q.matrix<-qmatrix(tess3.zczc123.pac, K = i)
  rownames(Q.matrix)<-rownames(zczc123.pac.geno)
  
  par(mar=c(10,5,5,0))
  barplot(Q.matrix, sort.by.Q = TRUE, 
          border = NA, space = 0,
          col.palette = my.palette2,
          las=3,
          cex.axis=0.8,
          ylab = "Ancestry coefficients", main=paste0("Admixture Proportions, Z. cavirostris, pacantic, k=",i,""))
  
  filename=paste0("zczc123.pac.Q.matrix.k=",i,".csv")
  write.table(Q.matrix, file = filename)
  dev.off()
}


par(mfrow=c(2,2))
plot(tess3.zczc123, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations, all Z. Cavirostris",
     ylab = "Cross-entropy score")

plot(tess3.zczc123.atl, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations, Atlantic",
     ylab = "Cross-entropy score")

plot(tess3.zczc123.med, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations, Mediterranean",
     ylab = "Cross-entropy score")

plot(tess3.zczc123.pac, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations, Indo-Pacific",
     ylab = "Cross-entropy score")
par(mfrow=c(1,1))

##for Mden
mden43.gl
mden43.coord


mden43.geno<-as.matrix(mden43.gl)
mden43.geno[1:5,1:5]
head(mden43.coord)
mden43.tess.coord<-mden43.coord[,c(3,6,7)]
rownames(mden43.tess.coord)<-mden43.tess.coord[,1]
mden43.tess.coord<-mden43.tess.coord[,2:3]

#to use:
mden43.geno
mden43.tess.coord<-as.matrix(mden43.tess.coord)

mden43.seppop.gl<-seppop(mden43.gl, drop=TRUE)
mden43.atl.gl<-mden43.seppop.gl$Atlantic
mden43.pac.gl<-mden43.seppop.gl$Indopacific
mden43.atl.geno<-as.matrix(mden43.atl.gl)
mden43.pac.geno<-as.matrix(mden43.pac.gl)
mden43.atl.coord<-subset(mden43.coord, mden43.coord$ocean=="Atlantic")
mden43.pac.coord<-subset(mden43.coord, mden43.coord$ocean=="Indo-Pacific")

mden43.atl.coord<-mden43.atl.coord[,c(3,6,7)]
rownames(mden43.atl.coord)<-mden43.atl.coord[,1]
mden43.atl.coord<-mden43.atl.coord[,2:3]
mden43.atl.coord<-as.matrix(mden43.atl.coord)

mden43.pac.coord<-mden43.pac.coord[,c(3,6,7)]
rownames(mden43.pac.coord)<-mden43.pac.coord[,1]
mden43.pac.coord<-mden43.pac.coord[,2:3]
mden43.pac.coord<-as.matrix(mden43.pac.coord)

k<-10
tess3.mden43<-tess3(X=mden43.geno, coord=mden43.tess.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)
tess3.mden43.atl<-tess3(X=mden43.atl.geno, coord=mden43.atl.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)
tess3.mden43.pac<-tess3(X=mden43.pac.geno, coord=mden43.pac.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)

par(mfrow=c(2,2))
plot(tess3.mden43, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations, all M. densirostris",
     ylab = "Cross-entropy score")
plot(tess3.mden43.atl, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations, Atlantic",
     ylab = "Cross-entropy score")
plot(tess3.mden43.pac, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations, Indo-Pacific",
     ylab = "Cross-entropy score")


#The plot function generates a plot for root mean-squared errors computed on a subset of loci used for cross-validation.
# specify crossentropy with error bars 
plot(tess3.zczc123, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-entropy score")
# zczc123 basic stats/fis/confidence intervals -----------------------------------------------------------------
#files
zczc123.gl
zczc123.pop.gl<-zczc123.gl
pop(zczc123.pop.gl)<-as.factor(c("Med_West",
                                 "Med_West",
                                 "Med_West",
                                 "Med_West",
                                 "Med_West",
                                 "Med_West",
                                 "Med_West",
                                 "Med_East",
                                 "Med_West",
                                 "Med_West",
                                 "Med_East",
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
                                 "Med_East",
                                 "Med_East",
                                 "Med_East",
                                 "Med_East",
                                 "Med_East",
                                 "Atl_CanIs",
                                 "Atl_CanIs",
                                 "Atl_CanIs",
                                 "Atl_CanIs",
                                 "Atl_CanIs",
                                 "Atl_France",
                                 "Atl_CanIs",
                                 "Atl_CanIs",
                                 "Atl_France",
                                 "Atl_France",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_CanIs",
                                 "Atl_NE",
                                 "Atl_Carib",
                                 "Atl_France",
                                 "Atl_CanIs",
                                 "Atl_Spain",
                                 "Atl_France",
                                 "Atl_Carib",
                                 "Atl_Spain",
                                 "Atl_CanIs",
                                 "Atl_CanIs",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_CanIs",
                                 "Atl_Carib",
                                 "Atl_CanIs",
                                 "Indo_NE",
                                 "Indo_NE",
                                 "Atl_Carib",
                                 "Med_East",
                                 "Indo_Mix",
                                 "Atl_Carib",
                                 "Indo_NE",
                                 "Atl_Carib",
                                 "Atl_NE",
                                 "Atl_CanIs",
                                 "Med_West",
                                 "Med_West",
                                 "Indo_Sou",
                                 "Med_West",
                                 "Med_East",
                                 "Med_West",
                                 "Atl_CanIs",
                                 "Atl_Carib",
                                 "Atl_Carib",
                                 "Atl_Carib",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_NE",
                                 "Atl_CanIs",
                                 "Atl_Carib",
                                 "Atl_Carib",
                                 "Atl_Carib",
                                 "Atl_Carib",
                                 "Atl_NE",
                                 "Indo_NE",
                                 "Indo_NE",
                                 "Med_West",
                                 "Med_West",
                                 "Indo_NE",
                                 "Indo_NE",
                                 "Atl_Carib",
                                 "Atl_Carib",
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
                                 "Indo_Mix",
                                 "Indo_Mix",
                                 "Indo_Sou",
                                 "Indo_Sou",
                                 "Indo_Sou",
                                 "Med_West",
                                 "Indo_NE",
                                 "Indo_Sou",
                                 "Indo_Sou",
                                 "Indo_Sou"))

pop(zczc123.pop.gl)
zczc120.pop.gl<-popsub(zczc123.pop.gl, exclude=c("Indo_Mix"))

zczc123.gi<-gl2gi(zczc123.gl)
zczc123.hfstat<-genind2hierfstat(zczc123.gi)

zczc123.pop.gi<-gl2gi(zczc123.pop.gl)
zczc123.pop.hfstat<-genind2hierfstat(zczc123.pop.gi)

zczc120.pop.gi<-gl2gi(zczc120.pop.gl)
zczc120.pop.hfstat<-genind2hierfstat(zczc120.pop.gi)

#calculate basic stats using hierfstat
zczc123.basicstats<-basic.stats(zczc123.hfstat, digits=3)
zczc123.pop.basicstats<-basic.stats(zczc123.pop.hfstat, digits=3)
zczc120.pop.basicstats<-basic.stats(zczc120.pop.hfstat, digits=3)

#calculate Fst and confidence intervals using dartr
zczc123.fst.ci<-gl.fst.pop(zczc123.gl, nboots=100, percent=95)
zczc123.pop.fst.ci<-gl.fst.pop(zczc123.pop.gl, nboots=100, percent=95)
zczc120.pop.fst.ci<-gl.fst.pop(zczc120.pop.gl, nboots=100, percent=95)
zczc123.fst.ci$Bootstraps[,c(1:2,103:106)]
zczc123.pop.fst.ci$Bootstraps[,c(1:2,103:106)]
zczc120.pop.fst.ci$Bootstraps[,c(1:2,103:106)]

#calculate Fst and p-values using stratag
zczc123.gt<-genind2gtypes(zczc123.gi)
zczc123.pop.gt<-genind2gtypes(zczc123.pop.gi)
zczc120.pop.gt<-genind2gtypes(zczc120.pop.gi)
zczc123.struct<-popStructTest(zczc123.gt, stats=c("fst"), nrep=100)
zczc123.pop.struct<-popStructTest(zczc123.pop.gt, stats=c("fst"), nrep=100)
zczc120.pop.struct<-popStructTest(zczc120.pop.gt, stats=c("fst", "d"), nrep=100, type="pairwise")

zczc123.pop.fst.hfstat<-pairwise.WCfst(zczc123.pop.hfstat)

#calculate dA and p-values using stratag
zczc123.d<-popStructTest(zczc123.gt, nrep=1000, stats="d", type="pairwise")
zczc123.pop.d<-popStructTest(zczc123.pop.gt, nrep=1000, stats="d", type="pairwise")

#calculate weir and cockerham fis point and confidence interval estimates with hierfstat
#first need to use seppop to generate separate genind files by population
zczc123.seppop.gi<-seppop(zczc123.gi)
zczc123.pop.seppop.gi<-seppop(zczc123.pop.gi)
zczc123.med.gi<-zczc123.seppop.gi$Mediterranean
zczc123.atl.gi<-zczc123.seppop.gi$Atlantic
zczc123.pac.gi<-zczc123.seppop.gi$`Indo-Pacific`

#calculate fis point estimate
zczc123.med.wc<-wc(zczc123.med.gi)
zczc123.atl.wc<-wc(zczc123.atl.gi)
zczc123.pac.wc<-wc(zczc123.pac.gi)

zczc123.med.wc #0.1219
zczc123.atl.wc #0.1208
zczc123.pac.wc #0.1176

#calculate fis 95% confidence intervals
zczc123.med.ci.fis<-boot.ppfis(zczc123.med.gi, nboot=1000)
zczc123.atl.ci.fis<-boot.ppfis(zczc123.atl.gi, nboot=1000)
zczc123.pac.ci.fis<-boot.ppfis(zczc123.pac.gi, nboot=1000)

zczc123.med.ci.fis #0.1156-0.1279
zczc123.atl.ci.fis #0.1169-0.1247
zczc123.pac.ci.fis #0.1136-0.1217

wc(zczc123.pop.seppop.gi$Atl_CanIs) #0.1151645
wc(zczc123.pop.seppop.gi$Atl_Carib) #0.1214143
wc(zczc123.pop.seppop.gi$Atl_France) #0.1164734
wc(zczc123.pop.seppop.gi$Atl_NE) #0.1097882
wc(zczc123.pop.seppop.gi$Atl_Spain) #-0.9530291
wc(zczc123.pop.seppop.gi$Indo_Cent) #0.1411067
wc(zczc123.pop.seppop.gi$Indo_Mix) #0.1494748
wc(zczc123.pop.seppop.gi$Indo_NE) #0.09693309
wc(zczc123.pop.seppop.gi$Indo_Sou) #0.1101989
wc(zczc123.pop.seppop.gi$Med_East) #0.07982834
wc(zczc123.pop.seppop.gi$Med_West) #0.07726678

zczc120.pop.seppop.gi<-seppop(zczc120.pop.gi)
wc(zczc120.pop.seppop.gi$Atl_CanIs) #0.1151645
wc(zczc120.pop.seppop.gi$Atl_Carib) #0.1214143
wc(zczc120.pop.seppop.gi$Atl_France) #0.1164734
wc(zczc120.pop.seppop.gi$Atl_Spain) #-0.9530291
wc(zczc120.pop.seppop.gi$Atl_NE) #0.1097882
wc(zczc120.pop.seppop.gi$Indo_Cent) #0.1411067
wc(zczc120.pop.seppop.gi$Indo_NE) #0.09693309
wc(zczc120.pop.seppop.gi$Indo_Sou) #0.1101989
wc(zczc120.pop.seppop.gi$Med_East) #0.07982834
wc(zczc120.pop.seppop.gi$Med_West) #0.07726678

boot.ppfis(zczc123.pop.seppop.gi$Atl_CanIs, nboot=1000) #0.1097 0.1205
boot.ppfis(zczc123.pop.seppop.gi$Atl_Carib, nboot=1000) #0.1154 0.1268
boot.ppfis(zczc123.pop.seppop.gi$Atl_France, nboot=1000) #0.1082 0.1254
boot.ppfis(zczc123.pop.seppop.gi$Atl_NE, nboot=1000) #0.1041 0.1154
test<-boot.ppfis(zczc123.pop.seppop.gi$Atl_Spain, nboot=1000) #-1.0352 -1.016
boot.ppfis(zczc123.pop.seppop.gi$Indo_Cent, nboot=1000) #0.1308 0.151
boot.ppfis(zczc123.pop.seppop.gi$Indo_Mix,  nboot=1000) #0.1355 0.1607
boot.ppfis(zczc123.pop.seppop.gi$Indo_NE, nboot=1000) #0.0923 0.1023
boot.ppfis(zczc123.pop.seppop.gi$Indo_Sou, nboot=1000) #0.1037 0.1178
boot.ppfis(zczc123.pop.seppop.gi$Med_East, nboot=1000) #0.0714 0.0883
boot.ppfis(zczc123.pop.seppop.gi$Med_West, nboot=1000) #0.0704 0.0842

boot.ppfis(zczc120.pop.seppop.gi$Atl_CanIs, nboot=1000) #0.1097 0.1205
boot.ppfis(zczc120.pop.seppop.gi$Atl_Carib, nboot=1000) #0.116 0.127
boot.ppfis(zczc120.pop.seppop.gi$Atl_France, nboot=1000) #0.1078 0.1248
boot.ppfis(zczc120.pop.seppop.gi$Atl_Spain, nboot=1000) #-1.0357 -1.0165
boot.ppfis(zczc120.pop.seppop.gi$Atl_NE, nboot=1000) #0.1044 0.1156
boot.ppfis(zczc120.pop.seppop.gi$Indo_Cent, nboot=1000) #0.1305 0.1502
boot.ppfis(zczc120.pop.seppop.gi$Indo_NE, nboot=1000) #0.0923 0.1022
boot.ppfis(zczc120.pop.seppop.gi$Indo_Sou, nboot=1000) #0.1035 0.1173
boot.ppfis(zczc120.pop.seppop.gi$Med_East, nboot=1000) #0.0714 0.0887
boot.ppfis(zczc120.pop.seppop.gi$Med_West, nboot=1000) #0.0708 0.084

#calculate confidence intervals around Hs and Ho
zczc123.basicstats
zczc123.pop.basicstats

zczc123.ho<-zczc123.basicstats$Ho
zczc123.hs<-zczc123.basicstats$Hs
zczc123.fis<-zczc123.basicstats$Fis
zczc123.pop.fis<-zczc123.pop.basicstats$Fis

zczc123.ho.melt<-melt(na.omit(zczc123.ho))
zczc123.ho.melt<-zczc123.ho.melt[,2:3]
colnames(zczc123.ho.melt)<-c("pop", "ho")
zczc123.ho.melt

Sum.zczc123.ho = groupwiseMean(ho ~ pop, data= zczc123.ho.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc123.ho

zczc123.hs.melt<-melt(na.omit(zczc123.hs))
zczc123.hs.melt<-zczc123.hs.melt[,2:3]
colnames(zczc123.hs.melt)<-c("pop", "hs")
zczc123.hs.melt

Sum.zczc123.hs = groupwiseMean(hs ~ pop, data= zczc123.hs.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc123.hs

zczc123.fis.melt<-melt(na.omit(zczc123.fis))
zczc123.fis.melt<-zczc123.fis.melt[,2:3]
colnames(zczc123.fis.melt)<-c("pop", "fis")
zczc123.fis.melt

Sum.zczc123.fis = groupwiseMean(fis ~ pop, data= zczc123.fis.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc123.fis

zczc123.pop.fis.melt<-melt(na.omit(zczc123.pop.fis))
zczc123.pop.fis.melt<-zczc123.pop.fis.melt[,2:3]
colnames(zczc123.pop.fis.melt)<-c("pop", "fis")
zczc123.pop.fis.melt

Sum.zczc123.pop.fis = groupwiseMean(fis ~ pop, data= zczc123.pop.fis.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc123.pop.fis

zczc123.pop.ho<-zczc123.pop.basicstats$Ho
zczc123.pop.hs<-zczc123.pop.basicstats$Hs
zczc123.pop.fis<-zczc123.pop.basicstats$Fis

zczc123.pop.ho.melt<-melt(na.omit(zczc123.pop.ho))
zczc123.pop.ho.melt<-zczc123.pop.ho.melt[,2:3]
colnames(zczc123.pop.ho.melt)<-c("pop", "ho")
zczc123.pop.ho.melt

Sum.zczc123.pop.ho = groupwiseMean(ho ~ pop, data= zczc123.pop.ho.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc123.pop.ho

zczc123.pop.hs.melt<-melt(na.omit(zczc123.pop.hs))
zczc123.pop.hs.melt<-zczc123.pop.hs.melt[,2:3]
colnames(zczc123.pop.hs.melt)<-c("pop", "hs")
zczc123.pop.hs.melt

Sum.zczc123.pop.hs = groupwiseMean(hs ~ pop, data= zczc123.pop.hs.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc123.pop.hs

zczc123.pop.fis.melt<-melt(na.omit(zczc123.pop.fis))
zczc123.pop.fis.melt<-zczc123.pop.fis.melt[,2:3]
colnames(zczc123.pop.fis.melt)<-c("pop", "fis")
zczc123.pop.fis.melt

Sum.zczc123.pop.fis = groupwiseMean(fis ~ pop, data= zczc123.pop.fis.melt, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc123.pop.fis

# isolation-by-distance ---------------------------------------------------
zczc123.gl
zczc123.gi
zczc123.coord
zczc123.coord2<-cbind(zczc123.coord[,2], zczc123.coord[,1])

mden43.gl
mden43.pop.gl
mden43.coord
mden43.loc<-as.data.frame(mden43.coord)

plot(mden43.loc)
map(add=T, interior=F) #x=longitude, y=latitude
colnames(mden43.loc)<-c("long", "lat")
head(mden43.loc)

mden43.coord2<-cbind(mden43.coord[,2], mden43.coord[,1])
head(mden43.coord2)
colnames(mden43.coord2)<-c("lat", "lon")
mden43.gl@other$latlong<-mden43.coord2
mden43.pop.gl@other$latlong<-mden43.coord2
cbind(sort(indNames(mden43.pop.gl)), rownames(mden43.coord2))

zczc123.gl@other$latlong<-zczc123.coord2
colnames(zczc123.gl@other$latlong)<-c("lat","lon")
zczc123.gl@other$latlong
zczc123.ibd<-gl.ibd(zczc123.gl)
zczc123.ibd

zczc123.pop.gl@other$latlong<-zczc123.coord2
colnames(zczc123.pop.gl@other$latlong)<-c("lat","lon")
zczc123.pop.ibd<-gl.ibd(zczc123.pop.gl)
zczc123.pop.ibd

zczc123.nopop.gl<-zczc123.gl
pop(zczc123.nopop.gl)<-c(1:123)
pop(zczc123.nopop.gl)
zczc123.nopop.ibd<-gl.ibd(zczc123.nopop.gl)

mden.pop.ibd<-gl.ibd(mden43.pop.gl)

#plotting ibd results showing the density of ploints. local density is measured using a 2-d kernel density estimate
library(MASS)
dens<-kde2d(zczc123.pop.ibd$Dgeo, zczc123.pop.ibd$Dgen, n=200)
myPal<-colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(zczc123.pop.ibd$Dgeo, zczc123.pop.ibd$Dgen, pch=20, cex=0.5)
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(zczc123.pop.ibd$Dgen~zczc123.pop.ibd$Dgeo))
title("Isolation by distance plot")

gl<-testset.gl

zczc123.gp<-genind2genpop(zczc123.gi)

zczc123.ibd.gen<-dist(zczc123.geno)
zczc123.ibd.geo<-dist(zczc123.coord)

zczc123.ibd<-mantel.randtest(zczc123.ibd.gen, zczc123.ibd.geo)

zczc123.nopop.gen<-dist(zczc123.nopop.gl)
zczc123.nopop.ibd<-mantel.randtest(zczc123.nopop.gen, zczc123.ibd.geo)
zczc123.nopop.ibd

library(MASS)
dens<-kde2d(zczc123.ibd.geo, zczc123.nopop.gen, n=500)
myPal<-colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(zczc123.ibd.geo, zczc123.nopop.gen, pch=20, cex=0.5)
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(zczc123.nopop.gen~zczc123.ibd.geo))
title("Isolation by distance plot")

#IBD within ocean basins

zczc123.atl.gl<-popsub(zczc123.pop.gl, sublist=c("Atl_CanIs", "Atl_Carib", "Atl_France", "Atl_NE", "Atl_Spain"))
zczc123.atl.gl$pop

zczc123.pac.gl<-popsub(zczc123.pop.gl, sublist=c("Indo_Cent", "Indo_Mix", "Indo_NE", "Indo_Sou"))
zczc123.pac.gl$pop

zczc123.med.gl<-popsub(zczc123.pop.gl, sublist=c("Med_East", "Med_West"))
zczc123.med.gl$pop

zczc123.nomed.gl<-popsub(zczc123.pop.gl, sublist=c("Atl_CanIs", "Atl_Carib", "Atl_France", "Atl_NE", "Atl_Spain","Indo_Cent", "Indo_Mix", "Indo_NE", "Indo_Sou"))

mden43.atl.gl<-popsub(mden43.pop.gl, sublist=c("Atl-Bahamas", "Atl-East", "Atl-Other"))
mden43.pac.gl<-popsub(mden43.pop.gl, sublist=c("Indo-Africa", "Indo-Hawaii", "Indo-South"))


atl.ibd<-gl.ibd(zczc123.atl.gl)
pac.ibd<-gl.ibd(zczc123.pac.gl)
med.ibd<-gl.ibd(zczc123.med.gl)
nomed.ibd<-gl.ibd(zczc123.nomed.gl)

md.ibd<-gl.ibd(mden43.pop.gl)
md.atl.ibd<-gl.ibd(mden43.atl.gl)
md.pac.ibd<-gl.ibd(mden43.pac.gl)


zczc123.atl.nopop.gl<-zczc123.atl.gl
pop(zczc123.atl.nopop.gl)<-c(1:54)
zczc123.pac.nopop.gl<-zczc123.pac.gl
pop(zczc123.pac.nopop.gl)<-c(1:36)
zczc123.med.nopop.gl<-zczc123.med.gl
pop(zczc123.med.nopop.gl)<-c(1:33)

atl.ibd<-gl.ibd(zczc123.atl.nopop.gl, permutations=99)
atl.geno<-dist(as.matrix(zczc123.atl.gl))
atl.coord<-dist(zczc123.atl.gl@other$latlong)
atl.ibd<-mantel.randtest(atl.geno, atl.coord)
atl.ibd

library(MASS)
dens<-kde2d(atl.coord, atl.geno, n=300)
myPal<-colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(atl.coord, atl.geno, pch=20, cex=0.5)
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(atl.geno~atl.coord))
title("Isolation by distance plot")

min(atl.geno)

library(plotly)
zc.ibd.plot<-plot_ly(x=atl.coord, y=atl.geno)
zc.ibd.plot

zczc123.ibd$plot
plot_ly(x=zczc123.pop.ibd$Dgeo, y=zczc123.pop.ibd$Dgen, type = 'scatter', mode = 'markers',
        text=rownames(zczc123.pop.ibd$Dgen))

# IBD and LC distance -----------------------------------------------------


#libraries
library(marmap)
library(tidyr)
library(ggmap)
library(ggplot2)
library(readxl)
library(fossil)

#starting coordinates
zczc123.coord
zczc125.coord<-zczc125_coord[,2:4]
zczc125.coord<-as.data.frame(zczc125.coord)
rownames(zczc125.coord)<-zczc125.coord[,1]
zczc125.coord<-zczc125.coord[,2:3]
head(zczc125.coord)

#remove these individuals
#Zcav20181_L3-30, row 10
#Zcav20181_L4-40, row 62

zczc123.coord<-zczc125.coord[c(1:9, 11:61,63:125),]
zczc123.coord
zczc123.coord2<-zczc123.coord[,c(2,1)]

#import blainville's coordinates
mden43.coord<-read.csv("./SampleLocationsPops.csv")
mden43.coord<-subset(mden43.coord, mden43.coord$species=="md")
mden43.coord2<-mden43.coord[,c(3,6,7)]
rownames(mden43.coord2)<-mden43.coord2[,1]
mden43.coord3<-mden43.coord2[,2:3]
mden43.coord4<-mden43.coord2[,c(3,2)]

#map data
wmap<-getNOAA.bathy(lon1 = -180, lon2 = 180, lat1 = -90, lat2 = 90, resolution = 10)
summary(wmap)
plot(wmap, image=TRUE, deep=-6000, shallow=0, step=1000)

wmap2<-getNOAA.bathy(lon1 = -180, lon2 = 180, lat1 = -90, lat2 = 90, resolution = 4)
summary(wmap2)
plot(wmap2, image=TRUE, deep=-6000, shallow=0, step=1000)
3
#centered on the antimeridian (middle of the pacific). ranges from 0-360
wmap3<-getNOAA.bathy(lon1 = 180, lon2 = -180, lat1 = -90, lat2 = 90, antimeridian = TRUE, keep=TRUE)
plot(wmap3, image=TRUE)

#check if coordinates are at sea, and adjust to nearest isobath
zc.dist.iso<-dist2isobath(wmap2, zczc123.coord2, isobath=-200)
zc.dist.iso
zc.dist.iso2<-dist2isobath(wmap3, zczc123.coord2, isobath=-200)
zc.dist.iso2
zc.dist<-zc.dist.iso[,4:5]
zc.dist
zc.dist2<-zc.dist.iso2[,4:5]
zc.dist2
zc.depth<-get.depth(wmap2, zc.dist, locator=FALSE)
zc.depth
zc.depth2<-get.depth(wmap3, zc.dist2, locator=FALSE)
zc.depth2

md.depth<-get.depth(wmap, mden43.coord4, locator=FALSE)
md.depth
md.depth2<-get.depth(wmap, md.dist, locator=FALSE)
md.depth2
md.dist.iso<-dist2isobath(wmap, mden43.coord4, isobath=-200)
md.dist<-md.dist.iso[,c(4,5)]

#trans1 <- trans.mat(NZMap)
trans <- trans.mat(wmap, min.depth =-10) #this ensured the path needed to have a minimum depth of 50m - if not put in place like 'trans1' then path cuts across land
trans2<-trans.mat(wmap, min.depth =-200)
trans3<-trans.mat(wmap, min.depth=-200)


#calculate the least cost distance between two locations through water, avoiding land. value in kilometers. 
out<-lc.dist(trans, zczc123.coord2, res="dist")
out2<-lc.dist(trans, zczc123.coord3, res="dist")
out3<-lc.dist(trans, zczc123.coord3, res="path")
out4<-lc.dist(trans, zc.dist, res="dist")
out.md<-lc.dist(trans, md.dist, res="dist")
out3.md<-as.matrix(out.md)
out3.md

out4
out4.zc<-as.matrix(out4)

#final LC distance matrices
out2.md
out3.md 

out4
out4.zc

plot(wmap2, image=TRUE)
points(zc.dist, pch = 21, col = "red", bg = col2alpha("red", .9),
       cex = 1.2)

plot(wmap2, image=TRUE)
points(md.dist, pch = 21, col = "red", bg = col2alpha("red", .9),
       cex = 1.2)

#new geographic distance matrix:
zc.geo.dist<-as.dist(out4.zc)
md.geo.dist<-as.dist(out3.md)
#genetic distance matrix:
zc.gen.dist<-dist(as.matrix(zczc123.gl))
md.gen.dist<-dist(as.matrix(mden43.gl))

#mantel test all ZC:
zczc123.mantel<-mantel.randtest(zc.geo.dist, zc.gen.dist)
zczc123.mantel
library(MASS)
dens<-kde2d(zc.geo.dist, zc.gen.dist, n=300)
myPal<-colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(zc.geo.dist, zc.gen.dist, pch=20, cex=0.5, xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(zc.gen.dist~zc.geo.dist))
title("Isolation by distance plot, All Cuvier's")
plot(zczc123.mantel)

#mantel test all MD:
mden43.mantel<-mantel.randtest(md.gen.dist, md.geo.dist)
mden43.mantel
library(MASS)
dens<-kde2d(md.geo.dist, md.gen.dist, n=300)
myPal<-colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(md.geo.dist, md.gen.dist, pch=20, cex=0.5, xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(md.gen.dist~md.geo.dist))
title("Isolation by distance plot, All Blainville's")

#mantel test,  within ocean basins
rownames(as.matrix(zc.gen.dist)) #both in the same order as the indNames in the genlight
rownames(as.matrix(zc.geo.dist))

zc.geo2.dist<-as.matrix(zc.geo.dist)
rownames(zc.geo2.dist)<-pop(zczc123.gl)
colnames(zc.geo2.dist)<-pop(zczc123.gl)
zc.atl.geo<-zc.geo2.dist[grepl("Atlantic", rownames(zc.geo2.dist)), grepl("Atlantic", colnames(zc.geo2.dist))]
zc.atl.geo
zc.pac.geo<-zc.geo2.dist[grepl("Indo-Pacific", rownames(zc.geo2.dist)), grepl("Indo-Pacific", colnames(zc.geo2.dist))]
zc.pac.geo
zc.med.geo<-zc.geo2.dist[grepl("Mediterranean", rownames(zc.geo2.dist)), grepl("Mediterranean", colnames(zc.geo2.dist))]
zc.med.geo
zc.nomed.geo<-zc.geo2.dist[grepl("Atlantic|Indo-Pacific", rownames(zc.geo2.dist)), grepl("Atlantic|Indo-Pacific", colnames(zc.geo2.dist))]
zc.nomed.geo
  
zc.gen2.dist<-as.matrix(zc.gen.dist)
rownames(zc.gen2.dist)<-pop(zczc123.gl)
colnames(zc.gen2.dist)<-pop(zczc123.gl)
zc.atl.gen<-zc.gen2.dist[grepl("Atlantic", rownames(zc.gen2.dist)), grepl("Atlantic", colnames(zc.gen2.dist))]
zc.atl.gen
zc.pac.gen<-zc.gen2.dist[grepl("Indo-Pacific", rownames(zc.gen2.dist)), grepl("Indo-Pacific", colnames(zc.gen2.dist))]
zc.pac.gen
zc.med.gen<-zc.gen2.dist[grepl("Mediterranean", rownames(zc.gen2.dist)), grepl("Mediterranean", colnames(zc.gen2.dist))]
zc.med.gen
zc.nomed.gen<-zc.gen2.dist[grepl("Atlantic|Indo-Pacific", rownames(zc.gen2.dist)), grepl("Atlantic|Indo-Pacific", colnames(zc.gen2.dist))]
zc.nomed.gen

#Atlantic only, ZC
zc.atl.geo<-as.dist(zc.atl.geo)
zc.atl.gen<-as.dist(zc.atl.gen)

zczc123.atl.mantel<-mantel.randtest(zc.atl.geo, zc.atl.gen)
zczc123.atl.mantel
library(MASS)
dens<-kde2d(zc.atl.geo, zc.atl.gen, n=300)
myPal<-colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(zc.atl.geo, zc.atl.gen, pch=20, cex=0.5,xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(zc.atl.gen~zc.atl.geo))
title("Isolation by distance plot, Atlantic Cuvier's")

#mantel test, Pacific ZC
zc.pac.geo<-as.dist(zc.pac.geo)
zc.pac.gen<-as.dist(zc.pac.gen)

zczc123.pac.mantel<-mantel.randtest(zc.pac.geo, zc.pac.gen)
zczc123.pac.mantel
library(MASS)
dens<-kde2d(zc.pac.geo, zc.pac.gen, n=300)
myPal<-colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(zc.pac.geo, zc.pac.gen, pch=20, cex=0.5,xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(zc.pac.gen~zc.pac.geo))
title("Isolation by distance plot, Indo-Pacific Cuvier's")

#IBD Med, ZC
zc.med.geo<-as.dist(zc.med.geo)
zc.med.gen<-as.dist(zc.med.gen)

zczc123.med.mantel<-mantel.randtest(zc.med.geo, zc.med.gen)
zczc123.med.mantel
library(MASS)
dens<-kde2d(zc.med.geo, zc.med.gen, n=300)
myPal<-colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(zc.med.geo, zc.med.gen, pch=20, cex=0.5,xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(zc.med.gen~zc.med.geo))
title("Isolation by distance plot, Mediterranean Cuvier's")

#IBD, All cuviers but mediterranean
zc.nomed.geo<-as.dist(zc.nomed.geo)
zc.nomed.gen<-as.dist(zc.nomed.gen)

zczc123.nomed.mantel<-mantel.randtest(zc.nomed.geo, zc.nomed.gen)
zczc123.nomed.mantel
library(MASS)
dens<-kde2d(zc.nomed.geo, zc.nomed.gen, n=300)
myPal<-colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(zc.nomed.geo, zc.nomed.gen, pch=20, cex=0.5,xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(zc.nomed.gen~zc.nomed.geo))
title("Isolation by distance plot, all Atlantic and Indo-Pacific Cuvier's")

#mantel test, Atl MD:

md.geo.dist
md.gen.dist

md.geo.dist2<-as.matrix(md.geo.dist)
rownames(md.geo.dist2)<-pop(mden43.gl)
colnames(md.geo.dist2)<-pop(mden43.gl)
mdatlgeo<-md.geo.dist2[grepl("Atlantic", rownames(md.geo.dist2)), grepl("Atlantic", colnames(md.geo.dist2))]
mdatlgeo

md.gen.dist2<-as.matrix(md.gen.dist)
rownames(md.gen.dist2)<-pop(mden43.gl)
colnames(md.gen.dist2)<-pop(mden43.gl)
mdatlgen<-md.gen.dist2[grepl("Atlantic", rownames(md.gen.dist2)), grepl("Atlantic", colnames(md.gen.dist2))]
mdatlgen

mdatlgeo2<-as.dist(mdatlgeo)
mdatlgen2<-as.dist(mdatlgen)

mden43.atl.mantel<-mantel.randtest(mdatlgen2, mdatlgeo2)
mden43.atl.mantel
library(MASS)
dens<-kde2d(mdatlgeo2, mdatlgen2, n=300)
myPal<-colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(mdatlgeo2, mdatlgen2, pch=20, cex=0.5, xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(mdatlgen2~mdatlgeo2))
title("Isolation by distance plot, Atlantic Blainville's")

#mantel test, Pacific MD:
md.geo.dist
md.gen.dist

mdpacgeo<-md.geo.dist2[grepl("Indopacific", rownames(md.geo.dist2)), grepl("Indopacific", colnames(md.geo.dist2))]
mdpacgeo

md.gen.dist2<-as.matrix(md.gen.dist)
rownames(md.gen.dist2)<-pop(mden43.gl)
colnames(md.gen.dist2)<-pop(mden43.gl)
mdpacgen<-md.gen.dist2[grepl("Indopacific", rownames(md.gen.dist2)), grepl("Indopacific", colnames(md.gen.dist2))]
mdpacgen

mdpacgeo2<-as.dist(mdpacgeo)
mdpacgen2<-as.dist(mdpacgen)

mden43.pac.mantel<-mantel.randtest(mdpacgen2, mdpacgeo2)
mden43.pac.mantel
library(MASS)
dens<-kde2d(mdpacgeo2, mdpacgen2, n=300)
myPal<-colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(mdpacgeo2, mdpacgen2, pch=20, cex=0.5, xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(mdpacgen2~mdpacgeo2))
title("Isolation by distance plot, Indo-Pacific Blainville's")

# DAPC --------------------------------------------------------------------
x<-zczc123.pop.gl
mat<-tab(x)
grp<-pop(x)
xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
xval
scatter(xval$DAPC)
assignplot(xval$DAPC) #shows both the membership probability (red rectangle) and the pre-assigned group membership (blue cross)

zczc123.nomed.pop.gl<-popsub(zczc123.pop.gl, exclude=c("Med_East", "Med_West"))
pop(zczc123.nomed.pop.gl)
mat<-tab(zczc123.nomed.pop.gl)
grp<-pop(zczc123.nomed.pop.gl)
nomed.val<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
nomed.val
scatter(nomed.val$DAPC)
assignplot(nomed.val$DAPC)

zczc123.nomedorspain.gl<-popsub(zczc123.pop.gl, exclude=c("Med_East", "Med_West", "Atl_Spain"))
mat<-tab(zczc123.nomedorspain.gl)
grp<-pop(zczc123.nomedorspain.gl)
nomedorspain.val<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
nomedorspain.val
scatter(nomedorspain.val$DAPC)
assignplot(nomedorspain.val$DAPC)

zczc121.pop.gl<-popsub(zczc123.pop.gl, exclude="Atl_Spain")
z<-zczc121.pop.gl
mat<-tab(z)
grp<-pop(z)
zval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
zval
scatter(zval$DAPC)
assignplot(zval$DAPC)

par(mar=c(5.1,1, 4.1, 2.1))
zczc118.pop.gl<-popsub(zczc123.pop.gl, exclude=c("Atl_Spain", "Indo_Mix"))
m<-zczc118.pop.gl
mat<-tab(m)
grp<-pop(m)
mval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
mval
scatter(mval$DAPC)
scatter(mval$DAPC, xax = 2, yax=3)
scatter(mval$DAPC, xax=1, yax=2)
scatter(mval$DAPC, xax=3, yax=4)
assignplot(mval$DAPC)
heatmap(mval$DAPC$posterior)
heatmap.bp(mval$DAPC$posterior)
mval$DAPC$posterior


a<-zczc123.atl.gl
mat<-tab(a)
grp<-pop(a)
aval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
aval
scatter(aval$DAPC)
assignplot(aval$DAPC)

p<-zczc123.pac.gl
mat<-tab(p)
grp<-pop(p)
pval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
pval
scatter(pval$DAPC)
assignplot(pval$DAPC)

m<-zczc123.med.gl
mat<-tab(m)
grp<-pop(m)
mval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
mval
scatter(mval$DAPC)
assignplot(mval$DAPC)

mden43.pop.gl
y<-mden43.pop.gl
mat<-tab(y)
grp<-pop(y)
yval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
yval
scatter(yval$DAPC)
assignplot(yval$DAPC)

mden43.atl.gi<-popsub(mden43.pop.gi, exclude=c("Indo-Africa", "Indo-Hawaii", "Indo-South", "NA"))
mden43.atl.gl<-gi2gl(mden43.atl.gi)
am<-mden43.atl.gl
mat<-tab(am)
grp<-pop(am)
amval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
amval
scatter(amval$DAPC, posi.da="topleft")
assignplot(amval$DAPC)

mden43.pac.gi<-popsub(mden43.pop.gi, exclude=c("Atl-Bahamas", "Atl-East", "Atl-Other"))
mden43.pac.gl<-gi2gl(mden43.pac.gi)
pm<-mden43.pac.gl
mat<-tab(pm)
grp<-pop(pm)
pmval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
pmval
scatter(pmval$DAPC, posi.da="topleft")
assignplot(pmval$DAPC)

par(mar = c(2, 2, 2, 2))
scatter (mval$DAPC)
scatter(yval$DAPC, xax = 2, yax = 3, xax="1")

scatter.dapc(yval$DAPC, var.contrib=TRUE)
summary.dapc(yval$DAPC)
summary.dapc(yval$DAPC$var)
scatter(yval$DAPC )

##DAPC using find.clusters()
zczc123.gl
zczc123.atl.gi
zczc123.med.gi
zczc123.pac.gi

levels(zczc123.pop.gl$pop)
zczc123.atl.nosp.gl<-popsub(zczc123.pop.gl, exclude=c("Atl_Spain","Indo_Cent","Indo_Mix", "Indo_NE","Indo_Sou", "Med_East","Med_West"))


mden43.gl
mden43.atl.gl
mden43.pac.gl

zczc123.grp<-find.clusters(zczc123.gl, max.n.clust=20) #2 clusters, clear elbow
zczc123.atl.grp<-find.clusters(zczc123.atl.gi, max.n.clust=20) #2 clusters, no dip, near linear increase
zczc123.atl.nosp.grp<-find.clusters(zczc123.atl.nosp.gl, max.n.clust=20) #2 clusters, no dip, linear increase
zczc123.med.grp<-find.clusters(zczc123.med.gi, max.n.clust=20) #2 clusters, clear elbow
zczc123.pac.grp<-find.clusters(zczc123.pac.gi, max.n.clust=20) #5 clusters, no dip, linear increase

mden43.grp<-find.clusters(mden43.gl, max.n.clust=20) #2 clusters, clear elbow
mden43.atl.grp<-find.clusters(mden43.atl.gl, max.n.clust=20) #2 clusters, slight decrease in slope
mden43.pac.grp<-find.clusters(mden43.pac.gl, max.n.clust=14) #9 clusters, no dip, increases until 9 and then falls
mden43.pac2.grp<-find.clusters(mden43.pac.gl, max.n.clust=14) #3 clusters

zczc123.dapc<-xvalDapc(tab(zczc123.gl), zczc123.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(zczc123.dapc$DAPC)
zczc123.dapc$DAPC$grp

zczc123.atl.dapc<-xvalDapc(tab(zczc123.atl.gl), zczc123.atl.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(zczc123.atl.dapc$DAPC)
zczc123.atl.dapc$DAPC$grp

zczc123.atl.nosp.dapc<-xvalDapc(tab(zczc123.atl.nosp.gl), zczc123.atl.nosp.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(zczc123.atl.nosp.dapc$DAPC)
zczc123.atl.nosp.dapc$DAPC$grp

zczc123.med.gl<-gi2gl(zczc123.med.gi)
zczc123.med.dapc<-xvalDapc(tab(zczc123.med.gl), zczc123.med.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(zczc123.med.dapc$DAPC)
zczc123.med.dapc$DAPC$grp

zczc123.pac.gl<-gi2gl(zczc123.pac.gi)
zczc123.pac.dapc<-xvalDapc(tab(zczc123.pac.gl), zczc123.pac.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(zczc123.pac.dapc$DAPC, posi.da = "bottomleft")
scatter(zczc123.pac.dapc$DAPC, xax = 2, yax = 3)
zczc123.pac.dapc$DAPC$grp

mden43.gl
mden43.dapc<-xvalDapc(tab(mden43.gl), mden43.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(mden43.dapc$DAPC)

mden43.atl.dapc<-xvalDapc(tab(mden43.atl.gl), mden43.atl.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(mden43.atl.dapc$DAPC)

mden43.pac.dapc<-xvalDapc(tab(mden43.pac.gl), mden43.pac.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(mden43.pac.dapc$DAPC, posi.da="topright")

mden43.pac2.dapc<-xvalDapc(tab(mden43.pac.gl), mden43.pac2.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(mden43.pac2.dapc$DAPC, posi.da="bottomleft")

zczc123.dapc$DAPC$grp
zczc123.atl.dapc$DAPC$grp
zczc123.atl.nosp.dapc$DAPC$grp
zczc123.med.dapc$DAPC$grp
zczc123.pac.dapc$DAPC$grp

mden43.dapc$DAPC$grp
mden43.atl.dapc$DAPC$grp
mden43.pac.dapc$DAPC$grp
mden43.pac2.dapc$DAPC$grp


plot(zczc123.grp$Kstat, xlab="No. Clusters", ylab="BIC Score", main="All Cuvier's")
plot(zczc123.atl.grp$Kstat, xlab="No. Clusters", ylab="BIC Score", main="Atlantic Cuvier's, All")
plot(zczc123.atl.nosp.grp$Kstat, xlab="No. Clusters", ylab="BIC Score", main="Atlantic Cuvier's, No Spain")
plot(zczc123.med.grp$Kstat, xlab="No. Clusters", ylab="BIC Score", main="Mediterranean Cuvier's")
plot(zczc123.pac.grp$Kstat, xlab="No. Clusters", ylab="BIC Score", main="Indo-Pacific Cuvier's")
plot(mden43.grp$Kstat, xlab="No. Clusters", ylab="BIC Score", main="All Blainville's")
plot(mden43.atl.grp$Kstat, xlab="No. Clusters", ylab="BIC Score", main="Atlantic Blainville's")
plot(mden43.pac.grp$Kstat, xlab="No. Clusters", ylab="BIC Score", main="Indo-Pacific Blainville's")

# BioNJ Trees with SRW --------------------------------
gl.install.vanilla.dartR()
BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)

zc_srw<-read.vcfR("~/Dropbox/Phd/Bioinformatics/bw_ddrad_global/global_paper/NEW_STUFF_OCT_2021/populations.snps.vcf")
zc_srw.gl<-vcfR2genlight(zc_srw)
zc_srw.gl<-gl.filter.monomorphs(zc_srw.gl)
gl.compliance.check(zc_srw.gl)

pop(zc_srw.gl)<-as.factor(c("Mediterranean",
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
                            "Indo-Pacific",
                            "SRW",
                            "SRW",
                            "SRW",
                            "SRW",
                            "SRW",
                            "SRW"))
zc_srw.pop.gl<-zc_srw.gl
pop(zc_srw.pop.gl)<-as.factor(c("Med_West",
                                "Med_West",
                                "Med_West",
                                "Med_West",
                                "Med_West",
                                "Med_West",
                                "Med_East",
                                "Med_West",
                                "Med_West",
                                "Med_East",
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
                                "Med_East",
                                "Med_East",
                                "Med_East",
                                "Atl_CanIs",
                                "Atl_CanIs",
                                "Atl_CanIs",
                                "Atl_CanIs",
                                "Atl_CanIs",
                                "Atl_France",
                                "Atl_CanIs",
                                "Atl_CanIs",
                                "Atl_France",
                                "Atl_France",
                                "Atl_NE",
                                "Atl_NE",
                                "Atl_NE",
                                "Atl_CanIs",
                                "Atl_NE",
                                "Atl_Carib",
                                "Atl_France",
                                "Atl_CanIs",
                                "Atl_Spain",
                                "Atl_France",
                                "Atl_Carib",
                                "Atl_Spain",
                                "Atl_CanIs",
                                "Atl_CanIs",
                                "Atl_NE",
                                "Atl_NE",
                                "Atl_NE",
                                "Atl_CanIs",
                                "Atl_Carib",
                                "Atl_CanIs",
                                "Indo_NE",
                                "Indo_NE",
                                "Atl_Carib",
                                "Med_East",
                                "Indo_Mix",
                                "Atl_Carib",
                                "Indo_NE",
                                "Atl_Carib",
                                "Atl_NE",
                                "Atl_CanIs",
                                "Med_West",
                                "Med_West",
                                "Indo_Sou",
                                "Med_West",
                                "Med_East",
                                "Med_West",
                                "Atl_CanIs",
                                "Atl_Carib",
                                "Atl_Carib",
                                "Atl_NE",
                                "Atl_NE",
                                "Atl_NE",
                                "Atl_CanIs",
                                "Atl_Carib",
                                "Atl_Carib",
                                "Atl_Carib",
                                "Atl_Carib",
                                "Indo_NE",
                                "Indo_NE",
                                "Med_West",
                                "Med_West",
                                "Indo_NE",
                                "Indo_NE",
                                "Atl_Carib",
                                "Atl_Carib",
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
                                "Indo_Mix",
                                "Indo_Mix",
                                "Indo_Sou",
                                "Indo_Sou",
                                "Indo_Sou",
                                "Med_West",
                                "Indo_NE",
                                "Indo_Sou",
                                "Indo_Sou",
                                "Indo_Sou",
                                "SRW",
                                "SRW",
                                "SRW",
                                "SRW",
                                "SRW",
                                "SRW"))

zc_srw.pop.gl
zc_srw.pop.gi<-gl2gi(zc_srw.pop.gl)

d<-as.data.frame(cbind(indNames(zc_srw.pop.gi), as.character(pop(zc_srw.pop.gi))))
colnames(d)<-c("label", "pop")
d

e<-as.data.frame(cbind(indNames(zc_srw.pop.gi), as.character(pop(zc_srw.gl))))
colnames(e)<-c("label", "ocean")
e


zc_srw.gi<-gl2gi(zc_srw.gl)

zc_srw.aboot<-aboot(zc_srw.gi, tree=bionj, tip.label=NULL, cutoff=50, sample=1000)
zc_srw2.aboot<-aboot(zc_srw.gi, tree=bionj, tip.label=NULL, sample=1000)
zc_srw.aboot
zc_srw2.aboot
# zc_srw.aboot2<-full_join(zc_srw.aboot, e, by='label' )
# zc_srw2.aboot2<-full_join(zc_srw2.aboot, e, by='label' )
# zc_srw.aboot2
# zc_srw2.aboot2
# ggtree(zc_srw.aboot2) + geom_tippoint(aes(colour=ocean))

zc_srw.root.aboot<-root(zc_srw.aboot, outgroup="NZ_Eau08AI032", edgelabel=TRUE)
zc_srw2.root.aboot<-root(zc_srw2.aboot, outgroup="NZ_Eau08AI032", edgelabel=TRUE)
class(zc_srw.root.aboot)
class(zc_srw2.root.aboot)
zc_srw.root.aboot<-full_join(zc_srw.root.aboot, e, by='label' )
zc_srw.root.aboot<-full_join(zc_srw.root.aboot, d, by="label")
ggtree(zc_srw.root.aboot) + geom_tippoint(aes(colour=ocean))
zc_srw2.root.aboot<-full_join(zc_srw2.root.aboot, e, by='label' )
zc_srw2.root.aboot<-full_join(zc_srw2.root.aboot, d, by="label")
ggtree(zc_srw2.root.aboot) + geom_tippoint(aes(colour=ocean))


zcsrw.root.aboot2<-drop.tip(zc_srw.root.aboot@phylo, c("NZ_Eau07AI006","NZ_Eau07AI008","NZ_Eau08AI006","NZ_Eau08AI013","NZ_Eau08AI026","NZ_Eau08AI032"))
zcsrw.root.aboot2
zcsrw2.root.aboot2<-drop.tip(zc_srw2.root.aboot@phylo, c("NZ_Eau07AI006","NZ_Eau07AI008","NZ_Eau08AI006","NZ_Eau08AI013","NZ_Eau08AI026","NZ_Eau08AI032"))
zcsrw2.root.aboot2
d2<-d[1:118,]
e2<-e[1:118,]
zcsrw.root.aboot2<-full_join(zcsrw.root.aboot2, e2, by='label' )
zcsrw.root.aboot2<-full_join(zcsrw.root.aboot2, d2, by="label")
ggtree(zcsrw.root.aboot2) +geom_tippoint(aes(colour=ocean))
ggtree(zcsrw.root.aboot2) +geom_tippoint(aes(colour=pop))
ggtree(zcsrw.root.aboot2) +geom_tiplab()
zcsrw2.root.aboot2<-full_join(zcsrw2.root.aboot2, e2, by='label' )
zcsrw2.root.aboot2<-full_join(zcsrw2.root.aboot2, d2, by="label")
ggtree(zcsrw2.root.aboot2) +geom_tippoint(aes(colour=ocean))
ggtree(zcsrw2.root.aboot2) +geom_tippoint(aes(colour=pop))
ggtree(zcsrw2.root.aboot2) +geom_tiplab()

##Final Cuvier's Tree for paper
png("./zczc.boot.phylo.png", width=4, height=9, units="in", res=1000)
p<-ggtree(zcsrw.root.aboot2) +
  geom_tippoint(aes(col=factor(ocean)), size=2, show.legend=NULL) +
  scale_colour_manual(values=c("#7b848f", "#c85200", "#5fa2ce")) +
  geom_label2(aes(label=label, subset=as.numeric(label)>50), fill='white', size=4, label.padding = unit(0.1, "lines"), label.r = unit(0, "lines"), nudge_x=-.001, show.legend=NULL) +
  theme_tree2() +
  theme(legend.position=c(.22,.15), legend.title = element_blank())
p
dev.off()

ggtree(zc_srw.root.aboot) + geom_text(aes(label=node))

png("./zczc.boot.phylo.srw.png", width=5, height=5, units="in", res=1000)
p<-ggtree(zc_srw.root.aboot) +
  geom_tippoint(aes(col=factor(ocean)), size=1.5) +
  geom_label2(aes(label=label, subset=as.numeric(label)>50), fill='white', size=2, label.padding = unit(0.1, "lines"), label.r = unit(0, "lines"), nudge_x=-.001, show.legend=NULL) +
  theme_tree2() +
  theme(legend.title = element_blank())
p2<-gridExtra::grid.arrange(ggtree::rotate(p, 144))

dev.off()

#generate bionj tree for blainvilles
md_srw.vcfR<-read.vcfR("~/Dropbox/Phd/Bioinformatics/bw_ddrad_global/SNP_files/no_populations/one_snp_per_locus/md+srw/populations.snps.vcf")
md_srw.gl<-vcfR2genlight(md_srw.vcfR)
md_srw.gi<-gl2gi(md_srw.gl)
md_srw.aboot<-aboot(md_srw.gi, tree=bionj, tip.label=NULL, cutoff=50, sample=1000)
md_srw.aboot2<-aboot(md_srw.gi, tree=bionj, tip.label=NULL, sample=1000) #re-running bootstrap with no cut-off at 50%
mdsrw.root<-root(md_srw.aboot, outgroup="NZ_Eau08AI032", edgelabel=TRUE)
mdsrw.root.pop<-full_join(mdsrw.root, md.meta, by="label")
mdsrw.root2<-drop.tip(mdsrw.root, c("NZ_Eau07AI006","NZ_Eau07AI008","NZ_Eau08AI006","NZ_Eau08AI013","NZ_Eau08AI026","NZ_Eau08AI032"))

#re-running with new bootstrap results that aren't cut off at 50%
mdsrw2.root<-root(md_srw.aboot2, outgroup="NZ_Eau08AI032", edgelabel=TRUE)
mdsrw2.root<-full_join(mdsrw2.root, md.meta, by="label")
mdsrw2.root2<-drop.tip(mdsrw2.root, c("NZ_Eau07AI006","NZ_Eau07AI008","NZ_Eau08AI006","NZ_Eau08AI013","NZ_Eau08AI026","NZ_Eau08AI032"))

md.meta<-data.frame(e,f,g)
colnames(md.meta)<-c("label", "ocean", "pop")
md.meta<-md.meta[1:42,]

mdsrw.root3<-full_join(mdsrw.root2, md.meta, by='label' )
mdsrw.root3
class(mdsrw.root3)
ggtree(mdsrw.root3)

mdsrw2.root2
class(mdsrw2.root2)
ggtree(mdsrw2.root2)

mdsrw2.root2@phylo$node.label<-as.numeric(mdsrw2.root2@phylo$node.label)


###Blainville's Tree from Paper
png("./md.boot.phylo.png", width=5, height=9, units="in", res=1000)
p<-ggtree(mdsrw.root3) +
  geom_tippoint(aes(col=factor(ocean)), size=4, show.legend=NULL) +
  scale_colour_manual(values=c("#7b848f", "#c85200")) +
  geom_label2(aes(label=label, subset=as.numeric(label)>50), fill='white', size=5, label.padding = unit(0.1, "lines"), label.r = unit(0, "lines"), nudge_x=-.001, show.legend=NULL) +
  theme_tree2() +
  theme(legend.position=c(.18,.2), legend.title = element_blank())
p
dev.off()



#This tree has the bootstrap values offset so they don't overlap
a<-ggtree(mdsrw.root3) +
  geom_tippoint(aes(col=factor(ocean)), size=4, show.legend=NULL) +
  scale_colour_manual(values=c("#7b848f", "#c85200")) +
  geom_label_repel(aes(label=as.numeric(label))) +
  theme_tree2() +
  theme(legend.position=c(.18,.2), legend.title = element_blank())

b<-ggtree(zcsrw.root.aboot2) +
  geom_tippoint(aes(col=factor(ocean)), size=2, show.legend=NULL) +
  scale_colour_manual(values=c("#7b848f", "#c85200", "#5fa2ce")) +
  geom_label_repel(aes(label=as.numeric(label))) +
  theme_tree2() +
  theme(legend.position=c(.22,.15), legend.title = element_blank())
  
#This tree has the bootstrap values converted into coloured dots on a scale
png("~/Dropbox/Phd/Bioinformatics/bw_ddrad_global/Results/trees/Oct2022/md.boot.phylo.png", width=5, height=9, units="in", res=1000)
ggtree(mdsrw2.root2) +
  geom_tippoint(aes(col=factor(ocean)), size=4, show.legend=FALSE) +
  scale_colour_manual(values=c("#7b848f", "#c85200")) +
  new_scale_color() +
  geom_nodepoint(aes(subset=!is.na(as.numeric(label)), colour=(as.numeric(label))), size=3) +
  scale_color_viridis_c(name="Bootstrap(%)") +
  theme_tree2(legend.position=c(0.20,0.2), legend.background = element_rect(fill='transparent'), legend.text=element_text(size=10), legend.key=element_blank()) 
dev.off()

png("~/Dropbox/Phd/Bioinformatics/bw_ddrad_global/Results/trees/Oct2022/zc.boot.phylo.png", width=5, height=11, units="in", res=1000)
ggtree(zcsrw2.root.aboot2) +
  geom_tippoint(aes(col=factor(ocean)), size=4, show.legend=FALSE) +
  scale_colour_manual(values=c("#7b848f", "#c85200", "#5fa2ce")) +
  new_scale_color() +
  geom_nodepoint(aes(subset=!is.na(as.numeric(label)), colour=(as.numeric(label))), size=3) +
  scale_color_viridis_c(name="Bootstrap(%)") +
  theme_tree2(legend.position=c(0.22, 0.15), legend.background = element_rect(fill='transparent'), legend.text=element_text(size=10), legend.key=element_blank())
dev.off()

ggarrange(a, b, c, d, ncol=2, nrow=2)


# png("./md.boot.phylo.srw.png", width=5, height=5, units="in", res=1000)
# p<-ggtree(mdsrw.root.pop) +
#   geom_tippoint(aes(col=factor(ocean)), size=1.5)  +
#   geom_label2(aes(label=label, subset=as.numeric(label)>50), fill='white', size=2, label.padding = unit(0.1, "lines"), label.r = unit(0, "lines"), nudge_x=-.001, show.legend=NULL) +
#   theme_tree2() +
#   theme(legend.title = element_blank())
# p
# dev.off()


# IQtree ------------------------------------------------------------------
setwd("~/Dropbox/Phd/Bioinformatics/bw_ddrad_global/global_paper/NEW_STUFF_OCT_2021/IQtree")
library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")

#read in treefile output from IQtree into R
md<-read.newick("md_srw/tree_finder/populations.fixed.phylip.treefile")
ggtree(md)
md$tip.label

md2<-read.newick("md_srw/tree_finder/populations.fixed.phylip.contree")
ggtree(md2)

#make some tables to add in extra information, including sample ID and population
d<-c("1",	"2",	"3",	"4",	"5",	"6",	"7",	"8",	"9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"20",	"21",	"22",	"23",	"24",	"25",	"26",	"27",	"28",	"29",	"30",	"31",	"32",	"33",	"34",	"35",	"36",	"37",	"38",	"39",	"40",	"41",	"42",	"43",	"44",	"45",	"46",	"47",	"48")
e<-c("BWLib_L4-1",	"BWLib_L4-2",	"BWLib_L4-4",	"BWLib_L4-5",	"BWLib_L4-6",	"BWLib_L4-7",	"BWLib_L4-8",	"Mde_L1-1",	"Mde_L1-3",	"Mde_L1-6",	"Mde_L1-7",	"Mde_L1-8",	"Mde_L2-10",	"Mde_L2-12",	"Mde_L2-13",	"Mde_L2-14",	"Mde_L3-19",	"Mde_L3-23",	"Mde_L3-24",	"Mde_L4-28",	"Mde_L4-30",	"Mde_L4-31",	"Mde_L5-33",	"Mde_L5-37",	"Mde_L5-38",	"Mde_L5-39",	"Mde_L5-40",	"Mde_L1-4",	"Mde_L2-9",	"Mde_L4-27",	"Mde_L4-29",	"Mde_L4-32",	"Mde_L5-34",	"Mde_L1-5",	"Mde_L2-15",	"Mde_L3-22",	"Mde_L4-26",	"Mde_L5-36",	"Mde_L1-2",	"Mde_L2-11",	"Mde_L3-21",	"Mde_L4-25",	"NZ_Eau07AI006",	"NZ_Eau07AI008",	"NZ_Eau08AI006",	"NZ_Eau08AI013",	"NZ_Eau08AI026",	"NZ_Eau08AI032")
f<-c("Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Atlantic",	"Indopacific",	"Indopacific",	"Indopacific",	"Indopacific",	"Indopacific",	"Indopacific",	"Indopacific",	"Indopacific",	"Indopacific",	"Indopacific",	"Indopacific",	"Indopacific",	"Indopacific",	"Indopacific",	"Indopacific",	"SRW",	"SRW",	"SRW",	"SRW",	"SRW",	"SRW")
g<-c("Atl_East",	"Atl_Oth",	"Atl_Bah",	"Atl_Bah",	"Atl_Bah",	"Atl_Bah",	"Atl_Bah",	"Atl_East",	"Atl_East",	"Atl_East",	"Atl_Oth",	"Atl_East",	"Atl_East",	"Atl_East",	"Atl_Oth",	"Atl_East",	"Atl_East",	"Atl_Bah",	"Atl_East",	"Atl_East",	"Atl_East",	"Atl_East",	"Atl_East",	"Atl_East",	"Atl_Bah",	"Atl_East",	"Atl_Oth",	"Indo_Haw",	"Indo_Haw",	"Indo_Haw",	"Indo_Haw",	"Indo_Haw",	"Indo_Haw",	"Indo_Afr",	"Indo_Afr",	"Indo_Afr",	"Indo_Afr",	"Indo_Afr",	"Indo_Sou",	"Indo_Sou",	"Indo_Sou",	"Indo_Sou",	"SRW",	"SRW",	"SRW",	"SRW",	"SRW",	"SRW")
cbind(e,g)

md.names<-as.data.frame(cbind(d, e))
colnames(md.names)<-c("label", "id")
md.names

md.ocean<-as.data.frame(cbind(d,f))
colnames(md.ocean)<-c("label", "ocean")
md.pop<-as.data.frame(cbind(d,g))
colnames(md.pop)<-c("label", "pop")

#plot tree using the sample IDs as tip labels
md<-full_join(md, md.names, by="label")
md

md@phylo$node.label
md.root<-root(md, outgroup="48", edgelabel=TRUE)
md.root
ggtree(md.root)
viewClade(md.root, "1", "42")

md.root2<-md.root
md.root2$edge.length
md.root2$edge.length<-(md.root2$edge.length*10)
ggtree(md.root2)

md.root2
md.root2<-full_join(md.root2, md.ocean[1:42,], by="label")
md.root3<-drop.tip(md.root2, c("43", "44", "45", "46", "47", "48"))
md.root4<-full_join(md.root3, md.names[1:42,], by="label")
md.root4<-full_join(md.root4, md.ocean[1:42,], by="label")
md.root4<-full_join(md.root4, md.pop[1:42,], by="label")
ggtree(md.root4)

ggtree(md.root4) + geom_tiplab(aes(label = id))
ggtree(md.root4) + geom_tiplab(aes(label = ocean))
ggtree(md.root4) + geom_tiplab(aes(label = pop))

ggtree(md.root2) +geom_tippoint(aes(colour=ocean))
ggtree(md.root4) +geom_tippoint(aes(colour=ocean))

#Zc tree
zc<-read.newick("zc_srw/zc4/populations.fixed.phylip.treefile")
ggtree(zc)
zc$tip.label

zc.pops<-read.csv("zc.csv")
zc.pops$label<-as.factor(zc.pops$label)
zc1<-full_join(zc, zc.pops, by="label")
ggtree(zc1) +geom_tippoint(aes(colour=ocean))

zc.root<-root(zc, outgroup="124", edgelabel=TRUE)
zc.root
ggtree(zc.root)
zc.root1b<-full_join(zc.root, zc.pops, by="label")
ggtree(zc.root1b) +geom_tippoint(aes(colour=ocean))

zc.root2<-drop.tip(zc.root, c("119", "120", "121", "122", "123", "124"))
ggtree(zc.root2)
zc.root3<-full_join(zc.root2, zc.pops[1:118,], by="label")
zc.root3
ggtree(zc.root3) +geom_tippoint(aes(colour=ocean))

zczc125.tree<-bionjs(dist(as.matrix(zczc125.gl)))
plot(zczc125.tree)
ggtree(zczc125.tree)



# Fis ---------------------------------------------------------------------
zczc123.basicstats
mden43.basicstats<-basic.stats(mden43.gi)

#How is fis distributed among snps? Are there snps that are consistently over multiple samples displaying particularly large heterozygote deficiencies?
zczc123.perloc.fis<-zczc123.basicstats$Fis
class(zczc123.perloc.fis)
dim(zczc123.perloc.fis)
zczc123.perloc.fis<-as.data.frame(zczc123.perloc.fis)
class(zczc123.perloc.fis)
dim(zczc123.perloc.fis)
hist.data.frame(zczc123.perloc.fis)
hist(zczc123.perloc.fis$Atlantic, col="red")
hist(zczc123.perloc.fis$`Indo-Pacific`, add=TRUE, col="blue")
hist(zczc123.perloc.fis$Mediterranean, add=TRUE, col="yellow")


fd<-gl.fixed.diff(zczc123.gl)
fd2<-gl.fixed.diff(mden43.gl)

zczc123.fis.gl<-filter_fis(zczc123.gi, approach="SNP", fis.min.threshold=-0.05, fis.max.threshold=0.05)



# dataset summary ---------------------------------------------------------
zczc123.gl
summary(NA.posi(zczc123.gl))
zczc123.gl
gl.report.rdepth(zczc123.gl)
zczc123.vcfR
gt <- extract.gt(zczc123.vcfR, element = 'DP', as.numeric = TRUE)
gt[1:5,1:5]
gt<-as.data.frame(gt)

gt %>% summarise_if(is.numeric, mean)
colMeans(gt[sapply(gt, is.numeric)], na.rm=TRUE) 

summary(NA.posi(mden43.gl))
gt <- extract.gt(mden43.vcfR, element = 'DP', as.numeric = TRUE)
gt<-as.data.frame(gt)
colMeans(gt[sapply(gt, is.numeric)], na.rm=TRUE) 
