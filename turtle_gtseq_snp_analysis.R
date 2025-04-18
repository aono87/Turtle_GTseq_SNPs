#Script to analyse vcf file of Leatherback turtle genotypes from Amy. 
#Mix of populations including W pacific, Mexico, Bycatch, St Croix and S. Africa
#Pretty basic population structure including a basic phylogenetic tree, DAPC, Fst and relatedness

#Required Libraries
library(adegenet)
library(vcfR)
library(remotes)
library(SNPRelate)
library(dartRverse)
library(poppr)
library(ape)
library(RColorBrewer)
library(dartR)
library(Demerelate)
library(tidyverse)
library(hierfstat)

##Import vcf into R for use in adegenet
dcor.vcfR<-read.vcfR(("./Dcor_DcPanel_205_maf.targetSNPs_012224.recode.vcf"), verbose = TRUE)
dcor.vcfR
head(getFIX(dcor.vcfR)) #View the fixed variables for the first six loci of the vcf file
is.biallelic(dcor.vcfR) #Tells you whether each locus is biallelic
getCHROM(dcor.vcfR) #give list of loci
getPOS(dcor.vcfR) #gives list of SNP positions within loci

#Prepare a file with the desired population info. 
#I like to make a table with the sample id in one column and another column with population ids.
#Import the population data into R
pop.data.dcor<-read.table("sample_pop.txt", sep="\t", header=TRUE)
pop.data.dcor
#Before you add the population data to the genotype data, the following code will let you know if the sample IDs of the population file are in the same order as the samples in the vcf file
all(colnames(dcor.vcfR@gt)[-1]==pop.data.dcor$ind) #check that sample IDs are the same and in the same order

#Before adding sample information, convert the vcfR file into a genlight file format
#this file format is extremely useful for downstreat QA/QC and analysis
#convert vcfR file to genlight files
dcor.gl<-vcfR2genlight(dcor.vcfR)
dcor.gl
#n=177 individuals, n=231 loci

#Now you can add the population info from the file
pop(dcor.gl)<-pop.data.dcor$NewLoc

#Check the features of the genlight file
dcor.gl$loc.names #unlike the vcfR file, the locus and positions are combined into one locus name
dcor.gl$ind.names #List of individuals in the genlight
dcor.gl$pop #population identity of each individual

##Remove duplicate individuals 
#use dartR to make a recode file from the individuals in the genlight to remove duplicates
#first sort genlight by individuals
dcor.gl<-gl.sort(dcor.gl, sort.by="ind")
gl.make.recode.ind(dcor.gl, out.recode.file="dcor.gl_recode_ind.csv", outpath="./")
#Open up csv file onto computer. 
#replace the sample name in the second column with "Delete" 
#save csv
#read back in the csv file and delete the chosen individuals
dcor.gl<-gl.recode.ind(dcor.gl, ind.recode="./dcor.gl_recode_ind.csv", recalc=TRUE, mono.rm=TRUE)
#check if individuals have been removed
dcor.gl

#If you want to export out a matrix of the sample ids in rows and snp loci in columns with genotypes
dcor.mat<-as.matrix(dcor.gl)
dcor.mat
write.csv(dcor.mat, file="./dcor.snp.genotype.matrix.csv")
#i used this file to make a list of SNP loci that I wanted to drop from the file. There are some loci with SNPs in multiple positions so for this first round of analysis I wanted to keep only a single SNP from each location. 

#remove desired loci from the file. 
#unlike individuals and populations, you can't rename or drop loci using a csv using dartR and must do it manually
dcor.gl<-gl.compliance.check(dcor.gl) #a funny feature of dartR package is it needs you to do these compliance checks. If things arent working, try this and it may help! 
dcor.snp.gl<-gl.drop.loc(dcor.gl, loc.list=c("locus022_63", "locus045_119", "locus053_99", "locus061_150", "locus100_115", "locus100_118", "locus129_60", "locus135_106", "locus188_138", "locus189_70", "locus237_111", "locus291_100", "locus334_109", "Dc01661_95", "Dc01707_114", "Dc03887_101", "Dc05864_137", "Dc07300_81", "Dc07937_134", "Dc08054_106", "Dc08189_82", "Dc08578_70", "Dc08865_124", "Dc10885_73", "Dc11062_41", "Dc11203_120", "Dc12329_111", "Dc13159_110", "Dc13592_122", "Dc16962_121", "Dc20232_74", "Dc21637_31", "Dc22013_123", "Dc22883_124", "Dc33714_113", "Dc36767_103", "Dc41234_117", "Dc46984_123", "Dc47572_144", "Dc47935_46", "Dc47935_56", "Dc48950_40", "Dc49184_94", "Dc52523_131", "Dc58248_125", "Dc62188_75"))
dcor.snp.gl
#n=148 individuals, #n=185 SNPs
#303 missing data, 1.11%
dcor.snp.gl<-gl.compliance.check(dcor.snp.gl)

#You can use dartR to check the call rate by locus or individual
dcor.snp2.gl<-gl.filter.callrate(dcor.snp.gl, method="ind", threshold=0.5, recalc=TRUE) #no samples with more than 50% missing data

#Create new subsets of the datasheet with different levels of missing data
dcor.snp100.gl<-gl.filter.callrate(dcor.snp.gl, method="loc", threshold=1, recalc=TRUE) #only kept loci with 0% missing data
dcor.snp95.gl<-gl.filter.callrate(dcor.snp.gl, method="loc", threshold=.95, recalc=TRUE) #only kept loci with <5% missing data
dcor.snp80.gl<-gl.filter.callrate(dcor.snp.gl, method="loc", threshold=.8, recalc=TRUE) #only kept loci with <20% missing data

dcor.snp100.gl #161 snps, 0% missing data
dcor.snp95.gl #179 snps, .16% missing data
dcor.snp80.gl #183 snps, .28% missing data

#sort gl files by population for the next plots
dcor.snp100.gl<-gl.sort(dcor.snp100.gl, sort.by="pop")
dcor.snp95.gl<-gl.sort(dcor.snp95.gl, sort.by="pop")
dcor.snp80.gl<-gl.sort(dcor.snp80.gl, sort.by="pop")

#these are heatmaps that show the genotype using a color for all individuals and snp loci
#they can be helpful in spotting linkage between SNPs and individuals (such as population structure)
glPlot(dcor.snp100.gl)
glPlot(dcor.snp95.gl)
glPlot(dcor.snp80.gl)
#there's no glaring linkage between snps or loci, with the exception of the last few samples from St Croix

###Create a phylogenetic tree of the samples based on the SNP genotypes
dcor.100.tree<-aboot(dcor.snp100.gl, tree="upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
cols <- brewer.pal(n = nPop(dcor.snp100.gl), name = "Dark2")
plot.phylo(dcor.100.tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(dcor.snp100.gl)])
nodelabels(dcor.100.tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend("topleft", legend=c("E. Pacific", "Indonesia", "Malaysia", "Papua New Guinea", "S. Africa", "Solomon Islands", "St Croix"), fill=c("#1B9E77", "#D95F02", "#7570B3" ,"#E7298A" ,"#66A61E", "#E6AB02" ,"#A6761D"), xpd=TRUE)


####Conduct a DAPC to visualise population structure and clustering
#Discriminant Analysis of Principle Components (DAPC) using the populations of origin for the samples
dcor.100.dapc<-dapc(dcor.snp100.gl)
#This will first plot the variance explained by PCA and you will need to select the number of PCs (Principle Components) to retain
##In this instance, retaining too many PCs may result in overfitting
##Best to limit the number of PCS based on where the cumulative variance starts to level out. 
##for this file, i will stop at 110
#You will then see a barplot of eigenvalues for the discriminant analysis and you will need to select the number of discriminant functions to retain
##Since there are not too many clusters, all can be retained
##WHen there are tens of clusters, it is likely that the first few dimensions will carry more informaton than the others
##I will retain all 6

scatter(dcor.100.dapc, posi.da="bottomleft") #this produces a scatter plot of the DAPC
#Each dot represents an individual and the elipses are the groups (populations)
#the presence and location of Eigenvalues and PC plots can be displayed in the scatterplot
scatter(dcor.100.dapc, posi.da="bottomleft", xax=1, yax=1, legend=TRUE) #this shows the same cluster data just on one axis
#Assignment of individuals to groups can help to idenfity how clear-cut the clusters are
#The assignment plot shows samples along the y axis and the populations along the x axis
#heat colors represent the membership probability (red=1, white=0)
#blue crosses represent the population of the individual in the orinal file
#Where blue crosses and red bars overlap, it means that the indivduals have re-assigned to their population of origin. 
#Individuals with less than 100% membership probability, they will show up in each population they have been assigned to. 
#In this example, theres one indidual that shares about 50% of its genetics with Indonesia and Solomon Islands populations (dcor102682_wpac)
assignplot(dcor.100.dapc)
round((dcor.100.dapc$posterior),3) #table of values that are used to make the assingment plot and compoplot below
#A compoplot is an alternate view to the same information
#each individual is presented as a bar with the value representing the membership probability to each population. 
compoplot(dcor.100.dapc)

###Adegenet also has a funtion whereby you can identify clusters based on the genotypes using k-means clustering and no prior population info
#Using "find clusters" to see how many clusters are in data
#this will help you find the natural number of clusters in the dataset, i have limited it to 40 clusters
grp<-find.clusters(dcor.snp100.gl, max.n.clust=40) 
##This will first plot the variance explained by the PCR and you must select the number to retain
#In this case, there's no reason not to include all PCs
#I am keeping all and will say 200 to make sure all are kept
##Then you will see a graph plotting the BIC vs number of clusters
#Typically, the optimal number of clusters will have the lowest BIC value, or is at least indicated by an elbow in the curve
#In this case, there is a clear deap at k=3 clusters

#the following table will display how individuals from the original populations cluster using k-means cluster
table(pop(dcor.snp100.gl), grp$grp)
##There are three clear clusters in the data: E. Pacific, W. Pacific and St Croix. 
##There are 2 E. pacific individuals that clustered with the W. Pacific group: 38952 and 30265, both bycaught animals from Peru

##It seems like structure is likely hierarchical within these samples, Could be worth running again just with the W. pacific samples

###NOTE: there is a function to do a DAPC with an optimisation step, helping to select the right number of PCs without overfitting the model. 
#when using xvalDAPC, the process is run 30 times with a subset of the data
#the number of PCs that has the highest proportion of successful predictions is selected
#in the example here, 60 PCs give the best model output
x<-dcor.snp100.gl
summary(x$pop)
mat<-tab(x)
grp<-pop(x)
xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
xval
scatter(xval$DAPC)
assignplot(xval$DAPC) 
#The scatter and assignment plots show generally the same information, however the divisions between populations are less clear, especially between Indonesia and the Solomon Islands. Some misclassifications when you do it this way. 

####Test differentiation between clusters using Fst
dcor.snp100.fst<-gl.fst.pop(dcor.snp100.gl, nboots=100)
dcor.snp100.fst
#All pairwise populations have Fst values that are statistically significant (<0.05)

####Demerelate
#Demerelate no longer on CRAN, have to install an old version
library(remotes)
install_github("cran/fts")
install_github("cran/demerelate")
library(Demerelate)

#Prepare Demerelate and run the analysis
#Demerelate creates a randomized offspring and non-relate individuals from the reference population (all data unless specificed)
#There can be no missing data in the file
#Thresholds for half and full sibling relationships are calculated 
#lots of different relatedness estimators can be used
dcor.snp100.gl
dcor.snp100.gl<-gl.filter.callrate(dcor.snp100.gl, method="loc", threshold=1, mono.rm=TRUE)#remove monomorphic loci and loci with missing values
dcor.snp100.gl #161 loci, none removed
dcor.snp100.dr<-gl2demerelate(dcor.snp100.gl, verbose=TRUE)#convert genlight to demerelate format
dcor.snp100.dr$`Sample-ID`<- dcor.snp100.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"=".")) #can't have "_" in sample ids, replacing with ".""
dcor.snp100.D <- Demerelate(dcor.snp100.dr, object = TRUE,
                     Fis = TRUE,
                     file.output = TRUE,
                     iteration = 10,  # increase to 1000 after testing
                     pairs = 10, 
                     p.correct=TRUE)

##Calculate genetic diversity parameters using Hierfstat package and dartR
#Calculates statistics per locus
dcor.snp100.stats<-gl.basic.stats(dcor.snp100.gl, digits=2)
dcor.snp100.stats
#Can summarise statistics per population
summary(dcor.snp100.stats$Ho)

#Calculate weir and cockerham Fis and 95% confidence intervals
#must first convert genlight to genind
dcor.snp100.gi<-gl2gi(dcor.snp100.gl)
#calculate 95% intervals around fis values per population
dcor.snp100.wcboot<-boot.ppfis(dcor.snp100.gi, nboot=1000)
dcor.snp100.wcboot
#add population IDs to the confidence intervals
dcor.snp100.wcfis<-cbind(levels(dcor.snp100.gi$pop), dcor.snp100.wcboot$fis.ci)
#convert genind to hierfstat format to calulate fis point estimates
dcor.snp100.hfstat<-genind2hierfstat(dcor.snp100.gi)
#Hierfstat will only calculate the Fis point estimates for all individuals together, so must first subset the data for each population
#This loop subsets the data into populations, calculates the Fis point estimates, and generates a table with the values
pop.wc.df<-data.frame()
pops<-unique(dcor.snp100.hfstat$pop)
for (i in pops){
  hfstat<-subset(dcor.snp100.hfstat, dcor.snp100.hfstat$pop==i) 
  wc<-wc(hfstat[,-1])
  pop.wc<-cbind(i, wc$FIS)
  pop.wc.df=rbind(pop.wc.df, pop.wc)
}
pop.wc.df #final table with population and wc Fis point estimates
#clean up the boot strap and point estimate tables so that they can be merged
colnames(pop.wc.df)<-c("pop", "wc_fis")
colnames(dcor.snp100.wcfis)<-c("pop", "lower_wc_fis", "upper_wc_fis")
#merge tables and reorder columns
dcor.snp100.wcfis<-merge(dcor.snp100.wcfis, pop.wc.df, by = "pop")
dcor.snp100.wcfis<-dcor.snp100.wcfis[,c(1,4,2,3)]
#reduce the number of decimal points to three
dcor.snp100.wcfis$wc_fis<-as.numeric(dcor.snp100.wcfis$wc_fis)
dcor.snp100.wcfis<-dcor.snp100.wcfis %>% mutate_if(is.numeric, ~round(.,3))
#Final table with population id, weir and cockerhap Fis point estimate and 95% confidence intervals
dcor.snp100.wcfis

                  