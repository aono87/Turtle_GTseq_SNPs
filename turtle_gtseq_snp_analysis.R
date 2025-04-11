#Script to analyse the 205 gtseq SNPs from leatherbacks

#Libraries
library(adegenet)
library(vcfR)
library(remotes)
library(SNPRelate)
library(dartRverse)
library(poppr)
library(ape)
library(RColorBrewer)

##Import vcf into R for use in adegenet
dcor.vcfR<-read.vcfR(("./Dcor_DcPanel_205_maf.targetSNPs_012224.recode.vcf"), verbose = TRUE)
dcor.vcfR
head(getFIX(dcor.vcfR))
is.biallelic(dcor.vcfR) #All are biallelic
getCHROM(dcor.vcfR)
getPOS(dcor.vcfR)
pop.data.dcor<-read.table("sample_pop.txt", sep="\t", header=TRUE)
pop.data.dcor
all(colnames(dcor.vcfR@gt)[-1]==pop.data.dcor$ind) #check that sample IDs are the same and in the same order

#convert vcfR file to genlight files
dcor.gl<-vcfR2genlight(dcor.vcfR)
dcor.gl
#n=177 individuals, n=231 loci

#Add pop IDs
pop(dcor.gl)<-pop.data.dcor$NewLoc

dcor.gl$loc.names
dcor.gl$loc.all
dcor.gl$ind.names
dcor.gl$pop
dcor.gl$gen

dcor.mat<-as.matrix(dcor.gl)
write.csv(dcor.mat, file="./dcor.snp.genotype.matrix.csv")

#Subset vcf file to only retain one snp per locus using gwscaR
drop_pop<-read.table("./loc_to_drop.txt")
drop_pop<-as.character(drop_pop)
print(drop_pop)

dcor.gl<-gl.compliance.check(dcor.gl)
dcor.snp.gl<-gl.drop.loc(dcor.gl, loc.list=c("locus022_63", "locus045_119", "locus053_99", "locus061_150", "locus100_115", "locus100_118", "locus129_60", "locus135_106", "locus188_138", "locus189_70", "locus237_111", "locus291_100", "locus334_109", "Dc01661_95", "Dc01707_114", "Dc03887_101", "Dc05864_137", "Dc07300_81", "Dc07937_134", "Dc08054_106", "Dc08189_82", "Dc08578_70", "Dc08865_124", "Dc10885_73", "Dc11062_41", "Dc11203_120", "Dc12329_111", "Dc13159_110", "Dc13592_122", "Dc16962_121", "Dc20232_74", "Dc21637_31", "Dc22013_123", "Dc22883_124", "Dc33714_113", "Dc36767_103", "Dc41234_117", "Dc46984_123", "Dc47572_144", "Dc47935_46", "Dc47935_56", "Dc48950_40", "Dc49184_94", "Dc52523_131", "Dc58248_125", "Dc62188_75"))
dcor.snp.gl
#n=177 individuals, #n=185 SNPs
#354 missing data, 1.08%

dcor.snp.gl<-gl.compliance.check(dcor.snp.gl)
dcor.snp2.gl<-gl.filter.callrate(dcor.snp.gl, method="ind", threshold=0.5, recalc=TRUE) #no samples with more than 50% missing data
dcor.snp100.gl<-gl.filter.callrate(dcor.snp.gl, method="loc", threshold=1, recalc=TRUE)
dcor.snp95.gl<-gl.filter.callrate(dcor.snp.gl, method="loc", threshold=.95, recalc=TRUE)
dcor.snp80.gl<-gl.filter.callrate(dcor.snp.gl, method="loc", threshold=.8, recalc=TRUE)

dcor.snp100.gl #161 snps, 0% missing data
dcor.snp95.gl #181 snps, .19% missing data
dcor.snp80.gl #183 snps, .25% missing data

#sort gl files by population
dcor.snp100.gl<-gl.sort(dcor.snp100.gl, sort.by="pop")
dcor.snp95.gl<-gl.sort(dcor.snp95.gl, sort.by="pop")
dcor.snp80.gl<-gl.sort(dcor.snp80.gl, sort.by="pop")

glPlot(dcor.snp100.gl)
glPlot(dcor.snp95.gl)
glPlot(dcor.snp80.gl)

#based on https://pmc.ncbi.nlm.nih.gov/articles/PMC9305793/
#remove individuals with >50% missing data
#Ho/He and Gis for each subpopulation in genodive
#structure k=1-10
#use k value to make dapc, retaining pcs to capture 80% of variance
#assingpop to assess assignment accuracy

###Tree
dcor.100.tree<-aboot(dcor.snp100.gl, tree="upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
cols <- brewer.pal(n = nPop(dcor.snp100.gl), name = "Dark2")
plot.phylo(dcor.100.tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(dcor.snp100.gl)])
nodelabels(dcor.100.tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend("topleft", legend=c("E. Pacific", "Indonesia", "Malaysia", "Papua New Guinea", "S. Africa", "Solomon Islands", "St Croix"), fill=c("#1B9E77", "#D95F02", "#7570B3" ,"#E7298A" ,"#66A61E", "#E6AB02" ,"#A6761D"), xpd=TRUE)


####DAPC
dcor.100.dapc<-dapc(dcor.snp100.gl)
scatter(dcor.100.dapc, posi.da="bottomleft")
scatter(dcor.100.dapc, posi.da="bottomleft", xax=2, yax=3)
assignplot(dcor.100.dapc)
compoplot(dcor.100.dapc)

###Using "find clusters" to see how many clusters are in data
grp<-find.clusters(dcor.snp100.gl, max.n.clust=40) #retain 150 PCs, #retain 3 clusters
names(grp)
grp$grp
table(pop(dcor.snp100.gl), grp$grp)
##There are three clear clusters in the data: E. Pacific, W. Pacific and St Croix. 
##There are 4 E. pacific individuals that clustered with the W. Pacific group: 38952, 38952b, 30265, 30265b
##38952: Peru Bycatch
##30265: Peru Bycatch

####Test differentiation between clusters
dcor.snp100.fst<-gl.fst.pop(dcor.snp100.gl, nboots=100)



