###Exploring vcf files in vcfR. Attempting to visualise

# libraries ---------------------------------------------------------------


install.packages("BiocManager")
install.packages("devtools")
BiocManager::install(c("SNPRelate", "qvalue"))
devtools::install_github("green-striped-gecko/dartR")
devtools::install_github('ramnathv/htmlwidgets')
devtools::install_github("ropensci/plotly")
install.packages("~/Downloads/related_1.0.tar.gz", repos=NULL, type="source")
devtools::install_github("gertstulp/ggplotgui")
devtools::install_github("thibautjombart/adegenet")
install.packages("leaflet.providers")
install.packages("lifecycle")
install.packages("farver")
install_github("green-striped-gecko/dartR", ref="dev")
install.packages("dartR")
library(dartR)
gl.install.vanilla.dartR()

library(leaflet.providers)
library(lifecycle)
library(farver)
library(BiocManager)
library(devtools)
library(dartR)
library(htmlwidgets)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap)
library(maps)
library(vcfR)
library(adegenet)
library(ggplot2)
library(devtools)
library(Rcpp)
library(plotly)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap)
library(tess3r)
library(R.filesets)
library(igraph)
library(poppr)
library(RColorBrewer)
library(Demerelate)
library(strataG)
library(tidyverse)
library(dplyr)
library(sfsmisc)
library(mlogit)
library(fts)
library(vegan)
library(stringr)
library(related)
library(ggplotgui)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(RColorBrewer)
library(hierfstat)

# Opt. uploading files and checking dataset stats ------------------------------


##Import vcf file into R from file
#convert vcf file to vcfR file
zczc125.vcfR<-read.vcfR(("~/Dropbox/Phd/Bioinformatics/bw_ddrad_global/SNP_files/no_populations/one_snp_per_locus/zczc_snp_files/zczc125.snps.vcf"), verbose = TRUE)
zczc125.vcfR


##Checking glplots

# md18_old_stacks<-read.vcfR("./vcf_files/mden_mden_18_old_stacks.vcf")
# md18_old_stacks.gl<-vcfR2genlight(md18_old_stacks)
# md18_old_stacks.gl
# glPlot(md18_old_stacks.gl)
# 
# md18_new_stacks<-read.vcfR("./vcf_files/mden_mden_18_new_stacks.vcf")
# md18_new_stacks.gl<-vcfR2genlight(md18_new_stacks)
# md18_new_stacks.gl
# glPlot(md18_new_stacks.gl)
# 
# md18_new_stacks_no_wl<-read.vcfR("./vcf_files/mden_mden_18_new_stacks_no_whitelist.vcf")
# md18_new_stacks_no_wl.gl<-vcfR2genlight(md18_new_stacks_no_wl)
# md18_new_stacks_no_wl.gl
# glPlot(md18_new_stacks_no_wl.gl)

md2018.vcfr<-read.vcfR("./md_md_2018_opt_final_sorted.vcf")
md2020.vcfr<-read.vcfR("./md_md_2020_opt_final_sorted.vcf")
mdcomb.vcfr<-read.vcfR("./md_md_combined_opt_final_sorted.vcf")
mdmerge.vcfr<-read.vcfR("./mden/vcf-files/md_md_merged.vcf")

zc2018.vcfR<-read.vcfR("./vcf_files/opt/zcav/zc_zc_2018_sorted.vcf")
zc2020.vcfR<-read.vcfR("./vcf_files/opt/zcav/zc_zc_2020_sorted.vcf")
zccomb.vcfR<-read.vcfR("./vcf_files/opt/zcav/zc_zc_comb_sorted.vcf")
zcmerge.vcfR<-read.vcfR("./zcav/vcf-files/zc_zc_merged.vcf")

md2018.vcfr


md2018.gl<-vcfR2genlight(md2018.vcfr)
md2020.gl<-vcfR2genlight(md2020.vcfr)
mdcomb.gl<-vcfR2genlight(mdcomb.vcfr)
mdmerge.gl<-vcfR2genlight(mdmerge.vcfr)

zc2018.gl<-vcfR2genlight(zc2018.vcfR)
zc2020.gl<-vcfR2genlight(zc2020.vcfR)
zccomb.gl<-vcfR2genlight(zccomb.vcfR)
zcmerge.gl<-vcfR2genlight(zcmerge.vcfR)

glPlot(md2018.gl, main="M. densirostris aligned to M. densirostris, optimisation subset, 2018")
glPlot(md2020.gl, main="M. densirostris aligned to M. densirostris, optimisation subset, 2020")
glPlot(mdcomb.gl, main="M. densirostris aligned to M. densirostris, optimisation subset, combined")

ggsave("./md2018_glplot.png", plot=glPlot(md2018.gl, main="M. densirostris aligned to M. densirostris, optimisation subset, 2018"))
ggsave("./md2020_glplot.png", plot=glPlot(md2020.gl, main="M. densirostris aligned to M. densirostris, optimisation subset, 2020"))
ggsave("./mdcombined_glplot.png", plot=glPlot(mdcomb.gl, main="M. densirostris aligned to M. densirostris, optimisation subset, combined"))
ggsave("./mdmerged_glplot.png", plot=glPlot(mdmerge.gl, main="M. densirostris aligned to M. densirostris, optimisation subset, merged"))


glPlot(zc2018.gl, main="Z. cavirostris aligned to Z. cavirostris, optimisation subset, 2018")
ggsave("./zc2018_glplot.png", plot=glPlot(zc2018.gl, main="Z. cavirostris aligned to Z. cavirostris, optimisation subset, 2018"))
glPlot(zc2020.gl, main="Z. cavirostris aligned to Z. cavirostris, optimisation subset, 2020")
ggsave("./zc2020_glplot.png", plot=glPlot(zc2020.gl, main="Z. cavirostris aligned to Z. cavirostris, optimisation subset, 2020"))
glPlot(zccomb.gl, main="Z. cavirostris aligned to Z. cavirostris, optimisation subset, combined")
ggsave("./zccombined_glplot.png", plot=glPlot(zccomb.gl, main="Z. cavirostris aligned to Z. cavirostris, optimisation subset, combined"))
ggsave("./zcmerged_glplot.png", plot=glPlot(zcmerge.gl, main="Z. cavirostris aligned to Z. cavirostris, optimisation subset, merged"))

md2018.gl #17 genotypes,  34,443 binary SNPs, size: 3.1 Mb, 17732 (3.03 %) missing data
md2020.gl #17 genotypes,  32,707 binary SNPs, size: 3.0 Mb, 28907 (5.2 %) missing data
mdcomb.gl #32 genotypes,  32,998 binary SNPs, size: 3.2 Mb, 44934 (4.26 %) missing data
mdmerge.gl #19 genotypes,  36,441 binary SNPs, size: 3.3 Mb, 27717 (4 %) missing data

zc2018.gl #24 genotypes,  81,211 binary SNPs, size: 8.7 Mb, 88123  (4.52 %) missing data
zc2020.gl #22 genotypes,  80,664 binary SNPs, size: 8.7 Mb, 100169 (5.64 %) missing data
zccomb.gl #46 genotypes,  64,093 binary SNPs, size: 7.9 Mb, 169275 (5.74 %) missing data
zcmerge.gl #24 genotypes, 75,375 binary SNPs, size: 8.5 Mb, 96396 (5.33 %) missing data


gl.report.callrate(md2018.gl)
# Loci with no missing values = 24096 [70%]
#        < 5% missing values =  24096 [70%]
#        < 10% missing values = 29227 [84.9%]
#        < 15% missing values = 32274 [93.7%]
#        < 20% missing values = 34443 [100%]
gl.report.callrate(md2020.gl)
# Loci with no missing values = 16303 [49.8%]
#         < 5% missing values = 16303 [49.8%]
#        < 10% missing values = 24091 [73.7%]
#        < 15% missing values = 28820 [88.1%]
#        < 20% missing values = 32707 [100%]
gl.report.callrate(mdcomb.gl)
# Loci with no missing values = 16871 [51.1%]
#         < 5% missing values = 22189 [67.2%]
#        < 10% missing values = 27559 [83.5%]
#        < 15% missing values = 29638 [89.8%]
#        < 20% missing values = 32998 [100%]
gl.report.callrate(mdmerge.gl)
#Loci with no missing values = 20024 [54.9%]
#       < 5% missing values =  20024 [54.9%]
#       < 10% missing values = 28587 [78.4%]
#       < 15% missing values = 32995 [90.5%]
#       < 20% missing values = 36441 [100%]
gl.report.callrate(zc2018.gl)
#Loci with no missing values = 40489 [49.9%]
#        < 5% missing values = 54625 [67.3%]
#       < 10% missing values = 66837 [82.3%]
#       < 15% missing values = 74770 [92.1%]
#       < 20% missing values = 81211 [100%]
gl.report.callrate(zc2020.gl)
#Loci with no missing values = 33210 [41.2%]
#        < 5% missing values = 52198 [64.7%]
#       < 10% missing values = 64110 [79.5%]
#       < 15% missing values = 72969 [90.5%]
#       < 20% missing values = 80664 [100%]
gl.report.callrate(zccomb.gl)
#Loci with no missing values = 19943 [31.1%]
#        < 5% missing values = 38370 [59.9%]
#       < 10% missing values = 48056 [75%]
#       < 15% missing values = 55101 [86%]
#       < 20% missing values = 64093 [100%]
gl.report.callrate(zcmerge.gl)
# Loci with no missing values = 27518 [36.5%]
#        < 5% missing values  = 48093 [63.8%]
#       < 10% missing values  = 60707 [80.5%]
#       < 15% missing values  = 68786 [91.3%]
#       < 20% missing values  = 75375 [100%]

gl.report.callrate(md2018.gl, method = "ind")
#  Individuals no missing values = 0  [0%]     across loci; all individuals would be filtered
#  with less than or equal to 5% = 16 [94.1%]; 1 individuals would be filtered
# with less than or equal to 10% = 16 [94.1%]; 1 individuals would be filtered
# with less than or equal to 15% = 17 [100%]  ; 0 individuals would be filtered
gl.report.callrate(md2020.gl, method = "ind")
#  Individuals no missing values = 0 [0%]      across loci; all individuals would be filtered
#  with less than or equal to 5% = 12 [70.6%]; 5 individuals would be filtered
# with less than or equal to 10% = 12 [70.6%]; 5 individuals would be filtered
# with less than or equal to 15% = 16 [94.1%]; 1 individuals would be filtered
# with less than or equal to 20% = 17 [100%] ; 0 individuals would be filtered
gl.report.callrate(mdcomb.gl, method = "ind")
# Individuals no missing values =  0 [0%]     across loci; all individuals would be filtered
# with less than or equal to 5% =  24 [75%];  8 individuals would be filtered
# with less than or equal to 10% = 30 [93.8%]; 2 individuals would be filtered
# with less than or equal to 15% = 31 [96.9%]; 1 individuals would be filtered
# with less than or equal to 20% = 32 [100%] ; 0 individuals would be filtered
gl.report.callrate(mdmerge.gl, method="ind")
# Individuals no missing values =   0 [0%]     across loci; all individuals would be filtered
# with less than or equal to 5% =  14 [73.7%]; 5 individuals would be filtered
# with less than or equal to 10% = 18 [94.7%]; 1 individuals would be filtered
# with less than or equal to 15% = 18 [94.7%]; 1 individuals would be filtered
# with less than or equal to 20% = 19 [100%] ;  0 individuals would be filtered

gl.report.callrate(zc2018.gl, method="ind")
#Individuals no missing values =  0 [0%]      across loci; all individuals would be filtered
#with less than or equal to 5% =  16 [66.7%]; 8 individuals would be filtered
#with less than or equal to 10% = 22 [91.7%]; 2 individuals would be filtered
#with less than or equal to 15% = 23 [95.8%]; 1 individuals would be filtered
#with less than or equal to 20% = 23 [95.8%]; 1 individuals would be filtered
#with less than or equal to 25% = 24 [100%] ; 0 individuals would be filtered
gl.report.callrate(zc2020.gl, method="ind")
#Individuals no missing values =  0 [0%]      across loci; all individuals would be filtered
#with less than or equal to 5% =  15 [68.2%]; 7 individuals would be filtered
#with less than or equal to 10% = 19 [86.4%]; 3 individuals would be filtered
#with less than or equal to 15% = 19 [86.4%]; 3 individuals would be filtered
#with less than or equal to 20% = 21 [95.5%]; 1 individuals would be filtered
#with less than or equal to 25% = 22 [100%]  ; 0 individuals would be filtered
gl.report.callrate(zccomb.gl, method="ind")
#Individuals no missing values =  0 [0%]      across loci; all individuals would be filtered
#with less than or equal to 5% =  28 [60.9%]; 18 individuals would be filtered
#with less than or equal to 10% = 40 [87%]   ; 6 individuals would be filtered
#with less than or equal to 15% = 42 [91.3%]; 4 individuals would be filtered
#with less than or equal to 20% = 45 [97.8%]; 1 individuals would be filtered
#with less than or equal to 25% = 46 [100%] ; 0 individuals would be filtered
gl.report.callrate(zcmerge.gl, method="ind")
# Individuals no missing values =  0 [0%]       across loci; all individuals would be filtered
# with less than or equal to 5% =  18 [75%];   6 individuals would be filtered
# with less than or equal to 10% = 20 [83.3%]; 4 individuals would be filtered
# with less than or equal to 15% = 22 [91.7%]; 2 individuals would be filtered
# with less than or equal to 20% = 23 [95.8%]; 1 individuals would be filtered
# with less than or equal to 25% = 24 [100%]  ; 0 individuals would be filtered
gl.report.monomorphs(md2018.gl) 
#Breakdown of       34443 loci
# Polymorphic loci: 34288 
# Monomorphic loci: 155 
# Loci with no scores (all NA): 0
# Loci with no scores (all NA): 0
gl.report.monomorphs(md2020.gl) 
# Breakdown of      32707 loci
# Polymorphic loci: 32591 
# Monomorphic loci: 116 
# Loci with no scores (all NA): 0 
gl.report.monomorphs(mdcomb.gl) 
# Breakdown of      32998 loci
# Polymorphic loci: 32973 
# Monomorphic loci: 25 
# Loci with no scores (all NA): 0 
gl.report.monomorphs(mdmerge.gl)
# Breakdown of      36441 loci
# Polymorphic loci: 36356 
# Monomorphic loci: 85 
# Loci with no scores (all NA): 0 
gl.report.monomorphs(zc2018.gl)
#Breakdown of      81211 loci
#Polymorphic loci: 81149 
#Monomorphic loci: 62 
#Loci with no scores (all NA): 0 
gl.report.monomorphs(zc2020.gl)
#Breakdown of      80664 loci
#Polymorphic loci: 80622 
#Monomorphic loci: 42 
#Loci with no scores (all NA): 0 
gl.report.monomorphs(zccomb.gl)
#Breakdown of      64093 loci
#Polymorphic loci: 64082 
#Monomorphic loci: 11 
#Loci with no scores (all NA): 0
gl.report.monomorphs(zcmerge.gl)
#Breakdown of      75375 loci
#Polymorphic loci: 75348 
#Monomorphic loci: 27 
#Loci with no scores (all NA): 0 



# Opt. stats with vcftools -----------------------------------------------------
#looking at stats with vcftools
library(readr)
library(scales)
library(reshape2)
library(tidyverse)
library(stringr)
library(ggplot2)

##calculate missing data percentage
# Produce a file with missingness per individual
system("vcftools --vcf md_md_2018_opt_final_sorted.vcf --out out.1 --missing-indv")

#Import missingness per individual into R from vcftools output file
out_1_imiss <- read_tsv("out.1.imiss")

#calculate average and SD of missingness per indidivual
mean.missingness<-(round(mean(out_1_imiss$F_MISS), digits=5)*100)
sd.missingness<-(round(sd(out_1_imiss$F_MISS), digits=5)*100)

q<-ggplot(out_1_imiss, aes(x = INDV, y=F_MISS)) +
  geom_bar(stat = "identity") + scale_y_continuous(labels=percent_format()) +
  theme(axis.text.x = element_text(size=6, angle = 90),
           axis.text.y = element_text(size=10),
           legend.text=element_text(size=10),
           legend.title=element_blank()) +
  labs(title="M. densirostris 2018 optimisation, missingness per individual") +
  annotate("text", x = 5, y = .2, label = (paste("mean missingness per individual=", (mean.missingness), "%"))) +
  annotate("text", x = 5, y = .19, label = (paste("s.d. missingness per individual=", (sd.missingness), "%")))
q

#save plot to file
ggsave("./mden18_missingness.png", plot=last_plot())

# Produce a file with missingness per individual
system("vcftools --vcf md_md_2020_opt_final_sorted.vcf --out out.1 --missing-indv")

#Import missingness per individual into R from vcftools output file
out_1_imiss <- read_tsv("out.1.imiss")

#calculate average and SD of missingness per indidivual
mean.missingness<-(round(mean(out_1_imiss$F_MISS), digits=5)*100)
sd.missingness<-(round(sd(out_1_imiss$F_MISS), digits=5)*100)

q<-ggplot(out_1_imiss, aes(x = INDV, y=F_MISS)) +
  geom_bar(stat = "identity") + scale_y_continuous(labels=percent_format()) +
  theme(axis.text.x = element_text(size=6, angle = 90),
        axis.text.y = element_text(size=10),
        legend.text=element_text(size=10),
        legend.title=element_blank()) +
  labs(title="M. densirostris 2020 optimisation, missingness per individual") +
  annotate("text", x = 5, y = .2, label = (paste("mean missingness per individual=", (mean.missingness), "%"))) +
  annotate("text", x = 5, y = .19, label = (paste("s.d. missingness per individual=", (sd.missingness), "%")))
q

#save plot to file
ggsave("./mden20_missingness.png", plot=last_plot())

# Produce a file with missingness per individual
system("vcftools --vcf md_md_combined_opt_final_sorted.vcf --out out.1 --missing-indv")

#Import missingness per individual into R from vcftools output file
out_1_imiss <- read_tsv("out.1.imiss")

#calculate average and SD of missingness per indidivual
mean.missingness<-(round(mean(out_1_imiss$F_MISS), digits=5)*100)
sd.missingness<-(round(sd(out_1_imiss$F_MISS), digits=5)*100)

q<-ggplot(out_1_imiss, aes(x = INDV, y=F_MISS)) +
  geom_bar(stat = "identity") + scale_y_continuous(labels=percent_format()) +
  theme(axis.text.x = element_text(size=6, angle = 90),
        axis.text.y = element_text(size=10),
        legend.text=element_text(size=10),
        legend.title=element_blank()) +
  labs(title="M. densirostris combined optimisation, missingness per individual") +
  annotate("text", x = 15, y = .2, label = (paste("mean missingness per individual=", (mean.missingness), "%"))) +
  annotate("text", x = 15, y = .19, label = (paste("s.d. missingness per individual=", (sd.missingness), "%")))
q

#save plot to file
ggsave("./mdencomb_missingness.png", plot=last_plot())

# Produce a file with missingness per individual
system("vcftools --vcf zc_zc_2018_sorted.vcf --out out.1 --missing-indv")

#Import missingness per individual into R from vcftools output file
out_1_imiss <- read_tsv("out.1.imiss")

#calculate average and SD of missingness per indidivual
mean.missingness<-(round(mean(out_1_imiss$F_MISS), digits=5)*100)
sd.missingness<-(round(sd(out_1_imiss$F_MISS), digits=5)*100)

q<-ggplot(out_1_imiss, aes(x = INDV, y=F_MISS)) +
  geom_bar(stat = "identity") + scale_y_continuous(labels=percent_format()) +
  theme(axis.text.x = element_text(size=6, angle = 90),
        axis.text.y = element_text(size=10),
        legend.text=element_text(size=10),
        legend.title=element_blank()) +
  labs(title="Z. cavirostris 2018 optimisation, missingness per individual") +
  annotate("text", x = 15, y = .2, label = (paste("mean missingness per individual=", (mean.missingness), "%"))) +
  annotate("text", x = 15, y = .19, label = (paste("s.d. missingness per individual=", (sd.missingness), "%")))
q

#save plot to file
ggsave("./zc18_missingness.png", plot=last_plot())

# Produce a file with missingness per individual
system("vcftools --vcf zc_zc_2020_sorted.vcf --out out.1 --missing-indv")

#Import missingness per individual into R from vcftools output file
out_1_imiss <- read_tsv("out.1.imiss")

#calculate average and SD of missingness per indidivual
mean.missingness<-(round(mean(out_1_imiss$F_MISS), digits=5)*100)
sd.missingness<-(round(sd(out_1_imiss$F_MISS), digits=5)*100)

q<-ggplot(out_1_imiss, aes(x = INDV, y=F_MISS)) +
  geom_bar(stat = "identity") + scale_y_continuous(labels=percent_format()) +
  theme(axis.text.x = element_text(size=6, angle = 90),
        axis.text.y = element_text(size=10),
        legend.text=element_text(size=10),
        legend.title=element_blank()) +
  labs(title="Z. cavirostris 2020 optimisation, missingness per individual") +
  annotate("text", x = 10, y = .2, label = (paste("mean missingness per individual=", (mean.missingness), "%"))) +
  annotate("text", x = 10, y = .19, label = (paste("s.d. missingness per individual=", (sd.missingness), "%")))
q

#save plot to file
ggsave("./zc20_missingness.png", plot=last_plot())

# Produce a file with missingness per individual
system("vcftools --vcf zc_zc_comb_sorted.vcf --out out.1 --missing-indv")

#Import missingness per individual into R from vcftools output file
out_1_imiss <- read_tsv("out.1.imiss")

#calculate average and SD of missingness per indidivual
mean.missingness<-(round(mean(out_1_imiss$F_MISS), digits=5)*100)
sd.missingness<-(round(sd(out_1_imiss$F_MISS), digits=5)*100)

q<-ggplot(out_1_imiss, aes(x = INDV, y=F_MISS)) +
  geom_bar(stat = "identity") + scale_y_continuous(labels=percent_format()) +
  theme(axis.text.x = element_text(size=6, angle = 90),
        axis.text.y = element_text(size=10),
        legend.text=element_text(size=10),
        legend.title=element_blank()) +
  labs(title="Z. cavirostris combined optimisation, missingness per individual") +
  annotate("text", x = 25, y = .21, label = (paste("mean missingness per individual=", (mean.missingness), "%"))) +
  annotate("text", x = 25, y = .20, label = (paste("s.d. missingness per individual=", (sd.missingness), "%")))
q

#save plot to file
ggsave("./zccomb_missingness.png", plot=last_plot())


# Opt. testing removing missing data -------------------------------------------
#testing using the two merged datasets
zcmerge.gl
mdmerge.gl

zcmerge.gl<-gl.compliance.check(zcmerge.gl)
mdmerge.gl<-gl.compliance.check(mdmerge.gl)
#dartR filters that can be used:
#gl.filter.callrate() calculate call rate (proportion with non-missing scores) for each locus or indivdiual and remove those loci or individual sbelow a specified threshold
#gl.filter.monomorphs() remove all monomorphic loci, including loci for which the scores are all missing (NA)

#filtering by callrate
#will try for 95% and 100%
#first will remove the monomorphic loci - DIDNT SEEM TO REMOVE ANY LOCI?

#filter loci by callrate
md95.gl<-gl.filter.callrate(mdmerge.gl, method="loc", threshold=0.95, plot=TRUE, v=3, mono.rm = TRUE)
md95.gl #19 genotypes,  19,988 binary SNPs, size: 3.2 Mb 0 (0 %) missing data
glPlot(md95.gl)

md100.gl<-gl.filter.callrate(mdmerge.gl, method="loc", threshold=1, plot=TRUE, v=3, mono.rm = TRUE)
md100.gl #19 genotypes,  19,988 binary SNPs, size: 3.2 Mb 0 (0 %) missing data
dev.off()
glPlot(md100.gl)

zc95.gl<-gl.filter.callrate(zcmerge.gl, method="loc", threshold=0.95, plot=TRUE, v=3, mono.rm=TRUE)
zc95.gl #24 genotypes,  48,079 binary SNPs, size: 9.1 Mb 20571 (1.78 %) missing data
dev.off() 
glPlot(zc95.gl)

zc100.gl<-gl.filter.callrate(zcmerge.gl, method="loc", threshold=1, plot=TRUE, v=3, mono.rm=TRUE)
zc100.gl #24 genotypes,  27,508 binary SNPs, size: 6 Mb 0 (0 %) missing data
dev.off()
glPlot(zc100.gl)

#add some population data to the two files
indNames(mdmerge.gl)
indNames(zcmerge.gl)

pop(mdmerge.gl)<-as.factor(c("Caribbean",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Caribbean",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Caribbean",
                             "Canaries",
                             "Caribbean",
                             "Caribbean"))
pop(zcmerge.gl)<-as.factor(c("Mediterranean",
                             "Canaries",
                             "Caribbean",
                             "Mediterranean",
                             "Canaries",
                             "Canaries",
                             "Mediterranean",
                             "Canaries",
                             "Caribbean",
                             "Caribbean",
                             "Caribbean",
                             "Mediterranean",
                             "Caribbean",
                             "Mediterranean",
                             "Canaries",
                             "Canaries",
                             "Mediterranean",
                             "Canaries",
                             "Mediterranean",
                             "Mediterranean",
                             "Caribbean",
                             "Mediterranean",
                             "Caribbean",
                             "Caribbean"))
pop(md95.gl)<-as.factor(c("Caribbean",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Caribbean",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Caribbean",
                             "Canaries",
                             "Caribbean",
                             "Caribbean"))
pop(zc95.gl)<-as.factor(c("Mediterranean",
                             "Canaries",
                             "Caribbean",
                             "Mediterranean",
                             "Canaries",
                             "Canaries",
                             "Mediterranean",
                             "Canaries",
                             "Caribbean",
                             "Caribbean",
                             "Caribbean",
                             "Mediterranean",
                             "Caribbean",
                             "Mediterranean",
                             "Canaries",
                             "Canaries",
                             "Mediterranean",
                             "Canaries",
                             "Mediterranean",
                             "Mediterranean",
                             "Caribbean",
                             "Mediterranean",
                             "Caribbean",
                             "Caribbean"))
pop(md100.gl)<-as.factor(c("Caribbean",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Caribbean",
                             "Canaries",
                             "Canaries",
                             "Canaries",
                             "Caribbean",
                             "Canaries",
                             "Caribbean",
                             "Caribbean"))
pop(zc100.gl)<-as.factor(c("Mediterranean",
                             "Canaries",
                             "Caribbean",
                             "Mediterranean",
                             "Canaries",
                             "Canaries",
                             "Mediterranean",
                             "Canaries",
                             "Caribbean",
                             "Caribbean",
                             "Caribbean",
                             "Mediterranean",
                             "Caribbean",
                             "Mediterranean",
                             "Canaries",
                             "Canaries",
                             "Mediterranean",
                             "Canaries",
                             "Mediterranean",
                             "Mediterranean",
                             "Caribbean",
                             "Mediterranean",
                             "Caribbean",
                             "Caribbean"))
mdmerge.gl
md95.gl
md100.gl
zcmerge.gl
zc95.gl
zc100.gl

mdmerge.pc<-gl.pcoa(mdmerge.gl, nfactors=5)
md95.pc<-gl.pcoa(md95.gl, nfactors=5)
md100.pc<-gl.pcoa(md100.gl, nfactors=5)
zcmerge.pc<-gl.pcoa(zcmerge.gl, nfactors=5, parallel=TRUE)
zc95.pc<-gl.pcoa(zc95.gl, nfactors=5, parallel=TRUE)
zc100.pc<-gl.pcoa(zc100.gl, nfactors=5, parallel=TRUE)

gl.pcoa.plot(pc, gl, labels="pop", xaxis=1, yaxis=2)

mdmerge.pcoa<-gl.pcoa.plot(mdmerge.pc, mdmerge.gl, labels="pop", xaxis=1, yaxis=2)
md95.pcoa   <-gl.pcoa.plot(md95.pc, md95.gl, labels="pop", xaxis=1, yaxis=2)
md100.pcoa  <-gl.pcoa.plot(md100.pc, md100.gl, labels="pop", xaxis=1, yaxis=2)
zcmerge.pcoa<-gl.pcoa.plot(zcmerge.pc, zcmerge.gl, labels="pop", xaxis=1, yaxis=2)
zc95.pcoa   <-gl.pcoa.plot(zc95.pc, zc95.gl, labels="pop", xaxis=1, yaxis=2)
zc100.pcoa  <-gl.pcoa.plot(zc100.pc, zc100.gl, labels="pop", xaxis=1, yaxis=2)

mdmerge.pcoa
md95.pcoa   
md100.pcoa  
zcmerge.pcoa
zc95.pcoa   
zc100.pcoa  

mdmerge.pc <-  glPca(mdmerge.gl, nf=20, parallel=TRUE)
md95.pc    <-  glPca(md95.gl  , nf=20, parallel=TRUE)
md100.pc   <-  glPca(md100.gl , nf=20, parallel=TRUE)
zcmerge.pc <-  glPca(zcmerge.gl, nf=20, parallel=TRUE)
zc95.pc    <-  glPca(zc95.gl  , nf=20, parallel=TRUE)
zc100.pc   <-  glPca(zc100.gl , nf=20, parallel=TRUE)

scatter(mdmerge.pc)
scatter(md95.pc   )
scatter(md100.pc  )
scatter(zcmerge.pc)
scatter(zc95.pc   )
scatter(zc100.pc  )

md.dapc1<-dapc(mdmerge.gl, n.pca=10, n.da=10)
scatter(md.dapc1)

zc.dapc1<-dapc(zcmerge.gl, n.pca=10, n.da=10)
scatter(zc.dapc1)

indNames(md2018.gl)
indNames(md2020.gl)
indNames(mdcomb.gl)

pop(md2018.gl)<-as.factor(c("Canaries-EH",
                            "Canaries-West",
                            "Canaries-EH",
                            "Florida",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Puerto Rico",
                            "Canaries-East",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-East",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Bahamas",
                            "Canaries-EH"))
pop(md2020.gl)<-as.factor(c("Bahamas",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-West",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Puerto Rico",
                            "Canaries-EH",
                            "Canaries-East",
                            "Canaries-EH",
                            "Puerto Rico",
                            "Canaries-EH",
                            "Bahamas"))
pop(mdcomb.gl)<-as.factor(c("Canaries-EH",
                            "Canaries-West",
                            "Canaries-EH",
                            "Florida",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Puerto Rico",
                            "Canaries-East",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-East",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Bahamas",
                            "Canaries-EH",
                            "Bahamas",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-West",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Canaries-EH",
                            "Puerto Rico",
                            "Canaries-EH",
                            "Canaries-East",
                            "Canaries-EH",
                            "Bahamas"))

md2018.dapc<-dapc(md2018.gl, n.pca=10, n.da=10)
md2020.dapc<-dapc(md2020.gl, n.pca=10, n.da=10)
mdcomb.dapc<-dapc(mdcomb.gl, n.pca=10, n.da=10)
md.dapc1<-dapc(mdmerge.gl, n.pca=10, n.da=5)


scatter(md2018.dapc)
scatter(md2020.dapc)
scatter(mdcomb.dapc)
scatter(md.dapc1)

indNames(mdmerge.gl)
pop(mdmerge.gl)<-as.factor(c("Bahamas",
                             "Canaries-East",
                             "Canaries-EH",
                             "Canaries-EH",
                             "Canaries-West",
                             "Canaries-EH",
                             "Canaries-EH",
                             "Canaries-EH",
                             "Canaries-EH",
                             "Canaries-EH",
                             "Canaries-EH",
                             "Puerto Rico",
                             "Canaries-EH",
                             "Canaries-East",
                             "Canaries-EH",
                             "Puerto Rico",
                             "Canaries-EH",
                             "Bahamas",
                             "Florida"))

#the one West Canary island sample seems to be skewing everying: Mde_L1-3		Mde_300713_LP				md-3_6
#will remove and re-try
mdmerge.gl

gl.make.recode.ind(mdmerge.gl, out.recode.file="./md_new_ind_assignments.csv", outpath="./")
mdmerge.trim.gl<-gl.recode.ind(mdmerge.gl, ind.recode="./md_new_ind_assignments.csv")
mdmerge.trim.gl<-gl.recalc.metrics(mdmerge.trim.gl)

md.trim.dapc<-dapc(mdmerge.trim.gl, n.pca=10, n.da=10)
scatter(md.trim.dapc)


# Final Files -------------------------------------------------------------

#upload vcf file into R and convert to a Genlight file 
#these vcf files were generated form stacks without population info
# mden.vcfR<-read.vcfR(("~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/final/mden_final.vcf"), verbose = TRUE)
# mden.gl<-vcfR2genlight(mden.vcfR)
# mden.gl #46 genotypes,  128,044 binary SNPs, size: 12.6 Mb, 269607 (4.58 %) missing data
# 
# zcav.vcfR<-read.vcfR(("~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/final/zcav_final.vcf"), verbose = TRUE)
# zcav.gl<-vcfR2genlight(zcav.vcfR)
# zcav.gl #92 genotypes,  171,163 binary SNPs, size: 22.7 Mb, 913998 (5.8 %) missing data

#set working directory for any results
setwd("~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/results/mden")
setwd("~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/results/zcav")

#these are the vcf files that were generated with pop information 
mden.vcfR<-read.vcfR("~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/final_pop/mden_pop_final.vcf")
zcav.vcfR<-read.vcfR("~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/final_pop/zcav_pop_final.vcf")

mden.dp<-extract.gt(mden.vcfR, element="DP", as.numeric=TRUE)
mden.dp[1:5,1:5]
colnames(mden.dp)


mden.fr<-mden.dp[,c("md-10_1_merged",
           "md-14_9_merged",
           "md-16_3",
           "BWLib_L4-7",
           "md-14_3_merged",
           "md-15_2",
           "md-9_5_merged",
           "md-11_9_merged",
           "md-8_4_merged",
           "md-9_4",
           "md-1_10",
           "md-10_10",
           "md-10_7",
           "md-12_3",
           "md-2_5",
           "md-3_2_merged",
           "md-3_6_merged",
           "md-3_7_merged",
           "md-4_4_merged",
           "md-4_5",
           "md-4_8_merged",
           "md-5_1_merged",
           "md-5_10_merged",
           "md-5_6_merged",
           "md-7_8_merged",
           "md-8_5_merged",
           "md-8_8_merged")]

mden.fr2<-mden.fr[c(mden.fr.gl$ind.names),]

dim(mden.fr)  #53862    27
mden.fr.gl$loc.names
mean(colMeans(mden.fr, na.rm = TRUE))
mean(colSds(mden.fr, na.rm = TRUE))

mden.fr2<-mden.fr[rowSums(is.na(mden.fr)) != ncol(mden.fr), ]
dim(mden.fr2) #53567    27

all_zero_mden <- function(r) all(row == 0)

zcav.vcfR

mden.gl<-vcfR2genlight(mden.vcfR)
zcav.gl<-vcfR2genlight(zcav.vcfR)

mden.gl #43 genotypes,  53,862 binary SNPs, size: 5.8 Mb 186168 (8.04 %) missing data
zcav.gl #88 genotypes,  149,616 binary SNPs, size: 23 Mb 1571087 (11.93 %) missing data

popNames(mden.gl)
popNames(zcav.gl)

mden.focal.gl<-gl.drop.pop(mden.gl, pop.list="Hawaii", recalc=TRUE)
zcav.focal.gl<-gl.drop.pop(zcav.gl, pop.list="Hawaii", recalc=TRUE)

mden.focal.gl #37 genotypes,  53,862 binary SNPs, size: 8.8 Mb, 142506 (7.15 %) missing data
zcav.focal.gl #84 genotypes,  149,616 binary SNPs, size: 31.6 Mb, 1472883 (11.72 %) missing data

#set wording directory for figures/results
setwd("~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/results/with_pop")
#generate and save a glplot to see if any individuals have strange genotypes
ggsave("./mden_pop_glplot.png", plot=glPlot(mden.gl, main="Blainville's, n=53,862 SNPs"))
#no weird looking individuals
ggsave("./zcav_pop_glplot.png", plot=glPlot(zcav.gl, main="Cuvier's, n=149,616 SNPs"))
#no weird looking individuals

# #removing bad glPlot individuals in dartR
# mden.gl<-gl.compliance.check(mden.gl)
# zcav.gl<-gl.compliance.check(zcav.gl)
# 
# gl.make.recode.ind(mden.gl, out.recode.file="./md_new_ind_assignments.csv", outpath="./")
# mden.remR.gl<-gl.recode.ind(mden.gl, ind.recode="./md_new_ind_assignments.csv")
# mden.remR.gl<-gl.recalc.metrics(mden.remR.gl)
# mden.remR.gl #44 genotypes,  128,044 binary SNPs, size: 20.2 Mb, 215094 (3.82 %) missing data
# mden.remR.gl<-gl.filter.monomorphs(mden.remR.gl)
# mden.remR.gl #44 genotypes,  77,273 binary SNPs, size: 12.7 Mb, 132784 (3.91 %) missing data
# ggsave("./mden_glplot_remR.png", plot=glPlot(mden.remR.gl, main="Blainville's, n=77,273 SNPs"))

#checking missing data
#filter loci by callrate

mden.gl<-gl.compliance.check(mden.gl)
zcav.gl<-gl.compliance.check(zcav.gl)

mden.focal.gl<-gl.compliance.check(mden.focal.gl)
zcav.focal.gl<-gl.compliance.check(zcav.focal.gl)

gl.filter.callrate(mden.gl, method="loc", threshold=0.95, plot=FALSE, v=3, mono.rm = TRUE) 
#43 genotypes,  37,824 binary SNPs, size: 6.2 Mb 19038 (1.17 %) missing data
gl.filter.callrate(mden.gl, method="loc", threshold=1, plot=FALSE, v=3, mono.rm = TRUE)
#43 genotypes,  23,253 binary SNPs, size: 3.9 Mb 0 (0 %) missing data

gl.filter.callrate(zcav.gl, method="loc", threshold=0.95, plot=FALSE, v=3, mono.rm=TRUE)
#88 genotypes,  91,745 binary SNPs, size: 17.8 Mb 121724 (1.51 %) missing data
gl.filter.callrate(zcav.gl, method="loc", threshold=1, plot=FALSE, v=3, mono.rm=TRUE)
#88 genotypes,  34,053 binary SNPs, size: 8 Mb 0 (0 %) missing data

#adding population data to the genlight file
#indNames(mden.remR.gl)
#pop(mden.remR.gl)<-as.factor(c("Caribbean",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Caribbean",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Caribbean",
 #                              "Canaries",
 #                              "Caribbean",
 #                              "Caribbean",
 #                              "Caribbean",
 #                              "Caribbean",
 #                              "Caribbean",
 #                              "Caribbean",
 #                              "Caribbean",
 #                              "Caribbean",
 #                              "Canaries",
 #                              "Caribbean",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Caribbean",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Canaries",
 #                              "Hawaii",
 #                              "Hawaii",
 #                              "Hawaii",
 #                              "Hawaii",
 #                              "Hawaii",
 #                              "Hawaii"))
#pop(mden.remR.gl)
indNames(mden.gl)
pop(mden.gl)<-as.factor(c("Caribbean",
                          "Canaries",
                          "Caribbean",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Caribbean",
                          "Caribbean",
                          "Caribbean",
                          "Caribbean",
                          "Caribbean",
                          "Caribbean",
                          "Caribbean",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Caribbean",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Caribbean",
                          "Canaries",
                          "Canaries",
                          "Caribbean",
                          "Caribbean",
                          "Hawaii",
                          "Hawaii",
                          "Hawaii",
                          "Hawaii",
                          "Hawaii",
                          "Hawaii"))
mden.gl
mden.pop<-cbind(indNames(mden.gl),as.character(pop(mden.gl)))
mden.pop

pop(zcav.gl)<-c("Canaries",
                          "Mediterranean-W",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Caribbean",
                          "Mediterranean-E",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Canaries",
                          "Mediterranean-E",
                          "Mediterranean-W",
                          "Canaries",
                          "Mediterranean-E",
                          "Mediterranean-E",
                          "Mediterranean-W",
                          "Mediterranean-E",
                          "Caribbean",
                          "Canaries",
                          "Caribbean",
                          "Mediterranean-W",
                          "Mediterranean-E",
                          "Caribbean",
                          "Caribbean",
                          "Caribbean",
                          "Canaries",
                          "Canaries",
                          "Mediterranean-E",
                          "Mediterranean-W",
                          "Canaries",
                          "Mediterranean-E",
                          "Canaries",
                          "Canaries",
                          "Caribbean",
                          "Mediterranean-W",
                          "Mediterranean-E",
                          "Caribbean",
                          "Mediterranean-W",
                          "Mediterranean-E",
                          "Caribbean",
                          "Caribbean",
                          "Mediterranean-W",
                          "Canaries",
                          "Mediterranean-E",
                          "Caribbean",
                          "Canaries",
                          "Mediterranean-W",
                          "Mediterranean-E",
                          "Mediterranean-W",
                          "Canaries",
                          "Mediterranean-W",
                          "Mediterranean-W",
                          "Canaries",
                          "Mediterranean-E",
                          "Mediterranean-W",
                          "Mediterranean-W",
                          "Canaries",
                          "Canaries",
                          "Caribbean",
                          "Mediterranean-W",
                          "Canaries",
                          "Mediterranean-W",
                          "Caribbean",
                          "Caribbean",
                          "Mediterranean-W",
                          "Mediterranean-W",
                          "Mediterranean-W",
                          "Caribbean",
                          "Caribbean",
                          "Mediterranean-W",
                          "Mediterranean-W",
                          "Caribbean",
                          "Caribbean",
                          "Mediterranean-E",
                          "Mediterranean-W",
                          "Caribbean",
                          "Canaries",
                          "Caribbean",
                          "Canaries",
                          "Mediterranean-W",
                          "Hawaii",
                          "Hawaii",
                          "Caribbean",
                          "Hawaii",
                          "Hawaii")
zcav.gl
zcav.pop<-cbind(indNames(zcav.gl),as.character(pop(zcav.gl)))
zcav.pop

#filter loci by callrate
# md95.gl<-gl.filter.callrate(mden.remR.gl, method="loc", threshold=0.95, plot=TRUE, v=3, mono.rm = TRUE)
# md95.gl #44 genotypes,  57,463 binary SNPs, size: 9.3 Mb, 30637 (1.21 %) missing data
# glPlot(md95.gl)
# 
# zc95.gl<-gl.filter.callrate(zcav.gl, method="loc", threshold=0.95, plot=TRUE, v=3, mono.rm=TRUE)
# zc95.gl #92 genotypes,  97,773 binary SNPs, size: 19 Mb, 166908 (1.86 %) missing data
# dev.off() 
# glPlot(zc95.gl)

#Make separate genlight file for Blainvilles in the atlantic
# mden.atl.gl<-gl.drop.pop(mden.gl, "Hawaii", mono.rm=TRUE, recalc=TRUE )
# mden.atl.gl<-gl.compliance.check(mden.atl.gl)
# mden.atl.gl$pop 

#Make separate genlight files for each of the focal regions
mden.sep.gl<-seppop(mden.gl)
zcav.sep.gl<-seppop(zcav.gl)

mden.can.gl<-mden.sep.gl$Canaries
mden.car.gl<-mden.sep.gl$Caribbean
mden.haw.gl<-mden.sep.gl$Hawaii
zcav.can.gl<-zcav.sep.gl$Canaries
zcav.car.gl<-zcav.sep.gl$Caribbean
zcav.haw.gl<-zcav.sep.gl$Hawaii
zcav.emed.gl<-zcav.sep.gl$`Mediterranean-E`
zcav.wmed.gl<-zcav.sep.gl$`Mediterranean-W`

mden.can.gl<-gl.compliance.check(mden.can.gl)
mden.car.gl<-gl.compliance.check(mden.car.gl)
mden.haw.gl<-gl.compliance.check(mden.haw.gl)
zcav.can.gl<-gl.compliance.check(zcav.can.gl)
zcav.car.gl<-gl.compliance.check(zcav.car.gl)
zcav.haw.gl<-gl.compliance.check(zcav.haw.gl)
zcav.emed.gl<-gl.compliance.check(zcav.emed.gl)
zcav.wmed.gl<-gl.compliance.check(zcav.wmed.gl)
zcav.med.gl<-gl.compliance.check(rbind(zcav.emed.gl, zcav.wmed.gl))

mden.can.gl #24 genotypes,  53,862 binary SNPs, size: 8.3 Mb 68096 (5.27 %) missing data
mden.car.gl #13 genotypes,  53,862 binary SNPs, size: 8.2 Mb 74410 (10.63 %) missing data
mden.haw.gl #6 genotypes,  53,862 binary SNPs, size: 8 Mb 43662 (13.51 %) missing data
zcav.can.gl #26 genotypes,  149,616 binary SNPs, size: 25.5 Mb 426315 (10.96 %) missing data
zcav.car.gl #21 genotypes,  149,616 binary SNPs, size: 24.6 Mb 259052 (8.24 %) missing data
zcav.haw.gl #4 genotypes,  149,616 binary SNPs, size: 23.4 Mb 98204 (16.41 %) missing data
zcav.emed.gl#14 genotypes,  149,616 binary SNPs, size: 24.6 Mb 313633 (14.97 %) missing data
zcav.wmed.gl#23 genotypes,  149,616 binary SNPs, size: 25.5 Mb 473883 (13.77 %) missing data
zcav.med.gl #37 genotypes,  149,616 binary SNPs, size: 23.8 Mb 787516 (14.23 %) missing data
popNames(zcav.med.gl)<-c("Mediterranean", "Mediterranean")

#filter datasets to not have any missing data
mden.100.can.gl <-gl.filter.callrate(mden.can.gl , method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
mden.100.car.gl <-gl.filter.callrate(mden.car.gl , method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
mden.100.haw.gl <-gl.filter.callrate(mden.haw.gl , method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.100.can.gl <-gl.filter.callrate(zcav.can.gl , method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.100.car.gl <-gl.filter.callrate(zcav.car.gl , method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.100.haw.gl <-gl.filter.callrate(zcav.haw.gl , method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.100.emed.gl<-gl.filter.callrate(zcav.emed.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.100.wmed.gl<-gl.filter.callrate(zcav.wmed.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.100.med.gl <-gl.filter.callrate(zcav.med.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)

zcav.100.gl<-gl.filter.callrate(zcav.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
mden.100.gl<-gl.filter.callrate(mden.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)

zcav.focal.100.gl<-gl.filter.callrate(zcav.focal.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
mden.focal.100.gl<-gl.filter.callrate(mden.focal.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)

mden.100.can.gl #24 genotypes,  28,314 binary SNPs, size: 4.5 Mb 0 (0 %) missing data
mden.100.car.gl #13 genotypes,  24,763 binary SNPs, size: 3.9 Mb 0 (0 %) missing data
mden.100.haw.gl #6 genotypes,  23,888 binary SNPs, size: 3.7 Mb 0 (0 %) missing data
zcav.100.can.gl #26 genotypes,  56,371 binary SNPs, size: 10.7 Mb 0 (0 %) missing data
zcav.100.car.gl #21 genotypes,  50,746 binary SNPs, size: 9.8 Mb 0 (0 %) missing data
zcav.100.haw.gl #4 genotypes,  48,479 binary SNPs, size: 9.2 Mb 0 (0 %) missing data
zcav.100.med.gl #37 genotypes,  27,097 binary SNPs, size: 3.9 Mb 0 (0 %) missing data
zcav.100.emed.gl#14 genotypes,  26,369 binary SNPs, size: 6.1 Mb 0 (0 %) missing data
zcav.100.wmed.gl#23 genotypes,  34,768 binary SNPs, size: 7.5 Mb 0 (0 %) missing data

zcav.100.gl # 88 genotypes,  34,053 binary SNPs, size: 8 Mb 0 (0 %) missing data
mden.100.gl #43 genotypes,  23,253 binary SNPs, size: 3.9 Mb 0 (0 %) missing data

zcav.focal.100.gl #84 genotypes,  34,637 binary SNPs, size: 8 Mb 0 (0 %) missing data
mden.focal.100.gl #37 genotypes,  23,167 binary SNPs, size: 3.9 Mb 0 (0 %) missing data

#Need to double check all monomorphic sites have been removed. First have to convert genlight files to genind files
 mden.100.can.gi<-gl2gi(mden.100.can.gl)
 mden.100.car.gi<-gl2gi(mden.100.car.gl)
 mden.100.haw.gi<-gl2gi(mden.100.haw.gl)
 zcav.100.can.gi<-gl2gi(zcav.100.can.gl)
 zcav.100.car.gi<-gl2gi(zcav.100.car.gl)
 zcav.100.haw.gi<-gl2gi(zcav.100.haw.gl)
 zcav.100.med.gi<-gl2gi(zcav.100.med.gl)
zcav.100.emed.gi<-gl2gi(zcav.100.emed.gl)
zcav.100.wmed.gi<-gl2gi(zcav.100.wmed.gl)
zcav.100.gi<-gl2gi(zcav.100.gl)
mden.100.gi<-gl2gi(mden.100.gl)
zcav.focal.100.gi<-gl2gi(zcav.focal.100.gl)
mden.focal.100.gi<-gl2gi(mden.focal.100.gl)

#check if there are any monomorphic sites using isPoly()
summary(isPoly(mden.100.can.gi)) #no
summary(isPoly(mden.100.car.gi)) #no
summary(isPoly(mden.100.haw.gi)) #yes
summary(isPoly(zcav.100.can.gi)) #no
summary(isPoly(zcav.100.car.gi)) #no
summary(isPoly(zcav.100.haw.gi)) #yes
summary(isPoly(zcav.100.med.gi)) #yes
summary(isPoly(zcav.100.emed.gi))#yes
summary(isPoly(zcav.100.wmed.gi))#yes
summary(isPoly(zcav.100.gi)) #yes
summary(isPoly(mden.100.gi)) #no
summary(isPoly(zcav.focal.100.gi))#yes still monomorphic sites
summary(isPoly(mden.focal.100.gi))#no more monomorphic sites

#use poppr to remove monomorphic sites since for some reason DartR didn't do it
mden.100.haw.gi<-informloci(mden.100.haw.gi, cutoff=1/nInd(mden.100.haw.gi))
zcav.100.haw.gi<-informloci(zcav.100.haw.gi, cutoff=1/nInd(zcav.100.haw.gi))
zcav.100.med.gi<-informloci(zcav.100.med.gi, cutoff=1/nInd(zcav.100.med.gi))
zcav.100.emed.gi<-informloci(zcav.100.emed.gi, cutoff=1/nInd(zcav.100.emed.gi))
zcav.100.wmed.gi<-informloci(zcav.100.wmed.gi, cutoff=1/nInd(zcav.100.wmed.gi))
zcav.100.gi<-informloci(zcav.100.gi, cutoff=1/nInd(zcav.100.gi))
mden.100.gi<-informloci(mden.100.gi, cutoff=1/nInd(mden.100.gi))
zcav.focal.100.gi<-informloci(zcav.focal.100.gi, cutoff=1/nInd(zcav.focal.100.gi))

#check that monomorphic sites have been removed
summary(isPoly(mden.100.haw.gi)) #all removed
summary(isPoly(zcav.100.haw.gi)) #all removed
summary(isPoly(zcav.100.med.gi)) #all removed
summary(isPoly(zcav.100.emed.gi))#all removed
summary(isPoly(zcav.100.wmed.gi))#all removed
summary(isPoly(zcav.100.gi)) #all removed
summary(isPoly(mden.100.gi)) #all removed
summary(isPoly(zcav.focal.100.gi))#

#convert back to genlight
mden.100.haw.gl<-gi2gl(mden.100.haw.gi)
zcav.100.haw.gl<-gi2gl(zcav.100.haw.gi)
zcav.100.med.gl<-gi2gl(zcav.100.med.gi)
zcav.100.emed.gl<-gi2gl(zcav.100.emed.gi)
zcav.100.wmed.gl<-gi2gl(zcav.100.wmed.gi)
zcav.100.gl<-gi2gl(zcav.100.gi)
mden.100.gl<-gi2gl(mden.100.gi)
zcav.focal.100.gl<-gi2gl(zcav.focal.100.gi)

mden.focal.100.gl #37 genotypes,  23,167 binary SNPs, size: 3.9 Mb 0 (0 %) missing data
zcav.focal.100.gl #84 genotypes,  34,135 binary SNPs, size: 5.4 Mb 0 (0 %) missing data


zcav.maf.gi<-informloci(zcav.100.gi, MAF=0.2)
zcav.maf.gi
zcav.maf.gl<-gi2gl(zcav.maf.gi)
zcav.maf.gl #88 genotypes,  3,426 binary SNPs

mden.maf.gi<-informloci(mden.100.gi, MAF=0.2)
mden.maf.gi
mden.maf.gl<-gi2gl(mden.maf.gi)
mden.maf.gl #43 genotypes,  7,134 binary SNPs

# exploring data with PCA/DAPC ----------------------------------------------------------
#PCA with data (following tutorial by Filip Kolář "SNP data analysis in R")
pca.1 <- glPcaFast(mden.remR.gl, nf=300)
pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis: 0.09277164
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis: 0.05615772
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis: 0.04661

col <- funky(3)
g1 <- s.class(pca.1$scores, pop(mden.remR.gl), xax=1, yax=2, col=transp(col,.6), 
              ellipseSize=0, starSize=0, ppoints.cex=4, 
              paxes.draw=T, pgrid.draw =F, plot = FALSE)
g1
g2 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels =
                 list(box = list(draw = FALSE),
                      optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
g2
ADEgS(c(g1, g2), layout = c(1, 2))

#zcav pca
pca.1 <- glPcaFast(zcav.gl, nf=300)
pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis: 0.1613515
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis: 0.03104519
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis: 0.02625066

col <- funky(3)
g1 <- s.class(pca.1$scores, pop(zcav.gl), xax=1, yax=2, col=transp(col,.6), 
              ellipseSize=0, starSize=0, ppoints.cex=4, 
              paxes.draw=T, pgrid.draw =F, plot = FALSE)
g1
g2 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels =
                 list(box = list(draw = FALSE),
                      optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
g2
ADEgS(c(g1, g2), layout = c(1, 2))

#DAPC with data (following same Tutorial)
grp <- find.clusters(mden.remR.gl, max.n.clust=20, glPca = pca.1, perc.pca = 100, n.iter=1e6, n.start=1000) #44, 2
write.table(grp$grp, file="grouping_Kmeans_all.txt", sep="\t", quote=F,
            col.names=F)
dapc1<-dapc(mden.remR.gl, grp$grp, glPca = pca.1) #44,2
scatter(dapc1)
compoplot(dapc1, cex.names = 0.4, legend=F, col=col)

grp2 <- find.clusters(mden.remR.gl, max.n.clust=20, glPca = pca.1, perc.pca = 100, n.iter=1e6, n.start=1000) #44, 3
write.table(grp2$grp2, file="grouping_Kmeans_all.txt", sep="\t", quote=F,
            col.names=F)
dapc2<-dapc(mden.remR.gl, grp2$grp2, glPca = pca.1) #44,3
scatter(dapc2)
compoplot(dapc2, cex.names = 0.4, legend=F, col=col, show.lab=TRUE, lab=pop)
pop<-pop(mden.remR.gl)

#dapc with zcav
pop<-pop(zcav.gl)
grp <- find.clusters(zcav.gl, max.n.clust=20, glPca = pca.1, perc.pca = 100, n.iter=1e6, n.start=1000) #92, 2
write.table(grp$grp, file="grouping_Kmeans_all.txt", sep="\t", quote=F,
            col.names=F)
dapc1<-dapc(zcav.gl, grp$grp, glPca = pca.1) #92,2
scatter(dapc1)
compoplot(dapc1, cex.names = 0.4, legend=F, col=col)

grp2 <- find.clusters(zcav.gl, max.n.clust=20, glPca = pca.1, perc.pca = 100, n.iter=1e6, n.start=1000) #92, 5
write.table(grp2$grp2, file="grouping_Kmeans_all.txt", sep="\t", quote=F,
            col.names=F)
dapc2<-dapc(zcav.gl, grp2$grp2, glPca = pca.1) #92,5
scatter(dapc2)
compoplot(dapc2, cex.names = 0.4, legend=F, col=col, show.lab=TRUE, lab=pop)


# exploring data with Tess3r ----------------------------------------------
mden.remR.gl
zcav.gl

#prepare coord file
mden.geo<-read.csv("./mden_coord.csv")
mden.geo<-mden.geo[1:44,]
mden.geo

zcav.geo<-read.csv("./zcav_metadata.csv")
dim(zcav.geo)
zcav.geo<-zcav.geo[1:92,]
tail(zcav.geo)

mden.coord<-mden.geo[,c(1,5,6)]
mden.coord
row.names(mden.coord)<-mden.coord$sample
mden.coord<-mden.coord[,c(2,3)]

zcav.coord<-zcav.geo[,c(1:3)]
zcav.coord
row.names(zcav.coord)<-zcav.coord$gl.filename
zcav.coord<-zcav.coord[,c(2,3)]

#check coordinates using a map
plot(x=mden.coord$long, y=mden.coord$lat, xlim=c(-180,180), ylim=c(-90,90))
map(add=T, interior=F)

mden.coord<-as.matrix(mden.coord) #must be a matrix to run in tess3r

plot(x=zcav.coord$Longitude, y=zcav.coord$Latitude, xlim=c(-180,180), ylim=c(-90,90))
map(add=T, interior=F)

zcav.coord<-as.matrix(zcav.coord)

#prepare geno file
mden.geno<-as.matrix(mden.remR.gl)
dim(mden.geno)
mden.geno[1:10,1:10]
row.names(mden.geno)

zcav.geno<-as.matrix(zcav.gl)
dim(zcav.geno)
zcav.geno[1:10,1:10]
row.names(zcav.geno)

#run Tess3r
k<-10
tess3.mden<-tess3(X=mden.geno, coord=mden.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)

k<-10
tess3.zcav<-tess3(X=zcav.geno, coord=zcav.coord, K=1:k, ploidy=2, openMP.core.num=2, rep=20, max.iteration=200, keep="best", mask=0, verbose=F)

#plot cross-entropy score
plot(tess3.mden, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations, Blainville's",
     ylab = "Cross-entropy score")

plot(tess3.zcav, crossvalid = FALSE, crossentropy = TRUE, pch = 19, col = "blue",
     xlab = "Number of ancestral populations, Cuvier's",
     ylab = "Cross-entropy score")

#generate colour pallete
library(wesanderson)
library(RColorBrewer)
my.col<-brewer.pal(10, "Paired")
my.palette<-CreatePalette(my.col)

#write q-matrix to file: Blainville's
for(i in 2:k) {
  Q.matrix <- qmatrix(tess3.mden, K = i)
  barplot(Q.matrix, sort.by.Q = TRUE, 
          border = NA, space = 0,
          col.palette = my.palette,
          xlab = "Blainville's Individuals", ylab = "Ancestry coefficients") -> bp
  axis(1, at = 1:nrow(Q.matrix), labels = row.names(mden.geno), las = 3, cex.axis = .3)
  write.table(Q.matrix, file=paste0("./tess3r/mden.qmatrix.sortedQ.k",i,".txt"))
}

#write q-matrix to file: Cuvier's
for(i in 2:k) {
  Q.matrix <- qmatrix(tess3.zcav, K = i)
  barplot(Q.matrix, sort.by.Q = TRUE, 
          border = NA, space = 0,
          col.palette = my.palette,
          xlab = "Cuvier's Individuals", ylab = "Ancestry coefficients") -> bp
  axis(1, at = 1:nrow(Q.matrix), labels = row.names(zcav.geno), las = 3, cex.axis = .3)
  write.table(Q.matrix, file=paste0("./tess3r/zcav.qmatrix.sortedQ.k",i,".txt"))
}

#generate barplots of q-matrices and write them to file: Blainville's
for(i in 2:k) {
  ppi<-300
  png(paste0("./tess3r/mden.barplots.sortedQ.k",i,".png"), width=8*ppi,height=6*ppi, res=ppi)
  Q.matrix.i <- qmatrix(tess3.mden, K = i)
  barplot(Q.matrix.i, sort.by.Q = TRUE, 
          border = NA, space = 0,
          col.palette = my.palette,
          xlab = "Individuals", ylab = "Ancestry coefficients", main=paste0("Admixture Proportions, Blainville's, k=",i,"")) -> bp.i
  axis(1, at = 1:nrow(Q.matrix.i), labels = row.names(mden.geno), las = 3, cex.axis = .3)
  dev.off()
}

#generate barplots of q-matrices and write them to file: Cuvier's
for(i in 2:k) {
  ppi<-300
  png(paste0("./tess3r/zcav.barplots.sortedQ.k",i,".png"), width=8*ppi,height=6*ppi, res=ppi)
  Q.matrix.i <- qmatrix(tess3.zcav, K = i)
  barplot(Q.matrix.i, sort.by.Q = TRUE, 
          border = NA, space = 0,
          col.palette = my.palette,
          xlab = "Individuals", ylab = "Ancestry coefficients", main=paste0("Admixture Proportions, Cuvier's, k=",i,"")) -> bp.i
  axis(1, at = 1:nrow(Q.matrix.i), labels = row.names(zcav.geno), las = 3, cex.axis = .3)
  dev.off()
}


# new code ----------------------------------------------------------------

# mden.pcoa<-gl.pcoa(mden.remR.gl, nfactors=5)
# gl.pcoa.scree(mden.pcoa)
# gl.pcoa.plot(mden.pcoa, mden.remR.gl,labels="ind" )
# gl.pcoa.plot.3d(mden.pcoa, mden.remR.gl)
# gl.dist.pop(mden.remR.gl)
# fd <- gl.fixed.diff(mden.remR.gl)
# Mygeo <- loadRDS("mden.remR.gl") 
# G <- gl.grm(mden.remR.gl) 
# 
# #md-15_10 is very much an outlier: NOAA_79840 (090602_Md1). In the global chapter, this individual was removed due to a bad glPlot. Perhaps try with just the md-15_10 (2020) sequence data. Would have to re-run gstacks
# 
# mden95.pcoa<-gl.pcoa(md95.gl, nfactors=5)
# gl.pcoa.scree(mden95.pcoa)
# gl.pcoa.plot(mden95.pcoa, md95.gl, labels="ind")
# gl.pcoa.plot.3d(mden95.pcoa, md95.gl)

#minimum spanning network: Blainville's
mden.dist <- bitwise.dist(mden.gl)
mden.msn <- poppr.msn(mden.gl, mden.dist, showplot = FALSE, include.ties = T)
node.size <- rep(2, times = nInd(mden.gl))
names(node.size) <- indNames(mden.gl)
vertex.attributes(mden.msn$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(mden.gl, mden.msn , palette = brewer.pal(n = nPop(mden.gl), name = "Dark2"), gadj = 70)

#minimum spanning network: Cuvier's
zcav.dist <- bitwise.dist(zcav.gl)
zcav.msn <- poppr.msn(zcav.gl, zcav.dist, showplot = FALSE, include.ties = T)
node.size <- rep(2, times = nInd(zcav.gl))
names(node.size) <- indNames(zcav.gl)
vertex.attributes(zcav.msn$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(zcav.gl, zcav.msn , palette = brewer.pal(n = nPop(zcav.gl), name = "Dark2"), gadj = 70)

#plot minor allele frequencies
mdenFreq<-glMean(mden.gl)
mdenFreq <- c(mdenFreq, 1-mdenFreq)
hist(mdenFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of Blainville's allele frequencies", nclass=20)
temp <- density(mdenFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)
#While a large number of loci are nearly fixed (frequencies close to 0 or 1), there is an
#appreciable number of alleles with intermediate frequencies and therefore susceptible to
#contain interesting biological signal. 

zcavFreq<-glMean(zcav.gl)
zcavFreq <- c(zcavFreq, 1-zcavFreq)
hist(zcavFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of Cuvier's allele frequencies", nclass=20)

# mden.random<-gl.subsample.loci(mden.remR.gl, 10000, method = "random", mono.rm = TRUE, verbose = 3)
# #443 genotypes,  10,000 binary SNPs, size: 2 Mb, 16581 (3.86 %) missing data
# mden.pic<-gl.subsample.loci(md.gl, 10000, method = "pic", mono.rm = TRUE, verbose = 3)
# #43 genotypes,  10,000 binary SNPs, size: 2 Mb, 0 (0 %) missing data
# 
# mden.random.pcoa<-gl.pcoa(mden.random, nfactors=5)
# gl.pcoa.plot(mden.random.pcoa, mden.random)
# 
# mden.pic.pcoa<-gl.pcoa(mden.pic, nfactors=5)
# gl.pcoa.plot(mden.pic.pcoa, mden.pic)
# 
# zcav.random<-gl.subsample.loci(zc95.gl, 10000, method = "random", mono.rm = TRUE, verbose = 3)
# #92 genotypes,  10,000 binary SNPs, size: 4.1 Mb, 17110 (1.86 %) missing data
# zcav.pic<-gl.subsample.loci(zc95.gl, 10000, method = "pic", mono.rm = TRUE, verbose = 3)
# #92 genotypes,  10,000 binary SNPs, size: 4.1 Mb, 17075 (1.86 %) missing data
# zcav.pic<-gl.subsample.loci(zc95.gl, 10000, method = "avgPIC")
# 
# 
# zcav.random.pcoa<-gl.pcoa(zcav.random, nfactors=5)
# gl.pcoa.plot(zcav.random.pcoa, zcav.random)
# 
# zcav.pic.pcoa<-gl.pcoa(zcav.pic, nfactors=5)
gl.pcoa.plot(zcav.pic.pcoa, zcav.pic)

# DAPC using xval optimisation GOOD ---------------------------------------
mden.gl
mden.atl.gl
zcav.gl

x<-mden.gl
mat<-tab(x)
grp<-pop(x)
xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
xval
scatter(xval$DAPC)
assignplot(xval$DAPC) #shows both the membership probability (red rectangle) and the pre-assigned group membership (blue cross)
#All good now!

y<-zcav.gl
mat<-tab(y)
grp<-pop(y)
yval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
yval
scatter(yval$DAPC)
compoplot(yval$DAPC, show.lab=TRUE, posi="top")
par(mar = c(4, 7, 2, 2))
assignplot(yval$DAPC, subset=45:88) #shows both the membership probability (heatmap: red rectangle is 100%) and the pre-assigned group membership (blue cross)
#all are perfect!

a<-mden.focal.100.gl
mat<-tab(a)
grp<-pop(a)
mden.xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
mden.xval
#$`Number of PCs Achieving Highest Mean Success`
#[1] "10"
#$`Number of PCs Achieving Lowest MSE`
#[1] "10"
#$call: dapc.data.frame(x = as.data.frame(x), grp = ..1, n.pca = ..2, n.da = ..3)
#$n.pca: 10 first PCs of PCA used
#$n.da: 1 discriminant functions saved
#$var (proportion of conserved variance): 0.398
scatter(mden.xval$DAPC)
assignplot(mden.xval$DAPC)

b<-zcav.focal.100.gl
mat<-tab(b)
grp<-pop(b)
zcav.xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
zcav.xval
#$`Number of PCs Achieving Highest Mean Success`
#[1] "20"
#$`Number of PCs Achieving Lowest MSE`
#[1] "10"
#$call: dapc.data.frame(x = as.data.frame(x), grp = ..1, n.pca = ..2, n.da = ..3)
#$n.pca: 10 first PCs of PCA used
#$n.da: 3 discriminant functions saved
#$var (proportion of conserved variance): 0.352
scatter(zcav.xval$DAPC)
scatter(zcav.xval$DAPC, leg=TRUE, clab=0, posi.da = "bottomleft", pch=20, cstar=0, cex=3, cell=0)
assignplot(zcav.xval$DAPC)

ggplot_shiny(dataset=yval$DAPC)
ggplot_gui(yval$DAPC)

#DAPCs with reduced datasets
#mden.random
#mden.pic
#zcav.random
#zcav.pic

#x<-mden.remR.gl
#mat<-tab(x)
#grp<-pop(x)
#xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
#xval
#scatter(xval$DAPC)
#
#x<-md95.gl
#mat<-tab(x)
#grp<-pop(x)
#xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
#xval
#scatter(xval$DAPC)
#
#x<-mden.random
#mat<-tab(x)
#grp<-pop(x)
#xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
#xval
#scatter(xval$DAPC)
#
#x<-mden.pic
#mat<-tab(x)
#grp<-pop(x)
#xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
#xval
#scatter(xval$DAPC)
#
#x<-zcav.random
#mat<-tab(x)
#grp<-pop(x)
#xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
#xval
#scatter(xval$DAPC)
#
#x<-zcav.pic
#mat<-tab(x)
#grp<-pop(x)
#xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
#xval
#scatter(xval$DAPC)
#
#x<-zcav.gl
#mat<-tab(x)
#grp<-pop(x)
#xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
#xval
#scatter(xval$DAPC)
#
#x<-zc95.gl
#mat<-tab(x)
#grp<-pop(x)
#xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
#xval
#scatter(xval$DAPC)

zcav.bah.dapc2<-xvalDapc(tab(zcav.fs.bah.gl), pop(zcav.fs.bah.gl), n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
zcav.can.dapc2<-xvalDapc(tab(zcav.fs.can.gl), pop(zcav.fs.can.gl), n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
zcav.med.dapc2<-xvalDapc(tab(zcav.fs.med.gl), pop(zcav.fs.med.gl), n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
mden.bah.dapc2<-xvalDapc(tab(mden.fs.bah.gl), pop(mden.fs.bah.gl), n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
mden.can.dapc2<-xvalDapc(tab(mden.fs.can.gl), pop(mden.fs.can.gl), n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)

scatter(zcav.bah.dapc2$DAPC)
scatter(zcav.can.dapc2$DAPC)
scatter(zcav.med.dapc2$DAPC)
scatter(mden.bah.dapc2$DAPC)
scatter(mden.can.dapc2$DAPC)

assignplot(zcav.bah.dapc2$DAPC)
assignplot(zcav.can.dapc2$DAPC)
assignplot(zcav.med.dapc2$DAPC)
assignplot(mden.bah.dapc2$DAPC)
assignplot(mden.can.dapc2$DAPC)

par(mar=c(4,5,2,2))
par(mfrow = c(5, 2))
dev.off()



#
# Demerelate Test Code --------------------------------------------------------------
#Trying out Phil Morin's Code (https://github.com/PAMorin/Demerelate_SNPs)

#sample data provided by package
data(demerelpop)
data(demerelref)
data(demereldist)
#Will Only work if you have NO NAs at all. So must filter the dataset first!
md.gl<-gl.filter.callrate(mden.remR.gl, method="loc", threshold=1, v=3, mono.rm=TRUE)
md.gl # 44 genotypes,  34,421 binary SNPs, size: 5.7 Mb, 0 (0 %) missing data
md.gl<-gl.recode.ind(md.gl, ind.recode="./md_new_ind_assignments.csv") #making sure wonky Caribbean sample is removed
md.gl<-gl.recalc.metrics(md.gl)
md.gl<-gl.filter.monomorphs(md.gl)
pop(md.gl)
#Filter 10k best SNPs based on PIC
mden.pic<-gl.subsample.loci(md.gl, 10000, method = "avgPIC", mono.rm = TRUE, verbose = 3)
#43 genotypes,  10,000 binary SNPs, size: 2 Mb, 0 (0 %) missing data
mden.pic<-gl.subsample.loci(md.gl, 10000, method = "avgPIC")
mden.pic
zcav.pic

#convert genlight files to demerelate format
mden.dr<-gl2demerelate(mden.pic, verbose=TRUE)
zcav.dr<-gl2demerelate(zcav.pic, verbose=TRUE)

dim(mden.dr) #43 20002
mden.dr[1:4,1:6]
#               Sample-ID       Population 10003:1:+_1 10003:1:+_2 100042:1:+_1 100042:1:+_2
#md-10_1_merged md-10_1_merged  Caribbean           1           1            1            1
#md-10_9_merged md-10_9_merged   Canaries          NA          NA            1            1
#md-11_9_merged md-11_9_merged   Canaries           1           1            2            2
#md-3_2_merged   md-3_2_merged   Canaries           1           1            1            1

dim(zcav.dr) #92 20002
zcav.dr[1:4,1:6]
#              Sample-ID      Population 100012:1:+_1 100012:1:+_2 100020:1:+_1 100020:1:+_2
#zc-4_1_merged zc-4_1_merged  Caribbean            1            1            1            1
#zc-4_6_merged zc-4_6_merged  Caribbean            1            1            1            1
#zc-5_4_merged zc-5_4_merged  Caribbean            1            1            1            2
#zc-1_1_merged zc-1_1_merged   Canaries            1            1            1            1

#Add lat/long data to files
mden.dr
mden.coord

mden.dr.coord<-merge(mden.coord, mden.dr, by="row.names")
dim(mden.dr.coord) #43 20005
mden.dr.coord[1:4,1:6]
##       Row.names    lat    long      Sample-ID Population 10003:1:+_1
##1     BWLib_L4-7 23.660 -76.780     BWLib_L4-7  Caribbean           1
##2        md-1_10 27.654 -18.098        md-1_10   Canaries           1
##3 md-10_1_merged 26.640 -76.910 md-10_1_merged  Caribbean           1
##4       md-10_10 27.654 -18.098       md-10_10   Canaries           1
mden.dr.coord

mden.dr.coord<-mden.dr.coord[,c(4,5,3,2)]
colnames(mden.dr.coord)<-c("Individual", "population", "X", "Y")
mden.dr.coord
colnames(mden.dr)[1:2]<-c("Individual", "population")

mden.dr<-mden.dr[order(mden.dr$`Sample-ID`),]
mden.dr.coord<-mden.dr.coord[order(mden.dr.coord$`Sample-ID`),]

is.data.frame(mden.dr.coord)
is.data.frame(mden.dr)

cbind(mden.dr$`Sample-ID`, mden.dr.coord$`Sample-ID`)

#remove the ":1:+" from each locus using stringr to see if this helps run the code without errors
loci<-colnames(mden.dr[-(1:2)])
loci
loci<-str_replace_all(loci, fixed(":1:+"), "")

colnames(mden.dr)<-c("Individual", "population", loci)
colnames(mden.dr)

is.character(mden.dr$Individual)
is.factor(mden.dr$population)
mden.dr$population<-as.factor(mden.dr$population)

#THE ISSUE IS THE SAMPLE NAMES! CAN HAVE - OR _, NEED TO REPLACE THEM ALL WITH A . AH!
ind<-mden.dr$Individual
ind
ind<-str_replace_all(ind, "-", ".")
ind<-str_replace_all(ind, "_", ".")
ind
mden.dr$Individual<-as.character(ind)
rownames(mden.dr)<-mden.dr$Individual

ind<-mden.dr.coord$Individual
ind
ind<-str_replace_all(ind, "-", ".")
ind<-str_replace_all(ind, "_", ".")
ind
mden.dr.coord$Individual<-as.character(ind)
rownames(mden.dr.coord)<-mden.dr.coord$Individual

#Run LociTest
# empirical data set # increase bootstrap (bt) to 1000 after testing; may take ≥several hours (took ~10-15 min for 10reps with 281 samples, 291 loci)
Loci.test(mden.dr, bt=1000, ref.pop = NA, object=TRUE, value = "Mxy", file.output=TRUE)

mden.D <- Demerelate(mden.dr, object = TRUE,
                Fis = TRUE,
                file.output = TRUE,
                iteration = 10,  # increase to 1000 after testing
                pairs = 10, 
                p.correct=TRUE)

Demerelate(mden.dr, tab.dist=mden.dr.coord, object=TRUE, 
           Fis=TRUE, 
           file.output=TRUE, 
           iteration=1000, 
           pairs=1000, 
           dis.data="decimal", 
           p.correct=TRUE,
           NA.rm=FALSE)
# RELATED-NOT GOOD WITH LOTS OF MARKERS -----------------------------------------------------------------
mden.pic
zcav.pic

#convert genlight files to the format most close to what is needed by RELATED
mden.df<-gl2demerelate(mden.maf.gl)
zcav.df<-gl2demerelate(zcav.pic)

#rename individuals with the population since you can't have a separate column in RELATED file
mden.df[,1:2]
mden.df<-subset(mden.df, select=-c(Population))
mden.df$`Sample-ID`<-c("CI_md-10_9_merged",
                       "CI_md-11_9_merged",
                       "CI_md-3_2_merged",
                       "CI_md-3_5_merged",
                       "CI_md-3_6_merged",
                       "CI_md-3_7_merged",
                       "CI_md-4_4_merged",
                       "CI_md-4_8_merged",
                       "CI_md-5_1_merged",
                       "CI_md-5_10_merged",
                       "CI_md-5_6_merged",
                       "CI_md-7_8_merged",
                       "CI_md-8_4_merged",
                       "CI_md-8_5_merged",
                       "CI_md-8_8_merged",
                       "CI_md-7_10_merged",
                       "CI_md-1_10",
                       "CI_md-10_10",
                       "CI_md-10_7",
                       "CI_md-11_7",
                       "CI_md-12_3",
                       "CI_md-2_5",
                       "CI_md-3_9",
                       "CI_md-4_5",
                       "CI_md-9_4",
                       "CB_md-10_1_merged",
                       "CB_md-7_2_merged",
                       "CB_md-8_7_merged",
                       "CB_md-9_5_merged",
                       "CB_md-9_9_merged",
                       "CB_md-14_3_merged",
                       "CB_md-14_7_merged",
                       "CB_md-14_9_merged",
                       "CB_md-15_8_merged",
                       "CB_md-16_3_merged",
                       "CB_BWLib_L4-7",
                       "CB_md-15_2",
                       "HW_Mde_L1-4",
                       "HW_Mde_L2-9",
                       "HW_Mde_L4-27",
                       "HW_Mde_L4-29",
                       "HW_Mde_L4-32",
                       "HW_Mde_L5-34")
mden.df$`Sample-ID`
mden.df[1:5,1:5]

mden.nonLik<-coancestry(mden.df, wang=1, lynchli = 1, quellergt = 1, ritland = 1, lynchrd = 1)
mden.nonLik$relatedness
mden.nonLik$inbreeding

#Simulating genotypes for 100 pairs. First need to input data and calculate the allele frequencies
##THIS RAN FOR 4 DAYS AND DIDNT DO ANYTHING
mden.rel<-readgenotypedata(mden.df)
mden.sim<-familysim(mden.rel$freqs, 5)
mden.sim

mden.rel$ninds
mden.rel$nloci
mden.rel$freqs

mden.compare<-compareestimators(mden.rel, 5)



# Demerelate Final Files --------------------------------------------------
mden.100.can.gl 
mden.100.car.gl 
mden.100.haw.gl 
zcav.100.can.gl 
zcav.100.car.gl 
zcav.100.haw.gl 
zcav.100.emed.gl
zcav.100.wmed.gl
zcav.100.med.gl
zcav.100.gl
mden.100.gl

popNames(zcav.100.med.gl)<-c("Mediterranean", "Mediterranean")
pop(zcav.100.med.gl)

#Convert to Demerelate files
#convert genlight files to demerelate format
mden.dr<-gl2demerelate(mden.100.gl, verbose=TRUE)
zcav.dr<-gl2demerelate(zcav.100.gl, verbose=TRUE)
zcav.maf.dr<-gl2demerelate(zcav.maf.gl, verbose=TRUE)


mden.can.dr <-gl2demerelate(mden.100.can.gl, verbose=TRUE)
mden.car.dr <-gl2demerelate(mden.100.car.gl, verbose=TRUE)
mden.haw.dr <-gl2demerelate(mden.100.haw.gl, verbose=TRUE)
zcav.can.dr <-gl2demerelate(zcav.100.can.gl, verbose=TRUE)
zcav.car.dr <-gl2demerelate(zcav.100.car.gl, verbose=TRUE)
zcav.haw.dr <-gl2demerelate(zcav.100.haw.gl, verbose=TRUE)
zcav.emed.dr<-gl2demerelate(zcav.100.emed.gl, verbose=TRUE)
zcav.wmed.dr<-gl2demerelate(zcav.100.wmed.gl, verbose=TRUE)
zcav.med.dr <-gl2demerelate(zcav.100.med.gl, verbose=TRUE)

mden.abah.dr<-gl2demerelate(mden.100.fs.abah.gl[,1:10000], verbose=TRUE)
mden.tbah.dr<-gl2demerelate(mden.100.fs.tbah.gl[,1:10000], verbose=TRUE)
mden.ecan.dr<-gl2demerelate(mden.100.fs.ecan.gl[,1:10000], verbose=TRUE)
mden.wcan.dr<-gl2demerelate(mden.100.fs.wcan.gl[,1:10000], verbose=TRUE)
zcav.abah.dr<-gl2demerelate(zcav.100.fs.abah.gl[,1:10000], verbose=TRUE)
zcav.tbah.dr<-gl2demerelate(zcav.100.fs.tbah.gl[,1:10000], verbose=TRUE)
zcav.ecan.dr<-gl2demerelate(zcav.100.fs.ecan.gl[,1:10000], verbose=TRUE)
zcav.wcan.dr<-gl2demerelate(zcav.100.fs.wcan.gl[,1:10000], verbose=TRUE)


#Remove dashes and underscores from sample names and replace with .
mden.can.dr$`Sample-ID`<- mden.can.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"="."))
mden.car.dr$`Sample-ID`<- mden.car.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"="."))
mden.haw.dr$`Sample-ID`<- mden.haw.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"="."))
zcav.can.dr$`Sample-ID`<- zcav.can.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"="."))
zcav.car.dr$`Sample-ID`<- zcav.car.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"="."))
zcav.haw.dr$`Sample-ID`<- zcav.haw.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"="."))
zcav.emed.dr$`Sample-ID`<- zcav.emed.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"="."))
zcav.wmed.dr$`Sample-ID`<- zcav.wmed.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"="."))
zcav.med.dr$`Sample-ID`<- zcav.med.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"="."))
zcav.dr$`Sample-ID`<- zcav.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"="."))
mden.dr$`Sample-ID`<- mden.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"="."))
zcav.maf.dr$`Sample-ID`<- zcav.maf.dr$`Sample-ID` %>% str_replace_all(c("-"=".", "_"="."))

mden.abah.dr$'Sample-ID'<-mden.abah.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))
mden.tbah.dr$'Sample-ID'<-mden.tbah.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))
mden.ecan.dr$'Sample-ID'<-mden.ecan.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))
mden.wcan.dr$'Sample-ID'<-mden.wcan.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))
zcav.abah.dr$'Sample-ID'<-zcav.abah.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))
zcav.tbah.dr$'Sample-ID'<-zcav.tbah.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))
zcav.ecan.dr$'Sample-ID'<-zcav.ecan.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))
zcav.wcan.dr$'Sample-ID'<-zcav.wcan.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))

mden.D <- Demerelate(mden.dr, object = TRUE,
                     Fis = TRUE,
                     file.output = TRUE,
                     iteration = 10,  # increase to 1000 after testing
                     pairs = 10, 
                     p.correct=TRUE)

#CANNOT HAVE ANY MISSING DATA OR MONOMORPHIC LOCI!

#calculate mxy
mden.can.D<-Demerelate(mden.can.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
mden.car.D<-Demerelate(mden.car.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
mden.haw.D<-Demerelate(mden.haw.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE) #WORKED!
zcav.can.D<-Demerelate(zcav.can.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
zcav.car.D<-Demerelate(zcav.car.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
zcav.haw.D<-Demerelate(zcav.haw.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE) #WORKED!
zcav.emed.D<-Demerelate(zcav.emed.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE) #WORKED!
zcav.wmed.D<-Demerelate(zcav.wmed.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE) #WORKED!
zcav.med.D<-Demerelate(zcav.med.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE) #WORKED!
zcav.D<-Demerelate(zcav.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE) #WORKED!
mden.D<-Demerelate(mden.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE) #WORKED!
zcav.maf.D<-Demerelate(zcav.maf.dr, object=TRUE, file.output=TRUE, pairs=100, p.correct=TRUE)

setwd("~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/results/with_pop/Demerelate/mxy_sepFocalSites_1000pairs")
mden.abah.D<-Demerelate(mden.abah.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
mden.tbah.D<-Demerelate(mden.tbah.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
mden.ecan.D<-Demerelate(mden.ecan.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
mden.wcan.D<-Demerelate(mden.wcan.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
zcav.abah.D<-Demerelate(zcav.abah.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
zcav.tbah.D<-Demerelate(zcav.tbah.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
zcav.ecan.D<-Demerelate(zcav.ecan.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
zcav.wcan.D<-Demerelate(zcav.wcan.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)


#calculate wang relatedness
setwd("~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/results/with_pop/Demerelate/wang_results_1000pairs/")
mden.can.wang<-Demerelate(mden.can.dr, object=TRUE, file.output=TRUE, value="wang", pairs=1000, p.correct=TRUE)
mden.car.wang<-Demerelate(mden.car.dr, object=TRUE, file.output=TRUE, value="wang", pairs=1000, p.correct=TRUE)
mden.haw.wang<-Demerelate(mden.haw.dr, object=TRUE, file.output=TRUE, value="wang", pairs=1000, p.correct=TRUE) #WORKED!
zcav.can.wang<-Demerelate(zcav.can.dr, object=TRUE, file.output=TRUE, value="wang", pairs=1000, p.correct=TRUE)
zcav.car.wang<-Demerelate(zcav.car.dr, object=TRUE, file.output=TRUE, value="wang", pairs=1000, p.correct=TRUE)
zcav.haw.wang<-Demerelate(zcav.haw.dr, object=TRUE, file.output=TRUE, value="wang", pairs=1000, p.correct=TRUE) #WORKED!
zcav.emed.wang<-Demerelate(zcav.emed.dr, object=TRUE, file.output=TRUE, value="wang", pairs=1000, p.correct=TRUE) #WORKED!
zcav.wmed.wang<-Demerelate(zcav.wmed.dr, object=TRUE, file.output=TRUE, value="wang", pairs=1000, p.correct=TRUE) #WORKED!

zcav.med.wang<-Demerelate(zcav.med.dr, object=TRUE, file.output=TRUE, value="wang", pairs=1000, p.correct=TRUE) #WORKED!
mden.wang<-Demerelate(mden.dr, object=TRUE, file.output=TRUE, value="wang", pairs=1000, p.correct=TRUE)
zcav.wang<-Demerelate(zcav.dr, object=TRUE, file.output=TRUE, value="wang", pairs=1000, p.correct=TRUE)  
zcav.maf.wang<-Demerelate(zcav.maf.dr, object=TRUE, file.output=TRUE, value="wang", pairs=100, p.correct=TRUE)

Loci.test(zcav.maf.dr, object=TRUE, bt=10, file.output=TRUE, value="Mxy")
Loci.test(zcav.maf.dr, object=TRUE, bt=10, file.output=TRUE, value="wang")


# Plotting Mxy Relatedness ------------------------------------------------
zcav.mxy<-read.csv("~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/results/with_pop/zcav_mxy.csv", header=TRUE)
head(zcav.mxy)

mden.mxy<-read.csv("~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/results/with_pop/mden_mxy.csv", header=TRUE)
head(mden.mxy)

zcav.mxy$Population<-factor(zcav.mxy$Population, levels=c("Canaries", "Caribbean", "Hawai'i","Mediterranean", "W_Mediterranean", "E_Mediterranean"))
summary(zcav.mxy$Population)

#Species	mden	mden	mden	zcav	zcav	zcav	zcav	zcav	zcav
#Pop	canaries	caribbean	hawaii	canaries	caribbean	hawaii	mediterranean	wmediterranean  emediterranean	
#HS_threshold	0.789	0.77	0.723	0.857	0.85	0.723	0.784	0.76  0.762	
#FS_threshold	0.857	0.845	0.815	0.902	0.897	0.814	0.855	0.839  0.84	
zcav.hs<-data.frame(value=c(0.857,0.85,0.723,0.784,0.76,0.762), boxplot.nr=c(1,2,3,4,5,6))
zcav.fs<-data.frame(value=c(0.902,0.897,0.814,0.855,0.839,0.84), boxplot.nr=c(1,2,3,4,5,6))
mden.hs<-data.frame(value=c(0.789,0.77,0.723), boxplot.nr=c(1,2,3))
mden.fs<-data.frame(value=c(0.857,0.845,0.815),boxplot.nr=c(1,2,3))

ggplot(zcav.mxy, aes(x=Population, y=Mxy)) +
  geom_boxplot() +
  theme_classic() +
  labs(x=NULL, y="Pairwise Relatedness (Mxy)") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  geom_segment(data=zcav.hs,aes(x=boxplot.nr-0.45,xend=boxplot.nr+0.45,
                                y=value,yend=value),inherit.aes=FALSE,size=1.5, linetype=2) +
  geom_segment(data=zcav.fs,aes(x=boxplot.nr-0.45,xend=boxplot.nr+0.45,
                                y=value,yend=value),inherit.aes=FALSE,size=1.5)

ggplot(mden.mxy, aes(x=Population, y=Mxy)) +
  geom_boxplot() +
  theme_classic() +
  labs(x=NULL, y="Pairwise Relatedness (Mxy)") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  geom_segment(data=mden.hs,aes(x=boxplot.nr-0.45,xend=boxplot.nr+0.45,
                                y=value,yend=value),inherit.aes=FALSE,size=1.5, linetype=2) +
  geom_segment(data=mden.fs,aes(x=boxplot.nr-0.45,xend=boxplot.nr+0.45,
                                y=value,yend=value),inherit.aes=FALSE,size=1.5)


hist(zcav.mxy$Mxy)
hist(mden.mxy$Mxy)

# Networks ----------------------------------------------------------------
head(zcav.mxy)
head(mden.mxy)

zcav.mxy<-zcav.mxy[,1:4]
mden.mxy<-mden.mxy[,1:4]
  
zcav.meta<-read.csv("zcav_mxy_meta.csv", header=TRUE)
mden.meta<-read.csv("mden_mxy_meta.csv", header=TRUE)

head(zcav.meta)
head(mden.meta)

#nodes table: list of individuals and any other attributes as columns (also vertices)
zc.nodes<-zcav.meta[,c(11,2,7,9)] #filename, sex, focal_site, focal_region
zc.nodes$Filename<- zc.nodes$Filename %>% str_replace_all(c("-"=".", "_"="."))
head(zc.nodes) 
tail(zc.nodes) 
zc.nodes<-zc.nodes[1:88,]
is.data.frame(zc.nodes)

head(mden.meta)
md.nodes<-mden.meta[,c(10,2,8,7)]
md.nodes$Filename<- md.nodes$Filename %>% str_replace_all(c("-"=".", "_"="."))
head(md.nodes)
tail(md.nodes)
is.data.frame(md.nodes)

#links: starting node, ending node, link attribute (also d)
zc.links<-zcav.mxy[,2:4]
head(zc.links)
tail(zc.links)
is.data.frame(zc.links)

zc.emed.links<-subset(zcav.mxy, zcav.mxy$Population=="E_Mediterranean")
zc.emed.links
zc.emed.links<-zc.emed.links[,2:4]
head(zc.emed.links)
tail(zc.emed.links)

zc.wmed.links<-subset(zcav.mxy, zcav.mxy$Population=="W_Mediterranean")
zc.wmed.links
zc.wmed.links<-zc.wmed.links[,2:4]
head(zc.wmed.links)
tail(zc.wmed.links)

md.links<-mden.mxy[,2:4]
head(md.links)
tail(md.links)
is.data.frame(md.links)

#Generating network for Zcav Canaries
zc.can.nodes<-subset(zc.nodes, zc.nodes$Focal_Site=="Canaries")
zc.can.nodes<-zc.can.nodes[order(zc.can.nodes$Filename),]
zc.can.nodes
zc.can.mxy<-subset(zcav.mxy, zcav.mxy$Population=="Canaries")
zc.can.links<-zc.can.mxy[,2:4]
#Build network
zc.can.net<-graph.data.frame(zc.can.links,directed=FALSE)
get.adjacency(zc.can.net, attr="Mxy", spars=FALSE)
#Add Focal Region as attribute
V(zc.can.net)$Focal_Region=as.character(zc.can.nodes$Focal_Region[match(V(zc.can.net)$name,zc.can.nodes$Filename)])
#Delete links that don't indicate sibship and color by relatedness
delete_edges(zc.can.net, E(zc.can.net)[Mxy<0.857])
E(zc.can.net)[Mxy>0.857]$color <- "green"
E(zc.can.net)[Mxy>0.902]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.can.net)$Focal_Region))]
plot(zc.can.net, vertex.color=my_color, layout=layout_with_graphopt(zc.can.net, niter=10000))
plot(zc.can.net, vertex.color=my_color, layout=layout_in_circle(zc.can.net))
legend("bottomleft", legend=levels(as.factor(V(zc.can.net)$Focal_Region)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

#Generating network for Zcav Caribbean
zc.car.nodes<-subset(zc.nodes, zc.nodes$Focal_Site=="Caribbean")
zc.car.nodes<-zc.car.nodes[order(zc.car.nodes$Filename),]
zc.car.nodes
zc.car.mxy<-subset(zcav.mxy, zcav.mxy$Population=="Caribbean")
zc.car.links<-zc.car.mxy[,2:4]
#Build network
zc.car.net<-graph.data.frame(zc.car.links,directed=FALSE)
get.adjacency(zc.car.net, attr="Mxy", spars=FALSE)
#Add Focal Region as attribute
V(zc.car.net)$Focal_Region=as.character(zc.car.nodes$Focal_Region[match(V(zc.car.net)$name,zc.car.nodes$Filename)])
#Delete links that don't indicate sibship and color by relatedness
delete_edges(zc.car.net, E(zc.car.net)[Mxy<0.85])
E(zc.car.net)[Mxy>0.85]$color <- "green"
E(zc.car.net)[Mxy>0.897]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.car.net)$Focal_Region))]
plot(zc.car.net, vertex.color=my_color, layout=layout_with_graphopt(zc.car.net, niter=1000))
plot(zc.car.net, vertex.color=my_color, layout=layout_in_circle(zc.car.net))
legend("bottomleft", legend=levels(as.factor(V(zc.car.net)$Focal_Region)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

#Generating network for Zcav Hawai'i
zc.haw.nodes<-subset(zc.nodes, zc.nodes$Focal_Site=="Hawai'i")
zc.haw.nodes<-zc.haw.nodes[order(zc.haw.nodes$Filename),]
zc.haw.nodes
zc.haw.mxy<-subset(zcav.mxy, zcav.mxy$Population=="Hawai'i")
zc.haw.links<-zc.haw.mxy[,2:4]
#Build network
zc.haw.net<-graph.data.frame(zc.haw.links,directed=FALSE)
get.adjacency(zc.haw.net, attr="Mxy", spars=FALSE)
#Add Focal Region as attribute
V(zc.haw.net)$Focal_Region=as.character(zc.haw.nodes$Focal_Region[match(V(zc.haw.net)$name,zc.haw.nodes$Filename)])
#Delete links that don't indicate sibship and color by relatedness
delete_edges(zc.haw.net, E(zc.haw.net)[Mxy<0.723])
E(zc.haw.net)[Mxy>0.723]$color <- "green"
E(zc.haw.net)[Mxy>0.814]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.haw.net)$Focal_Region))]
plot(zc.haw.net, vertex.color=my_color, layout=layout_in_circle(zc.haw.net))
legend("bottomleft", legend=levels(as.factor(V(zc.haw.net)$Focal_Region)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

#Generating network for Zcav Mediterranean
zc.med.nodes<-subset(zc.nodes, zc.nodes$Focal_Site=="Mediterranean")
zc.med.nodes<-zc.med.nodes[order(zc.med.nodes$Filename),]
zc.med.nodes
zc.med.mxy<-subset(zcav.mxy, zcav.mxy$Population=="Mediterranean")
zc.med.links<-zc.med.mxy[,2:4]
#Build network
zc.med.net<-graph.data.frame(zc.med.links,directed=FALSE)
get.adjacency(zc.med.net, attr="Mxy", spars=FALSE)
#Add Focal Region as attribute
V(zc.med.net)$Focal_Region=as.character(zc.med.nodes$Focal_Region[match(V(zc.med.net)$name,zc.med.nodes$Filename)])
#Delete links that don't indicate sibship and color by relatedness
delete_edges(zc.med.net, E(zc.med.net)[Mxy<0.784])
E(zc.med.net)[Mxy>0.784]$color <- "green"
E(zc.med.net)[Mxy>0.855]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.med.net)$Focal_Region))]
plot(zc.med.net, vertex.color=my_color, layout=layout_in_circle(zc.med.net))
legend("bottomleft", legend=levels(as.factor(V(zc.med.net)$Focal_Region)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

#Generating network for Zcav Mediterranean-East
zc.emed.nodes<-subset(zc.nodes, zc.nodes$Focal_Region=="Mediterranean-East" | zc.nodes$Focal_Region=="Corfu")
zc.emed.nodes
zc.emed.nodes<-zc.emed.nodes[order(zc.emed.nodes$Filename),]
zc.emed.nodes
zc.emed.mxy<-subset(zcav.mxy, zcav.mxy$Population=="E_Mediterranean")
zc.emed.mxy
zc.emed.links<-zc.emed.mxy[,2:4]
#Build network
zc.emed.net<-graph.data.frame(zc.emed.links,directed=FALSE)
get.adjacency(zc.emed.net, attr="Mxy", spars=FALSE)
#Add Focal Region as attribute
V(zc.emed.net)$Focal_Region=as.character(zc.emed.nodes$Focal_Region[match(V(zc.emed.net)$name,zc.emed.nodes$Filename)])
#Delete links that don't indicate sibship and color by relatedness
delete_edges(zc.emed.net, E(zc.emed.net)[Mxy<0.762])
E(zc.emed.net)[Mxy>0.762]$color <- "green"
E(zc.emed.net)[Mxy>0.840]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.emed.net)$Focal_Region))]
plot(zc.emed.net, vertex.color=my_color, layout=layout_in_circle(zc.emed.net))
legend("bottomleft", legend=levels(as.factor(V(zc.emed.net)$Focal_Region)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

#Generating network for Zcav Mediterranean-West
zc.wmed.nodes<-subset(zc.nodes, zc.nodes$Focal_Region=="Mediterranean-West")
zc.wmed.nodes
zc.wmed.nodes<-zc.wmed.nodes[order(zc.wmed.nodes$Filename),]
zc.wmed.nodes
zc.wmed.mxy<-subset(zcav.mxy, zcav.mxy$Population=="W_Mediterranean")
zc.wmed.mxy
zc.wmed.links<-zc.wmed.mxy[,2:4]
#Build network
zc.wmed.net<-graph.data.frame(zc.wmed.links,directed=FALSE)
get.adjacency(zc.wmed.net, attr="Mxy", spars=FALSE)
#Add Focal Region as attribute
V(zc.wmed.net)$Focal_Region=as.character(zc.wmed.nodes$Focal_Region[match(V(zc.wmed.net)$name,zc.wmed.nodes$Filename)])
#Delete links that don't indicate sibship and color by relatedness
delete_edges(zc.wmed.net, E(zc.wmed.net)[Mxy<0.76])
E(zc.wmed.net)[Mxy>0.76]$color <- "green"
E(zc.wmed.net)[Mxy>0.839]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.wmed.net)$Focal_Region))]
plot(zc.wmed.net, vertex.color=my_color, layout=layout_in_circle(zc.wmed.net))
legend("bottomleft", legend=levels(as.factor(V(zc.wmed.net)$Focal_Region)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")


#Generating network for Mden Canaries
md.can.nodes<-subset(md.nodes, md.nodes$Focal_Site=="Canaries")
md.can.nodes<-md.can.nodes[order(md.can.nodes$Filename),]
md.can.nodes
md.can.mxy<-subset(mden.mxy, mden.mxy$Population=="Canaries")
md.can.links<-md.can.mxy[,2:4]
#Build network
md.can.net<-graph.data.frame(md.can.links,directed=FALSE)
get.adjacency(md.can.net, attr="Mxy", spars=FALSE)
#Add Focal Region as attribute
V(md.can.net)$Focal_Region=as.character(md.can.nodes$Focal_Region[match(V(md.can.net)$name,md.can.nodes$Filename)])
V(md.can.net)$Focal_Region
#Delete links that don't indicate sibship and color by relatedness
delete_edges(md.can.net, E(md.can.net)[Mxy<0.789])
E(md.can.net)[Mxy>0.789]$color <- "green"
E(md.can.net)[Mxy>0.857]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(md.can.net)$Focal_Region))]
plot(md.can.net, vertex.color=my_color, layout=layout_in_circle(md.can.net))
legend("bottomleft", legend=levels(as.factor(V(md.can.net)$Focal_Region)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

#Generating network for Mden Caribbean
md.car.nodes<-subset(md.nodes, md.nodes$Focal_Site=="Caribbean")
md.car.nodes<-md.car.nodes[order(md.car.nodes$Filename),]
md.car.nodes
md.car.mxy<-subset(mden.mxy, mden.mxy$Population=="Caribbean")
md.car.links<-md.car.mxy[,2:4]
#Build network
md.car.net<-graph.data.frame(md.car.links,directed=FALSE)
get.adjacency(md.car.net, attr="Mxy", spars=FALSE)
#Add Focal Region as attribute
V(md.car.net)$Focal_Region=as.character(md.car.nodes$Focal_Region[match(V(md.car.net)$name,md.car.nodes$Filename)])
V(md.car.net)$name
#Delete links that don't indicate sibship and color by relatedness
delete_edges(md.car.net, E(md.car.net)[Mxy<0.77])
E(md.car.net)[Mxy>0.77]$color <- "green"
E(md.car.net)[Mxy>0.845]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(md.car.net)$Focal_Region))]
plot(md.car.net, vertex.color=my_color, layout=layout_in_circle(md.car.net))
legend("bottomleft", legend=levels(as.factor(V(md.car.net)$Focal_Region)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

#Generating network for Mden Hawaii
md.haw.nodes<-subset(md.nodes, md.nodes$Focal_Site=="Hawaii")
md.haw.nodes<-md.haw.nodes[order(md.haw.nodes$Filename),]
md.haw.nodes
md.haw.mxy<-subset(mden.mxy, mden.mxy$Population=="Hawai'i")
md.haw.links<-md.haw.mxy[,2:4]
#Build network
md.haw.net<-graph.data.frame(md.haw.links,directed=FALSE)
get.adjacency(md.haw.net, attr="Mxy", spars=FALSE)
#Add Focal Region as attribute
V(md.haw.net)$Focal_Region=as.character(md.haw.nodes$Focal_Region[match(V(md.haw.net)$name,md.car.nodes$Filename)])
V(md.haw.net)$Focal_Region
#Delete links that don't indicate sibship and color by relatedness
delete_edges(md.haw.net, E(md.haw.net)[Mxy<0.723])
E(md.haw.net)[Mxy>0.723]$color <- "green"
E(md.haw.net)[Mxy>0.815]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(md.haw.net)$Focal_Region))]
plot(md.haw.net, vertex.color=my_color, layout=layout_in_circle(md.haw.net))
legend("bottomleft", legend=levels(as.factor(V(md.haw.net)$Focal_Region)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

# Diversity stats ---------------------------------------------------------
mden.gl
mden.can.gl
mden.car.gl
mden.haw.gl

zcav.gl
zcav.can.gl
zcav.car.gl
zcav.haw.gl
zcav.wmed.gl
zcav.emed.gl
zcav.med.gl

zcav.100.fs.wmed.gi
zcav.100.fs.emed.gi
zcav.fs.bah.gi
zcav.fs.can.gi
mden.fs.bah.gi
mden.fs.can.gi

mden.gi<-gl2gi(mden.gl)
mden.sep.gi<-seppop(mden.gi)
mden.can.gi<-mden.sep.gi$Canaries
mden.car.gi<-mden.sep.gi$Caribbean
mden.haw.gi<-mden.sep.gi$Hawaii

zcav.gi<-gl2gi(zcav.gl)
zcav.sep.gi<-seppop(zcav.gi)
zcav.can.gi <-zcav.sep.gi$Canaries
zcav.car.gi <-zcav.sep.gi$Caribbean
zcav.haw.gi <-zcav.sep.gi$Hawaii
zcav.wmed.gi<-zcav.sep.gi$`Mediterranean-W`
zcav.emed.gi<-zcav.sep.gi$`Mediterranean-E`
zcav.med.gi<-repool(zcav.wmed.gi, zcav.emed.gi)
popNames(zcav.med.gi)<-c("Mediterranean","Mediterranean")
zcav.med.gi


mden.stats<-basic.stats(mden.gi, digits=3)
zcav.stats<-basic.stats(zcav.gi, digits=3)
zcav.med.stats<-basic.stats(zcav.med.gi, digits=3)

summary(mden.stats$Ho)
summary(mden.stats$Hs)
summary(mden.stats$Fis)

summary(zcav.stats$Ho)
summary(zcav.stats$Hs)
summary(zcav.stats$Fis)

summary(zcav.med.stats$Ho)
summary(zcav.med.stats$Hs)
summary(zcav.med.stats$Fis)

mden.can.fis<-wc(mden.can.gi)
mden.car.fis<-wc(mden.car.gi)
mden.haw.fis<-wc(mden.haw.gi)

mden.can.fis
mden.car.fis
mden.haw.fis

zcav.can.fis<-wc(zcav.can.gi)
zcav.car.fis<-wc(zcav.car.gi)
zcav.haw.fis<-wc(zcav.haw.gi)
zcav.wmed.fis<-wc(zcav.wmed.gi)
zcav.emed.fis<-wc(zcav.emed.gi)
zcav.med.fis<-wc(zcav.med.gi)

zcav.can.fis
zcav.car.fis
zcav.haw.fis
zcav.med.fis
zcav.emed.fis
zcav.wmed.fis

mden.boot.fis<-boot.ppfis(mden.gi, nboot=1000)
zcav.boot.fis<-boot.ppfis(zcav.gi, nboot=1000)
zcav.med.boot.fis<-boot.ppfis(zcav.med.gi, nboot=1000)

mden.boot.fis
zcav.boot.fis
zcav.med.boot.fis

zcav.100.gi
mden.100.gi

zcav.100.fs.wmed.gi
zcav.100.fs.emed.gi
zcav.fs.bah.gi
zcav.fs.can.gi
mden.fs.bah.gi
mden.fs.can.gi

zcav.wmed.stats<-basic.stats(zcav.100.fs.wmed.gi, digits=3)
zcav.emed.stats<-basic.stats(zcav.100.fs.emed.gi, digits=3)
zcav.bah.stats<-basic.stats(zcav.fs.bah.gi, digits=3)
zcav.can.stats<-basic.stats(zcav.fs.can.gi, digits=3)
mden.bah.stats<-basic.stats(mden.fs.bah.gi, digits=3)
mden.can.stats<-basic.stats(mden.fs.can.gi, digits=3)
zcav.med.stats<-basic.stats(zcav.100.med.gi, digits=3)

summary(zcav.wmed.stats$Ho)
summary(zcav.emed.stats$Ho)
summary(zcav.bah.stats$Ho)
summary(zcav.can.stats$Ho)
summary(mden.bah.stats$Ho)
summary(mden.can.stats$Ho)

summary(zcav.wmed.stats$Hs)
summary(zcav.emed.stats$Hs)
summary(zcav.bah.stats$Hs)
summary(zcav.can.stats$Hs)
summary(mden.bah.stats$Hs)
summary(mden.can.stats$Hs)

zcav.100.wmed.gi
zcav.100.emed.gi
zcav.fs.bah.gi$pop
zcav.fs.can.gi$pop
mden.fs.bah.gi$pop
mden.fs.can.gi$pop

zcav.sep.bah.gi<-seppop(zcav.fs.bah.gi)
zcav.sep.can.gi<-seppop(zcav.fs.can.gi)
mden.sep.bah.gi<-seppop(mden.fs.bah.gi)
mden.sep.can.gi<-seppop(mden.fs.can.gi)
zcav.sep.med.gi<-seppop(zcav.100.med.gi)

wc(zcav.sep.bah.gi$`Bahamas-Abaco`) #0.122886
wc(zcav.sep.bah.gi$`Bahamas-TOTO`) #0.1172785
wc(zcav.sep.can.gi$`Canary Islands-East`) #0.1144374
wc(zcav.sep.can.gi$`Canary Islands-West`) #0.0646678
wc(mden.sep.bah.gi$`Bahamas-Abaco`) #0.07898718
wc(mden.sep.bah.gi$`Bahamas-TOTO`) #0.02882785
wc(mden.sep.can.gi$`Canary Islands-East`) #0.04911995
wc(mden.sep.can.gi$`Canary Islands-West`) #0.03869325
wc(zcav.100.wmed.gi) #0.03425586
wc(zcav.100.emed.gi) #0.04493555
wc(zcav.sep.med.gi$`Mediterranean-E`) #0.03224651
wc(zcav.sep.med.gi$`Mediterranean-W`) #0.02726726

boot.ppfis(zcav.100.med.gi, nboot=1000) #1 0.0277 0.0364; 2 0.0239 0.0305
boot.ppfis(zcav.100.fs.wmed.gi, nboot=1000) #0.0313 0.0369
boot.ppfis(zcav.100.fs.emed.gi, nboot=1000) #0.0407 0.0488
boot.ppfis(zcav.fs.bah.gi, nboot=1000) #1 0.1189 0.1266; 2 0.1097 0.1250
boot.ppfis(zcav.fs.can.gi, nboot=1000) #1 0.1091 0.1194; 2 0.0618 0.0680
boot.ppfis(mden.fs.bah.gi, nboot=1000) #1 0.0709 0.0876; 2 0.0224 0.0359
boot.ppfis(mden.fs.can.gi, nboot=1000) #1 0.0404 0.0574; 2 0.0352 0.0425

#confidence intervals for ho and hs
library(rcompanion)

zczc125.ho.melt<-melt(zczc125.ho)
zczc125.ho.melt<-zczc125.ho.melt[,2:3]
colnames(zczc125.ho.melt)<-c("pop", "ho")
zczc125.ho.melt

zczc125.ho.melt<-na.omit(zczc125.ho.melt)

Sum.zczc125.pop1.ho = groupwiseMean(ho ~ pop, data= zczc125.pop1.ho, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
Sum.zczc125.pop1.ho

zcav.wmed<-zcav.wmed.stats$perloc
zcav.wmed<-zcav.wmed[,1:2]
zcav.wmed<-na.omit(zcav.wmed)
zcav.wmed<-melt(zcav.wmed)
head(zcav.wmed)

zcav.emed.stats
zcav.emed<-zcav.emed.stats$perloc
zcav.emed<-zcav.emed[,1:2]
zcav.emed<-na.omit(zcav.emed)
zcav.emed<-melt(zcav.emed)
head(zcav.emed)

zcav.med.stats
zcav.med.ho<-zcav.med.stats$Ho
zcav.med.ho<-melt(zcav.med.ho)
zcav.med.hs<-zcav.med.stats$Hs
zcav.med.hs<-melt(zcav.med.hs)
zcav.med.fis<-zcav.med.stats$Fis
zcav.med.fis<-melt(zcav.med.fis)

zcav.bah.stats
zcav.bah.ho<-zcav.bah.stats$Ho
zcav.bah.ho<-na.omit(zcav.bah.ho)
head(zcav.bah.ho)
zcav.bah.ho<-melt(zcav.bah.ho)
head(zcav.bah.ho)

zcav.bah.hs<-zcav.bah.stats$Hs
zcav.bah.hs<-na.omit(zcav.bah.hs)
head(zcav.bah.hs)
zcav.bah.hs<-melt(zcav.bah.hs)
head(zcav.bah.hs)

zcav.can.stats
mden.can.ho<-zcav.can.stats$Ho
zcav.can.ho<-na.omit(zcav.can.ho)
head(zcav.can.ho)
zcav.can.ho<-melt(zcav.can.ho)
head(zcav.can.ho)

zcav.can.hs<-zcav.can.stats$Hs
zcav.can.hs<-na.omit(zcav.can.hs)
head(zcav.can.hs)
zcav.can.hs<-melt(zcav.can.hs)
head(zcav.can.hs)


mden.bah.stats
mden.bah.ho<-mden.bah.stats$Ho
mden.bah.ho<-na.omit(mden.bah.ho)
head(mden.bah.ho)
mden.bah.ho<-melt(mden.bah.ho)
head(mden.bah.ho)

mden.bah.hs<-mden.bah.stats$Hs
mden.bah.hs<-na.omit(mden.bah.hs)
head(mden.bah.hs)
mden.bah.hs<-melt(mden.bah.hs)
head(mden.bah.hs)

mden.can.stats
mden.can.ho<-mden.can.stats$Ho
mden.can.ho<-na.omit(mden.can.ho)
head(mden.can.ho)
mden.can.ho<-melt(mden.can.ho)
head(mden.can.ho)

mden.can.hs<-mden.can.stats$Hs
mden.can.hs<-na.omit(mden.can.hs)
head(mden.can.hs)
mden.can.hs<-melt(mden.can.hs)
head(mden.can.hs)

groupwiseMean(value ~ Var2, data=zcav.med.ho, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, na.rm = TRUE)
#Var2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Mediterranean-E 27087 0.222       0.95      0.220      0.224            0.219            0.224
#2 Mediterranean-W 27087 0.259       0.95      0.257      0.261            0.257            0.261
groupwiseMean(value ~ Var2, data=zcav.med.hs, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, na.rm = TRUE)
#Var2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Mediterranean-E 27087 0.229       0.95      0.227      0.232            0.227            0.231
#2 Mediterranean-W 27087 0.266       0.95      0.264      0.268            0.264            0.268
groupwiseMean(value ~ Var2, data=zcav.med.fis, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, na.rm = TRUE)
#Var2     n   Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Mediterranean-E 20805 0.0250       0.95     0.0215     0.0284           0.0217           0.0285
#2 Mediterranean-W 25070 0.0262       0.95     0.0234     0.0290           0.0231           0.0288

groupwiseMean(value ~ variable, data=zcav.wmed, conf=0.95, R=1000, percentile=TRUE, bca=FALSE)
#variable     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1       Ho 34741 0.287       0.95      0.285      0.288            0.285            0.288
#2       Hs 34741 0.297       0.95      0.295      0.299            0.295            0.299
groupwiseMean(value ~ variable, data=zcav.emed, conf=0.95, R=1000, percentile=TRUE, bca=FALSE)
#variable     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1       Ho 26269 0.285       0.95      0.283      0.287            0.283            0.287
#2       Hs 26269 0.298       0.95      0.296      0.300            0.297            0.300
groupwiseMean(value ~ X2, data=zcav.bah.ho, conf=0.95, R=1000, percentile=TRUE, bca=FALSE)
#X2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Bahamas-Abaco 61643 0.241       0.95      0.239      0.242            0.239            0.242
#2  Bahamas-TOTO 61643 0.243       0.95      0.241      0.246            0.241            0.246
groupwiseMean(value ~ X2, data=zcav.bah.hs, conf=0.95, R=1000, percentile=TRUE, bca=FALSE)
#X2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Bahamas-Abaco 61643 0.275       0.95      0.273      0.276            0.273            0.276
#2  Bahamas-TOTO 61643 0.276       0.95      0.273      0.278            0.273            0.278
groupwiseMean(value ~ X2, data=zcav.can.ho, conf=0.95, R=1000, percentile=TRUE, bca=FALSE)
#X2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Canary Islands-East 56246 0.158       0.95      0.157       0.16            0.157             0.16
#2 Canary Islands-West 56246 0.169       0.95      0.168       0.17            0.168             0.17
groupwiseMean(value ~ X2, data=zcav.can.hs, conf=0.95, R=1000, percentile=TRUE, bca=FALSE)
#X2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Canary Islands-East 56246 0.179       0.95      0.177      0.180            0.177            0.180
#2 Canary Islands-West 56246 0.181       0.95      0.180      0.182            0.180            0.182
groupwiseMean(value ~ X2, data=mden.bah.ho, conf=0.95, R=1000, percentile=TRUE, bca=FALSE)
#X2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Bahamas-Abaco 26890 0.313       0.95      0.310      0.316            0.310            0.316
#2  Bahamas-TOTO 26890 0.312       0.95      0.309      0.315            0.309            0.315
groupwiseMean(value ~ X2, data=mden.bah.hs, conf=0.95, R=1000, percentile=TRUE, bca=FALSE)
#X2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Bahamas-Abaco 26890 0.340       0.95      0.337      0.343            0.337            0.343
#2  Bahamas-TOTO 26890 0.321       0.95      0.318      0.323            0.318            0.323
groupwiseMean(value ~ X2, data=mden.can.ho, conf=0.95, R=1000, percentile=TRUE, bca=FALSE)
#X2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Canary Islands-East 28911 0.256       0.95      0.253      0.259            0.253            0.259
#2 Canary Islands-West 28911 0.255       0.95      0.253      0.257            0.253            0.257
groupwiseMean(value ~ X2, data=mden.can.hs, conf=0.95, R=1000, percentile=TRUE, bca=FALSE)
#X2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Canary Islands-East 28911 0.269       0.95      0.266      0.272            0.267            0.272
#2 Canary Islands-West 28911 0.265       0.95      0.263      0.267            0.263            0.267
# PopGraph ----------------------------------------------------------------

install.packages( c("RgoogleMaps",
                    "geosphere",
                    "proto",
                    "sampling",
                    "seqinr",
                    "spacetime",
                    "spdep"), 
                  dependencies=TRUE )

library(devtools)
install_github("dyerlab/popgraph")
install_github("dyerlab/gstudio")

library(popgraph)
library(gstudio)

#data in package
data(arapat)
##   Population Latitude Longitude snp1 snp2 snp3 snp4
## 1    Olympia    47.15   -122.89  A:B  B:B  A:B  A:B
## 2 Bellingham    48.75   -122.49  A:B  A:A  A:B  B:B
## 3  St. Louis    38.81    -89.98  A:B  A:A          
## 4       Ames    42.26    -93.47  A:A  A:B          
## 5   Richmond    37.74    -77.16  A:A  B:B  A:B  B:B

file<-as.matrix(zcav.100.gl)
file<-as.data.frame(file)
dim(file)

colnames(zcav.meta)
file.meta<-cbind(zcav.meta$Filename, zcav.meta$Focal_Region, zcav.meta$Focal_Site, zcav.meta$Latitude, zcav.meta$Longitude)
file.meta<-file.meta[1:88,]
file.meta<-as.data.frame(file.meta)
file.meta

file.meta<-arrange(file.meta, file.meta$V1)

test<-cbind(file.meta$V1, rownames(file))
test
test<-merge(file.meta, file)
test[1:10,1:10]
colnames((test)[,1:5])<-c("sample", "focal_site", "focal_region", "lat", "long")


write.csv(test, file="./zcav.txt")

data <- read_population("./zcav.txt", type = "snp", header=TRUE, locus.columns = 6:ncol(test))
dim(data)
summary(data)



# Demerelate just focal sites ------------------------------------------
zcav.gl
mden.gl

zcav<-cbind(indNames(zcav.gl), as.character(pop(zcav.gl)))
zcav

zcav.fs.nohaw.gl<-zcav.gl
zcav.fs.nohaw.gl<-gl.drop.pop(zcav.fs.nohaw.gl, pop.list="Hawaii", recalc=TRUE)
zcav.fs.nohaw.gl<-gl.compliance.check(zcav.fs.nohaw.gl)
zcav.fs.nohaw.gl<-gl.filter.monomorphs(zcav.fs.nohaw.gl)
zcav.fs.nohaw.gl

pop(zcav.fs.gl)<-c("Canary Islands-East",
                   "Mediterranean-West",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Canary Islands-East",
                   "Remove",
                   "Mediterranean-East",
                   "Canary Islands-West",
                   "Remove",
                   "Canary Islands-East",
                   "Canary Islands-East",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Mediterranean-East",
                   "Mediterranean-West",
                   "Remove",
                   "Mediterranean-East",
                   "Mediterranean-East",
                   "Mediterranean-West",
                   "Mediterranean-East",
                   "Remove",
                   "Canary Islands-East",
                   "Bahamas-Abaco",
                   "Mediterranean-West",
                   "Mediterranean-East",
                   "Remove",
                   "Bahamas-Abaco",
                   "Remove",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Mediterranean-East",
                   "Mediterranean-West",
                   "Canary Islands-West",
                   "Mediterranean-East",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Remove",
                   "Mediterranean-West",
                   "Mediterranean-East",
                   "Bahamas-Abaco",
                   "Mediterranean-West",
                   "Mediterranean-East",
                   "Bahamas-Abaco",
                   "Bahamas-Abaco",
                   "Mediterranean-West",
                   "Canary Islands-West",
                   "Mediterranean-East",
                   "Remove",
                   "Canary Islands-West",
                   "Mediterranean-West",
                   "Mediterranean-East",
                   "Mediterranean-West",
                   "Remove",
                   "Mediterranean-West",
                   "Mediterranean-West",
                   "Canary Islands-West",
                   "Mediterranean-East",
                   "Mediterranean-West",
                   "Mediterranean-West",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Bahamas-Abaco",
                   "Mediterranean-West",
                   "Canary Islands-West",
                   "Mediterranean-West",
                   "Remove",
                   "Remove",
                   "Mediterranean-West",
                   "Mediterranean-West",
                   "Mediterranean-West",
                   "Remove",
                   "Remove",
                   "Mediterranean-West",
                   "Mediterranean-West",
                   "Bahamas-TOTO",
                   "Remove",
                   "Mediterranean-East",
                   "Mediterranean-West",
                   "Bahamas-TOTO",
                   "Canary Islands-West",
                   "Remove",
                   "Canary Islands-East",
                   "Mediterranean-West",
                   "Remove",
                   "Remove",
                   "Remove",
                   "Remove",
                   "Remove")
pop(zcav.fs.gl)
zcav.fs.gl<-gl.drop.pop(zcav.fs.gl, pop.list="Remove", recalc=TRUE)
zcav.fs.gl<-gl.compliance.check(zcav.fs.gl)
zcav.fs.gl<-gl.filter.monomorphs(zcav.fs.gl)
zcav.fs.gl #68 genotypes,  140,088 binary SNPs, size: 28.5 Mb 1070070 (11.23 %) missing data

mden<-cbind(indNames(mden.gl), as.character(pop(mden.gl)))
mden
mden.fs.nohaw.gl<-mden.gl
mden.fs.nohaw.gl<-gl.drop.pop(mden.fs.nohaw.gl, pop.list="Hawaii", mono.rm = TRUE)
mden.fs.gl<-gl.compliance.check(mden.fs.gl)
mden.fs.gl<-gl.filter.monomorphs(mden.fs.gl)

pop(mden.fs.gl)<-c("Bahamas-TOTO",
                   "Canary Islands-West",
                   "Bahamas-Abaco",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Remove",
                   "Canary Islands-East",
                   "Canary Islands-West",
                   "Bahamas-TOTO",
                   "Remove",
                   "Bahamas-Abaco",
                   "Remove",
                   "Bahamas-TOTO",
                   "Remove",
                   "Bahamas-Abaco",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Remove",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Remove",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Canary Islands-West",
                   "Remove",
                   "Remove",
                   "Canary Islands-West",
                   "Canary Islands-East",
                   "Canary Islands-West",
                   "Remove",
                   "Canary Islands-West",
                   "Canary Islands-East",
                   "Bahamas-TOTO",
                   "Remove",
                   "Remove",
                   "Remove",
                   "Remove",
                   "Remove",
                   "Remove",
                   "Remove")
pop(mden.fs.gl)
mden.fs.gl<-gl.drop.pop(mden.fs.gl, pop.list="Remove", recalc=TRUE)
mden.fs.gl<-gl.compliance.check(mden.fs.gl)
mden.fs.gl<-gl.filter.monomorphs(mden.fs.gl)
mden.fs.gl #27 genotypes,  47,127 binary SNPs, size: 7.6 Mb 71765 (5.64 %) missing data

#DAPC with focal site only snps and ind
zcav.fs.gl
mden.fs.gl

mden.sep.fs.gl<-seppop(mden.fs.gl)
popNames(mden.sep.fs.gl$`Bahamas-Abaco`)<-"Bahamas"
popNames(mden.sep.fs.gl$`Bahamas-TOTO`)<-"Bahamas"
mden.bah.fs.gl<-rbind(mden.sep.fs.gl$`Bahamas-Abaco`, mden.sep.fs.gl$`Bahamas-TOTO`)
pop(mden.bah.fs.gl)

popNames(mden.sep.fs.gl$`Canary Islands-East`)<-"Canaries"
popNames(mden.sep.fs.gl$`Canary Islands-West`)<-"Canaries"
mden.can.fs.gl<-rbind(mden.sep.fs.gl$`Canary Islands-East`, mden.sep.fs.gl$`Canary Islands-West`)
pop(mden.can.fs.gl)

mden.fr.gl<-rbind(mden.bah.fs.gl, mden.can.fs.gl)

zcav.fr.gl<-gl.filter.monomorphs(zcav.fr.gl)
mden.fr.gl<-gl.filter.monomorphs(mden.fr.gl)

zcav.fr.gl
mden.fr.gl

a<-mden.fr.gl
mat<-tab(a)
grp<-pop(a)
mden.xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
mden.xval
#$`Number of PCs Achieving Highest Mean Success`
#[1] "8"
#$`Number of PCs Achieving Lowest MSE`
#[1] "12"
#$call: dapc.data.frame(x = as.data.frame(x), grp = ..1, n.pca = ..2, n.da = ..3)
#$n.pca: 12 first PCs of PCA used
#$n.da: 1 discriminant functions saved
#$var (proportion of conserved variance): 0.568
par(mfrow=c(1,2))
par(mfrow=c(1,1))

scatter(mden.xval$DAPC, legend=TRUE, posi.leg="topleft")
#assignplot(mden.xval$DAPC)
#pheatmap(mden.xval[["DAPC"]][["posterior"]], treeheight_row = 0, treeheight_col = 0, labels_row=mden.xval[["DAPC"]][["grp"]])

x<-mden.xval[["DAPC"]]
n.grp <- ncol(x$posterior)
n.ind <- nrow(x$posterior)
Z <- t(x$posterior)
Z <- Z[,ncol(Z):1,drop=FALSE ]
image(x=1:n.grp, y=seq(.5, by=1, le=n.ind), Z, col=rev(heat.colors(100)), yaxt="n", ylab="", xaxt="n", xlab="Clusters")
axis(side=1, at=1:n.grp,tick=FALSE, labels=colnames(x$posterior))
axis(side=2, at=seq(.5, by=1, le=n.ind), labels=rev(rownames(x$posterior)), las=1, cex.axis=0.75)
abline(h=1:n.ind,  col="lightgrey")
abline(v=seq(0.5, by=1, le=2))
box()
newGrp <- colnames(x$posterior)
x.real.coord <- rev(match(x$grp, newGrp))
y.real.coord <- seq(.5, by=1, le=n.ind)
points(x.real.coord, y.real.coord, col="deepskyblue2", pch=3)

popNames(zcav.fs.gl)
zcav.sep.fs.gl<-seppop(zcav.fs.gl)
popNames(zcav.sep.fs.gl$`Bahamas-Abaco`)<-"Bahamas"
popNames(zcav.sep.fs.gl$`Bahamas-TOTO`)<-"Bahamas"
popNames(zcav.sep.fs.gl$`Canary Islands-East`)<-"Canaries"
popNames(zcav.sep.fs.gl$`Canary Islands-West`)<-"Canaries"
zcav.fr.gl<-rbind(zcav.sep.fs.gl$`Bahamas-Abaco`, zcav.sep.fs.gl$`Bahamas-TOTO`, zcav.sep.fs.gl$`Canary Islands-East`, zcav.sep.fs.gl$`Canary Islands-West`, zcav.sep.fs.gl$`Mediterranean-East`, zcav.sep.fs.gl$`Mediterranean-West`)
popNames(zcav.fr.gl)

b<-zcav.fr.gl
mat<-tab(b)
grp<-pop(b)
zcav.xval<-xvalDapc(mat, grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
zcav.xval
#$`Number of PCs Achieving Highest Mean Success`
#[1] "10"
#$`Number of PCs Achieving Lowest MSE`
#[1] "10"
#$call: dapc.data.frame(x = as.data.frame(x), grp = ..1, n.pca = ..2, n.da = ..3)
#$n.pca: 10 first PCs of PCA used
#$n.da: 3 discriminant functions saved
#$var (proportion of conserved variance): 0.376
scatter(zcav.xval$DAPC)
assignplot(zcav.xval$DAPC)
par(mar = c(10.1, 5.1, 1.1, 1.1), mfrow=c(1,2))
par(mar=c(5,5,5,5))

scatter(zcav.xval$DAPC, posi.da="Bottomleft", grid=TRUE)
scatter(zcav.xval$DAPC, leg=TRUE, clab=0, posi.da = "bottomleft", pch=20, cstar=0, cex=3, cell=0)

zcav.xval$DAPC$var.contr


#assignplot(mden.xval$DAPC)
#pheatmap(mden.xval[["DAPC"]][["posterior"]], treeheight_row = 0, treeheight_col = 0, labels_row=mden.xval[["DAPC"]][["grp"]])

x<-zcav.xval[["DAPC"]]
n.grp <- ncol(x$posterior)
n.ind <- nrow(x$posterior)
Z <- t(x$posterior)
Z <- Z[,ncol(Z):1,drop=FALSE ]
par(mar = c(10.1, 5.1, 1.1, 1.1))
image(x=1:n.grp, y=seq(.5, by=1, le=n.ind), Z, col=rev(heat.colors(100)), yaxt="n", ylab="", xaxt="n", xlab="")
axis(side=1, labels=FALSE)
text(x=1:length(grp), y=-10, labels=colnames(x$posterior), xpd=NA, srt=50, cex=1)
axis(side=2, at=seq(.5, by=1, le=n.ind), labels=rev(rownames(x$posterior)), las=1, cex.axis=0.5)
abline(h=1:n.ind,  col="lightgrey")
abline(v=seq(0.5, by=1, le=4))
box()
newGrp <- colnames(x$posterior)
x.real.coord <- rev(match(x$grp, newGrp))
y.real.coord <- seq(.5, by=1, le=n.ind)
points(x.real.coord, y.real.coord, col="deepskyblue2", pch=3)

zcav.fs.gl
mden.fs.gl
zcav.fs.seppop.gl<-seppop(zcav.fs.gl)
mden.fs.seppop.gl<-seppop(mden.fs.gl)

zcav.fs.ecan.gl<-zcav.fs.seppop.gl$`Canary Islands-East`
zcav.fs.wcan.gl<-zcav.fs.seppop.gl$`Canary Islands-West`
zcav.fs.abah.fs<-zcav.fs.seppop.gl$`Bahamas-Abaco`
zcav.fs.tbah.fs<-zcav.fs.seppop.gl$`Bahamas-TOTO`
zcav.fs.emed.gl<-zcav.fs.seppop.gl$`Mediterranean-East`
zcav.fs.wmed.gl<-zcav.fs.seppop.gl$`Mediterranean-West`
mden.fs.ecan.gl<-mden.fs.seppop.gl$`Canary Islands-East`
mden.fs.wcan.gl<-mden.fs.seppop.gl$`Canary Islands-West`
mden.fs.abah.fs<-mden.fs.seppop.gl$`Bahamas-Abaco`
mden.fs.tbah.fs<-mden.fs.seppop.gl$`Bahamas-TOTO`

zcav.fs.can.gl<-rbind(zcav.fs.ecan.gl, zcav.fs.wcan.gl)
zcav.fs.bah.gl<-rbind(zcav.fs.abah.fs, zcav.fs.tbah.fs)
zcav.fs.med.gl<-rbind(zcav.fs.emed.gl, zcav.fs.wmed.gl)
mden.fs.can.gl<-rbind(mden.fs.ecan.gl, mden.fs.wcan.gl)
mden.fs.bah.gl<-rbind(mden.fs.abah.fs, mden.fs.tbah.fs)

zcav.100.fs.gl<-gl.filter.callrate(zcav.fs.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
mden.100.fs.gl<-gl.filter.callrate(mden.fs.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.100.fs.ecan.gl<-gl.filter.callrate(zcav.fs.ecan.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.100.fs.wcan.gl<-gl.filter.callrate(zcav.fs.wcan.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.100.fs.abah.gl<-gl.filter.callrate(zcav.fs.abah.fs, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.100.fs.tbah.gl<-gl.filter.callrate(zcav.fs.tbah.fs, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.100.fs.emed.gl<-gl.filter.callrate(zcav.fs.emed.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.100.fs.wmed.gl<-gl.filter.callrate(zcav.fs.wmed.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
mden.100.fs.ecan.gl<-gl.filter.callrate(mden.fs.ecan.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
mden.100.fs.wcan.gl<-gl.filter.callrate(mden.fs.wcan.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
mden.100.fs.abah.gl<-gl.filter.callrate(mden.fs.abah.fs, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
mden.100.fs.tbah.gl<-gl.filter.callrate(mden.fs.tbah.fs, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)

zcav.fs.can.gl<-gl.compliance.check(zcav.fs.can.gl)
zcav.fs.bah.gl<-gl.compliance.check(zcav.fs.bah.gl)
zcav.fs.med.gl<-gl.compliance.check(zcav.fs.med.gl)
mden.fs.can.gl<-gl.compliance.check(mden.fs.can.gl)
mden.fs.bah.gl<-gl.compliance.check(mden.fs.bah.gl)

zcav.fs.can.gl<-gl.filter.callrate(zcav.fs.can.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.fs.bah.gl<-gl.filter.callrate(zcav.fs.bah.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
zcav.fs.med.gl<-gl.filter.callrate(zcav.fs.med.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
mden.fs.can.gl<-gl.filter.callrate(mden.fs.can.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)
mden.fs.bah.gl<-gl.filter.callrate(mden.fs.bah.gl, method="loc", threshold=1, mono.rm=TRUE, recalc=TRUE, plot=FALSE)

zcav.100.fs.gi<-gl2gi(zcav.100.fs.gl)
mden.100.fs.gi<-gl2gi(mden.100.fs.gl)
zcav.100.fs.ecan.gi<-gl2gi(zcav.100.fs.ecan.gl)
zcav.100.fs.wcan.gi<-gl2gi(zcav.100.fs.wcan.gl)
zcav.100.fs.abah.gi<-gl2gi(zcav.100.fs.abah.gl)
zcav.100.fs.tbah.gi<-gl2gi(zcav.100.fs.tbah.gl)
zcav.100.fs.emed.gi<-gl2gi(zcav.100.fs.emed.gl)
zcav.100.fs.wmed.gi<-gl2gi(zcav.100.fs.wmed.gl)
mden.100.fs.ecan.gi<-gl2gi(mden.100.fs.ecan.gl)
mden.100.fs.wcan.gi<-gl2gi(mden.100.fs.wcan.gl)
mden.100.fs.abah.gi<-gl2gi(mden.100.fs.abah.gl)
mden.100.fs.tbah.gi<-gl2gi(mden.100.fs.tbah.gl)

zcav.fs.can.gi<-gl2gi(zcav.fs.can.gl)
zcav.fs.bah.gi<-gl2gi(zcav.fs.bah.gl)
zcav.fs.med.gi<-gl2gi(zcav.fs.med.gl)
mden.fs.can.gi<-gl2gi(mden.fs.can.gl)
mden.fs.bah.gi<-gl2gi(mden.fs.bah.gl)

zcav.100.fs.gi
mden.100.fs.gi
zcav.100.fs.ecan.gi
zcav.100.fs.wcan.gi
zcav.100.fs.abah.gi
zcav.100.fs.tbah.gi
zcav.100.fs.emed.gi
zcav.100.fs.wmed.gi
mden.100.fs.ecan.gi
mden.100.fs.wcan.gi
mden.100.fs.abah.gi
mden.100.fs.tbah.gi

zcav.fs.can.gi
zcav.fs.bah.gi
zcav.fs.med.gi
mden.fs.can.gi
mden.fs.bah.gi

zcav.100.fs.gi<-informloci(zcav.100.fs.gi, cutoff=1/nInd(zcav.100.fs.gi))
mden.100.fs.gi<-informloci(mden.100.fs.gi, cutoff=1/nInd(mden.100.fs.gi))
zcav.100.fs.ecan.gi<-informloci(zcav.100.fs.ecan.gi, cutoff=1/nInd(zcav.100.fs.ecan.gi))
zcav.100.fs.wcan.gi<-informloci(zcav.100.fs.wcan.gi, cutoff=1/nInd(zcav.100.fs.wcan.gi))
zcav.100.fs.abah.gi<-informloci(zcav.100.fs.abah.gi, cutoff=1/nInd(zcav.100.fs.abah.gi))
zcav.100.fs.tbah.gi<-informloci(zcav.100.fs.tbah.gi, cutoff=1/nInd(zcav.100.fs.tbah.gi))
zcav.100.fs.emed.gi<-informloci(zcav.100.fs.emed.gi, cutoff=1/nInd(zcav.100.fs.emed.gi))
zcav.100.fs.wmed.gi<-informloci(zcav.100.fs.wmed.gi, cutoff=1/nInd(zcav.100.fs.wmed.gi))
mden.100.fs.ecan.gi<-informloci(mden.100.fs.ecan.gi, cutoff=1/nInd(mden.100.fs.ecan.gi))
mden.100.fs.wcan.gi<-informloci(mden.100.fs.wcan.gi, cutoff=1/nInd(mden.100.fs.wcan.gi))
mden.100.fs.abah.gi<-informloci(mden.100.fs.abah.gi, cutoff=1/nInd(mden.100.fs.abah.gi))
mden.100.fs.tbah.gi<-informloci(mden.100.fs.tbah.gi, cutoff=1/nInd(mden.100.fs.tbah.gi))

zcav.fs.can.gi<-informloci(zcav.fs.can.gi, cutoff=1/nInd(zcav.fs.can.gi))
zcav.fs.bah.gi<-informloci(zcav.fs.bah.gi, cutoff=1/nInd(zcav.fs.bah.gi))
zcav.fs.med.gi<-informloci(zcav.fs.med.gi, cutoff=1/nInd(zcav.fs.med.gi))
mden.fs.can.gi<-informloci(mden.fs.can.gi, cutoff=1/nInd(mden.fs.can.gi))
mden.fs.bah.gi<-informloci(mden.fs.bah.gi, cutoff=1/nInd(mden.fs.bah.gi))

summary(isPoly(zcav.100.fs.gi))
summary(isPoly(mden.100.fs.gi))
summary(isPoly(zcav.100.fs.ecan.gi))
summary(isPoly(zcav.100.fs.wcan.gi))
summary(isPoly(zcav.100.fs.abah.gi))
summary(isPoly(zcav.100.fs.tbah.gi))#go back
summary(isPoly(zcav.100.fs.emed.gi))
summary(isPoly(zcav.100.fs.wmed.gi))
summary(isPoly(mden.100.fs.ecan.gi))
summary(isPoly(mden.100.fs.wcan.gi))
summary(isPoly(mden.100.fs.abah.gi))
summary(isPoly(mden.100.fs.tbah.gi))

summary(isPoly(zcav.fs.can.gi))
summary(isPoly(zcav.fs.bah.gi))
summary(isPoly(zcav.fs.med.gi))
summary(isPoly(mden.fs.can.gi))
summary(isPoly(mden.fs.bah.gi))

zcav.fs.can.gl<-gi2gl(zcav.fs.can.gi)
zcav.fs.bah.gl<-gi2gl(zcav.fs.bah.gi)
zcav.fs.med.gl<-gi2gl(zcav.fs.med.gi)
mden.fs.can.gl<-gi2gl(mden.fs.can.gi)
mden.fs.bah.gl<-gi2gl(mden.fs.bah.gi)

zcav.fs.can.gl #23 genotypes,  56,246 binary SNPs, size: 7.9 Mb 0 (0 %) missing data
zcav.fs.bah.gl #8 genotypes,  61,643 binary SNPs, size: 8.4 Mb 0 (0 %) missing data
zcav.fs.med.gl #37 genotypes,  27,087 binary SNPs, size: 3.9 Mb 0 (0 %) missing data
mden.fs.can.gl #20 genotypes,  28,911 binary SNPs, size: 4 Mb 0 (0 %) missing data
mden.fs.bah.gl #7 genotypes,  26,890 binary SNPs, size: 3.7 Mb 0 (0 %) missing data

zcav.fs.can.dr<-gl2demerelate(zcav.fs.can.gl, verbose=TRUE)
zcav.fs.bah.dr<-gl2demerelate(zcav.fs.bah.gl, verbose=TRUE)
zcav.fs.med.dr<-gl2demerelate(zcav.fs.med.gl, verbose=TRUE)
mden.fs.can.dr<-gl2demerelate(mden.fs.can.gl, verbose=TRUE)
mden.fs.bah.dr<-gl2demerelate(mden.fs.bah.gl, verbose=TRUE)

#Remove dashes and underscores from sample names and replace with .
zcav.fs.can.dr$'Sample-ID'<-zcav.fs.can.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))
zcav.fs.bah.dr$'Sample-ID'<-zcav.fs.bah.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))
zcav.fs.med.dr$'Sample-ID'<-zcav.fs.med.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))
mden.fs.can.dr$'Sample-ID'<-mden.fs.can.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))
mden.fs.bah.dr$'Sample-ID'<-mden.fs.bah.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))

zcav.fs.can.D<-Demerelate(zcav.fs.can.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
zcav.fs.bah.D<-Demerelate(zcav.fs.bah.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
zcav.fs.med.D<-Demerelate(zcav.fs.med.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
mden.fs.can.D<-Demerelate(mden.fs.can.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
mden.fs.bah.D<-Demerelate(mden.fs.bah.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)

zcav.fs.emed.dr<-subset(zcav.fs.med.dr, zcav.fs.med.dr$Population=="Mediterranean-East")
zcav.fs.wmed.dr<-subset(zcav.fs.med.dr, zcav.fs.med.dr$Population=="Mediterranean-West")

zcav.fs.can.dr$Population<-"Canary Islands"
zcav.fs.bah.dr$Population<-"Bahamas"
zcav.fs.med.dr$Population<-"Mediterranean"
mden.fs.can.dr$Population<-"Canaries"
mden.fs.bah.dr$Population<-"Bahamas"
zcav.fs.emed.dr$Population
zcav.fs.wmed.dr

zcav.fs.can2.D<-Demerelate(zcav.fs.can.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
zcav.fs.bah2.D<-Demerelate(zcav.fs.bah.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
zcav.fs.med2.D<-Demerelate(zcav.fs.med.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
mden.fs.can2.D<-Demerelate(mden.fs.can.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
mden.fs.bah2.D<-Demerelate(mden.fs.bah.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
zcav.fs.emed2.D<-Demerelate(zcav.fs.emed.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)
zcav.fs.wmed2.D<-Demerelate(zcav.fs.wmed.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)

# Networks just focal sites -----------------------------------------------
#import mxy data into R
zcav.mxy<-read.csv("results/with_pop/Network_Files/22-4-21/zcav_fs_mxy.csv", header=TRUE)
head(zcav.mxy)
mden.mxy<-read.csv("./mden_mxy.csv", header=TRUE)
head(mden.mxy)
#import metadata into R
zcav.meta<-read.csv("./zcav_meta.csv", header=TRUE)
mden.meta<-read.csv("./mden_meta.csv", header=TRUE)
head(zcav.meta)
head(mden.meta)

#generate nodes table: list of individuals and any other attributes as columns (also vertices)
#make sure the sample IDs are in same format
zc.nodes<-zcav.meta[,c(2,7,3,5,4)] #filename, sex, focal_region, focal_site
head(zc.nodes)
colnames(zc.nodes)<-c("Filename", "sex", "focal_region", "focal_site", "fine-scale_site")
head(zc.nodes)
zc.nodes$Filename<- zc.nodes$Filename %>% str_replace_all(c("-"=".", "_"="."))
head(zc.nodes) 
tail(zc.nodes) 
is.data.frame(zc.nodes)

head(mden.meta)
md.nodes<-mden.meta[,c(2,8,3,5,4)]#filename, sex, focal_region, focal_site
head(md.nodes)
colnames(md.nodes)<-c("Filename", "sex", "focal_region", "focal_site", "fine-scale_site")
head(md.nodes)
md.nodes$Filename<- md.nodes$Filename %>% str_replace_all(c("-"=".", "_"="."))
head(md.nodes)
tail(md.nodes)
is.data.frame(md.nodes)

#Generate links (mxy values) between nodes
head(zcav.mxy)
zc.links<-zcav.mxy[,1:4]
head(zc.links)
tail(zc.links)
is.data.frame(zc.links)

head(mden.mxy)
md.links<-mden.mxy[,2:4]
head(md.links)
tail(md.links)
is.data.frame(md.links)

#Generating network for Zcav Canaries
zc.can.nodes<-subset(zc.nodes, zc.nodes$focal_region=="Canaries")
zc.can.nodes<-zc.can.nodes[order(zc.can.nodes$Filename),]
zc.can.nodes
zc.can.mxy<-subset(zc.links, zc.links$Region=="Canary Islands")
head(zc.can.mxy)
zc.can.links<-zc.can.mxy[,2:4]
head(zc.can.links)
#Build network
zc.can.net<-graph.data.frame(zc.can.links,directed=FALSE)
get.adjacency(zc.can.net, attr="Mxy", spars=FALSE)
V(zc.can.net)$name<-V(zc.can.net)$name[order(V(zc.can.net)$name)]
V(zc.can.net)$name

NEED TO RE DO THE FOCAL SITE ORDER TO MATCH THE SAMPLE ID ORDER


#Add Focal_Site attribute
length(V(zc.can.net)$name)
length(zc.can.nodes$Filename)
cbind(V(zc.can.net)$name,zc.can.nodes$Filename) #just to check that the sample names are in the same order
V(zc.can.net)$Focal_Site<-as.character(c("Canary Islands-East",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-East",
                                         "Canary Islands-East",
                                         "Canary Islands-East",
                                         "Canary Islands-East",
                                         "Canary Islands-East",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West"))
V(zc.can.net)$Focal_Site
# V(zc.can.net)$Focal_Site_Fine<-as.character(zc.can.nodes$`fine-scale_site`[match(V(zc.can.net)$name,zc.can.nodes$Filename)])
# V(zc.can.net)$Focal_Site_Fine
# cbind(V(zc.can.net)$Focal_Site,zc.can.nodes$focal_site) #just to check that the sample names are in the same order
#Delete links that don't indicate sibship and color by relatedness
delete_edges(zc.can.net, E(zc.can.net)[Mxy<0.852])
E(zc.can.net)[Mxy>=0.852]$color <- "green"
E(zc.can.net)[Mxy>=0.899]$color <- "blue"
get.adjacency(zc.can.net, attr="Mxy", spars=FALSE)

#plot network
library(RColorBrewer)
col  <- brewer.pal(5, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.can.net)$Focal_Site))]
plot(zc.can.net, vertex.color=my_color, vertex.label.degree=1, layout=layout_in_circle(zc.can.net, order=order(as.numeric(as.factor(V(zc.can.net)$Focal_Site)))))
legend("bottomleft", legend=levels(as.factor(V(zc.can.net)$Focal_Site)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

#Generating network for Zcav Bahamas
zc.bah.nodes<-subset(zc.nodes, zc.nodes$focal_region=="Caribbean")
head(zc.bah.nodes)
zc.bah.nodes<-zc.bah.nodes[order(zc.bah.nodes$Filename),]
zc.bah.nodes
zc.bah.mxy<-subset(zcav.mxy, zcav.mxy$Region=="Bahamas")
head(zc.bah.mxy)
zc.bah.links<-zc.bah.mxy[,2:4]
head(zc.bah.links)
#Build network
zc.bah.net<-graph.data.frame(zc.bah.links,directed=FALSE)
get.adjacency(zc.bah.net, attr="Mxy", spars=FALSE)
V(zc.bah.net)$name<-V(zc.bah.net)$name[order(V(zc.bah.net)$name)]
V(zc.bah.net)$name
#Add Focal Region as attribute
length(V(zc.bah.net)$name)
length(zc.bah.nodes$Filename)
cbind(V(zc.bah.net)$name,zc.bah.nodes$Filename) #just to check that the sample names are in the same order
V(zc.bah.net)$Focal_Site<-as.character(zc.bah.nodes$focal_site[match(V(zc.bah.net)$name,zc.bah.nodes$Filename)])
V(zc.bah.net)$Focal_Site
V(zc.bah.net)$Focal_Site_Fine<-as.character(zc.bah.nodes$`fine-scale_site`[match(V(zc.bah.net)$name,zc.bah.nodes$Filename)])
V(zc.bah.net)$Focal_Site_Fine
cbind(V(zc.bah.net)$Focal_Site,zc.bah.nodes$focal_site) #just to check that the sample names are in the same order
#Delete links that don't indicate sibship and color by relatedness
delete_edges(zc.bah.net, E(zc.bah.net)[Mxy<0.787])
E(zc.bah.net)[Mxy>0.787]$color <- "green"
E(zc.bah.net)[Mxy>0.856]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.bah.net)$Focal_Site))]
plot(zc.bah.net, vertex.color=my_color, layout=layout_in_circle(zc.bah.net, order=order(as.numeric(as.factor(V(zc.bah.net)$Focal_Site_Fine)))))
legend("bottomleft", legend=levels(as.factor(V(zc.bah.net)$Focal_Site_Fine)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

#Generating network for Zcav E-Mediterranean
zc.emed.nodes<-subset(zc.nodes, zc.nodes$focal_region==c("Mediterranean-E"))
head(zc.emed.nodes)
zc.emed.nodes<-zc.emed.nodes[order(zc.emed.nodes$Filename),]
zc.emed.nodes
zc.emed.mxy<-subset(zcav.mxy, zcav.mxy$Region=="Med-East")
head(zc.emed.mxy)
zc.emed.links<-zc.emed.mxy[,2:4]
head(zc.emed.links)
#Build network
zc.emed.net<-graph.data.frame(zc.emed.links,directed=FALSE)
get.adjacency(zc.emed.net, attr="Mxy", spars=FALSE)
V(zc.emed.net)$name<-V(zc.emed.net)$name[order(V(zc.emed.net)$name)]
V(zc.emed.net)$name
#Add Focal Region as attribute
length(V(zc.emed.net)$name)
length(zc.emed.nodes$Filename)
cbind(V(zc.emed.net)$name,zc.emed.nodes$Filename) #just to check that the sample names are in the same order
V(zc.emed.net)$Focal_Site<-as.character(zc.emed.nodes$focal_site[match(V(zc.emed.net)$name,zc.emed.nodes$Filename)])
V(zc.emed.net)$Focal_Site
V(zc.emed.net)$Focal_Site_Fine<-as.character(zc.emed.nodes$`fine-scale_site`[match(V(zc.emed.net)$name,zc.emed.nodes$Filename)])
V(zc.emed.net)$Focal_Site_Fine
cbind(V(zc.emed.net)$Focal_Site,zc.emed.nodes$Focal_Site_Fine) #just to check that the sample names are in the same order
#Delete links that don't indicate sibship and color by relatedness
delete_edges(zc.emed.net, E(zc.emed.net)[Mxy<0.762])
E(zc.emed.net)[Mxy>0.762]$color <- "green"
E(zc.emed.net)[Mxy>0.841]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(5, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.emed.net)$Focal_Site_Fine))]
plot(zc.emed.net, vertex.color=my_color, layout=layout_in_circle(zc.emed.net, order=order(as.numeric(as.factor(V(zc.emed.net)$Focal_Site_Fine)))))
legend("bottomleft", legend=levels(as.factor(V(zc.emed.net)$Focal_Site_Fine)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.bah.net)$Focal_Site))]
plot(zc.bah.net, vertex.color=my_color, layout=layout_in_circle(zc.bah.net, order=order(as.numeric(as.factor(V(zc.bah.net)$Focal_Site_Fine)))))
legend("bottomleft", legend=levels(as.factor(V(zc.bah.net)$Focal_Site_Fine)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

col  <- brewer.pal(5, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.can.net)$Focal_Site))]
plot(zc.can.net, vertex.color=my_color, vertex.label.degree=1, layout=layout_in_circle(zc.can.net, order=order(as.numeric(as.factor(V(zc.can.net)$Focal_Site)))))
legend("bottomleft", legend=levels(as.factor(V(zc.can.net)$Focal_Site)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

#Generating network for Zcav W-Mediterranean
zc.wmed.nodes<-subset(zc.nodes, zc.nodes$focal_region==c("Mediterranean-W"))
head(zc.wmed.nodes)
zc.wmed.nodes<-zc.wmed.nodes[order(zc.wmed.nodes$Filename),]
zc.wmed.nodes
zc.wmed.mxy<-subset(zcav.mxy, zcav.mxy$Region=="Med-West")
head(zc.wmed.mxy)
zc.wmed.links<-zc.wmed.mxy[,2:4]
head(zc.wmed.links)
#Build network
zc.wmed.net<-graph.data.frame(zc.wmed.links,directed=FALSE)
get.adjacency(zc.wmed.net, attr="Mxy", spars=FALSE)
V(zc.wmed.net)$name<-V(zc.wmed.net)$name[order(V(zc.wmed.net)$name)]
V(zc.wmed.net)$name
#Add Focal Region as attribute
length(V(zc.wmed.net)$name)
length(zc.wmed.nodes$Filename)
cbind(V(zc.wmed.net)$name,zc.wmed.nodes$Filename) #just to check that the sample names are in the same order
V(zc.wmed.net)$Focal_Site<-as.character(c("Corsica",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "CIMA (Ligurian)",
                                          "Corsica"))
V(zc.wmed.net)$Focal_Site
V(zc.wmed.net)$Focal_Site_Fine<-as.character(zc.wmed.nodes$`fine-scale_site`[match(V(zc.wmed.net)$name,zc.wmed.nodes$Filename)])
V(zc.wmed.net)$Focal_Site_Fine
cbind(V(zc.wmed.net)$Focal_Site,zc.wmed.nodes$focal_site) #just to check that the sample names are in the same order
#Delete links that don't indicate sibship and color by relatedness
delete_edges(zc.wmed.net, E(zc.wmed.net)[Mxy<0.759])
E(zc.wmed.net)[Mxy>0.759]$color <- "green"
E(zc.wmed.net)[Mxy>0.839]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.wmed.net)$Focal_Site))]
plot(zc.wmed.net, vertex.color=my_color, layout=layout_in_circle(zc.wmed.net, order=order(as.numeric(as.factor(V(zc.wmed.net)$Focal_Site)))))
legend("bottomleft", legend=levels(as.factor(V(zc.wmed.net)$Focal_Site)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

col  <- brewer.pal(5, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.emed.net)$Focal_Site_Fine))]
plot(zc.emed.net, vertex.color=my_color, layout=layout_in_circle(zc.emed.net, order=order(as.numeric(as.factor(V(zc.emed.net)$Focal_Site_Fine)))))
legend("bottomleft", legend=levels(as.factor(V(zc.emed.net)$Focal_Site_Fine)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.bah.net)$Focal_Site))]
plot(zc.bah.net, vertex.color=my_color, layout=layout_in_circle(zc.bah.net, order=order(as.numeric(as.factor(V(zc.bah.net)$Focal_Site_Fine)))))
legend("bottomleft", legend=levels(as.factor(V(zc.bah.net)$Focal_Site_Fine)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

col  <- brewer.pal(5, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.can.net)$Focal_Site_Fine))]
plot(zc.can.net, vertex.color=my_color, vertex.label.degree=1, layout=layout_in_circle(zc.can.net, order=order(as.numeric(as.factor(V(zc.can.net)$Focal_Site_Fine)))))
legend("bottomleft", legend=levels(as.factor(V(zc.can.net)$Focal_Site_Fine)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

#Generating network for Mden Bahamas
md.bah.nodes<-subset(md.nodes, md.nodes$focal_region=="Caribbean")
head(md.bah.nodes)
md.bah.nodes<-md.bah.nodes[order(md.bah.nodes$Filename),]
md.bah.nodes
md.bah.mxy<-subset(mden.mxy, mden.mxy$Region=="Bahamas")
head(md.bah.mxy)
md.bah.links<-md.bah.mxy[,2:4]
head(md.bah.links)
md.bah.links

#Build network
md.bah.net<-graph.data.frame(md.bah.links,directed=FALSE)
get.adjacency(md.bah.net, attr="Mxy", spars=FALSE)
V(md.bah.net)$name<-V(md.bah.net)$name[order(V(md.bah.net)$name)]
V(md.bah.net)$name
#Add Focal Region as attribute
length(V(md.bah.net)$name)
length(md.bah.nodes$Filename)
cbind(V(md.bah.net)$name,md.bah.nodes$Filename) #just to check that the sample names are in the same order
V(md.bah.net)$Focal_Site<-as.character(c("Abaco","Abaco", "Abaco", "TOTO", "TOTO", "TOTO","TOTO"))
V(md.bah.net)$Focal_Site
# V(md.bah.net)$Focal_Site_Fine<-as.character(md.bah.nodes$`fine-scale_site`[match(V(md.bah.net)$name,md.bah.nodes$Filename)])
# V(md.bah.net)$Focal_Site_Fine
# cbind(V(md.bah.net)$Focal_Site,md.bah.nodes$focal_site) #just to check that the sample names are in the same order
#Delete links that don't indicate sibship and color by relatedness
delete_edges(md.bah.net, E(md.bah.net)[Mxy<0.744])
E(md.bah.net)[Mxy>0.744]$color <- "green"
E(md.bah.net)[Mxy>0.829]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(4, "Set1") 
my_color <- col[as.numeric(as.factor(V(md.bah.net)$Focal_Site))]
plot(md.bah.net, vertex.color=my_color, layout=layout_in_circle(md.bah.net, order=order(as.numeric(as.factor(V(md.bah.net)$Focal_Site)))))
legend("bottomleft", legend=levels(as.factor(V(md.bah.net)$Focal_Site)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

#Generating network for Mden Canaries
md.can.nodes<-subset(md.nodes, md.nodes$focal_region=="Canaries")
md.can.nodes<-md.can.nodes[order(md.can.nodes$Filename),]
md.can.nodes
md.can.mxy<-subset(mden.mxy, mden.mxy$Region=="Canaries")
head(md.can.mxy)
md.can.links<-md.can.mxy[,2:4]
head(md.can.links)
#Build network
md.can.net<-graph.data.frame(md.can.links,directed=FALSE)
get.adjacency(md.can.net, attr="Mxy", spars=FALSE)
# V(md.can.net)$name<-V(md.can.net)$name[order(V(md.can.net)$name)]
# V(md.can.net)$name
#Add Focal_Site attribute
length(V(md.can.net)$name)
length(md.can.nodes$Filename)
cbind(V(md.can.net)$name,md.can.nodes$Filename) #just to check that the sample names are in the same order
V(md.can.net)$Focal_Site<-as.character(c("Canary Islands-East",
                                         "Canary Islands-East",
                                         "Canary Islands-East",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West",
                                         "Canary Islands-West"))

V(md.can.net)$Focal_Site
# V(md.can.net)$Focal_Site_Fine<-as.character(md.can.nodes$`fine-scale_site`[match(V(md.can.net)$name,md.can.nodes$Filename)])
# V(md.can.net)$Focal_Site_Fine
# cbind(V(md.can.net)$Focal_Site,md.can.nodes$focal_site) #just to check that the sample names are in the same order
#Delete links that don't indicate sibship and color by relatedness
delete_edges(md.can.net, E(md.can.net)[Mxy<0.784])
E(md.can.net)[Mxy>=0.784]$color <- "green"
E(md.can.net)[Mxy>=0.855]$color <- "blue"

#plot network
library(RColorBrewer)
col  <- brewer.pal(5, "Set1") 
my_color <- col[as.numeric(as.factor(V(md.can.net)$Focal_Site))]
plot(md.can.net, vertex.color=my_color, layout=layout_in_circle(md.can.net,order=order(as.numeric(as.factor(V(md.can.net)$Focal_Site))) ))
plot(md.can.net, vertex.color=my_color, layout=layout_in_circle(md.can.net))

legend("bottomleft", legend=levels(as.factor(V(md.can.net)$Focal_Site)), col=col, bty="n", pch=20, pt.cex=3)
legend("topleft", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

##Plot Mden Networks together
par(mfrow=c(1,2))
par(mar=c(2,2,2,2))
col  <- brewer.pal(5, "Set1")  
my_color <- col[as.numeric(as.factor(V(md.bah.net)$Focal_Site))]
plot(md.bah.net, main="Bahamas", vertex.color=my_color, layout=layout_in_circle(md.bah.net, order=order(as.numeric(as.factor(V(md.bah.net)$Focal_Site)))))
legend("bottomleft", legend=levels(as.factor(V(md.bah.net)$Focal_Site)), col=col, bty="n", pch=20, pt.cex=3)
legend("bottomright", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

my_color <- col[as.numeric(as.factor(V(md.can.net)$Focal_Site))]
plot(md.can.net, main="Canary Islands", vertex.color=my_color, layout=layout_in_circle(md.can.net,order=order(as.numeric(as.factor(V(md.can.net)$Focal_Site))) ))
legend("bottomleft", legend=levels(as.factor(V(md.can.net)$Focal_Site)), col=col, bty="n", pch=20, pt.cex=3)
legend("bottomright", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

##Plot Zcav Networks together
par(mar=c(2,2,2,2))
par(mfrow=c(2,2))
col  <- brewer.pal(5, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.bah.net)$Focal_Site))]
plot(zc.bah.net, main="Bahamas", vertex.color=my_color, layout=layout_in_circle(zc.bah.net, order=order(as.numeric(as.factor(V(zc.bah.net)$Focal_Site)))))
legend("bottomleft", legend=levels(as.factor(V(zc.bah.net)$Focal_Site)), col=col, bty="n", pch=20, pt.cex=3)
legend("bottomright", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

col  <- brewer.pal(5, "Set1") 
my_color <- col[as.numeric(as.factor(V(zc.can.net)$Focal_Site))]
plot(zc.can.net, main="Canary Islands", vertex.color=my_color, vertex.label.degree=1, layout=layout_in_circle(zc.can.net, order=order(as.numeric(as.factor(V(zc.can.net)$Focal_Site)))))
legend("bottomleft", legend=levels(as.factor(V(zc.can.net)$Focal_Site)), col=col, bty="n", pch=20, pt.cex=3)
legend("bottomright", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

my_color <- col[as.numeric(as.factor(V(zc.wmed.net)$Focal_Site))]
plot(zc.wmed.net, main="Mediterranean-West", vertex.color=my_color, layout=layout_in_circle(zc.wmed.net, order=order(as.numeric(as.factor(V(zc.wmed.net)$Focal_Site)))))
legend("bottomleft", legend=levels(as.factor(V(zc.wmed.net)$Focal_Site)), col=col, bty="n", pch=20, pt.cex=3)
legend("bottomright", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

my_color <- col[as.numeric(as.factor(V(zc.emed.net)$Focal_Site_Fine))]
plot(zc.emed.net, main="Mediterranean-East", vertex.color=my_color, layout=layout_in_circle(zc.emed.net, order=order(as.numeric(as.factor(V(zc.emed.net)$Focal_Site_Fine)))))
legend("bottomleft", legend=levels(as.factor(V(zc.emed.net)$Focal_Site_Fine)), col=col, bty="n", pch=20, pt.cex=3)
legend("bottomright", legend=c("Half Sibling", "Full Sibling"), col=c("green", "blue"), lty=1, bty="n")

dev.off()

# Boxplots of relatedness, focal only -------------------------------------
zcav.mxy
mden.mxy

boxplot(zcav.mxy$Mxy~zcav.mxy$Region)
boxplot(mden.mxy$Mxy~mden.mxy$Region)

#Species	zcav	zcav	zcav	zcav
#Pop	Bahamas	Canaries	Emediterranean  Wmediterranean	
#HS_threshold	0.787 0.852 0.762 0.759	
#FS_threshold	0.856 0.899 0.841 0.839	
#NS_thrshold  0.744 0.824 0.713 0.710
zcav.hs<-data.frame(value=c(0.787, 0.852, 0.762, 0.759), boxplot.nr=c(1,2,3,4))
zcav.fs<-data.frame(value=c(.856, 0.899, 0.841, 0.839), boxplot.nr=c(1,2,3,4))
zcav.ns<-data.frame(value=c(0.744, 0.824, 0.713, 0.710), boxplot.nr=c(1,2,3,4))

mean(mden.fs.bah.D[[1]][["Randomized_Populations_for_Relatedness_Statistics"]][[3]][["Randomized_Non"]])
mean(mden.fs.can.D[[1]][["Randomized_Populations_for_Relatedness_Statistics"]][[3]][["Randomized_Non"]])
mean(zcav.fs.bah2.D[[1]][["Randomized_Populations_for_Relatedness_Statistics"]][[3]][["Randomized_Non"]])
mean(zcav.fs.can2.D[[1]][["Randomized_Populations_for_Relatedness_Statistics"]][[3]][["Randomized_Non"]])
mean(zcav.fs.bah2.D[[1]][["Randomized_Populations_for_Relatedness_Statistics"]][[3]][["Randomized_Non"]])
mean(zcav.emed.D[[1]][["Randomized_Populations_for_Relatedness_Statistics"]][[3]][["Randomized_Non"]])
mean(zcav.wmed.D[[1]][["Randomized_Populations_for_Relatedness_Statistics"]][[3]][["Randomized_Non"]])

#Species	mden  mden
#Pop	Bahamas	Canaries	
#HS_threshold	0.744 0.784 
#FS_threshold	0.829 0.855 
#NS_threshold 0.692 0.741
mden.hs<-data.frame(value=c(0.744, 0.784), boxplot.nr=c(1,2))
mden.fs<-data.frame(value=c(0.829, 0.855),boxplot.nr=c(1,2))
mden.ns<-data.frame(value=c(0.692, 0.741), boxplot.nr=c(1,2))

ggplot(zcav.mxy, aes(x=Region, y=Mxy)) +
  geom_boxplot() +
  theme_classic() +
  labs(x=NULL, y="Pairwise Relatedness (Mxy)") +
  scale_x_discrete(labels=c("Bahamas", "Canary Islands", "Mediterranean-East", "Mediterranean-West")) +
  theme(axis.text.x=element_text(angle = 45, hjust = .5, vjust=0.5)) +
  geom_segment(data=zcav.hs,aes(x=boxplot.nr-0.45,xend=boxplot.nr+0.45,
                                y=value,yend=value),inherit.aes=FALSE,size=1.5, linetype=2) +
  geom_segment(data=zcav.fs,aes(x=boxplot.nr-0.45,xend=boxplot.nr+0.45,
                                y=value,yend=value),inherit.aes=FALSE,size=1.5) +
  geom_segment(data=zcav.ns,aes(x=boxplot.nr-0.45,xend=boxplot.nr+0.45,
                                y=value,yend=value),inherit.aes=FALSE,size=1.5, col="grey", linetype=2)



ggplot(mden.mxy, aes(x=Region, y=Mxy)) +
  geom_boxplot() +
  theme_classic() +
  labs(x=NULL, y="Pairwise Relatedness (Mxy)") +
  scale_x_discrete(labels=c("Bahamas", "Canary Islands")) +
  theme(axis.text.x=element_text(angle = 45, hjust = .5, vjust=0.5)) +
  geom_segment(data=mden.hs,aes(x=boxplot.nr-0.45,xend=boxplot.nr+0.45,
                                y=value,yend=value),inherit.aes=FALSE,size=1.5, linetype=2) +
  geom_segment(data=mden.fs,aes(x=boxplot.nr-0.45,xend=boxplot.nr+0.45,
                                y=value,yend=value),inherit.aes=FALSE,size=1.5) +
  geom_segment(data=mden.ns,aes(x=boxplot.nr-0.45,xend=boxplot.nr+0.45,
                                y=value,yend=value),inherit.aes=FALSE,size=1.5, col="grey", linetype=2)
#merging databases
zcav.mxy
zc.nodes
colnames(zcav.mxy)<-c("Region", "Filename", "Ind2", "Mxy")
colnames(zcav.mxy)
colnames(zc.nodes)
zc.join1<-dplyr::full_join(zcav.mxy, zc.nodes, by="Filename")
head(zc.join1)
colnames(zc.join1)<-c("Region", "Ind1", "Filename", "Mxy", "sex", "focal_region", "focal_site")
colnames(zc.join1)
zc.join2<-dplyr::full_join(zc.join1, zc.nodes, by="Filename")
head(zc.join2)

mden.mxy
md.nodes
colnames(mden.mxy)<-c("Region", "Filename", "Ind2", "Mxy")
colnames(mden.mxy)
colnames(md.nodes)
md.join1<-dplyr::full_join(mden.mxy, md.nodes, by="Filename")
head(md.join1)
colnames(md.join1)<-c("Region", "Ind1", "Filename", "Mxy", "sex", "focal_region", "focal_site")
colnames(md.join1)
md.join2<-dplyr::full_join(md.join1, md.nodes, by="Filename")
head(md.join2)

zc.join2
md.join2

zc.join2$keep<-ifelse(zc.join2$focal_site.x==zc.join2$focal_site.y, "Keep", "Remove")
head(zc.join2)
md.join2$keep<-ifelse(md.join2$focal_site.x==md.join2$focal_site.y, "Keep", "Remove")
head(md.join2)

zc.fs.mxy<-subset(zc.join2, zc.join2$keep=="Keep")
head(zc.fs.mxy)
zc.fs.mxy<-zc.fs.mxy[,c(1,4,7,8)]
colnames(zc.fs.mxy)<-c("focal_region", "mxy", "focal_site", "site")
head(zc.fs.mxy)

md.fs.mxy<-subset(md.join2, md.join2$keep=="Keep")
head(md.fs.mxy)
md.fs.mxy<-md.fs.mxy[,c(1,4,7,8)]
colnames(md.fs.mxy)<-c("focal_region", "mxy", "focal_site", "site")
head(md.fs.mxy)

ggplot(zc.fs.mxy, aes(x=focal_site, y=mxy)) +
  geom_boxplot() +
  theme_classic() +
  labs(x=NULL, y="Pairwise Relatedness (Mxy)") +
  theme(axis.text.x=element_text(angle = 45, hjust = .5, vjust=0.5)) +
  geom_segment(data=zcav.hs,aes(x=c(0.5,2.55,4.55,5.6),xend=c(2.5,4.5,5.5,6.5),
                                y=value,yend=value),inherit.aes=FALSE,size=1.5, linetype=2) +
  geom_segment(data=zcav.fs,aes(x=c(0.5,2.55,4.55,5.6),xend=c(2.5,4.5,5.5,6.5),
                                y=value,yend=value),inherit.aes=FALSE,size=1.5) +
  geom_segment(data=zcav.ns,aes(x=c(0.5,2.55,4.55,5.6),xend=c(2.5,4.5,5.5,6.5),
                                y=value,yend=value),inherit.aes=FALSE,size=1.5, col="grey", linetype=2)

ggplot(md.fs.mxy, aes(x=focal_site, y=mxy)) +
  geom_boxplot() +
  theme_classic() +
  labs(x=NULL, y="Pairwise Relatedness (Mxy)") +
  theme(axis.text.x=element_text(angle = 45, hjust = .5, vjust=0.5)) +
  geom_segment(data=mden.hs,aes(x=c(0.5,2.55),xend=c(2.45,4.45),
                                y=value,yend=value),inherit.aes=FALSE,size=1.5, linetype=2) +
  geom_segment(data=mden.fs,aes(x=c(0.5,2.55),xend=c(2.5,4.5),
                                y=value,yend=value),inherit.aes=FALSE,size=1.5) +
  geom_segment(data=mden.ns,aes(x=c(0.5,2.55),xend=c(2.5,4.5),
                                y=value,yend=value),inherit.aes=FALSE,size=1.5, col="grey", linetype=2)

# Internal Relatedness ----------------------------------------------------


zcav.fr.gl$pop
mden.fr.gl

zcav.fs.gl<-rbind(zcav.fs.seppop.gl$`Bahamas-Abaco`, zcav.fs.seppop.gl$`Bahamas-TOTO`, zcav.fs.seppop.gl$`Canary Islands-East`, zcav.fs.seppop.gl$`Canary Islands-West`, zcav.fs.seppop.gl$`Mediterranean-East`, zcav.fs.seppop.gl$`Mediterranean-West`)
zcav.fs.gl #this is made from combining the seppop files (so it will be all SNPs, not just region-specificc ones)

mden.fs.gl<-rbind(mden.fs.seppop.gl$`Bahamas-Abaco`, mden.fs.seppop.gl$`Bahamas-TOTO`, mden.fs.seppop.gl$`Canary Islands-East`, mden.fs.seppop.gl$`Canary Islands-West`)
mden.fs.gl


zcav.fr.gl#68 genotypes,  140,088 binary SNPs, size: 16.1 Mb, 1070070 (11.23 %) missing data; Levels: Bahamas Canaries Mediterranean-East Mediterranean-West
zcav.fs.gl#68 genotypes,  140,088 binary SNPs, size: 16.1 Mb, 1070070 (11.23 %) missing data; Levels: Bahamas-Abaco Bahamas-TOTO Canary Islands-East Canary Islands-West Mediterranean-East Mediterranean-West
mden.fr.gl#27 genotypes,  47,127 binary SNPs, size: 3.9 Mb, 71765 (5.64 %) missing data; Levels: Bahamas Canaries
mden.fs.gl#27 genotypes,  47,127 binary SNPs, size: 3.9 Mb, 71765 (5.64 %) missing data; Levels: Bahamas-Abaco Bahamas-TOTO Canary Islands-East Canary Islands-West


library(haplo.stats)
zcav.fr.2col<-geno1to2(as.matrix(zcav.fr.gl), locus.label=NULL)
zcav.fr.2col[1:5,1:5]
zcav.fs.2col<-geno1to2(as.matrix(zcav.fs.gl), locus.label = NULL)
zcav.fs.2ccol[1:5,1:5]

mden.fr.2col<-geno1to2(as.matrix(mden.fr.gl), locus.label=NULL)
mden.fr.2col[1:5,1:5]
mden.fs.2col<-geno1to2(as.matrix(mden.fs.gl), locus.label=NULL)
mden.fs.2col[1:5,1:5]

library(Rhh)
mden.all.fr.ir<-ir(mden.fr.2col)
mden.all.fr.ir
zcav.all.fr.ir<-ir(zcav.fr.2col)
zcav.all.fr.ir

mden.all.fs.ir<-ir(mden.fs.2col)
mden.all.fs.ir
zcav.all.fs.ir<-ir(zcav.fs.2col)
zcav.all.fs.ir

#these are exactly the same
cbind(mden.all.fr.ir, mden.all.fs.ir)
cbind(zcav.all.fr.ir, zcav.all.fs.ir)

zcav.pop.fr.ir<-as.data.frame(cbind(as.character(zcav.fr.gl$ind.names), as.character(zcav.fr.gl$pop), as.numeric(zcav.all.fr.ir)))
zcav.pop.fr.ir
zcav.pop.fr.ir[,1]
rownames(zcav.pop.fr.ir)<-zcav.pop.fr.ir[,1]
zcav.pop.fr.ir<-zcav.pop.fr.ir[,2:3]
colnames(zcav.pop.fr.ir)<-c("population", "IR")
zcav.pop.fr.ir

zcav.pop.fs.ir<-as.data.frame(cbind(as.character(zcav.fs.gl$ind.names), as.character(zcav.fs.gl$pop), as.numeric(zcav.all.fs.ir)))
zcav.pop.fs.ir
zcav.pop.fs.ir[,1]
rownames(zcav.pop.fs.ir)<-zcav.pop.fs.ir[,1]
zcav.pop.fs.ir<-zcav.pop.fs.ir[,2:3]
colnames(zcav.pop.fs.ir)<-c("population", "IR")
zcav.pop.fs.ir

mden.pop.fr.ir<-as.data.frame(cbind(as.character(mden.fr.gl$ind.names), as.character(mden.fr.gl$pop), as.numeric(mden.all.fr.ir)))
mden.pop.fr.ir
mden.pop.fr.ir[,1]
rownames(mden.pop.fr.ir)<-mden.pop.fr.ir[,1]
mden.pop.fr.ir<-mden.pop.fr.ir[,2:3]
colnames(mden.pop.fr.ir)<-c("population", "IR")
mden.pop.fr.ir

mden.pop.fs.ir<-as.data.frame(cbind(as.character(mden.fs.gl$ind.names), as.character(mden.fs.gl$pop), as.numeric(mden.all.fs.ir)))
mden.pop.fs.ir
mden.pop.fs.ir[,1]
rownames(mden.pop.fs.ir)<-mden.pop.fs.ir[,1]
mden.pop.fs.ir<-mden.pop.fs.ir[,2:3]
colnames(mden.pop.fs.ir)<-c("population", "IR")
mden.pop.fs.ir

boxplot(as.numeric(as.character(zcav.pop.fr.ir$IR)) ~ zcav.pop.fr.ir$population)
boxplot(as.numeric(as.character(mden.pop.fr.ir$IR)) ~ mden.pop.fr.ir$population)
boxplot(as.numeric(as.character(zcav.pop.fs.ir$IR)) ~ zcav.pop.fs.ir$population)
boxplot(as.numeric(as.character(mden.pop.fs.ir$IR)) ~ mden.pop.fs.ir$population)

##Calculate IR within populations

zcav.fs.gl$pop

zcav.fs.bah.gl
zcav.fs.can.gl
zcav.100.emed.gl
zcav.100.wmed.gl
mden.fs.bah.gl
mden.fs.can.gl

zcav.fs.bah.2col<- geno1to2(as.matrix(zcav.fs.bah.gl), locus.label=NULL)
zcav.fs.can.2col<- geno1to2(as.matrix(zcav.fs.can.gl), locus.label=NULL)
zcav.fs.emed.2col<-geno1to2(as.matrix(zcav.100.emed.gl), locus.label=NULL)
zcav.fs.wmed.2col<-geno1to2(as.matrix(zcav.100.wmed.gl), locus.label=NULL)
mden.fs.bah.2col<- geno1to2(as.matrix(mden.fs.bah.gl), locus.label=NULL)
mden.fs.can.2col<- geno1to2(as.matrix(mden.fs.can.gl), locus.label=NULL)
zcav.fs.med.2col<-geno1to2(as.matrix(zcav.100.med.gi), locus.label=NULL)

zcav.fs.bah.ir<- ir(zcav.fs.bah.2col)
zcav.fs.can.ir<- ir(zcav.fs.can.2col)
zcav.fs.emed.ir<-ir(zcav.fs.emed.2col)
zcav.fs.wmed.ir<-ir(zcav.fs.wmed.2col)
mden.fs.bah.ir<- ir(mden.fs.bah.2col)
mden.fs.can.ir<- ir(mden.fs.can.2col)
zcav.fs.med.ir<- ir(zcav.fs.med.2col)

zcav.fs.bah.ir<-as.data.frame(cbind(as.character(zcav.fs.bah.gl$ind.names), as.character(zcav.fs.bah.gl$pop), as.numeric(zcav.fs.bah.ir)))
zcav.fs.can.ir<-as.data.frame(cbind(as.character(zcav.fs.can.gl$ind.names), as.character(zcav.fs.can.gl$pop), as.numeric(zcav.fs.can.ir)))
zcav.fs.emed.ir<-as.data.frame(cbind(as.character(zcav.fs.emed.gl$ind.names), as.character(zcav.fs.emed.gl$pop), as.numeric(zcav.fs.emed.ir)))
zcav.fs.wmed.ir<-as.data.frame(cbind(as.character(zcav.fs.wmed.gl$ind.names), as.character(zcav.fs.wmed.gl$pop), as.numeric(zcav.fs.wmed.ir)))
mden.fs.bah.ir<-as.data.frame(cbind(as.character(mden.fs.bah.gl$ind.names), as.character(mden.fs.bah.gl$pop), as.numeric(mden.fs.bah.ir)))
mden.fs.can.ir<-as.data.frame(cbind(as.character(mden.fs.can.gl$ind.names), as.character(mden.fs.can.gl$pop), as.numeric(mden.fs.can.ir)))
zcav.fs.med.ir<-as.data.frame(cbind(as.character(indNames(zcav.100.med.gi)), as.character(zcav.fs.med.gi$pop), as.numeric(zcav.fs.med.ir)))

zcav.pops.ir<-rbind(zcav.fs.bah.ir, zcav.fs.can.ir, zcav.fs.emed.ir, zcav.fs.wmed.ir)
mden.pops.ir<-rbind(mden.fs.can.ir, mden.fs.bah.ir)

zcav.pops.ir<-as.data.frame(zcav.pops.ir)
mden.pops.ir<-as.data.frame(mden.pops.ir)

zcav.pops.ir$V3<-as.numeric(zcav.pops.ir$V3)
mden.pops.ir$V3<-as.numeric(mden.pops.ir$V3)

ggplot(zcav.pops.ir, aes(x=V2, y=V3)) + 
  geom_boxplot() +
  ylab("Internal Relatedness") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=0.5), axis.text.y=element_text(size=14), axis.title.y = element_text(size=14))
ggplot(mden.pops.ir, aes(x=V2, y=V3)) +
  geom_boxplot() +
  ylab("Internal Relatedness") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=0.5), axis.text.y=element_text(size=14), axis.title.y = element_text(size=14))

#filges to calculate mean and 95% confidence intervals from
zcav.pop.fr.ir #calculated with all individuals together and using the SNP dataset that includes all animals
zcav.pops.ir #calculated using region-specific SNP files

mden.pop.fr.ir
mden.pops.ir

zcav.pop.fr.ir$IR<-as.numeric(zcav.pop.fr.ir$IR)
zcav.pops.ir$V3<-as.numeric(zcav.pops.ir$V3)
mden.pop.fr.ir$IR<-as.numeric(mden.pop.fr.ir$IR)
mden.pops.ir$V3<-as.numeric(mden.pops.ir$V3)
zcav.fs.med.ir$V3<-as.numeric(zcav.fs.med.ir$V3)

library(rcompanion)
zcav.fr.ir.sum<-groupwiseMean(IR ~ population, data=zcav.pop.fr.ir, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
zcav.fs.ir.sum<-groupwiseMean(V3 ~ V2, data=zcav.pops.ir, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
mden.fr.ir.sum<-groupwiseMean(IR ~ population, data=mden.pop.fr.ir, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
mden.fs.ir.sum<-groupwiseMean(V3 ~ V2, data=mden.pops.ir, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)
zcav.med.sum<-groupwiseMean(V3 ~ V2, data=zcav.fs.med.ir, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3)

zcav.med.sum
#V2  n     Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Mediterranean-East 14  0.09700       0.95     0.0627    0.13100           0.0630          0.12100
#2 Mediterranean-West 23 -0.00384       0.95    -0.0127    0.00502          -0.0117          0.00388

zcav.fr.ir.sum
#population  n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#Bahamas  8 0.144       0.95     0.1260      0.163           0.1310            0.159
#Canaries 23 0.107       0.95     0.0828      0.132           0.0833            0.127
#Mediterranean-East 14 0.279       0.95     0.2460      0.312           0.2480            0.302
#Mediterranean-West 23 0.206       0.95     0.1980      0.215           0.1980            0.214
zcav.fs.ir.sum
#V2  n    Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#Bahamas-Abaco  6 -0.0915       0.95    -0.1200   -0.06340          -0.1120          -0.0724
#Bahamas-TOTO  2 -0.0986       0.95    -0.1560   -0.04130          -0.1030          -0.0941
#Canary Islands-East  6  0.0192       0.95    -0.0448    0.08310          -0.0118           0.0698
#Canary Islands-West 17 -0.0458       0.95    -0.0845   -0.00713          -0.0832          -0.0186
#Mediterranean-East 14 -0.0506       0.95    -0.0819   -0.01930          -0.0805          -0.0276
#Mediterranean-West 23 -0.0215       0.95    -0.0315   -0.01160          -0.0309          -0.0133
mden.fr.ir.sum
#population  n   Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#Bahamas  7 0.0403       0.95    0.01760     0.0631          0.02540           0.0585
#Canaries 20 0.0118       0.95   -0.00508     0.0286         -0.00345           0.0293
mden.fs.ir.sum
#V2  n    Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#Bahamas-Abaco  3 -0.1130       0.95    -0.2280    0.00208          -0.1480          -0.0604
#Bahamas-TOTO  4 -0.1470       0.95    -0.1990   -0.09450          -0.1740          -0.1190
#Canary Islands-East  3 -0.0317       0.95    -0.0517   -0.01170          -0.0381          -0.0228
#Canary Islands-West 17 -0.0367       0.95    -0.0591   -0.01430          -0.0561          -0.0146

# focal region div.  stats -----------------------------------------------
# focal region data files
zcav.fr.gl
mden.fr.gl

zcav.fr.gi<-gl2gi(zcav.fr.gl)
mden.fr.gi<-gl2gi(mden.fr.gl)

gl2genepop(zcav.fr.gl, outfile = "./genepop/zcav.fr.genepop")

zcav.fr.stats<-basic.stats(zcav.fr.gi)
mden.fr.stats<-basic.stats(mden.fr.gi)

zcav.100.fs.stats<-basic.stats(zcav.focal.100.gi)

zcav.med.gi
zcav.medonly.stats<-basic.stats(zcav.med.gi)
zcav.med.het<-cbind(zcav.medonly.stats$Ho[,1],zcav.medonly.stats$Hs[,1])
head(zcav.med.het)
colnames(zcav.med.het)<-c("Ho", "Hs")
melt(zcav.med.het)
groupwiseMean(value~Var2, data=melt(zcav.med.het), conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#Var2      n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1   Ho 132534 0.108       0.95      0.107      0.109            0.107            0.109
#2   Hs 132534 0.122       0.95      0.121      0.123            0.121            0.123
wc(zcav.med.gi) #0.1133961
boot.ppfis(zcav.med.gi) #0.1108 0.1154

zcav.fr.gl$pop
zcav.med.gl<-zcav.fr.gl
zcav.med.gl<-gl.recalc.metrics(zcav.med.gl)
zcav.med.gl<-gl.compliance.check(zcav.med.gl)
zcav.med.gl<-gl.drop.pop(zcav.med.gl, pop.list = c("Bahamas", "Canaries"), recalc = TRUE)
zcav.med.gl$pop
zcav.med.gl
zcav.med.gl$other$loc.metrics <- as.data.frame(zcav.med.gl$other$loc.metrics)
zcav.med.gl<-gl.merge.pop(zcav.med.gl, old=c("Mediterranean-East", "Mediterranean-West"), new="Mediterranean")
zcav.med.gl$pop
zcav.med.gl<-gl.compliance.check(zcav.med.gl)
zcav.med.gl<-gl.filter.allna(zcav.med.gl)
zcav.medonly.stats2<-gl.basic.stats(zcav.med.gl)
zcav.medonly.stats2
zcav.med.het2<-cbind(zcav.medonly.stats2$Ho[,1],zcav.medonly.stats2$Hs[,1])
head(zcav.med.het2)
colnames(zcav.med.het2)<-c("Ho", "Hs")
melt(zcav.med.het2)
groupwiseMean(value~Var2, data=melt(zcav.med.het2), conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#Var2      n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1   Ho 125769 0.114       0.95      0.113      0.115            0.113            0.115
#2   Hs 125769 0.128       0.95      0.127      0.129            0.127            0.129
zcav.med.gi2<-gl2gi(zcav.med.gl)
wc(zcav.med.gi2)#0.1138104
boot.ppfis(zcav.med.gi2) #0.1118 0.1164

zcav.fs.med.gl$pop
zcav.fs.med.gl2<-gl.merge.pop(zcav.fs.med.gl, old=c("Mediterranean-East", "Mediterranean-West"), new="Mediterranean")
zcav.fs.med.gl2$pop
zcav.fs.medonly.stats2<-gl.basic.stats(zcav.fs.med.gl2)
zcav.fs.med.het<-cbind(zcav.fs.medonly.stats2$Ho[,1],zcav.fs.medonly.stats2$Hs[,1])
head(zcav.fs.med.het)
colnames(zcav.fs.med.het)<-c("Ho", "Hs")
melt(zcav.fs.med.het)
groupwiseMean(value~Var2, data=melt(zcav.fs.med.het), conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#Var2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1   Ho 27087 0.245       0.95      0.243      0.247            0.243            0.247
#2   Hs 27087 0.264       0.95      0.262      0.266            0.262            0.266

zcav.fs.med.gi2<-gl2gi(zcav.fs.med.gl2)
wc(zcav.fs.med.gi2) #0.07117613
boot.ppfis(zcav.fs.med.gi2) #0.0684 0.0734


zcav.fr.boot.fis<-boot.ppfis(zcav.fr.gi, nboot = 1000)
#ll     hl
#1 0.1517 0.1580
#2 0.1161 0.1202
#3 0.0676 0.0740
#4 0.0711 0.0766
mden.fr.boot.fis<-boot.ppfis(mden.fr.gi, nboot = 1000)
#1 0.0703 0.0803
#2 0.0657 0.0717

zcav.seppop.fr.gi<-seppop(zcav.fr.gi)
zcav.fr.bah.gi<-zcav.seppop.fr.gi$Bahamas
zcav.fr.zcav.seppop.fr.gi$Canaries
zcav.fr.zcav.seppop.fr.gi$`Mediterranean-East`
zcav.fr.zcav.seppop.fr.gi$`Mediterranean-West`

mden.seppop.fr.gi<-seppop(mden.fr.gi)
mden.seppop.fr.gi$Bahamas
mden.seppop.fr.gi$Canaries

wc(zcav.seppop.fr.gi$Bahamas) #0.1548455
wc(zcav.seppop.fr.gi$Canaries) #0.1180368
wc(zcav.seppop.fr.gi$`Mediterranean-East`) #0.07096432
wc(zcav.seppop.fr.gi$`Mediterranean-West`) #0.07389351
wc(mden.seppop.fr.gi$Bahamas) #0.0754514
wc(mden.seppop.fr.gi$Canaries) #0.06862256

zcav.fr.stats.ho<-zcav.fr.stats$Ho
zcav.fr.stats.hs<-zcav.fr.stats$Hs
mden.fr.stats.ho<-mden.fr.stats$Ho
mden.fr.stats.hs<-mden.fr.stats$Hs
zcav.fr.stats.fis<-zcav.fr.stats$Fis
boxplot(zcav.fr.stats.fis)
boxplot(zcav.fr.stats$Ho)
boxplot(zcav.fr.stats$Hs)

zcav.fr.stats.ho<-melt(zcav.fr.stats.ho)
zcav.fr.stats.hs<-melt(zcav.fr.stats.hs)
mden.fr.stats.ho<-melt(mden.fr.stats.ho)
mden.fr.stats.hs<-melt(mden.fr.stats.hs)

zcav.100.fs.stats.ho<-melt(zcav.100.fs.stats$Ho)
zcav.100.fs.stats.hs<-melt(zcav.100.fs.stats$Hs)

zcav.100.ho.sum<-groupwiseMean(value ~ Var2, data=zcav.100.fs.stats.ho, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
zcav.100.hs.sum<-groupwiseMean(value ~ Var2, data=zcav.100.fs.stats.hs, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)

zcav.fr.ho.sum<-groupwiseMean(value ~ Var2, data=zcav.fr.stats.ho, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#Var2      n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1            Bahamas 135301 0.152       0.95      0.151      0.153            0.151            0.153
#2           Canaries 131180 0.155       0.95      0.155      0.156            0.155            0.156
#3 Mediterranean-East 125769 0.106       0.95      0.105      0.107            0.105            0.107
#4 Mediterranean-West 125769 0.118       0.95      0.117      0.119            0.117            0.119
zcav.fr.hs.sum<-groupwiseMean(value ~ Var2, data=zcav.fr.stats.hs, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#Var2      n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1            Bahamas 135301 0.180       0.95      0.179      0.180            0.178            0.180
#2           Canaries 131180 0.176       0.95      0.175      0.177            0.175            0.177
#3 Mediterranean-East 125769 0.114       0.95      0.113      0.115            0.113            0.115
#4 Mediterranean-West 125769 0.128       0.95      0.127      0.129            0.127            0.129
mden.fr.ho.sum<-groupwiseMean(value ~ Var2, data=mden.fr.stats.ho, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#Var2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1  Bahamas 44046 0.232       0.95      0.230      0.234            0.230            0.234
#2 Canaries 46499 0.237       0.95      0.235      0.238            0.235            0.239
mden.fr.hs.sum<-groupwiseMean(value ~ Var2, data=mden.fr.stats.hs, conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#Var2     n  Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1  Bahamas 44046 0.251       0.95      0.249      0.252            0.249            0.253
#2 Canaries 46499 0.254       0.95      0.253      0.256            0.253            0.256



# Allelic Richness --------------------------------------------------------
#calculating alleleic richness
mden.fr.ar<-allelic.richness(mden.fr.gi)
zcav.fr.ar<-allelic.richness(zcav.fr.gi)

zcav.fs.wmed.ar<-allelic.richness(as.matrix(zcav.100.fs.wmed.gl))
zcav.fs.emed.ar<-allelic.richness(zcav.100.fs.emed.gi)
zcav.fs.bah.ar <-allelic.richness(zcav.fs.bah.gi)
zcav.fs.can.ar <-allelic.richness(zcav.fs.can.gi)
zcav.fs.med.ar <-allelic.richness(zcav.fs.med.gi)
mden.fs.bah.ar <-allelic.richness(mden.fs.bah.gi)
mden.fs.can.ar <-allelic.richness(mden.fs.can.gi)

zcav.100.fs.wmed.gi

mden.fs.gi

boxplot(mden.fr.ar$Ar)
#Bahamas         Canaries    
#Min.   :1.000   Min.   :1.000  
#1st Qu.:1.000   1st Qu.:1.462  
#Median :1.934   Median :1.783  
#Mean   :1.682   Mean   :1.703  
#3rd Qu.:1.999   3rd Qu.:1.976  
#Max.   :2.000   Max.   :2.000  
#NA's   :3081    NA's   :628    
boxplot(zcav.fr.ar$Ar)
#Bahamas         Canaries     Mediterranean-East Mediterranean-West
#Min.   :1.000   Min.   :1.000   Min.   :1.000      Min.   :1.000     
#1st Qu.:1.000   1st Qu.:1.182   1st Qu.:1.000      1st Qu.:1.000     
#Median :1.500   Median :1.444   Median :1.000      Median :1.000     
#Mean   :1.483   Mean   :1.482   Mean   :1.276      Mean   :1.311     
#3rd Qu.:1.900   3rd Qu.:1.782   3rd Qu.:1.652      3rd Qu.:1.764     
#Max.   :2.000   Max.   :1.997   Max.   :2.000      Max.   :1.997     
#NA's   :4787    NA's   :8908    NA's   :14319      NA's   :14319     

summary(mden.fr.ar$Ar)
#Bahamas         Canaries    
#Min.   :1.000   Min.   :1.000  
#1st Qu.:1.000   1st Qu.:1.462  
#Median :1.934   Median :1.783  
#Mean   :1.682   Mean   :1.703  
#3rd Qu.:1.999   3rd Qu.:1.976  
#Max.   :2.000   Max.   :2.000  
#NA's   :3081    NA's   :628    

summary(zcav.fr.ar$Ar)
#Bahamas         Canaries     Mediterranean-East Mediterranean-West
#Min.   :1.000   Min.   :1.000   Min.   :1.000      Min.   :1.000     
#1st Qu.:1.000   1st Qu.:1.182   1st Qu.:1.000      1st Qu.:1.000     
#Median :1.500   Median :1.444   Median :1.000      Median :1.000     
#Mean   :1.483   Mean   :1.482   Mean   :1.276      Mean   :1.311     
#3rd Qu.:1.900   3rd Qu.:1.782   3rd Qu.:1.652      3rd Qu.:1.764     
#Max.   :2.000   Max.   :1.997   Max.   :2.000      Max.   :1.997     
#NA's   :4787    NA's   :8908    NA's   :14319      NA's   :14319   

mden.fr.ar.sum<-groupwiseMean(value ~ variable, data=melt(mden.fr.ar$Ar), conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#variable     n Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1  Bahamas 44046 1.68       0.95       1.68       1.69             1.68             1.69
#2 Canaries 46499 1.70       0.95       1.70       1.71             1.70             1.71
zcav.fr.ar.sum<-groupwiseMean(value ~ variable, data=melt(zcav.fr.ar$Ar), conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#variable      n Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1            Bahamas 135301 1.48       0.95       1.48       1.48             1.48             1.48
#2           Canaries 131180 1.48       0.95       1.48       1.48             1.48             1.48
#3 Mediterranean-East 125769 1.28       0.95       1.27       1.28             1.27             1.28
#4 Mediterranean-West 125769 1.31       0.95       1.31       1.31             1.31             1.31

groupwiseMean(value~variable, data=melt(zcav.fs.wmed.ar$Ar), conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
groupwiseMean(value~variable, data=melt(zcav.fs.emed.ar$Ar), conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
groupwiseMean(value~variable, data=melt (zcav.fs.bah.ar$Ar ), conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#variable     n Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Bahamas-Abaco 61643 1.50       0.95       1.50       1.51             1.50             1.51
#2  Bahamas-TOTO 61643 1.49       0.95       1.48       1.49             1.48             1.49
groupwiseMean(value~variable, data=melt (zcav.fs.can.ar$Ar ), conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#variable     n Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Canary Islands-East 56246 1.63       0.95       1.62       1.63             1.62             1.63
#2 Canary Islands-West 56246 1.63       0.95       1.63       1.63             1.63             1.63
groupwiseMean(value~variable, data=melt (zcav.fs.med.ar$Ar ), conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#variable     n Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Mediterranean-East 27087 1.77       0.95       1.76       1.77             1.76             1.77
#2 Mediterranean-West 27087 1.88       0.95       1.87       1.88             1.87             1.88
groupwiseMean(value~variable, data=melt (mden.fs.bah.ar$Ar ), conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#variable     n Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Bahamas-Abaco 26890 1.75       0.95       1.75       1.76             1.75             1.76
#2  Bahamas-TOTO 26890 1.72       0.95       1.71       1.72             1.71             1.72
groupwiseMean(value~variable, data=melt (mden.fs.can.ar$Ar ), conf=0.95, R=1000, percentile=TRUE, bca=FALSE, digits=3, na.rm=TRUE)
#variable     n Mean Conf.level Trad.lower Trad.upper Percentile.lower Percentile.upper
#1 Canary Islands-East 28911  1.6       0.95        1.6       1.61              1.6             1.61
#2 Canary Islands-West 28911  1.6       0.95        1.6       1.60              1.6             1.60

# Effective population size -----------------------------------------------

#Calculating Ne
devtools::install_github(repo="zakrobinson/RLDNe")
library(RLDNe)

#workflow to calculate ne
data("wgp_example")
wgp_example
#write the LDNe Parameter file for the LD method in NeEstimator V2.1 using the genepop file. 
param_files<- NeV2_LDNe_create(input_file = gp_file$Output_File ,param_file = "Ne_params.txt" ,NE_out_file = "Ne_out.txt")
#function to run LDNe from NeEstimator 2.1. 
run_LDNe(LDNe_params = param_files$param_file)
Ne_estimates<-readLDNe_tab(path = param_files$Ne_out_tab)

gl2genepop(mden.fr.gl, outfile = "./mden.fr.genepop")
mden_param_files<- NeV2_LDNe_create(input_file = "./mden.fr.genepop" ,param_file = "mden_Ne_params.txt" ,NE_out_file = "mden_Ne_out.txt")
#manually opened and changed the mating system to random
run_LDNe(LDNe_params = mden_param_files$param_file)
mden_Ne_estimates<-readLDNe_tab(path = mden_param_files$Ne_out_tab)
mden_Ne_estimates

zcav_param_files<-NeV2_LDNe_create(input_file = "./zcav.fr.genepop" ,param_file = "zcav_Ne_params.txt" ,NE_out_file = "zcav_Ne_out.txt")
run_LDNe(LDNe_params = zcav_param_files$param_file)


gl2genepop(zcav.fr.gl, outfile="./zcav.fr.genepop")
gl2genepop(zcav.fs.bah.gl,   outfile=  "./zcav.fs.bah.genepop")
gl2genepop(zcav.fs.can.gl,   outfile=  "./zcav.fs.can.genepop")
gl2genepop(zcav.100.emed.gl, outfile=  "./zcav.100.emed.genepop")
gl2genepop(zcav.100.wmed.gl, outfile=  "./zcav.100.wmed.genepop")
gl2genepop(zcav.fs.med.gl,   outfile=  "./zcav.fs.med.genepop")
gl2genepop(mden.fs.bah.gl,   outfile=  "./mden.fs.bah.genepop")
gl2genepop(mden.fs.can.gl,   outfile=  "./mden.fs.can.genepop")

zc.bah
zc.can
zc.emd
zc.wmd
zc.med
md.bah
md.can

zc.bah_param_files<- NeV2_LDNe_create(input_file = "./zcav.fs.bah.genepop", param_file = "zc.bah_Ne_params.txt" ,NE_out_file = "zc.bah_Ne_out.txt")
zc.can_param_files<- NeV2_LDNe_create(input_file = "./zcav.fs.can.genepop", param_file = "zc.can_Ne_params.txt" ,NE_out_file = "zc.can_Ne_out.txt")
zc.emd_param_files<- NeV2_LDNe_create(input_file = "./zcav.100.emed.genepop", param_file = "zc.emd_Ne_params.txt" ,NE_out_file = "zc.emd_Ne_out.txt")
zc.wmd_param_files<- NeV2_LDNe_create(input_file = "./zcav.100.wmed.genepop" ,param_file = "zc.wmd_Ne_params.txt" ,NE_out_file = "zc.wmd_Ne_out.txt")
zc.med_param_files<- NeV2_LDNe_create(input_file = "./Genepop files/zcav.fs.med.genepop", param_file = "zc.med_Ne_params.txt" ,NE_out_file = "zc.med_Ne_out.txt")
md.bah_param_files<- NeV2_LDNe_create(input_file = "./mden.fs.bah.genepop", param_file = "md.bah_Ne_params.txt" ,NE_out_file = "md.bah_Ne_out.txt")
md.can_param_files<- NeV2_LDNe_create(input_file = "./mden.fs.can.genepop", param_file = "md.can_Ne_params.txt" ,NE_out_file = "md.can_Ne_out.txt")

run_LDNe(LDNe_params = zc.bah_param_files$param_file)
run_LDNe(LDNe_params = zc.can_param_files$param_file)
run_LDNe(LDNe_params = zc.emd_param_files$param_file)
run_LDNe(LDNe_params = zc.wmd_param_files$param_file)
run_LDNe(LDNe_params = zc.med_param_files$param_file)
run_LDNe(LDNe_params = md.bah_param_files$param_file)
run_LDNe(LDNe_params = md.can_param_files$param_file)

zc.bah_Ne_estimates<-readLDNe_tab(path = zc.bah_param_files$Ne_out_tab)

zc.bah_Ne_estimates

zcav.ne.fr.gl<-zcav.fr.gl
zcav.ne.fr.gl$pop<-as.factor(c("Bahamas",
                                "Bahamas",
                                "Bahamas",
                                "Bahamas",
                                "Bahamas",
                                "Bahamas",
                                "Bahamas",
                                "Bahamas",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "Canaries",
                                "EMedtierranean",
                                "EMedtierranean",
                                "EMedtierranean",
                                "EMedtierranean",
                                "EMedtierranean",
                                "EMedtierranean",
                                "EMedtierranean",
                                "EMedtierranean",
                                "EMedtierranean",
                                "EMedtierranean",
                                "EMedtierranean",
                                "EMedtierranean",
                                "EMedtierranean",
                                "EMedtierranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean",
                                "Wmediterranean"))
zcav.ne.fr.gl$pop
zcav.ne.fr.gl$ind.names<-as.character(c(1:68))
zcav.ne.fr.gl$ind.names
zcav.LDNe<-gl.LDNe(zcav.ne.fr.gl[,1:1000], outpath=getwd(), critical=c(0,0.05), neest.path = "./Neestimator/",mating="random")
zcav.LDNe

zcav.LDNe.10k<-gl.LDNe(zcav.ne.fr.gl[,1:10000], outpath=getwd(), critical=c(0,0.05), neest.path = "./Neestimator/",mating="random")
zcav.LDNe.10k

mden.LDNe.10k<-gl.LDNe(mden.fr.gl[,1:10000],outpath=getwd(), critical=c(0,0.05), neest.path = "./Neestimator/",mating="random")
mden.LDNe.10k
mden.fs.LDNe.10k<-gl.LDNe(mden.ne.fs.gl[,1:10000], outpath=getwd(), critical=c(0,0.05), neest.path = "./Neestimator/",mating="random")
mden.fs.LDNe.20k<-gl.LDNe(mden.ne.fs.gl[,1:20000], outpath=getwd(), critical=c(0,0.02, 0.05), neest.path = "./Neestimator/",mating="random")

mden.bah.LDne
mden.can.LDne

mden.ne.fs.gl<-mden.fs.gl
mden.ne.fs.gl$pop<-as.factor(c("Abaco",
                               "Abaco",
                               "Abaco",
                               "TOTO",
                               "TOTO",
                               "TOTO",
                               "TOTO",
                               "Ecan",
                               "Ecan",
                               "Ecan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan",
                               "Wcan"))
mden.ne.fs.gl$ind.names<-as.character(c(1:27))

zcav.ne.fs.med.gl<-zcav.fs.med.gl
zcav.ne.fs.med.gl$pop<-as.factor(c("EMedtierranean",
                         "EMedtierranean",
                         "EMedtierranean",
                         "EMedtierranean",
                         "EMedtierranean",
                         "EMedtierranean",
                         "EMedtierranean",
                         "EMedtierranean",
                         "EMedtierranean",
                         "EMedtierranean",
                         "EMedtierranean",
                         "EMedtierranean",
                         "EMedtierranean",
                         "EMedtierranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean",
                         "Wmediterranean"))
zcav.ne.fs.med.gl$pop
zcav.ne.fs.med.gl$ind.names<-as.character(c(1:37))
zcav.med.LDNe.10k<-gl.LDNe(zcav.ne.fs.med.gl[,1:10000],outpath=getwd(), critical=c(0,0.05), neest.path = "./Neestimator/",mating="random")

mden.100.fs.gl$pop
zcav.100.fs.gl$pop

mden.100.fs.10k.gl<-gl.subsample.loci(mden.100.fs.gl, n=10000, method="random")
zcav.100.fs.10k.gl<-gl.subsample.loci(zcav.100.fs.gl, n=10000, method="random")
mden.fr.gl<-gl.compliance.check(mden.fr.gl)
mden.100.fr.gl<-gl.recalc.metrics(mden.fr.gl)
mden.100.fr.gl@other$loc.metrics<-as.data.frame(mden.100.fr.gl@other$loc.metrics)
mden.100.fr.gl<-gl.compliance.check(mden.100.fr.gl)
mden.100.fr.gl<-gl.filter.callrate(mden.100.fr.gl, method="loc", threshold=1, mono.rm = TRUE)
mden.100.fr.10k.gl<-gl.subsample.loci(mden.100.fr.gl, n=10000, )

zcav.100.fs.10k.gl$ind.names<-as.character(c(1:68))
zcav.100.fs.10k.gl$pop<-as.factor(c("Ecan",
                                    "Wmed",
                                    "Wcan",
                                    "Wcan",
                                    "Ecan",
                                    "Emed",
                                    "Wcan",
                                    "Ecan",
                                    "Ecan",
                                    "Wcan",
                                    "Wcan",
                                    "Emed",
                                    "Wmed",
                                    "Emed",
                                    "Emed",
                                    "Wmed",
                                    "Emed",
                                    "Ecan",
                                    "AbacoBAH",
                                    "Wmed",
                                    "Emed",
                                    "AbacoBAH",
                                    "Wcan",
                                    "Wcan",
                                    "Emed",
                                    "Wmed",
                                    "Wcan",
                                    "Emed",
                                    "Wcan",
                                    "Wcan",
                                    "Wmed",
                                    "Emed",
                                    "AbacoBAH",
                                    "Wmed",
                                    "Emed",
                                    "AbacoBAH",
                                    "AbacoBAH",
                                    "Wmed",
                                    "Wcan",
                                    "Emed",
                                    "Wcan",
                                    "Wmed",
                                    "Emed",
                                    "Wmed",
                                    "Wmed",
                                    "Wmed",
                                    "Wcan",
                                    "Emed",
                                    "Wmed",
                                    "Wmed",
                                    "Wcan",
                                    "Wcan",
                                    "AbacoBAH",
                                    "Wmed",
                                    "Wcan",
                                    "Wmed",
                                    "Wmed",
                                    "Wmed",
                                    "Wmed",
                                    "Wmed",
                                    "Wmed",
                                    "TOTOBAH",
                                    "Emed",
                                    "Wmed",
                                    "TOTOBAH",
                                    "Wcan",
                                    "Ecan",
                                    "Wmed"))

mden.100.fr.LDNe.10k<-gl.LDNe(mden.100.fr.10k.gl,outpath=getwd(), critical=c(0,0.01,0.02,0.05), neest.path = "./Neestimator/",mating="random")
mden.100.fs.LDNe.10k<-gl.LDNe(mden.100.fs.10k.gl,outpath=getwd(), critical=c(0,0.01,0.02,0.05), neest.path = "./Neestimator/",mating="random")
zcav.100.fs.LDNe.10k<-gl.LDNe(zcav.100.fs.10k.gl,outpath=getwd(), critical=c(0,0.01,0.02,0.05), neest.path = "./Neestimator/",mating="random")
zcav.100.fr.LDNe.10k<-gl.LDNe(zcav.100.fr.10k.gl,outpath=getwd(), critical=c(0,0.01,0.02,0.05), neest.path = "./Neestimator/",mating="random")

zcav.100.med.10k.gl<-gl.subsample.loci(zcav.100.med.gl, n=10000, method="random")
zcav.100.med.LDNe.10k<-gl.LDNe(zcav.100.med.10k.gl,outpath=getwd(), critical=c(0,0.01,0.02,0.05), neest.path = "./Neestimator/",mating="random")

zcav.100.fr.10k.gl<-zcav.100.fs.10k.gl
zcav.100.fr.10k.gl$pop<-as.factor(c("CAN",
                                    "Wmed",
                                    "CAN",
                                    "CAN",
                                    "CAN",
                                    "Emed",
                                    "CAN",
                                    "CAN",
                                    "CAN",
                                    "CAN",
                                    "CAN",
                                    "Emed",
                                    "Wmed",
                                    "Emed",
                                    "Emed",
                                    "Wmed",
                                    "Emed",
                                    "CAN",
                                    "BAH",
                                    "Wmed",
                                    "Emed",
                                    "BAH",
                                    "CAN",
                                    "CAN",
                                    "Emed",
                                    "Wmed",
                                    "CAN",
                                    "Emed",
                                    "CAN",
                                    "CAN",
                                    "Wmed",
                                    "Emed",
                                    "BAH",
                                    "Wmed",
                                    "Emed",
                                    "BAH",
                                    "BAH",
                                    "Wmed",
                                    "CAN",
                                    "Emed",
                                    "CAN",
                                    "Wmed",
                                    "Emed",
                                    "Wmed",
                                    "Wmed",
                                    "Wmed",
                                    "CAN",
                                    "Emed",
                                    "Wmed",
                                    "Wmed",
                                    "CAN",
                                    "CAN",
                                    "BAH",
                                    "Wmed",
                                    "CAN",
                                    "Wmed",
                                    "Wmed",
                                    "Wmed",
                                    "Wmed",
                                    "Wmed",
                                    "Wmed",
                                    "BAH",
                                    "Emed",
                                    "Wmed",
                                    "BAH",
                                    "CAN",
                                    "CAN",
                                    "Wmed"))

# IR and div stats plots --------------------------------------------------

par(mfrow=c(2,2))
boxplot(zcav.fr.stats.ho$value ~ zcav.fr.stats.ho$Var2)
boxplot(zcav.fr.stats.hs$value ~ zcav.fr.stats.hs$Var2)
boxplot(zcav.pops.ir$V3 ~ zcav.pops.ir$V2)
boxplot(zcav.pop.fr.ir$IR ~ zcav.pop.fr.ir$population)

boxplot(mden.fr.stats.ho$value ~ mden.fr.stats.ho$Var2)
boxplot(mden.fr.stats.hs$value ~ mden.fr.stats.hs$Var2)
boxplot(mden.pops.ir$V3 ~ mden.pops.ir$V2)
boxplot(mden.pop.fr.ir$IR ~ mden.pop.fr.ir$population)


# DAPC a priori and within focal regions ----------------------------------

#a priori
zcav.fs.gl$pop
zcav.fs.gi<-gl2gi(zcav.fs.gl)

zcav.grp<-find.clusters(zcav.fs.gi, max.n.clust=20) #80 PCs, 6 clusters
zcav.dapc<-xvalDapc(tab(zcav.fs.gl), zcav.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(zcav.dapc$DAPC)
cbind(zcav.grp$grp, zcav.fs.gl$ind.names, as.character(zcav.fs.gl$pop))

zcav.bah.grp<-find.clusters(zcav.fs.bah.gi) 
zcav.bah.dapc<-xvalDapc(tab(zcav.fs.bah.gl), zcav.bah.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(zcav.bah.dapc$DAPC)
cbind(zcav.bah.grp$grp, zcav.fs.bah.gl$ind.names, as.character(zcav.fs.bah.gl$pop))

zcav.can.grp<-find.clusters(zcav.fs.can.gi) 
zcav.can.dapc<-xvalDapc(tab(zcav.fs.can.gl), zcav.can.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(zcav.can.dapc$DAPC)
cbind(zcav.can.grp$grp, zcav.fs.can.gl$ind.names, as.character(zcav.fs.can.gl$pop))

zcav.med.grp<-find.clusters(zcav.fs.med.gi)
zcav.med.dapc<-xvalDapc(tab(zcav.fs.med.gl), zcav.med.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(zcav.med.dapc$DAPC)
cbind(zcav.med.grp$grp, zcav.fs.med.gl$ind.names, as.character(zcav.fs.med.gl$pop))

mden.bah.grp<-find.clusters(mden.fs.bah.gi) 
mden.bah.dapc<-xvalDapc(tab(mden.fs.bah.gl), mden.bah.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(mden.bah.dapc$DAPC)
cbind(mden.bah.grp$grp, mden.fs.bah.gl$ind.names, as.character(mden.fs.bah.gl$pop))

mden.can.grp<-find.clusters(mden.fs.can.gi) 
mden.can.dapc<-xvalDapc(tab(mden.fs.can.gl), mden.can.grp$grp, n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
scatter(mden.can.dapc$DAPC)
cbind(mden.can.grp$grp, mden.fs.can.gl$ind.names, as.character(mden.fs.can.gl$pop))


#de Novo populations DAPC in focal regions
zcav.bah.dapc2<-xvalDapc(tab(zcav.fs.bah.gl), pop(zcav.fs.bah.gl), n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
zcav.can.dapc2<-xvalDapc(tab(zcav.fs.can.gl), pop(zcav.fs.can.gl), n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
zcav.med.dapc2<-xvalDapc(tab(zcav.fs.med.gl), pop(zcav.fs.med.gl), n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
mden.bah.dapc2<-xvalDapc(tab(mden.fs.bah.gl), pop(mden.fs.bah.gl), n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)
mden.can.dapc2<-xvalDapc(tab(mden.fs.can.gl), pop(mden.fs.can.gl), n.pca.max=300, training.set=0.9, result="groupMean", center=TRUE, scale=FALSE, n.pca=NULL, n.rep=30, xval.plot=TRUE)

zcav.fs.can.gl
mden.fs.can.gl

scatter(zcav.bah.dapc2$DAPC)
scatter(zcav.can.dapc2$DAPC)
scatter(zcav.med.dapc2$DAPC)
scatter(mden.bah.dapc2$DAPC)
scatter(mden.can.dapc2$DAPC)

assignplot(zcav.bah.dapc2$DAPC)
assignplot(zcav.can.dapc2$DAPC)
assignplot(zcav.med.dapc2$DAPC)
assignplot(mden.bah.dapc2$DAPC)
assignplot(mden.can.dapc2$DAPC)

par(mar=c(4,7,2,2))
par(mfrow = c(1, 2))
dev.off()

# filtering vcf files based on final genlight files ---------------------------------------------------------
mden.fr.gl #27 genotypes,  47,127 binary SNPs
zcav.fr.gl #68 genotypes,  140,088 binary SNPs

cbind(mden.fr.gl$ind.names, as.character(mden.fr.gl$pop))
cbind(zcav.fr.gl$ind.names, as.character(zcav.fr.gl$pop))

./populations -V ~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/final_pop/mden_pop_final.vcf -M ~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/final_pop/final_mden_popmap.txt -O ~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/final_pop/mden/ --whitelist ~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/final_pop/mden.fr.loci.txt -p 1 --vcf
./populations -V ~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/final_pop/zcav_pop_final.vcf -M ~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/final_pop/final_zcav_popmap2.txt -O ~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/final_pop/zcav2/ --whitelist ~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/final_pop/zcav.fr.loci.txt -p 1 --vcf
  
cut -f 3 zcav_pop_final.vcf > whitelist
sed -e 's/_/\t/g' whitelist > whitelistv2
sed '1,15d' whitelistv2 > whitelistv3

write.csv(as.character(zcav.fr.gl$loc.names), file="zcav.fr.loci.txt")
write.csv(as.character(mden.fr.gl$loc.names), file="mden.fr.loci.txt")

vcftools --vcf mden_pop_final.p.snps.vcf --missing-indv
vcftools --vcf mden_pop_final.p.snps.vcf --site-depth

vcftools --vcf zcav_pop_final.p.snps.vcf --missing-indv
vcftools --vcf zcav_pop_final.p.snps.vcf --site-depth



# gstudio -----------------------------------------------------------------
remotes::install_github("dyerlab/popgraph")
remotes::install_github("dyerlab/gstudio")

library(popgraph)
library(gstudio)

data("arapat")
genetic_diversity(arapat, mode="Ho", small.N=TRUE)

mden.fr.df[1:5,1:5]

mden.fr.df<-tab(mden.seppop.fr.gi$Bahamas)
mden.fr.df<-as.data.frame(mden.fr.df)
genetic_diversity(mden.fr.df, mode="Ho", small.N=TRUE)

# exploring indv with different loci --------------------------------------

zc.opt.vcfR<-read.vcfR("~/Dropbox/Phd/Bioinformatics/bw_ddrad_focal/ddrad_focal_R/vcf_files/opt/zcav/vcf-files/zc_zc_comb_sorted.vcf", verbose = TRUE)
zc.opt.gl<-vcfR2genlight(zc.opt.vcfR)
zc.opt.gl #46 genotypes,  64,093 binary SNPs, size: 7.9 Mb, 169275 (5.74 %) missing data
zc.opt.gl<-gl.compliance.check(zc.opt.gl)
zc.opt.100.gl<-gl.filter.callrate(zc.opt.gl, method="loc", threshold=1, v=3, mono.rm = TRUE)
zc.opt.100.gl#46 genotypes,  19,943 binary SNPs, size: 5.3 Mb, 0 (0 %) missing data

zc.opt.dr<-gl2demerelate(zc.opt.100.gl)
zc.opt.dr$'Sample-ID'<-zc.opt.dr$'Sample-ID' %>% str_replace_all(c("-"=".", "_"="."))

#install demerelate from source code downloaded from CRAN
zc.opt.D<-Demerelate(zc.opt.dr, object=TRUE, file.output=TRUE, pairs=1000, p.correct=TRUE)

#genetic distance matrix
zc.opt.dist<-gl.dist.ind(zc.opt.100.gl)
zc.opt.dist.mat<-as.matrix(zc.opt.dist)
write.csv(x = zc.opt.dist.mat, "./matrix")

rownames(zc.opt.dist.mat)
heatmap(zc.opt.dist.mat)

ggplot(zc.opt.dist.mat) + 
  geom_tile()

# Maps --------------------------------------------------------------------
devtools::install_github("MikkoVihtakari/ggOceanMapsData")
devtools::install_github("MikkoVihtakari/ggOceanMaps")

library(ggOceanMapsData)
library(ggOceanMaps)

basemap(limits = c(-80, -75, 23, 28), bathymetry = TRUE)
basemap(limits = c(-80, -75, 23, 28), bathymetry = TRUE, bathy.style = "poly_greys")
basemap(limits = c(0, 46, 70, 81), bathymetry = TRUE, bathy.style = "contour_blues")
