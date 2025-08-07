##Script to compare individuals for duplicates

library(tidyverse)
library(strataG)
library(dplyr)
library(swfscMisc)

#gtypes file
g

#compare all individuals with eachother with 80% similarity
all.dups<-dupGenotypes(g, num.shared = 0.8)
all.dups
#     ids.1     ids.2         strata.1         strata.2     mismatch.loci   num.loci.genotyped  num.loci.shared   prop.loci.shared
#1  z0012504  z0012514  SOLOMON_ISLANDS  SOLOMON_ISLANDS          <NA>                158             158                1
#2  z0026457  z0040945 PAPUA_NEW_GUINEA PAPUA_NEW_GUINEA          <NA>                157             157                1
#3 z0012509b z0012509c  SOLOMON_ISLANDS  SOLOMON_ISLANDS          <NA>                147             147                1
#4  z0012496  z0061464  SOLOMON_ISLANDS  SOLOMON_ISLANDS          <NA>                144             144                1
#5  z0012509  z0012520  SOLOMON_ISLANDS  SOLOMON_ISLANDS          <NA>                140             140                1
#6  z0026450  z0040943 PAPUA_NEW_GUINEA PAPUA_NEW_GUINEA          <NA>                139             139                1
write.csv(all.dups, "results-raw/all.duplicate.genotypes.csv")

#Forgot to remove 40943 and 40945 after merging with samples. will remove and re-run
#ids.1    ids.2         strata.1         strata.2 mismatch.loci num.loci.genotyped num.loci.shared prop.loci.shared
#1 z0026457 z0040945 PAPUA_NEW_GUINEA PAPUA_NEW_GUINEA          <NA>                160             160                1
#2 z0026450 z0040943 PAPUA_NEW_GUINEA PAPUA_NEW_GUINEA          <NA>                139             139                1

dups.50<-dupGenotypes(g, num.shared = 0.5)
write.csv(dups.50, "results-raw/duplicate.genotypes.50.csv")
