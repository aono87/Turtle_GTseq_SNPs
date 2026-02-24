# Merging microsatellite and gtseq genotypes
library(dplyr)
library(strataG)

#gt.seq genotypes
geno.table


#microsat genotypes
msat.table<-read_excel("data-raw/msats/Dcor_Msat_Genotypes.xlsx")

#Add column and zpad the values in Lab_ID
msat.table<-msat.table %>% 
  mutate(Indiv=Lab_ID) %>%
  relocate(Indiv)

msat.table$Indiv <- sprintf("z%07d", msat.table$Indiv)

str(msat.table)
str(geno.table)

#merge tables, keeping only msat values for individuals that match the gtseq genotypes sheet
all.genotypes<-left_join(geno.table, msat.table, by="Indiv")

