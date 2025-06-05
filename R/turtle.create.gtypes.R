library(tidyverse)
library(strataG)
library(dplyr)
library(swfscMisc)
library(readxl)

project <- "dcor.wpac.dupe-check"
min.reads <- 20

#Create a strata file for the geno.table to make a gtypes file
#export final sample ID list
Dcor.strata<-data.frame(Indiv=geno.table$Indiv)
Dcor.strata$Genotype<-"Yes"
head(Dcor.strata)

#read in sample info from qa.qc file
Dcor.sample.info<-read_xlsx("results-raw/turtle.qa.qc.xlsx", sheet="SampleData")

#Select desired columns from sample info
Dcor.sample.info<-Dcor.sample.info %>% select(mplot.id,	id,	lab.id,	sampling.location,	Turtle_ID,	Lab_ID,	Field_ID,	Year_collected,	Location,	Country,	Sex,	Collection_Method)
head(Dcor.sample.info)

#Merge final sample list with sample info sheet
Dcor.strata$Indiv #z coded sample id with b and c for 12509
Dcor.sample.info$mplot.id #z coded sample id with bs and cs
Dcor.sample.info<-Dcor.sample.info %>% rename(Indiv=mplot.id)
head(Dcor.sample.info)    

Dcor.strata.all <- full_join(Dcor.strata, Dcor.sample.info, by="Indiv")
head(Dcor.strata.all)

#Filter strata sheet to only include samples with final genotypes
Dcor.strata.all<- filter(Dcor.strata.all, Genotype=="Yes")
head(Dcor.strata.all)

genos <- column_to_rownames(geno.table, var = "Indiv")

split.genos <- alleleSplit(genos, sep= "/") %>% 
  data.frame() %>%
  rownames_to_column(var = "Indiv")

df <- right_join(Dcor.strata.all, split.genos, by = "Indiv")
head(df)

df.strata <- select(df, c(Indiv, sampling.location)) %>%
  column_to_rownames(var = "Indiv")

g <- df2gtypes(df, ploidy = 2, id.col = 1, strata.col = 5, loc.col = 14, schemes=df.strata)
g

save(g, df, file = paste0("data/gtypes_", project, "_minReads.", min.reads, ".rda"))
