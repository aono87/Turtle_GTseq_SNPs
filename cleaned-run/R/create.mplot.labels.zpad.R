#Prepare sample IDs with z padding
library(swfscMisc)
inds<-read.csv("data-raw/mplot_labels/w-leatherback-inds.txt", sep = "\t", header=FALSE)
head(inds)
dim(inds)
ind.name<-inds$V2
ind.name
ind.name.norep<-ind.name %<>%
  gsub("b", "", .) %>%
  gsub("c", "", .) 
ind.name.norep.pad<-zero.pad(as.numeric(ind.name.norep))
ind.name.norep.pad.z<-paste("z0", ind.name.norep.pad, sep="")
ind.name.norep.pad.z
new.inds<-cbind(inds, ind.name.norep.pad.z)
write.table(new.inds, "data-raw/mplot_labels/w-leatherback-inds-zpad.txt", sep="\t", col.names =FALSE, row.names = FALSE)
#Manually go into file and add the "b" and "c" for duplicates. 
#Manually reformat into required format for microhaploty: filename, ID, pop (NA here)