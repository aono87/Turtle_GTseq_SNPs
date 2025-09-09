#Notes and things to change

#Questionable haps
#Need to figure out a place to evaluate and address questionable haplotypes

# Karen's code to change genotypes for questionable haps
genos_to_change <- read.csv('data-raw/genos_to_change.csv')
for (i in 1:nrow(genos_to_change)){
  idx <- which(tgt$locus == genos_to_change$locus[i] & tgt$Indiv == genos_to_change$Indiv[i])
  tgt$gt[idx] <- genos_to_change$gt[i]
}


#Replicate samples
#Need to figure out a place and way to evaluate and address replicated samples

#Writing a new vcf file
#When and where?

#Allelic balance vs Allelic ratio

#Ideas
#Run through Part 1 and Part 2 with all replicates. 
#After part 2, evaluate replicates, merge replicates
#Run through Part 1 and Part 2
#After part 2, evaluate questionable haplotypes
#Run part 3 and 4

#Change all "messages" to just commented lines. not really possible to automate this. 

#zpadding samples
