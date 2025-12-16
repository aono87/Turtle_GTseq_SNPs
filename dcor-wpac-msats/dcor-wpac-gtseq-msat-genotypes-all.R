##Import genotypes GTseq generated data (50 read minimum to call genotype)
#generate rplots using terminal
#rscript Mplot.R Output_all Output_all/


#Importing data
gtseq_genotypes_wide <- read_table("megasat-wreps-minreads50/Output_all/Genotype.txt")
head(gtseq_genotypes_wide)
#remove last column, which is empty
gtseq_genotypes_wide<-gtseq_genotypes_wide[,-68]

#Fix the sample ID column , should read "LabID"
gtseq_genotypes_wide<- gtseq_genotypes_wide %>% rename(LabID=Sample_idx1_idx2)
head(gtseq_genotypes_wide)

#Updated wide to long format
gtseq_genos<-gtseq_genotypes_wide
gtseq_genos <- gtseq_genos %>%
  pivot_longer(
    cols = -c(LabID),
    names_to = "raw_locus_name",
    values_to = "value"
  ) %>% mutate(
    allele_num = if_else(str_ends(raw_locus_name, "-b"), 2, 1), # Create the allele number: if the name ends with "-b", it's 2, otherwise 1.
    Locus = str_remove(raw_locus_name, "GT_"), Locus = str_remove(Locus, "-b") # Create the clean Locus name: remove the "GT_" prefix and then the "-b" suffix
  ) %>%
  pivot_wider(
    id_cols = c(LabID,Locus),
    names_from = allele_num,
    values_from = value
  )

#Cleanup dataset by replacing blanks, x and 0s with NA
#We want to be sure that missing data is not classified as a "0" later on for calculations
gtseq_genos<-replace(gtseq_genos, gtseq_genos==0, NA)
gtseq_genos<-replace(gtseq_genos, gtseq_genos=='X', NA)
gtseq_genos<-replace(gtseq_genos, gtseq_genos=='Unscored', NA)
head(gtseq_genos)

#Look at distribution of alleles in gtseq dataset
gtseq_genos_longer <- gtseq_genos %>%
  pivot_longer(
    cols = c(`1`, `2`), 
    names_to = "Allele_Type", 
    values_to = "Allele_Value"
  ) %>%
  filter(!is.na(Allele_Value)) %>% # Remove rows where Allele_Value is NA
  mutate(Allele_Value = as.numeric(Allele_Value))# Remove rows where Allele_Value is NA

# histogram plot
# We map the allele values to the x-axis.
# facet_wrap(~ Locus) creates a separate plot for each unique Locus.
# `scales = "free"` allows each plot to have its own x and y axis scales, which is useful
# if the allele value ranges differ significantly between loci.
ggplot(gtseq_genos_longer, aes(x = Allele_Value)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black", alpha = 0.8) +
  facet_wrap(~ Locus, scales = "free") +
  labs(
    title = "Distribution of GTseq Allele Values by Locus",
    subtitle = "Alleles 1 and 2 are combined for frequency count",
    x = "Allele Value",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "lightgray", color = "black"),
    strip.text = element_text(face = "bold")
  )

#Convert genotypes >200bp to NA (value chosen based on histogram of allele lengths)
#Megasat will genotype any allele if the forward primer and forward flanking region are present
#If the sequence cuts off in the middle of the reverse flanking region, the program will 
#automatically add the rest of the reverse flanking region.
#But, if the microsatellite is longer than the sequencing chemistry allows, and it 
#cuts off before the reverse flanking region, it will give an arbitrarily large size (>300bp) fragment
#and it is impossible to know what the true genotype would be in this circumstance
#Need to turn any genotype with a >200bp allele to NA for both alleles
gtseq_genos <- gtseq_genos %>%
  mutate(across(c(`1`, `2`), as.numeric)) %>%
  mutate(
    `1` = ifelse(`1` > 200 | `2` > 200, NA, `1`),
    `2` = ifelse(`1` > 200 | `2` > 200, NA, `2`)
  )

#use code above to make new histograms to be sure that the desired alleles have been removed


##Check for mismatches in replicated individuals
source("./R/Compare.replicates.msats.R")
##This function reports mismatched genotypes including NAs
na_mismatches_all_locs <- na_mismatches(
  df_long = gtseq_genos,
  ind_col = "LabID",
  locus_col = "Locus",
  allele_cols = c("1", "2")
)
##This function reports only mismatched genotypes that don't include NAs
no_na_mismatches_all_locs<- no_na_mismatches(
  df_long = gtseq_genos,
  ind_col = "LabID",
  locus_col = "Locus",
  allele_cols = c("1", "2")
)

write.csv(na_mismatches_all_locs, file = "./results/na_gtseq_mismatches_after-200-filter-all-locs.csv")
write.csv(no_na_mismatches_all_locs, file = "./results/no_na_gtseq_mismatches_after-200-filter-all-locs.csv")

##After updating genotypes for mismatches, update the correct genotypes and 
#amend the gtseq_genotypes file
new.genos<-read.csv("results/genos.to.change.csv")

#Merge new genotypes with original genotype file
genos <- gtseq_genos %>%
  left_join(new.genos, by = c("LabID", "Locus")) %>%
  mutate(
    `1` = coalesce(`New.Allele.1`, `1`),
    `2` = coalesce(`New.Allele.2`, `2`)
  ) %>%
  select(-`New.Allele.1`, -`New.Allele.2`)

dim(genos)
dim(gtseq_genos)

gtseq_genotypes<-genos

#Merge replicates
#first, need to drop 12509 as its wrongly labelled
genos <- gtseq_genos %>%
  filter(LabID != "12509")
dim(genos)

genos <- genos %>%
  mutate(base_id = sub("[bc]$", "", LabID)) %>%
  group_by(base_id, Locus) %>%
  summarise(
    `1` = first(na.omit(`1`)),
    `2` = first(na.omit(`2`)),
    .groups = 'drop') %>%
  rename(LabID = base_id)
dim(genos)

gtseq_genos<-genos

