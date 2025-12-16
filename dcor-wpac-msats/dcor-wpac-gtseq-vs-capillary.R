######Assessing how well microsatellites amplify and genotype from GTseq sequence data
####Compare genotypes from GTseq data to those generated from capillary electrophoresis
####Two methods to generate genotypes ("msat", "gtseq")

####Load required libraries
library(scales) 
library(ggplot2)
library(psych)  
library(readr)  
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl) 
library(stats)
library(graphics)
library(grDevices)
library(utils)
library(datasets)
library(methods)
library(base)
library(strataG)

####Import Genotype files from capillary generated data
#There are some differences between the original msat genotype files compared to the gtseq file formats 
#Different steps required to format them for comparisons later. 

##Capillary generated microsatellite genotypes
msat_genotypes_wide <- read_excel("Dc_pac_sample_info.xlsx", sheet = "GenotypeDc_pac_GTseq_022525.") #combined genotypes for Isla and Micho
head(msat_genotypes_wide) 

##Sort columns so that the loci are in alphabetical order
##Note: genotypes are presented in a single column per locus, with a space between the two allele sizes
##This step wasn't necessary in these datasets since the loci were already in alphabetical order. But may be usefull in future. 
#a<-sort(colnames(Cm_msat_genotypes_wide[,3:12])) #creates the alphabetical ordered list of loci
#Cm_msat_genotypes_wide[,3:12]<-Cm_msat_genotypes_wide[,a]   
#head(Cm_msat_genotypes_wide)

##NOT USING FOR DCOR AS ALREADY IN SEPARATE COLUMNS
##Separate alleles into different columns
##Note: first convert data frame into tibble to use tidyr/dplyr
#Cm_msat_genotypes_wide<-as_tibble(Cm_msat_genotypes_wide)
#head(Cm_msat_genotypes_wide)
##Separate the specified columns (3-12) based on the specified delimiter (" ") and rename column headings with an _1 or _2
#Cm_msat_genotypes_wide<-Cm_msat_genotypes_wide %>% separate_wider_delim(cols=c(3:12), delim=" ", names_sep="_")
#head(Cm_msat_genotypes_wide)

#Standardize dataframe 
#Add a column that describes method used  and move it to sit after the sample location
#Delete columns that won't be used
msat_genotypes_wide <- msat_genotypes_wide %>% mutate(Method="Capillary")
msat_genotypes_wide<- msat_genotypes_wide %>% relocate(Method, .after=Comment)
msat_genotypes_wide<- msat_genotypes_wide %>% select(-c(ID, Discard, Other_ID, Turtle_ID, Turtle_Archive_ID, Capture_Type, Location, Location_Code, Comment, EditDate, EditUserID, RecordCreationDate))
head(msat_genotypes_wide)

#Convert dataframe from wide to long format
#Each line will represent one sample, one locus, two alleles. 
msat_genotypes <- msat_genotypes_wide %>%
  pivot_longer(
    cols = -c(Lab_ID, Method),
    names_to = c("Locus", ".value"),
    names_sep = "_"
  ) 

#Extract only the loci that appear in both datasets:
loci_to_keep <- c("14-5", "C102", "D1", "LB123", "LB133", "LB141", "LB142", "LB143", "LB145", "LB158")

# Filter the dataframe
msat_genotypes <- msat_genotypes %>%
  filter(Locus %in% loci_to_keep)

#Change ID name to match
msat_genotypes <- msat_genotypes %>%
  rename(LabID=Lab_ID)

#Look at distribution of alleles in gtseq dataset
msat_genotypes_longer <- msat_genotypes %>%
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
ggplot(msat_genotypes_longer, aes(x = Allele_Value)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black", alpha = 0.8) +
  facet_wrap(~ Locus, scales = "free") +
  labs(
    title = "Distribution of Capillary Allele Values by Locus",
    subtitle = "Alleles 1 and 2 are combined for frequency count",
    x = "Allele Value",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "lightgray", color = "black"),
    strip.text = element_text(face = "bold")
  )

##Import genotypes GTseq generated data
#Importing data
gtseq_genotypes_wide <- read_table("megasat/Output_all/Genotype.txt")
head(gtseq_genotypes_wide)
#remove last column, which is empty
gtseq_genotypes_wide<-gtseq_genotypes_wide[,-68]

#Fix the sample ID column to match the format of the Capillary datasheet, should read "LabID"
gtseq_genotypes_wide<- gtseq_genotypes_wide %>% rename(LabID=Sample_idx1_idx2)
head(gtseq_genotypes_wide)

#Add column with method like the other dataset
#Note: need column names to be the same for when the datasets are merged later
gtseq_genotypes_wide<- gtseq_genotypes_wide %>% mutate(Method="GTseq") %>% relocate("Method", .after=LabID)
head(gtseq_genotypes_wide)
#checking individual genotypes
gtseq_genotypes_wide %>%
  select(c("LabID","GT_LB142", "GT_LB142-b")) %>%
  filter(LabID=="102605b")

#Check genotypes from replicates and merge
gtseq_genotypes_wide
gt<-df2gtypes(gtseq_genotypes_wide, id.col=1, strata.col=NULL, loc.col=3, ploidy=2 )
gt@data
gt.dupes<-dupGenotypes(gt)

colnames(msat_genotypes_wide)
colnames(gtseq_genotypes_wide)

##Copy locus names from the capillary dataset since they are in the right format
##Need to have each locus allele labelled with "_1" or "_2" with identical locus names first
#oldcol<-colnames(Cm_gtseq_genotypes_micho) #locus/allele names we want to replace
#newcol<-colnames(Cm_msat_genotypes_wide) #correctly formatted locus/allele names
#Cm_gtseq_genotypes_micho<-Cm_gtseq_genotypes_micho %>% rename_at(vars(oldcol), ~ newcol)
#head(Cm_gtseq_genotypes_micho)
#
##Convert each column to characters before converting dataframe from wide to long format
#Cm_gtseq_genotypes_micho<- Cm_gtseq_genotypes_micho %>% mutate_all(as.character)

#Convert dataset from wide to long format
#Each line will represent one sample, one locus, two alleles. 
#Cm_gtseq_genotypes_micho<- Cm_gtseq_genotypes_micho %>% pivot_longer(cols=-c(1:3), names_to=c("Locus", ".value"), names_sep="_")
#head(Cm_gtseq_genotypes_micho)

#Updated wide to long format
gtseq_genotypes<-gtseq_genotypes_wide
gtseq_genotypes <- gtseq_genotypes %>%
  pivot_longer(
    cols = -c(LabID, Method),
    names_to = "raw_locus_name",
    values_to = "value"
  ) %>% mutate(
    allele_num = if_else(str_ends(raw_locus_name, "-b"), 2, 1), # Create the allele number: if the name ends with "-b", it's 2, otherwise 1.
    Locus = str_remove(raw_locus_name, "GT_"), Locus = str_remove(Locus, "-b") # Create the clean Locus name: remove the "GT_" prefix and then the "-b" suffix
    ) %>%
  pivot_wider(
    id_cols = c(LabID, Method, Locus),
    names_from = allele_num,
    values_from = value
  )

#Extract only the loci that appear in both datasets:
loci_to_keep <- c("14-5", "C102", "D1", "LB123", "LB133", "LB141", "LB142", "LB143", "LB145", "LB158")
# Filter the dataframe
gtseq_genotypes <- gtseq_genotypes %>%
  filter(Locus %in% loci_to_keep)

#Cleanup dataset by replacing blanks, x and 0s with NA
#We want to be sure that missing data is not classified as a "0" later on for calculations
gtseq_genotypes<-replace(gtseq_genotypes, gtseq_genotypes==0, NA)
gtseq_genotypes<-replace(gtseq_genotypes, gtseq_genotypes=='X', NA)
gtseq_genotypes<-replace(gtseq_genotypes, gtseq_genotypes=='Unscored', NA)
head(gtseq_genotypes)

#Look at distribution of alleles in gtseq dataset
gtseq_genotypes_longer <- gtseq_genotypes %>%
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
ggplot(gtseq_genotypes_longer, aes(x = Allele_Value)) +
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

##Look at sequence lengths from fastq files since some of the megasat genotypes are too long
source("R/plot_all_sequence_lengths.R")
source("R/plot_max_sequence_lengths.R")

all_sequence_lengths<-plot_all_sequence_lengths(directory="data-raw/all/", output_file = "results/all_sequence_lengths.png")
max_sequence_lengths<-plot_max_sequence_lengths(directory="data-raw/all/", output_file = "results/max_sequence_lengths.png")

#####The max. sequence length across all files is 151 bp

#Convert genotypes >300bp to NA (value chosen based on histogram of allele lengths)
#Megasat will genotype any allele if the forward primer and forward flanking region are present
#If the sequence cuts off in the middle of the reverse flanking region, the program will 
#automatically add the rest of the reverse flanking region.
#But, if the microsatellite is longer than the sequencing chemistry allows, and it 
#cuts off before the reverse flanking region, it will give an arbitrarily large size (>300bp) fragment
#and it is impossible to know what the true genotype would be in this circumstance
#Need to turn any genotype with a 300bp allele to NA for both alleles
gtseq_genotypes <- gtseq_genotypes %>%
  mutate(across(c(`1`, `2`), as.numeric)) %>%
  mutate(
    `1` = ifelse(`1` > 300 | `2` > 300, NA, `1`),
    `2` = ifelse(`1` > 300 | `2` > 300, NA, `2`)
  )

#use code above to make new histograms to be sure that the desired alleles have been removed

##Check for mismatches in replicated individuals
source("./R/Compare.replicates.msats.R")
##This function reports mismatched genotypes including NAs
na_mismatches <- na_mismatches(
  df_long = gtseq_genotypes,
  ind_col = "LabID",
  locus_col = "Locus",
  allele_cols = c("1", "2")
)
##This function reports only mismatched genotypes that don't include NAs
no_na_mismatches<- no_na_mismatches(
  df_long = gtseq_genotypes,
  ind_col = "LabID",
  locus_col = "Locus",
  allele_cols = c("1", "2")
)

write.csv(na_mismatches, file = "./results/na_gtseq_mismatches_after-large-allele-filter.csv")
write.csv(no_na_mismatches, file = "./results/no_na_gtseq_mismatches_after-large-allele-filter.csv")

## Ensure alleles are numeric
gtseq_genotypes_numeric <- gtseq_genotypes %>%
  mutate(across(c(`1`, `2`), as.numeric))

## Find all replicate sets that have NO missing data
valid_comparisons <- gtseq_genotypes_numeric %>%
  mutate(base_id = gsub("[a-zA-Z]+$", "", LabID)) %>%
  group_by(base_id) %>%
  filter(n_distinct(LabID) > 1) %>% # Keep only replicate sets
  ungroup() %>%
  group_by(base_id, Locus) %>%
  summarise(
    has_na = any(is.na(`1`) | is.na(`2`)), # Check for NAs in the group
    .groups = 'drop'
  ) %>%
  filter(!has_na) # Keep only groups with NO NAs

## Count the total valid comparisons for each locus
total_comparisons_per_locus <- valid_comparisons %>%
  group_by(Locus) %>%
  summarise(total_comparisons = n())

## Count the number of mismatches per locus from your dataframe
total_mismatches_per_locus <- no_na_mismatches %>%
  group_by(Locus) %>%
  summarise(total_mismatches = n_distinct(base_id))

## Join mismatches and comparisons
gtseq_error_summary <- left_join(total_comparisons_per_locus, 
                                 total_mismatches_per_locus, 
                                 by = "Locus") %>%
## Replace NA mismatches with 0 (for loci with 0 errors)
mutate(total_mismatches = ifelse(is.na(total_mismatches), 0, total_mismatches)) %>%
## Calculate the final rate
mutate(error_rate = total_mismatches / total_comparisons)

## Print the summary table
cat("--- GTseq Per-Locus Error Rates ---\n")
print(gtseq_error_summary)

## Calculate and print the OVERALL error rate
overall_gtseq_error <- sum(gtseq_error_summary$total_mismatches) / sum(gtseq_error_summary$total_comparisons)
cat("\nOverall GTseq Error Rate:", scales::percent(overall_gtseq_error, accuracy = 0.1), "\n")



##Generate list of replicates
replicates <- gtseq_genotypes %>%
  mutate(base_id = str_remove(LabID, "[a-z]$")) %>%
  group_by(base_id) %>%
  summarise(replicate_count = n()) %>%
  filter(replicate_count > 1) %>%
  pull(base_id)

write.csv(replicates, "data-raw/replicate_individuals.csv")

##After updating genotypes for mismatches, update the correct genotypes and 
#amend the gtseq_genotypes file
new.genos<-read.csv("results/genos.to.change.csv")

#Merge new genotypes with original genotype file
genos <- gtseq_genotypes %>%
  left_join(new.genos, by = c("LabID", "Locus")) %>%
  mutate(
    `1` = coalesce(`New.allele.1`, `1`),
    `2` = coalesce(`New.allele.2`, `2`)
  ) %>%
  select(-`New.allele.1`, -`New.allele.2`)

dim(genos)
dim(gtseq_genotypes)

gtseq_genotypes<-genos

#Merge replicates
#first, need to drop 12509 as its wrongly labelled
genos <- gtseq_genotypes %>%
  filter(LabID != "12509")
dim(genos)

genos <- genos %>%
  mutate(base_id = sub("[bc]$", "", LabID)) %>%
  group_by(base_id, Locus) %>%
  summarise(
    Method = first(Method),
    `1` = first(na.omit(`1`)),
    `2` = first(na.omit(`2`)),
    .groups = 'drop') %>%
  rename(LabID = base_id)
dim(genos)

gtseq_genotypes<-genos %>%
  select(LabID, Method, Locus, `1`, `2`)

####Compare capillary and GTseq generated genotypes
#Final longform datasets
msat_genotypes #genotypes gnerated from capillary electrophoresis
gtseq_genotypes #genotypes generated from GTseq

#Merge the datasets for each genotyping method
merged_genotypes<-merge(msat_genotypes, gtseq_genotypes, by.x=c("LabID", "Locus"), by.y=c("LabID", "Locus"))

#Drop individuals/loci that weren't genotyped across both methods
merged_genotypes <- merged_genotypes %>%
  filter( !if_all(c(`1.x`, `2.x`, `1.y`, `2.y`), is.na)) %>%
  drop_na(`1.x`, `2.x`, `1.y`, `2.y`)

#Genotype columns to must be treated as numeric for the calculations
#Note, the column names have changed since unique values are requred for column names
merged_genotypes$`1.x`<-as.numeric(merged_genotypes$`1.x`) #Capillary generated allele 1
merged_genotypes$`2.x`<-as.numeric(merged_genotypes$`2.x`) #Capillary generated allele 2
merged_genotypes$`1.y`<-as.numeric(merged_genotypes$`1.y`) #GTseq generated allele 1
merged_genotypes$`2.y`<-as.numeric(merged_genotypes$`2.y`) #GTseq generated allele 2


#Zygosity comparison
#This calculation tells whether the genotypes from the two methods share the same zygosity
#This also tells whether the magnitude of difference between heterozygous alleles is the same between methods
#A zero indicates that the two datasets share the same zygosity and the same magnitude of difference
#Positive values indicate that the magnitude of difference between allele 1 and 2 is greater in GTseq dataset
#This could mean one allele dropped out entirely in the capillary method, or it was estimated to be a different size 
#Negative values indicate that the magnitude of difference between allele 1 and 2 is smaller in GTseq dataset
#This could mean one allele dropped out entirely in GTseq method, or it was estimated to be a different size
#Equation: (GTseqAllele2-MsatAllele2)-(GTseqAllele1-MsatAllele1)
merged_genotypes$zygtest<-(merged_genotypes$`2.y`-merged_genotypes$`2.x`)-(merged_genotypes$`1.y`-merged_genotypes$`1.x`)
#summary statistics of zygosity value per locus
describeBy(merged_genotypes$zygtest, merged_genotypes$Locus, mat=TRUE) 
#histogram of the zygosity values per locus
ggplot(merged_genotypes, 
       aes(x= zygtest, fill = Locus))+ geom_histogram (binwidth = 2, alpha = 0.5, position = "identity")+
  facet_wrap(~Locus, ncol=4) +
  ggtitle("Zygosity test")

#Offset calculation
#This calculation tells you how many bases the are different between the genotypes generated using capillary and GTseq.
#If consistent, GTseq and capillary genotypes could be combined with an offset value pplied
#Equation: GTseq Allele 1 - Capillary Allele 1
merged_genotypes$offset<-merged_genotypes$`1.y`-merged_genotypes$`1.x`
head(merged_genotypes)
#Summary statistics for offset by locus
describeBy(merged_genotypes$offset, merged_genotypes$Locus, mat=TRUE)
#Histograms of offset value by locus
ggplot(merged_genotypes, 
       aes(x= offset, fill = Locus))+ geom_histogram (binwidth = 1, alpha = 0.5, position = "identity")+
  facet_wrap(~Locus, ncol=4, scales="free") +
  scale_x_continuous(breaks=breaks_pretty())+
  ggtitle("Offset")

#####Plotting zygtest and offest results against eachother.
## The size of the dot represents the sample size
p<-ggplot(merged_genotypes, aes(x= zygtest, y=offset)) + geom_count()  +
  facet_wrap(~Locus, ncol=4, scales="free") + scale_x_continuous(label = scales::label_comma(accuracy = 1)) + scale_y_continuous(label = scales::label_comma(accuracy = 1)) +
  ggtitle("Zygtest vs Offset") +
  theme_bw()
q<-layer_data(p) #this pulls out the count data for each combination of zygtest and offset combinations
#pull out the associated zygtest and offset values for the point size with the biggest n
q<-q %>% group_by(PANEL) %>% filter(n==max(n))
q<-q[,c(1:3)] #only keep the panel number, x (zygtest) and y (offset) values
colnames(q)<-c("PANEL", "zygtest", "offset")

#Modify this dataframe with Locus names that correspond with panel number
#make a data frame with the corresponding 
locus_panels<-data.frame(panel=c(1:8), locus=c("14-5", "D1", "LB123", "LB133", "LB142", "LB143", "LB145", "LB158"))
colnames(locus_panels)<-c("PANEL", "Locus")

#combine two dataframes to change panel id to locus name
r<-merge(q, locus_panels, by="PANEL", all.x = TRUE, all.y=TRUE)
r<-r[order(r$Locus),]

#join up new dataframe with old one to make a new figure
#old=merged_genotypes
#new=q
s<-merge(merged_genotypes, r, by="Locus", all.x=TRUE, all.y=TRUE )
s$facet_label<-paste(s$Locus, "(", s$zygtest.y, ",", s$offset.y, ")")

##Make a new plot with the zygtest and offset values printed in the facets
zygtest.plot<-ggplot(s, aes(x= zygtest.x, y=offset.x)) + geom_count()  +
  facet_wrap(~facet_label, ncol=4, scales="free") + scale_x_continuous(label = scales::label_comma(accuracy = 1)) + scale_y_continuous(label = scales::label_comma(accuracy = 1)) +
  ggtitle("Zygtest vs Offset") +
  theme_bw()
dim(merged_genotypes)
length(unique(merged_genotypes$LabID)) #135
ggsave("./results/zygtest-v-offset.png", plot = zygtest.plot)


#Estimating concordance
#Join the systematic offset from 'r' back into the main 'merged_genotypes' df
#'r' has the 'Locus' and the 'offset.y' (the main offset)
comparison_df <- left_join(merged_genotypes, r, by = "Locus")

# 2. Normalize alleles and calculate match
comparison_df <- comparison_df %>%
  rowwise() %>%
  mutate(
    # Sort capillary alleles
    cap_a1 = min(`1.x`, `2.x`),
    cap_a2 = max(`1.x`, `2.x`),
    
    # Sort GTseq alleles
    gt_a1 = min(`1.y`, `2.y`),
    gt_a2 = max(`1.y`, `2.y`),
    
    # Apply systematic offset (from offset.y) to capillary alleles
    cap_a1_adj = cap_a1 + offset.y,
    cap_a2_adj = cap_a2 + offset.y,
    
    # 4. Check for an exact match
    is_match = (cap_a1_adj == gt_a1) & (cap_a2_adj == gt_a2)
  ) %>%
  ungroup()

# 5. Summarize the results
# Overall Concordance
total_comparisons <- nrow(comparison_df)
total_matches <- sum(comparison_df$is_match)
overall_concordance_rate <- total_matches / total_comparisons
cat("Overall Concordance:", overall_concordance_rate, "\n")

# Per-Locus Concordance
locus_concordance <- comparison_df %>%
  group_by(Locus) %>%
  summarise(
    concordance = mean(is_match, na.rm = TRUE),
    n = n()
  )

# Create the bar plot
geno_concordance_plot<-ggplot(locus_concordance, aes(x = reorder(Locus, concordance), y = concordance)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_text(aes(label = paste0(round(concordance * 100), "%\n(n=", n, ")")), 
            vjust = -0.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) + # Increased limit
  labs(
    title = "Genotype Concordance (GTseq vs. Capillary)",
    subtitle = "After applying systematic locus offset",
    x = "Locus",
    y = "Concordance Rate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./results/Genotype_concordance.png", geno_concordance_plot)

# Create status columns (using the normalized alleles from above)
comparison_df <- comparison_df %>%
  mutate(
    cap_status = ifelse(cap_a1 == cap_a2, "Homozygous", "Heterozygous"),
    gt_status = ifelse(gt_a1 == gt_a2, "Homozygous", "Heterozygous")
  )

# Create the contingency table
mismatch_summary <- table(
  Capillary = comparison_df$cap_status,
  GTseq = comparison_df$gt_status
)
print(mismatch_summary)

library(ggplot2)
library(dplyr)

#Create zygosity status columns
zygosity_summary <- comparison_df %>%
  filter(!is.na(is_match)) %>% # Only use complete comparisons
  mutate(
    cap_status = ifelse(cap_a1 == cap_a2, "Homozygous", "Heterozygous"),
    gt_status = ifelse(gt_a1 == gt_a2, "Homozygous", "Heterozygous"),
    comparison = case_when(
      cap_status == "Homozygous" & gt_status == "Homozygous"   ~ "Match (Hom)",
      cap_status == "Heterozygous" & gt_status == "Heterozygous" ~ "Match (Het)",
      cap_status == "Heterozygous" & gt_status == "Homozygous"   ~ "Mismatch (Cap-Het / GT-Hom)",
      cap_status == "Homozygous" & gt_status == "Heterozygous"   ~ "Mismatch (Cap-Hom / GT-Het)"
    )
  )

#Create a stacked bar plot
zygosity_plot<-ggplot(zygosity_summary, aes(x = Locus, fill = comparison)) +
  geom_bar(position = "fill") +
  geom_text(data = locus_concordance, 
            aes(x = Locus, y = 1.05, label = paste("n=", n)), 
            inherit.aes = FALSE, 
            size = 3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.15)) + # Increased limit
  labs(
    title = "Zygosity Agreement by Locus",
    x = "Locus",
    y = "Proportion of Samples",
    fill = "Comparison Status"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("./results/zygosity_barplot.png", zygosity_plot)

#Pivot the data into a "long" format for easier plotting
alleles_long <- comparison_df %>%
  select(LabID, Locus, cap_a1_adj, gt_a1, cap_a2_adj, gt_a2) %>%
  pivot_longer(
    cols = c(cap_a1_adj, cap_a2_adj),
    names_to = "cap_allele_num",
    values_to = "cap_adj_size"
  ) %>%
  mutate(
    gt_size = ifelse(cap_allele_num == "cap_a1_adj", gt_a1, gt_a2)
  )

#Create the scatter plot
allele_size_scatter<-ggplot(alleles_long, aes(x = cap_adj_size, y = gt_size)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  
  # ADDED: This layer adds the 'n' value to the top-left corner of each facet
  geom_text(data = locus_concordance, 
            aes(x = -Inf, y = Inf, label = paste("n=", n)), 
            inherit.aes = FALSE, 
            hjust = -0.2, 
            vjust = 1.5) +
  
  facet_wrap(~ Locus, scales = "free") +
  labs(
    title = "GTseq vs. Capillary Allele Size Correlation",
    subtitle = "Capillary alleles adjusted for systematic offset. Red line = perfect match.",
    x = "Adjusted Capillary Allele Size (bp)",
    y = "GTseq Allele Size (bp)"
  ) +
  theme_bw()
ggsave("./results/allele_size_scatter.png", allele_size_scatter)


