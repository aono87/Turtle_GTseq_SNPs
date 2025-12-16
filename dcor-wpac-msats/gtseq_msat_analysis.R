library(adegenet)
library(dplyr)
library(strataG)
library(PopGenReport)
library(pegas)
library(poppr)
library(ggplot2)
library(ggrepel)

#Final msat genotypes
genos<-gtseq_genos #long format

#Final gtseq snp gtypes file
load(file = "~/Documents/GitHub/Turtle_GTseq_SNPs/dcor-wpac-new/data/gtypes_dcor_wpac_new_final_nodups_minReads20.rda")
g

#Convert from long format to wide format
gtseq_genos
msats_wide<-gtseq_genos %>%
  pivot_wider(
    names_from = Locus,
    values_from = c(`1`,`2`)) %>%
  rename_with(
    ~ sub("^([0-9]+)_(.*)$", "\\2.\\1", .),
    matches("^[0-9]+_")
  )
str(msats_wide)

#sort columns
locus_cols <- setdiff(colnames(msats_wide), "LabID")
sorted_locus_cols <- locus_cols[order(
  str_remove(locus_cols, "\\.\\d$"), # Sorts by base name (e.g., "A1")
  str_extract(locus_cols, "\\d$")     # Then sorts by allele number (e.g., "1")
)]
msats_wide<- msats_wide %>%
  select(LabID, all_of(sorted_locus_cols))

#Final msat dataframe
msats_wide

#Get sample final sample ID and stratum
strata<-g@data %>%
  select(c(id, stratum)) %>%
  distinct() %>%
  rename(LabID = id)
str(strata)

#add zpadding to sample names
msats_wide <- msats_wide %>%
  mutate(LabID = sprintf("z%07d", as.numeric(LabID)))
str(msats_wide)

#Merge dataframes
msat_genos<-left_join(
  strata, msats_wide, by="LabID")
str(msat_genos)
dim(msat_genos) 

#check for duplicates using stratag
msats.g<-df2gtypes(msat_genos, ploidy=2, strata.col=2, loc.col=3)
msat.dups<-dupGenotypes(msats.g, num.shared=0.9)
#No duplicates found

#allele frequencies using stratag
msat.af<-alleleFreqs(msats.g, by.strata = TRUE)
msat.af

#Summaries
msat.loc.sum<-summarizeLoci(msats.g)
msat.loc.sum
hist(msat.loc.sum$prop.genotyped, breaks = 15)

msat.ind.sum<-summarizeInds(msats.g)
msat.ind.sum
hist(msat.ind.sum$pct.loci.missing.genotypes, breaks = 15)

#Remove loci with <50% genotyped
locs.to.keep<- msat.loc.sum %>%
  filter(prop.genotyped > 0.8) %>%
  pull(locus) 

msats_filtered<-msat_genos %>%
  select(c(LabID, stratum, matches(locs.to.keep)))

msats.g.filtered<-df2gtypes(msats_filtered, ploidy=2, strata.col=2, loc.col=3)

#Summaries
hist(summarizeInds(msats.g.filtered)$pct.loci.missing.genotypes, breaks=20)
hist(summarizeLoci(msats.g.filtered)$prop.genotyped, breaks=20)

#Remove inds with >30% missing data
inds.to.keep<-summarizeInds(msats.g.filtered) %>%
  filter(pct.loci.missing.genotypes <0.3) %>%
  pull(id)

msats_filtered<-msats_filtered %>%
  filter(LabID %in% inds.to.keep)

msats.g.filtered<-df2gtypes(msats_filtered, ploidy=2, strata.col=2, loc.col=3)

#Summaries
hist(summarizeInds(msats.g.filtered)$pct.loci.missing.genotypes)
hist(summarizeLoci(msats.g.filtered)$prop.genotyped)

msat.loc.sum<-summarizeLoci(msats.g.filtered)
msat.loc.sum

msat.ind.sum<-summarizeInds(msats.g.filtered)
msat.ind.sum

##Worked through this tutorial: https://bookdown.org/hhwagner1/LandGenCourse_book/Week5.html
##Check if al loci are polymorphic
msat.gi<-gtypes2genind(msats.g.filtered)
summary(msat.gi)
#NO, two loci only have one allele
msat.gi@loc.n.all
#Need to drop A1 and B116 as they are monomorphic across amples
nLoc(msat.gi)
allLocs<-locNames(msat.gi)
loci_to_remove <- c("A1", "B116")
loci_to_keep <- allLocs[!allLocs %in% loci_to_remove]
msat.gi.poly <- msat.gi[loc = loci_to_keep]
summary(msat.gi.poly)

## Check for null alleles
msat.null.all<-null.all(msat.gi.poly)
msat.null.all$homozygotes
msat.null.all$null.allele.freq

#the 95% CI spans 0 for all loci in both methods. The frequency of null alleles is not different than zero. 
msat.null.all$null.allele.freq$summary1
msat.null.all$null.allele.freq$summary2

#Check for HWE-all inds
round(pegas::hw.test(msat.gi.poly, B = 1000), digits = 3)
#No loci out of HWE with exact test result

#Check for HWE-pops
HWE.test <- data.frame(sapply(seppop(msat.gi.poly), 
                              function(ls) pegas::hw.test(ls, B=0)[,3]))
HWE.test.chisq <- t(data.matrix(HWE.test))
{cat("Chi-squared test (p-values):", "\n")
  round(HWE.test.chisq,3)}

HWE.test <- data.frame(sapply(seppop(msat.gi.poly), 
                              function(ls) pegas::hw.test(ls, B=1000)[,4]))
HWE.test.MC <- t(data.matrix(HWE.test))
{cat("MC permuation test (p-values):", "\n")
  round(HWE.test.MC,3)}

alpha=0.05/90
Prop.loci.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 2, mean), 
                                   MC=apply(HWE.test.MC<alpha, 2, mean))
Prop.loci.out.of.HWE 
#non out out HWE

Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
                                   MC=apply(HWE.test.MC<alpha, 1, mean))
Prop.pops.out.of.HWE   

Chisq.fdr <- matrix(p.adjust(HWE.test.chisq,method="fdr"), 
                    nrow=nrow(HWE.test.chisq))
MC.fdr <- matrix(p.adjust(HWE.test.MC, method="fdr"), 
                 nrow=nrow(HWE.test.MC))

Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
                                   MC=apply(HWE.test.MC<alpha, 1, mean),
                                   Chisq.fdr=apply(Chisq.fdr<alpha, 1, mean),
                                   MC.fdr=apply(MC.fdr<alpha, 1, mean))
Prop.pops.out.of.HWE


#Checking for linkage disequilibrium
poppr::ia(msat.gi.poly, sample=199)

LD.pair <- poppr::pair.ia(msat.gi.poly)
LD.pair

#Genetic diversity
Sum <- adegenet::summary(msat.gi.poly)
names(Sum)
par(mar=c(10, 5.5,1,1))
barplot(Sum$pop.n.all, las=3, 
        xlab = "", ylab = "Number of alleles")

plot(Sum$n.by.pop, Sum$pop.n.all, 
     xlab = "Sample size", ylab = "Number of alleles")
abline(lm(Sum$pop.n.all ~ Sum$n.by.pop), col = "red")

Richness <- PopGenReport::allel.rich(msat.gi.poly, min.alleles = NULL)
Richness$alleles.sampled #2

par(mar=c(10, 4.5,1,1))
barplot(Richness$mean.richness, las=3, ylab="Rarefied allelic richness (Ar)")
plot(colMeans(Richness$pop.sizes), Richness$mean.richness,
     xlab="Valid sample size", 
     ylab="Rarefied allelic richness (Ar)")
abline(lm(Richness$mean.richness ~ colMeans(Richness$pop.sizes)), col="red")

#Heterozygosity
par(mar=c(3, 4.5,1,1))
barplot(Sum$Hexp, ylim=c(0,1), ylab="Expected heterozygosity")
barplot(Sum$Hobs, ylim=c(0,1), ylab="Observed heterozygosity")

Hobs <- t(sapply(seppop(msat.gi.poly), function(ls) summary(ls)$Hobs))
Hexp <- t(sapply(seppop(msat.gi.poly), function(ls) summary(ls)$Hexp))
{cat("Expected heterozygosity (Hexp):", "\n")
  round(Hexp, 2)}
{cat("\n", "Observed heterozygosity (Hobs):", "\n")
  round(Hobs, 2)}

par(mar=c(5.5, 4.5, 1, 1))
Hobs.pop <- apply(Hobs, MARGIN = 1, FUN = mean)
Hexp.pop <- apply(Hexp, 1, mean) 
barplot(Hexp.pop, ylim=c(0,1), las=3, ylab="Expected heterozygosity")
barplot(Hobs.pop, ylim=c(0,1), las=3, ylab="Observed heterozygosity")

msat.diversity <- data.frame(Pop = names(Hobs.pop),
                              n = Sum$n.by.pop,
                              Hobs = Hobs.pop,
                              Hexp = Hexp.pop,
                              Ar = Richness$mean.richness)
msat.diversity

#Population level
msat.genpop <- adegenet::genind2genpop(msat.gi.poly)
Freq <- adegenet::makefreq(msat.genpop)
round(Freq, 2)
apply(Freq, MARGIN = 1, FUN = sum)


#Combining SNP and microsatellite data
msat.gi.poly
g

g.snps<-g
g.msat<-genind2gtypes(msat.gi.poly)

#gi.snps<-gtypes2genind(g.snps)
#gi.msat<-msat.gi.poly
#
#snp.df<-as.data.frame(gi.snps@tab)
#msat.df<-as.data.frame(gi.msat@tab)
#
#snp.df <- snp.df %>%
#  rownames_to_column(var = "Indiv")
#msat.df<-msat.df%>%
#  rownames_to_column(var="Indiv")
#gi.msat@pop
#gi.snps@pop
#
#snp.df<-snp.df %>%
#  mutate(stratum=gi.snps@pop) %>%
#  relocate(stratum, .after=Indiv)
#
#msat.df<-msat.df %>%
#  mutate(stratum=gi.msat@pop) %>%
#  relocate(stratum, .after=Indiv)
#
#all.df<-full_join(snp.df, msat.df, by=c("Indiv", "stratum"))

snp.df<-as.data.frame(g.snps@data)
msat.df<-as.data.frame(g.msat@data)

str(snp.df)
str(msat.df)

snp.df.wide <- snp.df %>%
  group_by(id, stratum, locus) %>%
  summarise(genotype = paste(sort(allele), collapse = "/"), .groups = 'drop') %>%
  pivot_wider(
    names_from = locus,
    values_from = genotype,
    values_fill = NA_character_
  )

msat.df.wide <- msat.df %>%
  group_by(id, stratum, locus) %>%
  summarise(genotype = paste(sort(allele), collapse = "/"), .groups = 'drop') %>%
  pivot_wider(
    names_from = locus,
    values_from = genotype,
    values_fill = NA_character_
  )

all.df.wide <- full_join(snp.df.wide, msat.df.wide, by=c("id", "stratum"))
all.gi<-df2genind(
  all.df.wide[, -c(1, 2)],    # Genotype data (all columns except id and stratum)
  sep = "/",                  # The character separating alleles
  ind.names = all.df.wide$id,     # A vector of individual names
  pop = all.df.wide$stratum,      # A vector of population assignments
  ploidy = 2                  # Assuming diploid data
)

#DAPC
mat_ref <- tab(all.gi, NA.method = 'mean')
grp_ref <- pop(all.gi)

# Run the cross-validation. This can be time-consuming.
# Cross-validation parameters.
xval_n_pca_max <- 200
xval_training_set <- 0.9
xval_reps <- 30

# DAPC parameters for prediction (if using the custom predictAllIDsDAPC function).
loov_n_da <- 8
loov_n_pca <- 100
xval_results <- xvalDapc(mat_ref, grp_ref, n.pca.max = xval_n_pca_max,
                         training.set = xval_training_set, result = 'groupMean',
                         n.rep = xval_reps, xval.plot = TRUE)

dapc_results <- xval_results$DAPC
message(paste("   - Optimal number of PCs identified:", dapc_results$n.pca))
message(paste("   - Number of discriminant axes retained:", dapc_results$n.da))


# ====================================================================
# STEP 3: VISUALIZE AND SUMMARIZE DAPC RESULTS
# ====================================================================
message("\nStep 3: Generating plots and summarizing DAPC results...")

# Define a color palette for plotting.
num_groups <- nlevels(pop(all.gi))
plot_colors <- hcl.colors(num_groups, palette = "Zissou 1")

# --- Visualise DAPC Scatter Plot ---
scatter.dapc(dapc_results, col = plot_colors, scree.da = FALSE, legend = TRUE, posi.leg = "topright")

message("\nStep 3: Generating advanced plot with non-overlapping labels...")

# --- 1. Prepare Data for ggplot ---
# Extract individual coordinates and add population info
ind_coords_df <- as.data.frame(dapc_results$ind.coord)
ind_coords_df$population <- dapc_results$grp

population_counts <- ind_coords_df %>% 
  count(population)
print("Number of individuals per population:")
print(population_counts)

# Extract centroid coordinates
cent_coords_df <- as.data.frame(dapc_results$grp.coord)
cent_coords_df$population <- rownames(cent_coords_df)

# Calculate percentage of variance explained by the axes
percent_explained <- (dapc_results$eig[1:2] / sum(dapc_results$eig)) * 100
percent_explained

# --- 2. Create the ggplot object ---
dapc_plot <- ggplot() +
  # Add individual points
  geom_point(
    data = ind_coords_df,
    aes(x = LD1, y = LD2, color = population),
    size = 2, alpha = 0.6
  ) +
  # Add ellipses to group points by population
  stat_ellipse(
    data = ind_coords_df,
    aes(x = LD1, y = LD2, color = population, fill = population),
    geom = "polygon", alpha = 0.1, linetype = "dashed"
  ) +
  # Add the non-overlapping population labels at the centroids
  geom_label_repel(
    data = cent_coords_df,
    aes(x = LD1, y = LD2, label = population, fill = population),
    color = "white", fontface = "bold", point.padding = 1, box.padding = 0.5
  ) +
  # CHANGED: Apply the custom manual color scheme
  scale_color_manual(values = plot_colors, guide = "none") +
  scale_fill_manual(values = plot_colors, guide = "none") +
  # Add informative labels and a clean theme
  labs(
    title = "Discriminant Analysis of Principal Components (DAPC)",
    x = paste0("LD1 (", round(percent_explained[1], 1), "%)"),
    y = paste0("LD2 (", round(percent_explained[2], 1), "%)")
  ) +
  theme_classic()
dapc_plot

#DAPC without in-water samples
all.gi
pop.gi<-all.gi

pops_to_remove <- c("In-water, CA", "In-water, Indonesia")
all_pops <- levels(pop(all.gi))
pops_to_keep <- all_pops[!all_pops %in% pops_to_remove]
pop.gi <- pop.gi[pop = pops_to_keep]
levels(pop(pop.gi))
pop.gi

#re-run dapc
mat_ref_pop <- tab(pop.gi, NA.method = 'mean')
grp_ref_pop <- pop(pop.gi)
xval_results_pop <- xvalDapc(mat_ref_pop, grp_ref_pop, n.pca.max = xval_n_pca_max,
                         training.set = xval_training_set, result = 'groupMean',
                         n.rep = xval_reps, xval.plot = TRUE)

dapc_results_pop <- xval_results_pop$DAPC
message(paste("   - Optimal number of PCs identified:", dapc_results$n.pca))
message(paste("   - Number of discriminant axes retained:", dapc_results$n.da))

num_groups_pop <- nlevels(pop(pop.gi))
plot_colors_pop <- hcl.colors(num_groups_pop, palette = "Zissou 1")
scatter.dapc(dapc_results_pop, col = plot_colors_pop, scree.da = FALSE, legend = TRUE, posi.leg = "topright")

#Try without Malaysia too
all.gi
pop.gi<-all.gi

pops_to_remove <- c("In-water, CA", "In-water, Indonesia", "Malaysia", "PNG", "Bird's Head-Winter")
all_pops <- levels(pop(all.gi))
pops_to_keep <- all_pops[!all_pops %in% pops_to_remove]
pop.gi <- pop.gi[pop = pops_to_keep]
levels(pop(pop.gi))
pop.gi

#re-run dapc
mat_ref_pop <- tab(pop.gi, NA.method = 'mean')
grp_ref_pop <- pop(pop.gi)
xval_results_pop <- xvalDapc(mat_ref_pop, grp_ref_pop, n.pca.max = xval_n_pca_max,
                             training.set = xval_training_set, result = 'groupMean',
                             n.rep = xval_reps, xval.plot = TRUE)

dapc_results_pop <- xval_results_pop$DAPC
message(paste("   - Optimal number of PCs identified:", dapc_results$n.pca))
message(paste("   - Number of discriminant axes retained:", dapc_results$n.da))

num_groups_pop <- nlevels(pop(pop.gi))
plot_colors_pop <- hcl.colors(num_groups_pop, palette = "Zissou 1")
scatter.dapc(dapc_results_pop, col = plot_colors_pop, scree.da = FALSE, legend = TRUE, posi.leg = "topright")

####Re-doing DAPC for SNPs only
gi.snps

#DAPC
mat_ref <- tab(gi.snps, NA.method = 'mean')
grp_ref <- pop(gi.snps)

# Run the cross-validation. This can be time-consuming.
# Cross-validation parameters.
xval_n_pca_max <- 200
xval_training_set <- 0.9
xval_reps <- 30

# DAPC parameters for prediction (if using the custom predictAllIDsDAPC function).
loov_n_da <- 8
loov_n_pca <- 100
xval_results <- xvalDapc(mat_ref, grp_ref, n.pca.max = xval_n_pca_max,
                         training.set = xval_training_set, result = 'groupMean',
                         n.rep = xval_reps, xval.plot = TRUE)

dapc_results <- xval_results$DAPC
message(paste("   - Optimal number of PCs identified:", dapc_results$n.pca))
message(paste("   - Number of discriminant axes retained:", dapc_results$n.da))


# ====================================================================
# STEP 3: VISUALIZE AND SUMMARIZE DAPC RESULTS
# ====================================================================
message("\nStep 3: Generating plots and summarizing DAPC results...")

# Define a color palette for plotting.
num_groups <- nlevels(pop(gi.snps))
plot_colors <- hcl.colors(num_groups, palette = "Zissou 1")

# --- Visualise DAPC Scatter Plot ---
scatter.dapc(dapc_results, col = plot_colors, scree.da = FALSE, legend = TRUE, posi.leg = "topright")


#DAPC without in-water samples
gi.snps

pops_to_remove <- c("In-water, CA", "In-water, Indonesia", "Malaysia", "PNG", "Bird's Head-Winter")
all_pops <- levels(pop(gi.snps))
pops_to_keep <- all_pops[!all_pops %in% pops_to_remove]
pop.gi <- gi.snps[pop = pops_to_keep]
levels(pop(pop.gi))
pop.gi

#re-run dapc
mat_ref_pop <- tab(pop.gi, NA.method = 'mean')
grp_ref_pop <- pop(pop.gi)
xval_results_pop <- xvalDapc(mat_ref_pop, grp_ref_pop, n.pca.max = xval_n_pca_max,
                             training.set = xval_training_set, result = 'groupMean',
                             n.rep = xval_reps, xval.plot = TRUE)

dapc_results_pop <- xval_results_pop$DAPC
message(paste("   - Optimal number of PCs identified:", dapc_results$n.pca))
message(paste("   - Number of discriminant axes retained:", dapc_results$n.da))

num_groups_pop <- nlevels(pop(pop.gi))
plot_colors_pop <- hcl.colors(num_groups_pop, palette = "Zissou 1")
scatter.dapc(dapc_results_pop, col = plot_colors_pop, scree.da = FALSE, legend = TRUE, posi.leg = "topright")


pops_to_remove <- c("In-water, CA", "In-water, Indonesia")
all_pops <- levels(pop(all.gi))
pops_to_keep <- all_pops[!all_pops %in% pops_to_remove]
pop.gi.7 <- all.gi[pop = pops_to_keep]
msat.gi.7<-gi.msat[pop = pops_to_keep]
popNames(pop.gi.7)
popNames(msat.gi.7)

pops_to_remove <- c("In-water, CA", "In-water, Indonesia", "Malaysia")
all_pops <- levels(pop(all.gi))
pops_to_keep <- all_pops[!all_pops %in% pops_to_remove]
pop.gi.6 <- all.gi[pop = pops_to_keep]
msat.gi.6 <- gi.msat[pop = pops_to_keep]
popNames(pop.gi.6)
popNames(msat.gi.6)

pops_to_remove <- c("In-water, CA", "In-water, Indonesia", "Malaysia", "PNG")
pops_to_remove2 <- c("In-water, CA", "In-water, Indonesia", "Malaysia", "Bird's Head-Winter")
all_pops <- levels(pop(all.gi))
pops_to_keep <- all_pops[!all_pops %in% pops_to_remove]
pops_to_keep2<- all_pops[!all_pops %in% pops_to_remove2]
pop.gi.5 <- all.gi[pop = pops_to_keep]
msat.gi.5 <- gi.msat[pop = pops_to_keep2]
popNames(pop.gi.5)
popNames(msat.gi.5)

pops_to_remove <- c("In-water, CA", "In-water, Indonesia", "Malaysia", "PNG", "Bird's Head-Winter")
all_pops <- levels(pop(all.gi))
pops_to_keep <- all_pops[!all_pops %in% pops_to_remove]
pop.gi.4 <- all.gi[pop = pops_to_keep]
msat.gi.4 <- gi.msat[pop = pops_to_keep]
popNames(pop.gi.4)
popNames(msat.gi.4)

mat_ref_msat <- tab(gi.msat, NA.method = 'mean')
grp_ref_msat <- pop(gi.msat)
xval_results_msat <- xvalDapc(mat_ref_msat, grp_ref_msat, n.pca.max = xval_n_pca_max,
                             training.set = xval_training_set, result = 'groupMean',
                             n.rep = xval_reps, xval.plot = TRUE)
scatter.dapc(xval_results_msat$DAPC, legend = TRUE, posi.leg = "topright")
table.value(table(xval_results_msat$DAPC$assign, grp_ref_msat), 
            col.lab = levels(grp_ref_msat))

mat_ref_msat <- tab(msat.gi.7, NA.method = 'mean')
grp_ref_msat <- pop(msat.gi.7)
xval_results_msat <- xvalDapc(mat_ref_msat, grp_ref_msat, n.pca.max = xval_n_pca_max,
                              training.set = xval_training_set, result = 'groupMean',
                              n.rep = xval_reps, xval.plot = TRUE)
scatter.dapc(xval_results_msat$DAPC, legend = TRUE, posi.leg = "topright")
table.value(table(xval_results_msat$DAPC$assign, grp_ref_msat), 
            col.lab = levels(grp_ref_msat))

mat_ref_msat <- tab(msat.gi.6, NA.method = 'mean')
grp_ref_msat <- pop(msat.gi.6)
xval_results_msat <- xvalDapc(mat_ref_msat, grp_ref_msat, n.pca.max = xval_n_pca_max,
                              training.set = xval_training_set, result = 'groupMean',
                              n.rep = xval_reps, xval.plot = TRUE)
scatter.dapc(xval_results_msat$DAPC, legend = TRUE, posi.leg = "topright")
table.value(table(xval_results_msat$DAPC$assign, grp_ref_msat), 
            col.lab = levels(grp_ref_msat))

mat_ref_msat <- tab(msat.gi.5, NA.method = 'mean')
grp_ref_msat <- pop(msat.gi.5)
xval_results_msat <- xvalDapc(mat_ref_msat, grp_ref_msat, n.pca.max = xval_n_pca_max,
                              training.set = xval_training_set, result = 'groupMean',
                              n.rep = xval_reps, xval.plot = TRUE)
scatter.dapc(xval_results_msat$DAPC, legend = TRUE, posi.leg = "topright")
table.value(table(xval_results_msat$DAPC$assign, grp_ref_msat), 
            col.lab = levels(grp_ref_msat))

###Let's try Rubias
all.gi #all samples, all loci
#read in rubias metadata file
rub.meta<-read_excel("./results/gtseq-all-locs.qa.qc.xlsx", sheet = "rubias")

genotype_cols <- all.df.wide %>%
  select(-id, -stratum) %>%
  names()

all.df.rubias <- all.df.wide %>%
  # Gather all genotype columns into a 'long' format.
  pivot_longer(
    cols = all_of(genotype_cols),
    names_to = "locus",
    values_to = "genotype"
  ) %>%
  # Split the 'genotype' column into two new allele columns.
  separate_wider_delim(
    genotype,
    delim = "/",
    names = c("allele1", "allele2"),
    too_few = "align_end"
  ) %>%
  # Pivot back to wide format, combining locus and allele names.
  pivot_wider(
    names_from = locus,
    values_from = c(allele1, allele2),
    names_sep = "." # Separates allele and locus with a dot (e.g., allele1.Dc00544)
  ) %>%
  # --- Step 3: Refine column names for allele1 and allele2 ---
  # Use a single rename_with() to handle both allele prefixes.
  rename_with(~ case_when(
    # For allele1 columns, remove the prefix entirely
    startsWith(.x, "allele1.") ~ sub("^allele1\\.", "", .x),
    # For allele2 columns, capture the locus name and add the .1 suffix
    startsWith(.x, "allele2.") ~ sub("^allele2\\.(.*)", "\\1.1", .x),
    # For all other columns (id, stratum), keep as is
    .default = .x
  )) %>%
  # Final cleaning: Move the 'id' and 'stratum' columns to the front.
  relocate(id, stratum, everything()) %>%
  relocate(id, stratum) %>%
  select(id, stratum, sort(names(.)[!(names(.) %in% c("id", "stratum"))])) %>%
  rename(LabID=id)

#merge metadata and genotypes
all.df.rubias
rub.meta
rubias.dat<-full_join(rub.meta, all.df.rubias, by="LabID") %>%
  rename(indiv=LabID) %>%
  mutate(across(7:ncol(.), ~na_if(., "")))
str(rubias.dat)

#Self assign reference populations
self_assignment<-self_assign(reference=filter(rubias.dat, sample_type=='reference'), gen_start_col = 7)
assignment_summary <- self_assignment %>%
  group_by(indiv) %>%
  top_n(1, scaled_likelihood) %>%
  ungroup() %>%
  mutate(across(where(is.list), as.character))

assignment_matrix_collection <- table(
  Observed = assignment_summary$collection,
  Inferred = assignment_summary$inferred_collection
)

assignment_matrix_repunit <- table(
  Observed = assignment_summary$repunit,
  Inferred = assignment_summary$inferred_repunit
)
print(assignment_matrix)
print(assignment_matrix_repunit)

###Rubias with just microsatellites
#read in rubias metadata file
rub.meta<-read_excel("./results/gtseq-all-locs.qa.qc.xlsx", sheet = "rubias")

genotype_cols <- msat.df.wide %>%
  select(-id, -stratum) %>%
  names()

msat.df.rubias <- msat.df.wide %>%
  # Gather all genotype columns into a 'long' format.
  pivot_longer(
    cols = all_of(genotype_cols),
    names_to = "locus",
    values_to = "genotype"
  ) %>%
  # Split the 'genotype' column into two new allele columns.
  separate_wider_delim(
    genotype,
    delim = "/",
    names = c("allele1", "allele2"),
    too_few = "align_end"
  ) %>%
  # Pivot back to wide format, combining locus and allele names.
  pivot_wider(
    names_from = locus,
    values_from = c(allele1, allele2),
    names_sep = "." # Separates allele and locus with a dot (e.g., allele1.Dc00544)
  ) %>%
  # --- Step 3: Refine column names for allele1 and allele2 ---
  # Use a single rename_with() to handle both allele prefixes.
  rename_with(~ case_when(
    # For allele1 columns, remove the prefix entirely
    startsWith(.x, "allele1.") ~ sub("^allele1\\.", "", .x),
    # For allele2 columns, capture the locus name and add the .1 suffix
    startsWith(.x, "allele2.") ~ sub("^allele2\\.(.*)", "\\1.1", .x),
    # For all other columns (id, stratum), keep as is
    .default = .x
  )) %>%
  # Final cleaning: Move the 'id' and 'stratum' columns to the front.
  relocate(id, stratum, everything()) %>%
  relocate(id, stratum) %>%
  select(id, stratum, sort(names(.)[!(names(.) %in% c("id", "stratum"))])) %>%
  rename(LabID=id)

#merge metadata and genotypes
msat.df.rubias
rub.meta
msat.rubias.dat<-full_join(rub.meta, msat.df.rubias, by="LabID") %>%
  rename(indiv=LabID) %>%
  mutate(across(7:ncol(.), ~na_if(., ""))) %>%
str(rubias.dat)

#Self assign reference populations
msat_self_assignment<-self_assign(reference=filter(msat.rubias.dat, sample_type=='reference'), gen_start_col = 7)
msat_assignment_summary <- msat_self_assignment %>%
  group_by(indiv) %>%
  top_n(1, scaled_likelihood) %>%
  ungroup() %>%
  mutate(across(where(is.list), as.character))

msat_assignment_matrix_collection <- table(
  Observed = msat_assignment_summary$collection,
  Inferred = msat_assignment_summary$inferred_collection
)

msat_assignment_matrix_repunit <- table(
  Observed = msat_assignment_summary$repunit,
  Inferred = msat_assignment_summary$inferred_repunit
)
print(msat_assignment_matrix_collection)
print(msat_assignment_matrix_repunit)
