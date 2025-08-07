# Identify samples with Ns or Xs in their haplotypes or more than 2 haplotypes
questionable.hap <- sapply(1:nrow(tgt), function(i){
  ifelse(length(grep("N", tgt$gt[i],)) > 0 || length(grep("X", tgt$gt[i],)) > 0 
         || tgt$num.haps[i] > 2, TRUE, FALSE) 
})
genos.to.check <- filter(tgt, questionable.hap == TRUE)
table(genos.to.check$locus)
#locus016 locus019 locus069 
#3       44        1

tgt <- filter(tgt, questionable.hap == FALSE)
if(nrow(genos.to.check > 0)) {
  print("Some samples with Ns or Xs in their haplotypes or more than 2 haplotypes")
  print(paste0("Questionable genotypes saved to results-R/", project, ".genos.to.check.rda"))
  save(genos.to.check, file = paste0("results-R/", project, ".genos.to.check.rda"))
}

# summarize individual data
missing.data.ind <- data.frame(table(tgt$Indiv[!is.na(tgt$gt)])) %>%
  mutate(missing = num.locs-Freq)
names(missing.data.ind) <- c("labID", "genos", "missing")
length(which(missing.data.ind$genos >= min.genos.per.ind))
inds.2.keep <- filter(missing.data.ind, genos >= min.genos.per.ind) |> 
  pull(labID)
tgt <- filter(tgt, Indiv %in% inds.2.keep) #filters out individuals that have fewer genotypes than "min.genos.per.ind"
num.inds <- length(unique(tgt$Indiv))
write.csv(missing.data.ind, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.num.genos.per.ind.csv"))
dim(missing.data.ind)

# summarize locus data
tgt_long <- tgt |> 
  select(locus, Indiv, gt, depth.1, depth.2) |> 
  separate_wider_delim(
    cols = gt, 
    delim = '/', 
    names = c('haplo.1', 'haplo.2'), 
    cols_remove = FALSE
  ) |> 
  mutate(depth.2 = ifelse(haplo.2 == haplo.1, NA, depth.2))

loc.sum <- tgt_long %>%
  pivot_longer(cols = c(haplo.1, haplo.2), names_to = 'hap') |> 
  mutate(tmp = strsplit(as.character(value), "")) %>%
  unnest(tmp) %>%
  group_by(locus, Indiv, hap) %>%
  mutate(name = 1:n()) %>%
  pivot_wider(id_cols = c(locus, Indiv, gt, hap), values_from = tmp, names_prefix = 'snp') |> 
  ungroup() |> 
  filter(!is.na(gt)) |> 
  group_by(locus) |> 
  summarise(
    inds.genoed = n() / 2,
    num.unique.genos = length(unique(gt)),
    num.alleles.pos1 = length(unique(snp1)),
    num.alleles.pos2 = length(unique(snp2)),
    num.alleles.pos3 = length(unique(snp3))
  )
write.csv(loc.sum, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.locus.summary.csv"))