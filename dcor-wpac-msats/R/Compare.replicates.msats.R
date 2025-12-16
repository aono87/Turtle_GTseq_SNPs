# install.packages("tidyverse")
library(tidyverse)

#' Find Mismatches Caused by NA Values
#' This function identifies replicate sets where the genotypes differ AND
#' at least one of the genotypes in the comparison contains an NA value.
#'
#' @param df_long A dataframe in long format.
#' @param ind_col The name of the column with individual IDs.
#' @param locus_col The name of the column with locus names.
#' @param allele_cols A character vector of the allele column names.
#'
#' @return A dataframe showing only the mismatches that involve NA values.
#'
na_mismatches <- function(df_long, ind_col, locus_col, allele_cols) {
  
  df_long %>%
    # Step 1: Directly check for NAs in any allele column.
    mutate(
      has_na = if_any(all_of(allele_cols), is.na)
    ) %>%
    # Step 2: Create the sorted genotype string.
    mutate(
      sorted_genotype = pmap_chr(
        select(., all_of(allele_cols)),
        ~ paste(sort(c(...)), collapse = "/")
      )
    ) %>%
    # Step 3: Create the base ID for grouping.
    mutate(base_id = gsub("[a-zA-Z]+$", "", .data[[ind_col]])) %>%
    # Step 4: Group by replicate set and locus.
    group_by(base_id, .data[[locus_col]]) %>%
    # Step 5: Apply the new filtering logic.
    filter(
      # Condition 1: Must be a replicate group.
      n() > 1 &
        # Condition 2 (NA Check): The group MUST contain at least one NA.
        any(has_na) &
        # Condition 3: The genotypes must be different.
        n_distinct(sorted_genotype) > 1
    ) %>%
    # Step 6: Clean up and return the results.
    ungroup() %>%
    arrange(base_id, .data[[locus_col]], .data[[ind_col]]) %>%
    select(base_id, all_of(ind_col), all_of(locus_col), all_of(allele_cols))
}

no_na_mismatches <- function(df_long, ind_col, locus_col, allele_cols) {
  df_long %>%
    # Step 1: Directly check for NAs in any allele column for each row.
    mutate(
      has_na = if_any(all_of(allele_cols), is.na)
    ) %>%
    # Step 2: Create the sorted genotype string.
    mutate(
      sorted_genotype = pmap_chr(
        select(., all_of(allele_cols)),
        ~ paste(sort(c(...)), collapse = "/")
      )
    ) %>%
    # Step 3: Create the base ID for grouping.
    mutate(base_id = gsub("[a-zA-Z]+$", "", .data[[ind_col]])) %>%
    # Step 4: Group by replicate set and locus.
    group_by(base_id, .data[[locus_col]]) %>%
    # Step 5: Apply the robust filtering logic.
    filter(
      # Condition 1: Must be a replicate group (more than one sample).
      n() > 1 &
        # Condition 2 (Robust NA Check): The entire group must have zero NAs.
        !any(has_na) &
        # Condition 3: Genotypes must be different (a true mismatch).
        n_distinct(sorted_genotype) > 1
    ) %>%
    # Step 6: Clean up and return the results.
    ungroup() %>%
    arrange(base_id, .data[[locus_col]], .data[[ind_col]]) %>%
    select(base_id, all_of(ind_col), all_of(locus_col), all_of(allele_cols))
}