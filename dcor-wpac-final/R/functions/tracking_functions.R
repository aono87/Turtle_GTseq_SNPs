# R/functions/tracking_functions.R
# This script contains consolidated functions for tracking the genotyping workflow.

#' Update summary, individual, and locus trackers after a filtering step.
#'
#' @return A list containing the three updated tracker data frames.
#'
update_all_trackers <- function(progress_tracker, individual_tracker, locus_tracker,
                                tgt_before, tgt_after,
                                ind_summary_before, loc_summary_before,
                                stage_prefix, indiv_threshold, locus_threshold,
                                manual_inds_remove, manual_loci_remove,
                                previous_ind_status_col = NULL,
                                previous_locus_status_col = NULL) {
  
  # --- 1. Update the DETAILED `individual_tracker` ---
  
  # a. Determine the status ("Kept" or "Removed") based on the PRE-filter data
  individual_tracker <- individual_tracker %>%
    dplyr::left_join(ind_summary_before %>% dplyr::select(Indiv, loci.genoed, prop.genoed), by = "Indiv") %>%
    dplyr::rename(
      !!paste0(stage_prefix, "_loci_genoed_pre") := loci.genoed,
      !!paste0(stage_prefix, "_prop_genoed_pre") := prop.genoed
    )
  
  ind_status_col <- paste0(stage_prefix, "_status")
  ind_prop_col_pre <- paste0(stage_prefix, "_prop_genoed_pre")
  
  if (is.null(previous_ind_status_col)) {
    individual_tracker <- individual_tracker %>%
      dplyr::mutate(!!ind_status_col := dplyr::case_when(
        Indiv %in% manual_inds_remove ~ "Removed (Manual)",
        !!rlang::sym(ind_prop_col_pre) < indiv_threshold ~ "Removed (Low Call Rate)",
        is.na(!!rlang::sym(ind_prop_col_pre)) ~ "Removed (No Data)",
        TRUE ~ "Kept"
      ))
  } else {
    individual_tracker <- individual_tracker %>%
      dplyr::mutate(!!ind_status_col := dplyr::case_when(
        stringr::str_starts(!!rlang::sym(previous_ind_status_col), "Removed") ~ !!rlang::sym(previous_ind_status_col),
        Indiv %in% manual_inds_remove ~ "Removed (Manual)",
        !!rlang::sym(ind_prop_col_pre) < indiv_threshold ~ "Removed (Low Call Rate)",
        is.na(!!rlang::sym(ind_prop_col_pre)) ~ "Removed (No Data)",
        TRUE ~ "Kept"
      ))
  }
  
  # b. Calculate and add the POST-filter stats for the samples that were kept
  num_locs_post <- dplyr::n_distinct(tgt_after$locus)
  ind_sum_post_filter <- tgt_after %>%
    dplyr::filter(!is.na(gt)) %>%
    dplyr::group_by(Indiv) %>%
    dplyr::summarise(
      post_loci_count = dplyr::n_distinct(locus),
      post_prop_genoed = dplyr::n_distinct(locus) / num_locs_post,
      .groups = 'drop'
    ) %>%
    dplyr::rename_with(~ paste0(stage_prefix, "_", .), .cols = dplyr::starts_with("post_"))
  
  individual_tracker <- individual_tracker %>%
    dplyr::left_join(ind_sum_post_filter, by = "Indiv") %>%
    dplyr::select(-dplyr::starts_with(paste0(stage_prefix, "_loci_genoed_pre")),
                  -dplyr::starts_with(paste0(stage_prefix, "_prop_genoed_pre")))
  
  # --- 2. Update the DETAILED `locus_tracker` ---
  
  # a. Determine status ("Kept" or "Removed") based on PRE-filter data
  locus_tracker <- locus_tracker %>%
    dplyr::left_join(loc_summary_before %>% dplyr::select(locus, inds.genoed, prop.genoed), by = "locus") %>%
    dplyr::rename(
      !!paste0(stage_prefix, "_inds_genoed_pre") := inds.genoed,
      !!paste0(stage_prefix, "_prop_genoed_pre") := prop.genoed
    )
  
  locus_status_col <- paste0(stage_prefix, "_status")
  locus_prop_col_pre <- paste0(stage_prefix, "_prop_genoed_pre")
  
  if (is.null(previous_locus_status_col)) {
    locus_tracker <- locus_tracker %>%
      dplyr::mutate(!!locus_status_col := dplyr::case_when(
        locus %in% manual_loci_remove ~ "Removed (Manual)",
        !!rlang::sym(locus_prop_col_pre) < locus_threshold ~ "Removed (Low Call Rate)",
        is.na(!!rlang::sym(locus_prop_col_pre)) ~ "Removed (No Data)",
        TRUE ~ "Kept"
      ))
  } else {
    locus_tracker <- locus_tracker %>%
      dplyr::mutate(!!locus_status_col := dplyr::case_when(
        stringr::str_starts(!!rlang::sym(previous_locus_status_col), "Removed") ~ !!rlang::sym(previous_locus_status_col),
        locus %in% manual_loci_remove ~ "Removed (Manual)",
        !!rlang::sym(locus_prop_col_pre) < locus_threshold ~ "Removed (Low Call Rate)",
        is.na(!!rlang::sym(locus_prop_col_pre)) ~ "Removed (No Data)",
        TRUE ~ "Kept"
      ))
  }
  
  # b. Calculate and add POST-filter stats for kept loci
  num_inds_post <- dplyr::n_distinct(tgt_after$Indiv)
  loc_sum_post_filter <- tgt_after %>%
    dplyr::filter(!is.na(gt)) %>%
    dplyr::group_by(locus) %>%
    dplyr::summarise(
      post_inds_count = dplyr::n_distinct(Indiv),
      post_prop_genoed = dplyr::n_distinct(Indiv) / num_inds_post,
      .groups = 'drop'
    ) %>%
    dplyr::rename_with(~ paste0(stage_prefix, "_", .), .cols = dplyr::starts_with("post_"))
  
  locus_tracker <- locus_tracker %>%
    dplyr::left_join(loc_sum_post_filter, by = "locus") %>%
    dplyr::select(-dplyr::starts_with(paste0(stage_prefix, "_inds_genoed_pre")),
                  -dplyr::starts_with(paste0(stage_prefix, "_prop_genoed_pre")))
  
  
  # --- 3. Update the SUMMARY `progress_tracker` ---
  
  notes <- paste0("Locus threshold: ", locus_threshold, "; Indiv threshold: ", indiv_threshold)
  loci_removed_str <- if (length(manual_loci_remove) > 0 && manual_loci_remove[1] != "") { paste(manual_loci_remove, collapse = ", ") } else { NA_character_ }
  inds_removed_str <- if (length(manual_inds_remove) > 0 && manual_inds_remove[1] != "") { paste(manual_inds_remove, collapse = ", ") } else { NA_character_ }
  
  new_summary_row <- tibble::tibble(
    Stage = paste0("0", nrow(progress_tracker) + 1, "_", stage_prefix, "_Filter"),
    Loci_Count = dplyr::n_distinct(tgt_after$locus),
    Loci_Removed = dplyr::n_distinct(tgt_before$locus) - dplyr::n_distinct(tgt_after$locus),
    Manual_Loci_Removed = loci_removed_str,
    Individuals_Count = dplyr::n_distinct(tgt_after$Indiv),
    Individuals_Removed = dplyr::n_distinct(tgt_before$Indiv) - dplyr::n_distinct(tgt_after$Indiv),
    Manual_Individuals_Removed = inds_removed_str,
    Notes = notes
  )
  
  progress_tracker <- dplyr::bind_rows(progress_tracker, new_summary_row)
  
  # --- 4. Return the updated trackers ---
  
  return(list(
    progress_tracker = progress_tracker,
    individual_tracker = individual_tracker,
    locus_tracker = locus_tracker
  ))
}
