
## NEW FEATURE: Workflow Summary Functions

### workflow_summary_functions.R

```r
# workflow_summary_functions.R
# Functions to create and maintain a progressive workflow summary report

#' Initialize Workflow Summary
#' Creates the initial summary data structure
#' @param project_name Character string of project name
#' @param analysis_stage Character string of current stage
#' @param results_path Path where summary will be saved
initialize_workflow_summary <- function(project_name, analysis_stage, results_path) {
  
  summary_data <- list(
    project_info = list(
      project_name = project_name,
      analysis_stage = analysis_stage,
      start_time = Sys.time(),
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      user = Sys.info()["user"]
    ),
    
    workflow_steps = data.frame(
      step = character(0),
      timestamp = as.POSIXct(character(0)),
      n_individuals = integer(0),
      n_loci = integer(0),
      individuals_dropped = integer(0),
      loci_dropped = integer(0),
      reason = character(0),
      notes = character(0),
      stringsAsFactors = FALSE
    ),
    
    filtering_history = list(),
    replicate_issues = data.frame(
      individual1 = character(0),
      individual2 = character(0),
      similarity = numeric(0),
      action_taken = character(0),
      stringsAsFactors = FALSE
    ),
    
    problematic_loci = data.frame(
      locus = character(0),
      issue_type = character(0),
      description = character(0),
      step_identified = character(0),
      action_taken = character(0),
      stringsAsFactors = FALSE
    ),
    
    problematic_individuals = data.frame(
      individual = character(0),
      issue_type = character(0),
      value = numeric(0),
      threshold = numeric(0),
      step_identified = character(0),
      action_taken = character(0),
      stringsAsFactors = FALSE
    ),
    
    parameters_used = list()
  )
  
  # Save initial summary
  summary_file <- file.path(results_path, paste0(project_name, "_", analysis_stage, "_workflow_summary.rds"))
  saveRDS(summary_data, summary_file)
  
  message("✅ Workflow summary initialized: ", summary_file)
  return(summary_data)
}

#' Add Step to Workflow Summary
#' Records each major step in the workflow
add_workflow_step <- function(summary_data, step_name, n_individuals, n_loci, 
                             individuals_dropped = 0, loci_dropped = 0, 
                             reason = "", notes = "", results_path, 
                             project_name, analysis_stage) {
  
  new_step <- data.frame(
    step = step_name,
    timestamp = Sys.time(),
    n_individuals = n_individuals,
    n_loci = n_loci,
    individuals_dropped = individuals_dropped,
    loci_dropped = loci_dropped,
    reason = reason,
    notes = notes,
    stringsAsFactors = FALSE
  )
  
  summary_data$workflow_steps <- rbind(summary_data$workflow_steps, new_step)
  
  # Save updated summary
  summary_file <- file.path(results_path, paste0(project_name, "_", analysis_stage, "_workflow_summary.rds"))
  saveRDS(summary_data, summary_file)
  
  message("📝 Step recorded: ", step_name, " (", n_individuals, " individuals, ", n_loci, " loci)")
  return(summary_data)
}

#' Add Filtering Parameters
#' Records the filtering thresholds used
add_filtering_parameters <- function(summary_data, step_name, parameters, 
                                   results_path, project_name, analysis_stage) {
  
  summary_data$parameters_used[[step_name]] <- c(
    timestamp = as.character(Sys.time()),
    parameters
  )
  
  # Save updated summary
  summary_file <- file.path(results_path, paste0(project_name, "_", analysis_stage, "_workflow_summary.rds"))
  saveRDS(summary_data, summary_file)
  
  return(summary_data)
}

#' Add Problematic Loci
#' Records loci that had issues
add_problematic_loci <- function(summary_data, loci_df, results_path, 
                               project_name, analysis_stage) {
  
  summary_data$problematic_loci <- rbind(summary_data$problematic_loci, loci_df)
  
  # Save updated summary
  summary_file <- file.path(results_path, paste0(project_name, "_", analysis_stage, "_workflow_summary.rds"))
  saveRDS(summary_data, summary_file)
  
  return(summary_data)
}

#' Add Problematic Individuals
#' Records individuals that had issues
add_problematic_individuals <- function(summary_data, individuals_df, results_path, 
                                      project_name, analysis_stage) {
  
  summary_data$problematic_individuals <- rbind(summary_data$problematic_individuals, individuals_df)
  
  # Save updated summary
  summary_file <- file.path(results_path, paste0(project_name, "_", analysis_stage, "_workflow_summary.rds"))
  saveRDS(summary_data, summary_file)
  
  return(summary_data)
}

#' Add Replicate Information
#' Records replicate comparison results
add_replicate_info <- function(summary_data, replicate_df, results_path, 
                             project_name, analysis_stage) {
  
  summary_data$replicate_issues <- rbind(summary_data$replicate_issues, replicate_df)
  
  # Save updated summary
  summary_file <- file.path(results_path, paste0(project_name, "_", analysis_stage, "_workflow_summary.rds"))
  saveRDS(summary_data, summary_file)
  
  return(summary_data)
}

#' Generate Summary Report
#' Creates a formatted markdown/HTML report
generate_summary_report <- function(summary_data, results_path, project_name, analysis_stage) {
  
  # Create report filename
  report_file <- file.path(results_path, paste0(project_name, "_", analysis_stage, "_workflow_summary_report.html"))
  
  # Create the report content
  report_content <- paste0(
    "# GT-seq Workflow Summary Report\n\n",
    "**Project:** ", summary_data$project_info$project_name, "\n\n",
    "**Analysis Stage:** ", summary_data$project_info$analysis_stage, "\n\n",
    "**Generated:** ", Sys.time(), "\n\n",
    "**User:** ", summary_data$project_info$user, "\n\n",
    "**R Version:** ", summary_data$project_info$r_version, "\n\n",
    "---\n\n"
  )
  
  # Workflow Steps Table
  if (nrow(summary_data$workflow_steps) > 0) {
    report_content <- paste0(report_content, "## Workflow Steps\n\n")
    
    # Create a nice summary table
    steps_summary <- summary_data$workflow_steps %>%
      mutate(
        timestamp = format(timestamp, "%Y-%m-%d %H:%M:%S"),
        change_individuals = ifelse(individuals_dropped == 0, "", 
                                  paste0("-", individuals_dropped)),
        change_loci = ifelse(loci_dropped == 0, "", 
                           paste0("-", loci_dropped))
      ) %>%
      select(step, timestamp, n_individuals, change_individuals, 
             n_loci, change_loci, reason)
    
    report_content <- paste0(report_content, 
                           knitr::kable(steps_summary, format = "html"), 
                           "\n\n")
  }
  
  # Parameters Used
  if (length(summary_data$parameters_used) > 0) {
    report_content <- paste0(report_content, "## Filtering Parameters Used\n\n")
    for (step_name in names(summary_data$parameters_used)) {
      params <- summary_data$parameters_used[[step_name]]
      report_content <- paste0(report_content, "### ", step_name, "\n")
      report_content <- paste0(report_content, "**Timestamp:** ", params["timestamp"], "\n\n")
      for (param_name in names(params)[names(params) != "timestamp"]) {
        report_content <- paste0(report_content, "- **", param_name, ":** ", params[param_name], "\n")
      }
      report_content <- paste0(report_content, "\n")
    }
  }
  
  # Problematic Loci
  if (nrow(summary_data$problematic_loci) > 0) {
    report_content <- paste0(report_content, "## Problematic Loci\n\n")
    report_content <- paste0(report_content, 
                           knitr::kable(summary_data$problematic_loci, format = "html"), 
                           "\n\n")
  }
  
  # Problematic Individuals
  if (nrow(summary_data$problematic_individuals) > 0) {
    report_content <- paste0(report_content, "## Problematic Individuals\n\n")
    report_content <- paste0(report_content, 
                           knitr::kable(summary_data$problematic_individuals, format = "html"), 
                           "\n\n")
  }
  
  # Replicate Issues
  if (nrow(summary_data$replicate_issues) > 0) {
    report_content <- paste0(report_content, "## Replicate Comparison Results\n\n")
    report_content <- paste0(report_content, 
                           knitr::kable(summary_data$replicate_issues, format = "html"), 
                           "\n\n")
  }
  
  # Final Summary
  if (nrow(summary_data$workflow_steps) > 0) {
    final_step <- summary_data$workflow_steps[nrow(summary_data$workflow_steps), ]
    initial_step <- summary_data$workflow_steps[1, ]
    
    report_content <- paste0(report_content, "## Final Summary\n\n",
                           "- **Initial dataset:** ", initial_step$n_individuals, " individuals, ", initial_step$n_loci, " loci\n",
                           "- **Final dataset:** ", final_step$n_individuals, " individuals, ", final_step$n_loci, " loci\n",
                           "- **Total individuals removed:** ", initial_step$n_individuals - final_step$n_individuals, 
                           " (", round((initial_step$n_individuals - final_step$n_individuals) / initial_step$n_individuals * 100, 1), "%)\n",
                           "- **Total loci removed:** ", initial_step$n_loci - final_step$n_loci,
                           " (", round((initial_step$n_loci - final_step$n_loci) / initial_step$n_loci * 100, 1), "%)\n\n")
  }
  
  # Convert markdown to HTML and save
  html_content <- markdown::markdownToHTML(text = report_content, 
                                         options = c("use_xhtml", "smartypants", "base64_images", "mathjax", "highlight_code"))
  
  writeLines(html_content, report_file)
  
  message("📊 Summary report generated: ", report_file)
  
  # Also save a simple CSV version of the steps
  csv_file <- file.path(results_path, paste0(project_name, "_", analysis_stage, "_workflow_steps.csv"))
  write.csv(summary_data$workflow_steps, csv_file, row.names = FALSE)
  
  return(report_file)
}
```

### Integration Example for Script 01

```r
# Example integration into 01_genotyping_and_qc_refactored_v2.R
# Add this after loading custom functions:

source("R/functions/workflow_summary_functions.R")

# Add this after the MASTER CONFIGURATION section:

# Initialize workflow summary
workflow_summary <- initialize_workflow_summary(
  project_name = project_name,
  analysis_stage = analysis_stage, 
  results_path = results_raw_path
)

# PART 1 modifications:
# After creating initial tgt object:

# Record initial data load
workflow_summary <- add_workflow_step(
  summary_data = workflow_summary,
  step_name = "01_initial_data_load",
  n_individuals = n_distinct(tgt$Indiv),
  n_loci = n_distinct(tgt$locus),
  reason = "Initial genotype calling from SAM files",
  notes = paste0("Used microhaplot with min_read_depth=", min_read_depth, 
                ", ab_min_het=", ab_min_het, ", ab_max_homo=", ab_max_homo),
  results_path = results_raw_path,
  project_name = project_name,
  analysis_stage = analysis_stage
)

# PART 2 modifications:
# Before applying first filter, record parameters:

workflow_summary <- add_filtering_parameters(
  summary_data = workflow_summary,
  step_name = "02_loose_filtering_parameters",
  parameters = list(
    locus_threshold = pass1.locus.threshold,
    individual_threshold = pass1.indiv.threshold,
    manual_loci_removed = paste(pass1.manual.loci.to.remove, collapse = ", "),
    manual_individuals_removed = paste(pass1.manual.inds.to.remove, collapse = ", ")
  ),
  results_path = results_raw_path,
  project_name = project_name,
  analysis_stage = analysis_stage
)

# After applying first filter:
workflow_summary <- add_workflow_step(
  summary_data = workflow_summary,
  step_name = "02_loose_filtering_applied",
  n_individuals = n_distinct(tgt$Indiv),
  n_loci = n_distinct(tgt$locus),
  individuals_dropped = n_inds_before - n_distinct(tgt$Indiv),
  loci_dropped = n_loci_before - n_distinct(tgt$locus),
  reason = "First pass (loose) filtering",
  notes = paste0("Removed loci with <", pass1.locus.threshold*100, "% genotyping and individuals with <", pass1.indiv.threshold*100, "% genotyping"),
  results_path = results_raw_path,
  project_name = project_name,
  analysis_stage = analysis_stage
)

# Add problematic loci identified in loose filtering
problematic_loci_loose <- loc_sum_initial %>%
  filter(prop.genoed < pass1.locus.threshold | locus %in% pass1.manual.loci.to.remove) %>%
  mutate(
    issue_type = case_when(
      prop.genoed < pass1.locus.threshold ~ "Low genotyping rate",
      locus %in% pass1.manual.loci.to.remove ~ "Manual removal",
      TRUE ~ "Other"
    ),
    description = case_when(
      prop.genoed < pass1.locus.threshold ~ paste0("Genotyped in only ", round(prop.genoed*100, 1), "% of individuals"),
      locus %in% pass1.manual.loci.to.remove ~ "Manually flagged for removal",
      TRUE ~ ""
    ),
    step_identified = "02_loose_filtering",
    action_taken = "Removed"
  ) %>%
  select(locus, issue_type, description, step_identified, action_taken)

if (nrow(problematic_loci_loose) > 0) {
  workflow_summary <- add_problematic_loci(
    summary_data = workflow_summary,
    loci_df = problematic_loci_loose,
    results_path = results_raw_path,
    project_name = project_name,
    analysis_stage = analysis_stage
  )
}

# PART 3 modifications (similar pattern):
# Record final filtering parameters and results

# At the very end, generate the final report:
message("\n📊 Generating workflow summary report...")
report_file <- generate_summary_report(
  summary_data = workflow_summary,
  results_path = results_raw_path,
  project_name = project_name,
  analysis_stage = analysis_stage
)

message("✅ Workflow summary report saved: ", report_file)
```

## Additional Integration Points:
### For Script 03 (Duplicate Check):
```r
# After finding duplicates:
if (!is.null(duplicate.pairs)) {
  replicate_info <- duplicate.pairs %>%
    mutate(
      action_taken = "User review required",
      similarity = round(similarity, 3)
    ) %>%
    rename(individual1 = id.1, individual2 = id.2)
    
  workflow_summary <- add_replicate_info(
    summary_data = workflow_summary,
    replicate_df = replicate_info,
    results_path = results_raw_path,
    project_name = project_name,
    analysis_stage = analysis_stage
  )
}
```

### For Script 04 (Population Genetics):
```r
# Record high homozygosity individuals:
if (length(high.homo.samps) > 0) {
  high_homo_df <- ind.summary %>%
    filter(id %in% high.homo.samps) %>%
    mutate(
      issue_type = "High homozygosity",
      value = pct.loci.homozygous,
      threshold = high.homo.threshold,
      step_identified = "04_homozygosity_check",
      action_taken = "Flagged for review"
    ) %>%
    select(individual = id, issue_type, value, threshold, step_identified, action_taken)
    
  workflow_summary <- add_problematic_individuals(...)
}
```
---

## Workflow Assessment Summary

### Strengths
- **Excellent Structure**: Clear 4-step progression with logical dependencies
- **Consistent Naming**: Automated path generation reduces user error
- **Great Documentation**: Comprehensive README with step-by-step instructions
- **Robust Design**: Interactive approach with mandatory examination points
- **Two-pass Filtering**: Methodologically sound loose → strict approach

### Areas for Enhancement
1. **Prerequisites**: Add R version requirements and installation verification
2. **Example Data**: Small test dataset to verify installation
3. **Parameter Guidance**: Recommended threshold ranges with rationale
4. **Troubleshooting**: Common error messages and solutions
5. **Output Interpretation**: Brief explanations of result files

### New Summary Report Features
- **Progressive Tracking**: Summary builds up as workflow runs
- **Automatic Persistence**: Saves after each step
- **Rich Information**: Captures parameters, problematic samples/loci, replicate issues
- **Multiple Formats**: RDS, CSV, and HTML outputs
- **Publication Ready**: Formatted tables perfect for methods sections

### Overall Rating: Very Good (8/10)

The workflow is well-designed and user-friendly. The addition of the summary reporting functionality would make it excellent for publication and audit trail purposes.

---

