# 06_STRUCTURE_ANALYSIS.R
#
# Workflow Step 6:
# This script performs a population structure analysis using the program STRUCTURE,
# run via the `strataG` package. It estimates the most likely number of genetic
# clusters (K) in a subset of the data, summarizes results using Evanno's method,
# and visualizes the cluster assignments for each value of K.
#
# INPUTS:
#   - A `gtypes` object (`.rda` file) from `02_create_gtypes.R`.
#
# OUTPUTS:
#   - An RDA file containing the raw results from the `structureRun`.
#   - A CSV file and a plot of the Evanno metrics (Delta K).
#   - A PDF file containing all the STRUCTURE bar plots for each K value.

# --- LOAD LIBRARIES ---
library(strataG)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(pophelper)

# ----------------------------------------------------------------------
# --- CONFIGURATION ---
# ----------------------------------------------------------------------

# --- Project and File Naming ---
# These settings MUST match the previous scripts.
project_name <- "dcor_wpac"
analysis_stage <- "final_nodups"
min_reads <- 20
run_label <- paste(project_name, analysis_stage, sep = "_")

# --- Auto-Generated Paths (Do not change) ---
#run_label <- paste(project_name, analysis_stage, sep = "_")
gtypes_path <- file.path("data", paste0("gtypes_", run_label, "_minReads", min_reads, ".rda"))

# --- Output Paths ---
results_raw_path <- file.path("results-raw", analysis_stage)
results_r_path <- file.path("results-R", analysis_stage)
dir.create(results_raw_path, showWarnings = FALSE, recursive = TRUE)
dir.create(results_r_path, showWarnings = FALSE, recursive = TRUE)

# --- Analysis Parameters ---
# Stratification scheme and specific strata to include in the analysis.
stratification_scheme <- "Stratum_ABO"
#strata_to_include <- c('CentAm-CA.OR.WA', 'HI-SEAK.NBC') #if you only want to run some strata

# STRUCTURE parameters
k_range <- 1:10       # The range of K (number of clusters) to test.
num_k_reps <- 3      # Number of replicate runs for each value of K.
burnin <- 100000      # Burn-in period for the MCMC.
numreps <- 500000     # Number of MCMC repetitions after burn-in.
num_cores <- 12       # Number of CPU cores to use.
struct_path <- "~/Software/structure/"
#clumpp_path <- "~/Software/CLUMPP/"
Sys.setenv(PATH=paste(Sys.getenv("PATH"), struct_path, clumpp_path, sep=":"))

# -------------------------
# --- END CONFIGURATION ---
# -------------------------


# ====================================================================
# STEP 1: LOAD AND PREPARE DATA
# ====================================================================
message("Step 1: Loading gtypes object and preparing data...")

# Check if the gtypes file exists before proceeding.
if (!file.exists(gtypes_path)) {
  stop(paste("Gtypes file not found at:", gtypes_path,
             "\nPlease ensure script 02 has been run successfully with matching configuration."))
}
load(gtypes_path)

# Stratify data based on the chosen scheme.
g_stratified <- stratify(g, scheme = stratification_scheme, drop = TRUE)

# Filter to include only the specified strata for this analysis.
#g_subset <- g_stratified[,,which(getStrataNames(g_stratified) %in% strata_to_include)]
#message(paste("Data subset to include only:", paste(strata_to_include, collapse = ", ")))


# ====================================================================
# STEP 2: RUN STRUCTURE ANALYSIS
# ====================================================================
message("\nStep 2: Running STRUCTURE. This may take a very long time...")

# Run the STRUCTURE analysis across the specified range of K.
structure_run_label <- paste(run_label, "noloc", sep = "_")
sr_noloc <- structureRun(
  g_stratified, 
  k.range = k_range, 
  num.k.rep = num_k_reps, 
  label = structure_run_label, 
  delete.files = FALSE, # Set to FALSE for debugging if needed
  num.cores = num_cores,
  burnin = burnin, 
  numreps = numreps,
  noadmix = FALSE, 
  freqscorr = TRUE
)

structure_run_label <- paste(run_label, "wloc", sep = "_")
sr_wloc <- structureRun(
  g_stratified, 
  k.range = k_range, 
  num.k.rep = num_k_reps, 
  label = structure_run_label, 
  delete.files = FALSE, # Set to FALSE for debugging if needed
  num.cores = num_cores,
  burnin = burnin, 
  numreps = numreps,
  noadmix = FALSE, 
  freqscorr = TRUE, 
  pop.prior = "locprior"
)

##Using popflag
strata <- as.data.frame(getStrata(g_stratified))
head(strata)
strata<-rownames_to_column(strata)
head(strata)
colnames(strata)<-c("LabID", "Stratum_ABO")
head(strata)
dim(strata)
unique(strata$Stratum_ABO)

population_key <- strata %>%
  distinct(Stratum_ABO) %>% # Get unique populations in their order of appearance
  mutate(Number = row_number()) %>% # Assign a number to each (1, 2, 3...)
  select(Number, Population_Name = Stratum_ABO) # R
population_key

popflag<-data.frame(LabID=unique(getIndNames(g_stratified))) %>%
  left_join(strata, by="LabID") %>%
  transmute(flag = ifelse(Stratum_ABO %in% c("In-water, CA", "In-water, Indonesia"), 0, 1))

structure_run_label <- paste(run_label, "popflag", sep = "_")
sr_popflag_nogensback<-structureRun(
  g_stratified, 
  k.range = k_range, 
  num.k.rep = num_k_reps, 
  label = structure_run_label, 
  delete.files = FALSE, # Set to FALSE for debugging if needed
  num.cores = num_cores,
  burnin = burnin, 
  numreps = numreps,
  noadmix = FALSE, 
  freqscorr = TRUE, 
  pop.prior="usepopinfo",
  popflag=popflag$flag)
#This works but throws up errors if you are running a k value that is less than the number of actual populations you have.


# Save the raw STRUCTURE results object immediately.
sr_noloc_path <- file.path(results_r_path, paste0(run_label, "_structure_sr_noloc.rda"))
save(sr_noloc, file = sr_noloc_path)

sr_wloc_path <- file.path(results_r_path, paste0(run_label, "_structure_sr_wloc.rda"))
save(sr_wloc, file = sr_wloc_path)

sr_popflag_path <- file.path(results_r_path, paste0(run_label, "_structure_sr_popflag.rda"))
save(sr_popflag_nogensback, file=sr_popflag_path)

# ====================================================================
# STEP 3: SUMMARIZE RESULTS WITH EVANNO METHOD
# ====================================================================
message("\nStep 3: Calculating and plotting Evanno metrics to find the optimal K...")

#sr<-sr_noloc
#sr<-sr_wloc
sr<-sr_popflag_nogensback

# Calculate Evanno metrics (requires K > 1 and multiple reps).
evno <- evanno(sr)
message("Evanno Metrics (summary of Delta K):")
print(evno$df)
evno$plots

# Save the Evanno metrics table.
evanno_path <- file.path(results_raw_path, paste0(run_label, "_evanno_metrics_sr_popflag.csv"))
write.csv(evno$df, file = evanno_path, row.names = FALSE)
message(paste("Evanno metrics table saved to:", evanno_path))

# Create and save a plot of Delta K.
delta_k_plot <- ggplot(evno$df, aes(x = k, y = delta.k)) +
  geom_point(size = 3) +
  geom_line() +
  labs(
    title = "Evanno's Delta K",
    subtitle = "Peak indicates the most likely number of clusters (K)",
    x = "Number of Clusters (K)",
    y = "Delta K"
  ) +
  theme_minimal()
delta_k_plot
plot_path <- file.path(results_raw_path, paste0(run_label, "_delta_k_plot_sr_popflag.pdf"))
ggsave(plot_path, plot = delta_k_plot, device = "pdf", width = 7, height = 5)
message(paste("Delta K plot saved to:", plot_path))

 #Plotting qmatrix for one run
q_matrix_4<-sr$
head(q_matrix_7)
colnames(q_matrix_7)<-c("Indiv", "%Missing", "Stratum", "K=1", "K=2", "K=3", "K=4", "K=5", "K=6", "K=7")
pop_order<-c("Bird's Head-Summer", "Bird's Head-Winter", "HaevoSI-Summer", "HaevoSI-Winter", "IsabelSI-South", "PNG", "Malaysia", "In-water, CA", "In-water, Indonesia")
q_matrix_7<-q_matrix_7 %>%
  mutate(Stratum=factor(Stratum, levels=pop_order)) %>%
  arrange(Stratum)

q_long<-q_matrix_7 %>%
  pivot_longer(
    cols=starts_with("K"),
    names_to = "cluster",
    values_to = "ancestry"
  )
head(q_long)
indiv_order <- unique(q_long$Indiv)

pop_boundaries <- q_matrix_7 %>%
  group_by(Stratum) %>%
  summarise(start_id = min(Indiv), .groups = 'drop')
pop_boundaries

vline_positions <- q_long %>%
  group_by(Stratum) %>%
  summarise(start_pos = match(first(Indiv), indiv_order), .groups = 'drop')

label_positions <- vline_positions %>%
  mutate(
    # The end position is the start of the next stratum minus 1
    end_pos = lead(start_pos) - 1,
    # For the very last stratum, the end is the total number of individuals
    end_pos = ifelse(is.na(end_pos), length(indiv_order), end_pos),
    # Calculate the midpoint
    mid_pos = start_pos + (end_pos - start_pos) / 2
  )


ggplot(q_long, aes(x = factor(Indiv, levels = indiv_order), y = ancestry, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  geom_vline(
    # Use the numeric positions you just calculated
    xintercept = vline_positions$start_pos[-1] - 0.5,
    color = "black",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_text(
    data = label_positions, # Use the new data frame for labels
    aes(x = mid_pos, y = 1.05, label = Stratum), # Set aesthetics
    angle = 60,         # Angle the labels for readability
    hjust = 0,          # Adjust horizontal alignment
    inherit.aes = FALSE # Don't inherit aesthetics from the main plot
  ) +
  theme_minimal() +
  labs(x = "Individual", y = "Ancestry Proportion") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing.x = unit(0, "cm"),
    plot.margin = unit(c(3.5, 1.7, 0.5, 0.5), "cm"), # Increased top margin to 3cm
    legend.position = "none"         
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(clip="off")

# ====================================================================
# STEP 4: PROCESS AND PLOT CLUMPP RESULTS FOR EACH K
# Following this vignette: https://www.royfrancis.com/pophelper/articles/index.html
# ====================================================================
sfiles<-list.files(path="dcor.wpac.final.nodups.noloc.structureRun/", pattern="_f$", full.names=T)
slist<-readQ(files=sfiles, indlabfromfile = TRUE)
head(slist[[1]])
tr1<-tabulateQ(qlist=slist)
head(tr1)
sr1<-summariseQ(tr1)
head(sr1)
em<-evannoMethodStructure(data=sr1)
em_plot<-evannoMethodStructure(data=sr1, exportplot=F, returnplot=T, returndata=F, basesize=12,linesize=0.7)
grid.arrange(em_plot)
#plotting
slist1<-alignK(qlist=slist)
slist1<-mergeQ(qlist=slist1)
p1 <- plotQ(slist1,imgoutput="join",returnplot=T,exportplot=F,basesize=11)
grid.arrange(p1$plot[[1]])

slist2<-sortQ2(qlist=slist)
slist2<-mergeQ2(qlist=slist2)
p2<-plotQ(slist2,imgoutput="join",returnplot=T,exportplot=F,basesize=11, 
          grplab=strata_only, grplabangle = 45,
          showindlab = FALSE, ordergrp = TRUE,
          grplabpos=1, grplabspacer = -.1, linesize = 0, 
          showdiv = TRUE, pointsize = 0) 
ggsave("results-raw/structure_plots/structure_plots_k1.10_wloc.png", grid.arrange(p2$plot[[1]]), width=7, height=5)
grid.arrange(p2$plot[[1]])

p3<-plotQ(slist2[7],returnplot=T,exportplot=F,basesize=10, 
          grplab=strata_only, grplabangle = 45,
          showindlab = FALSE, ordergrp = TRUE, sortind="all",
          grplabpos=1, grplabspacer = 0, linesize = 0, grplabsize = 4, 
          showdiv = TRUE, pointsize = 0) 
grid.arrange(p3$plot[[1]])
ggsave("results-raw/structure_plots/structure_plots_k7_wloc.png", grid.arrange(p3$plot[[1]]), width=7.5, height=5)

strata_only<-strata$Stratum_ABO
strata_only<-as.data.frame(strata_only)
colnames(strata_only)<-"Stratum"

#Workaround functions from https://github.com/royfrancis/pophelper/issues/91
mergeQ2 <- function(qlist) {
  
  is.qlist(qlist)
  if(diff(range(as.integer(tabulateQ(qlist)$ind)))!=0) stop("mergeQ: Number of individuals differ between runs.")
  
  # Computes mean cell-wise across dataframes
  # x A list of numeric dataframes
  # 
  mergy <- function(x) {
    return(list(Reduce(`+`, x)/length(x)))
  }
  
  # if all runs have same K, merge as is
  if(diff(range(as.integer(tabulateQ(qlist)$k)))==0) {
    labels <- summariseQ(tabulateQ(qlist))$k
    x <- mergy(qlist)
    names(x) <- labels
  }else{
    # if runs have different K, split them and merge within sublists
    qlist <- sortQ2(qlist)
    labels <- summariseQ(tabulateQ(qlist,sorttable=FALSE))$k
    x <- unlist(lapply(splitQ(qlist),mergy),recursive=FALSE)
    names(x) <- labels
  }
  
  return(as.qlist(x))
}
sortQ2 <- function(qlist,by="k",decreasing=FALSE,debug=FALSE) {
  
  is.qlist(qlist)
  if(length(by)==0) stop("sortQ: Argument 'by' must not be length zero.")
  if(!is.character(by)) stop("sortQ: Argument 'by' must be a character.")
  if(!is.logical(decreasing)) stop("sortQ: Argument 'decreasing' must be a logical datatype.")
  
  fun1 <- function(x) as.matrix(unlist(attributes(x)))
  a <- lapply(qlist,fun1)
  if(debug) print(a)
  if(any(!sapply(a,function(x) any(grepl(paste0(by,collapse="|"),rownames(x)))))) {
    stop(paste0("One or more of the attributes provided in by (",by,") is missing in one or more runs. If 'ind' or 'k' is missing, use 'as.qlist()' to add them."))
  }
  
  # get df of attributes
  b <- as.data.frame(t(as.data.frame(lapply(a,function(x,y) x[y,],by),stringAsFactors=FALSE)),stringsAsFactors=FALSE)
  fun2 <- function(x) if(all(!is.na(as.numeric(as.character(x))))) {return(as.numeric(as.character(x)))}else{return(x)}
  b <- as.data.frame(sapply(b,fun2),stringAsFactors=FALSE)
  
  if(debug) {print(str(b)); print(b)}
  
  # order
  ord <- do.call(order,b[,by,drop=FALSE])
  if(decreasing) ord <- rev(ord)
  # sort qlist
  return(qlist[ord])
}

##Script for processing the "popflag" results that are in a different format
# install.packages("readr") # Run this line if you don't have the 'readr' package
library(readr)
library(dplyr)
library(stringr)

#' Extracts the "Inferred ancestry of individuals" table from STRUCTURE output files.
#'
#' @param input_directory The path to the folder containing your '_f.txt' files.
#' @param output_directory The path to the folder where CSV files will be saved.

process_structure_files(input_directory = "dcor.wpac.final.nodups.popflag.structureRun/", output_directory = "dcor.wpac.final.nodups.popflag.structureRun/")
qlist<-import_ancestry_csvs(input_directory = "dcor.wpac.final.nodups.popflag.structureRun/")

qlist<-as.qlist(qlist)
head(slist[[1]])
tr1<-tabulateQ(qlist=qlist)
head(tr1)
sr1<-summariseQ(tr1)
head(sr1)
em<-evannoMethodStructure(data=sr1)
em_plot<-evannoMethodStructure(data=sr1, exportplot=F, returnplot=T, returndata=F, basesize=12,linesize=0.7)
grid.arrange(em_plot)

qlist2<-sortQ2(qlist=qlist)
qlist2<-mergeQ2(qlist=qlist2)
qp<-plotQ(qlist2,imgoutput="join",returnplot=T,exportplot=F,basesize=11, 
          grplab=strata_only, grplabangle = 45,
          showindlab = FALSE, ordergrp = TRUE,
          grplabpos=1, grplabspacer = -.1, linesize = 0, 
          showdiv = TRUE, pointsize = 0) 
grid.arrange(qp$plot[[1]])
ggsave("results-raw/structure_plots/structure_plots_k1.10_popflag.png", grid.arrange(qp$plot[[1]]), width=7, height=5)

process_structure_files <- function(input_directory = ".", output_directory = ".") {
  
  # Find all files in the directory that end with "_f"
  files_to_process <- list.files(
    path = input_directory,
    pattern = "_f$",
    full.names = TRUE
  )
  
  if (length(files_to_process) == 0) {
    stop("No files ending with '_f' found in the specified directory.")
  }
  
  cat("Found", length(files_to_process), "files to process.\n\n")
  
  # Loop through each file found
  for (file_path in files_to_process) {
    cat("Processing:", basename(file_path), "\n")
    
    # Read all lines from the file
    lines <- read_lines(file_path)
    
    # --- Dynamically find the number of clusters (K) ---
    k_line <- lines[str_detect(lines, "populations assumed")]
    num_clusters <- as.numeric(str_extract(k_line, "\\d+"))
    
    if (is.na(num_clusters) || length(num_clusters) == 0) {
      warning(paste("Could not determine K (number of clusters) in", basename(file_path), "- Skipping."))
      next
    }
    
    # Find the starting line of the target table
    start_index <- which(str_detect(lines, "^Inferred ancestry of individuals:"))
    
    if (length(start_index) == 0) {
      warning(paste("Could not find the ancestry table in", basename(file_path)))
      next
    }
    
    # The actual data starts 2 lines after the header
    data_lines <- lines[(start_index + 2):length(lines)]
    
    # Find the end of the table (it ends when a blank line is encountered)
    end_index <- which(data_lines == "")[1]
    if (!is.na(end_index)) {
      data_lines <- data_lines[1:(end_index - 1)]
    }
    
    # Process each line to extract the data
    all_rows <- lapply(data_lines, function(line) {
      
      clean_line <- str_trim(line)
      # This pattern matches "1 z0006796 (9) 7 :"
      main_info_part <- str_extract(clean_line, "^\\d+\\s+\\S+\\s+\\([:alnum:]+\\)\\s+\\d+\\s+:")
      
      if (is.na(main_info_part)) {
        warning(paste("Could not parse line in", basename(file_path), ":\n  '", line, "'\n  Skipping line."))
        return(NULL)
      }
      
      # Extract "Label", "(%Miss)", and "Pop"
      line_parts <- str_split(str_trim(main_info_part), "\\s+")[[1]]
      sample_id <- line_parts[2]
      population <- as.numeric(line_parts[4])
      
      q_values_part <- str_remove(clean_line, fixed(main_info_part))
      
      if (is.na(q_values_part)) return(NULL)
      
      # --- Updated Logic to handle both data formats ---
      if (str_detect(q_values_part, "\\|")) {
        # Format 1: "0.961 | Pop 1: 0.000 0.001 0.004 | ..."
        parts <- str_split(q_values_part, "\\|")[[1]]
        
        q_vector <- numeric(num_clusters)
        
        # The first part is the main ancestry for the assumed population
        main_ancestry <- as.numeric(str_trim(parts[1]))
        if (!is.na(main_ancestry)) {
          q_vector[population] <- main_ancestry
        }
        
        # Process the remaining parts for other populations
        if (length(parts) > 1) {
          for (i in 2:length(parts)) {
            other_pop_str <- parts[i]
            other_pop_num <- as.numeric(str_extract(other_pop_str, "(?<=Pop )\\d+"))
            
            # --- CORRECTED LOGIC ---
            # Extract all numbers and take the FIRST one.
            other_pop_values <- as.numeric(str_extract_all(other_pop_str, "\\d+\\.\\d+")[[1]])
            
            if (!is.na(other_pop_num) && length(other_pop_values) >= 1) {
              point_estimate <- other_pop_values[1] # Always take the first value
              q_vector[other_pop_num] <- point_estimate
            }
          }
        }
        
      } else {
        # Format 2: "0.171 0.153 0.157 ..."
        q_values_str <- str_trim(q_values_part)
        if (nchar(q_values_str) == 0) {
          q_vector <- c(1.0) # Case for K=1
        } else {
          q_vector <- as.numeric(str_split(q_values_str, "\\s+")[[1]])
        }
      }
      
      if(length(q_vector) < num_clusters) {
        length(q_vector) <- num_clusters
      }
      
      return(c(SampleID = sample_id, AssumedPop = population, q_vector))
    })
    
    all_rows <- all_rows[!sapply(all_rows, is.null)]
    
    if (length(all_rows) == 0) {
      cat(" -> No valid data rows found. Skipping file.\n\n")
      next
    }
    
    result_df <- as.data.frame(do.call(rbind, all_rows))
    
    num_q_cols <- ncol(result_df) - 2
    colnames(result_df) <- c("SampleID", "AssumedPop", paste0("Cluster", 1:num_q_cols))
    
    result_df <- result_df %>%
      mutate(across(-SampleID, as.numeric))
    
    output_filename <- str_replace(basename(file_path), "_f$", "_ancestry.csv")
    output_path <- file.path(output_directory, output_filename)
    
    write.csv(result_df, output_path, row.names = FALSE)
    cat(" -> Successfully saved to:", output_path, "\n\n")
  }
}

#' Imports and formats ancestry CSV files into a list of q-matrices.
#'
#' This function searches a directory for CSV files ending in "_ancestry.csv",
#' reads each one, and formats it into a q-matrix. The SampleID becomes the
#' row names, and only the cluster membership columns are kept.
#'
#' @param input_directory The path to the folder containing your '_ancestry.csv' files.
#' @return A named list where each element is a q-matrix (data.frame). The list
#'   is also given the class 'qlist' for compatibility with the pophelper package.

import_ancestry_csvs <- function(input_directory = ".") {
  
  # Find all files in the directory that end with "_ancestry.csv"
  files_to_import <- list.files(
    path = input_directory,
    pattern = "_ancestry\\.csv$",
    full.names = TRUE
  )
  
  if (length(files_to_import) == 0) {
    stop("No files ending with '_ancestry.csv' found in the specified directory.")
  }
  
  cat("Found", length(files_to_import), "ancestry CSV files to import.\n\n")
  
  # Use lapply to loop through each file path, read, and process it
  q_matrix_list <- lapply(files_to_import, function(file_path) {
    
    # Read the CSV file
    ancestry_df <- read_csv(file_path, show_col_types = FALSE)
    
    # --- Format into a q-matrix ---
    q_matrix <- as.data.frame(ancestry_df %>% select(starts_with("Cluster")))
    rownames(q_matrix) <- ancestry_df$SampleID
    
    # --- Find and parse the original _f file for metadata ---
    original_f_file <- str_replace(file_path, "_ancestry\\.csv$", "_f")
    
    if (file.exists(original_f_file)) {
      cat(" -> Reading metadata from:", basename(original_f_file), "\n")
      f_lines <- read_lines(original_f_file)
      
      # Helper function to safely extract a numeric value from a line
      extract_value <- function(pattern) {
        line <- f_lines[str_detect(f_lines, pattern)]
        # Extracts floating point or integer numbers, including negative signs
        as.numeric(str_extract(line, "-?\\d+(\\.\\d+)?"))
      }
      
      # Attach attributes to the q-matrix object
      attr(q_matrix, "ind") <- extract_value("individuals")
      attr(q_matrix, "k") <- extract_value("populations assumed")
      attr(q_matrix, "loci") <- extract_value("loci")
      attr(q_matrix, "burnin") <- extract_value("Burn-in period")
      attr(q_matrix, "reps") <- extract_value("Reps")
      attr(q_matrix, "elpd") <- extract_value("Estimated Ln Prob of Data")
      attr(q_matrix, "mvll") <- extract_value("Mean value of ln likelihood")
      attr(q_matrix, "vll") <- extract_value("Variance of ln likelihood")
      
    } else {
      warning(paste("Could not find corresponding _f file for:", basename(file_path)))
    }
    
    return(q_matrix)
  })
  
  # --- Set the names for the list elements ---
  list_names <- basename(files_to_import) %>%
    str_remove("_ancestry\\.csv$")
  
  names(q_matrix_list) <- list_names
  
  # Add the 'qlist' class for pophelper compatibility
  class(q_matrix_list) <- c("qlist", "list")
  
  cat("\n✅ Successfully imported and formatted", length(q_matrix_list), "files into a qlist.\n")
  
  return(q_matrix_list)
}

