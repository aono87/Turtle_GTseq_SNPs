# Load necessary libraries
# Ensure you have these installed: install.packages(c("dplyr", "tidyr", "ggplot2"))
library(dplyr)
library(tidyr)
library(ggplot2)

# --- 1. SETUP YOUR PARAMETERS HERE ---

# Your main results object from strataG
# Make sure this object is loaded in your R environment
structure_results <- sr_popflag_nogensback 

# Define the range of K values you want to plot
k_range <- 2:7 

# The base name of your runs (everything before the '.k' part)
run_basename <- "dcor.wpac.final.nodups.structureRun"

# Define the specific order for your populations in the plot
pop_order <- c(
  "Bird's Head-Summer", "Bird's Head-Winter", "HaevoSI-Summer", 
  "HaevoSI-Winter", "IsabelSI-South", "PNG", "Malaysia", 
  "In-water, CA", "In-water, Indonesia" # CORRECTED: Changed "In-water" to "In-water, CA"
)


# --- 2. FUNCTION TO CREATE THE STRUCTURE PLOT ---

create_structure_plot <- function(q_matrix, k_value, pop_order) {
  
  # Convert matrix to a data frame and set column names
  q_df <- as.data.frame(q_matrix)
  colnames(q_df) <- c("Indiv", "%Missing", "Stratum", paste0("Cluster", 1:k_value))
  
  # Re-order data and pivot to long format for ggplot
  q_long <- q_df %>%
    mutate(Stratum = factor(Stratum, levels = pop_order)) %>%
    # Drop any individuals not in the specified pop_order
    filter(!is.na(Stratum)) %>% 
    arrange(Stratum) %>%
    pivot_longer(
      cols = starts_with("Cluster"),
      names_to = "cluster",
      values_to = "ancestry"
    )
  
  # Get the final, sorted order of individuals for the x-axis
  indiv_order <- unique(q_long$Indiv)
  
  # Calculate positions for labels and boundary lines
  label_positions <- q_long %>%
    group_by(Stratum) %>%
    summarise(start_pos = match(first(Indiv), indiv_order), .groups = 'drop') %>%
    mutate(
      end_pos = lead(start_pos) - 1,
      end_pos = ifelse(is.na(end_pos), length(indiv_order), end_pos),
      mid_pos = start_pos + (end_pos - start_pos) / 2
    )
  
  # Generate the plot
  ggplot(q_long, aes(x = factor(Indiv, levels = indiv_order), y = ancestry, fill = cluster)) +
    geom_bar(stat = "identity", width = 1) +
    
    # CORRECTED: Replaced geom_vline with geom_segment to control line height
    geom_segment(
      data = label_positions[-1, ],
      aes(x = start_pos - 0.5, xend = start_pos - 0.5, y = 1, yend = 0),
      color = "black", linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE
    ) +
    
    geom_text(
      data = label_positions, 
      aes(x = mid_pos, y = -0.05, label = Stratum),
      # CORRECTED: Changed angle and alignment for readability
      angle = 45,         
      hjust = 1,          
      vjust = 0.5,        
      inherit.aes = FALSE
    ) +
    theme_minimal() +
    labs(
      title = paste("STRUCTURE Plot: K =", k_value),
      x = NULL, # Remove the x-axis label
      y = "Ancestry Proportion"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.spacing.x = unit(0, "cm"),
      # Adjust margin for bottom labels
      plot.margin = unit(c(1, 1, 4, 1), "cm"), 
      legend.position = "none"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(clip = "off")
}


# --- 3. LOOP TO GENERATE AND PRINT PLOTS ---

for (k in k_range) {
  result_name <- paste0(run_basename, ".k", k, ".r1")
  q_mat <- try(structure_results[[result_name]]$q.mat, silent = TRUE)
  
  if (inherits(q_mat, "try-error") || is.null(q_mat)) {
    warning(paste("Skipping K =", k, "- results not found."))
    next
  }
  
  # Create and print the plot using the function
  p <- create_structure_plot(q_mat, k, pop_order)
  print(p)
}
