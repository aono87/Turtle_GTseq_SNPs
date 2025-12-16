# This script defines a function to merge FASTQ files and create a histogram
# of ALL sequence lengths.

#' Plot a histogram of all sequence lengths from FASTQ files in a directory.
#'
#' @param directory A string path to the directory containing FASTQ files.
#' @param output_file A string path for the output histogram image (e.g., "all_lengths.png").
#' @param plot_title A string for the main title of the plot.
#' @return Invisible NULL. The function saves a plot to a file.
#' @examples
#' # To use this function, you would call it like this:
#' # plot_all_sequence_lengths(directory = "/path/to/fastq", output_file = "histogram.png")

plot_all_sequence_lengths <- function(directory, output_file, plot_title = "Distribution of All Sequence Lengths in FASTQ Files") {
  
  # --- Main Logic ---
  
  # Find all .fastq or .fastq.gz files in the specified directory
  cat("Searching for FASTQ files in:", directory, "\n")
  files <- list.files(path = directory, pattern = "\\.fastq(\\.gz)?$", full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No FASTQ files found in the specified directory.")
  }
  cat("Found", length(files), "files to process.\n")
  
  # Function to extract all sequence lengths from a single FASTQ file
  get_all_lengths <- function(filepath) {
    cat("Processing:", basename(filepath), "\n")
    tryCatch({
      # Open a connection that can handle both plain and gzipped files
      con <- if (grepl("\\.gz$", filepath)) gzfile(filepath, "r") else file(filepath, "r")
      # Reading all lines can be memory intensive for huge files.
      lines <- readLines(con)
      close(con)
      
      # The sequence is on the 2nd line of every 4-line record
      sequence_lines <- lines[seq(2, length(lines), by = 4)]
      return(nchar(sequence_lines))
    }, error = function(e) {
      warning("Could not process file: ", basename(filepath), " | Error: ", e$message)
      return(NULL) # Return NULL on error
    })
  }
  
  # Apply the function to all files and combine the lengths into one large vector
  all_lengths <- unlist(lapply(files, get_all_lengths))
  
  if (length(all_lengths) == 0) {
    stop("No sequences were found in any of the files.")
  }
  
  cat("\nTotal sequences found:", length(all_lengths), "\n")
  
  # --- Plotting ---
  cat("Generating histogram and saving to:", output_file, "\n")
  
  # Get the minimum length from the data to start the axis
  min_len <- min(all_lengths)
  
  # Set a fixed upper limit for the plot as requested, e.g., 175bp
  max_bp_limit <- 175
  
  # Create breaks for every base pair up to the limit.
  # This gives a very granular view of the length distribution.
  plot_breaks <- seq(floor(min_len), max_bp_limit, by = 1)
  
  png(output_file, width = 1200, height = 700, res = 100)
  hist(all_lengths,
       breaks = plot_breaks,
       col = "skyblue",
       border = "black",
       main = plot_title,
       xlab = "Sequence Length (bases)",
       ylab = "Frequency",
       xlim = c(floor(min_len), max_bp_limit)) # Explicitly set x-axis limits
  grid(ny = NULL, nx = NA, col = "gray", lty = "dotted")
  dev.off()
  
  cat("Histogram saved successfully!\n")
}

# --- Example Usage ---
#
# # 1. Source this file in your R session:
# #    source("plot_all_lengths_function.R")
# #
# # 2. Then call the function with your parameters:
# #    plot_all_sequence_lengths(directory = "path/to/your/fastq_files", 
# #                              output_file = "all_sequence_lengths.png")
#

