# This script defines a function to find the longest read in each FASTQ file
# and plot a histogram of these maximum lengths.

#' Plot a histogram of the longest sequence length from each FASTQ file in a directory.
#'
#' @param directory A string path to the directory containing FASTQ files.
#' @param output_file A string path for the output histogram image (e.g., "max_lengths.png").
#' @param plot_title A string for the main title of the plot.
#' @return Invisible NULL. The function saves a plot to a file.
#' @examples
#' # To use this function, you would call it like this:
#' # plot_max_sequence_lengths(directory = "/path/to/fastq", output_file = "max_lengths_hist.png")

plot_max_sequence_lengths <- function(directory, output_file, plot_title = "Distribution of Longest Read Lengths per File") {
  
  # --- Main Logic ---
  
  # Find all .fastq or .fastq.gz files in the specified directory
  cat("Searching for FASTQ files in:", directory, "\n")
  files <- list.files(path = directory, pattern = "\\.fastq(\\.gz)?$", full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No FASTQ files found in the specified directory.")
  }
  cat("Found", length(files), "files to process.\n")
  
  # Function to find the maximum sequence length from a single FASTQ file
  get_max_length <- function(filepath) {
    cat("Scanning:", basename(filepath), "\n")
    max_len <- 0
    tryCatch({
      # Open a connection that can handle both plain and gzipped files
      con <- if (grepl("\\.gz$", filepath)) gzfile(filepath, "r") else file(filepath, "r")
      lines <- readLines(con)
      close(con)
      
      # The sequence is on the 2nd line of every 4-line record
      sequence_lines <- lines[seq(2, length(lines), by = 4)]
      
      if (length(sequence_lines) > 0) {
        max_len <- max(nchar(sequence_lines))
      }
      cat("  -> Longest read found:", max_len, "bases\n")
      return(max_len)
    }, error = function(e) {
      warning("Could not process file: ", basename(filepath), " | Error: ", e$message)
      return(0) # Return 0 on error
    })
  }
  
  # Apply the function to all files to get a vector of max lengths
  max_lengths <- sapply(files, get_max_length)
  
  # Filter out files where no sequences were found
  max_lengths <- max_lengths[max_lengths > 0]
  
  if (length(max_lengths) == 0) {
    stop("No sequences were found in any of the files.")
  }
  
  # --- Plotting ---
  cat("\nGenerating histogram and saving to:", output_file, "\n")
  png(output_file, width = 1200, height = 700, res = 100)
  
  # Adjust bins based on the number of files, but cap it to avoid overly thin bars
  num_bins <- min(length(max_lengths), 50)
  
  hist(max_lengths,
       breaks = num_bins,
       col = "coral",
       border = "black",
       main = plot_title,
       xlab = "Longest Sequence Length (bases)",
       ylab = "Number of Files")
  grid(ny = NULL, nx = NA, col = "gray", lty = "dotted")
  dev.off()
  
  cat("Histogram saved successfully!\n")
}

# --- Example Usage ---
#
# # 1. Source this file in your R session:
# #    source("plot_max_lengths_function.R")
# #
# # 2. Then call the function with your parameters:
# #    plot_max_sequence_lengths(directory = "path/to/your/fastq_files",
# #                              output_file = "max_sequence_lengths.png")
#
