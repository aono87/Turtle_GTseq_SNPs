#!/bin/bash

# =============================================================================
#  Standalone Bash Script to Merge Replicate FASTQ Files
# =============================================================================
#
#  Description:
#  This script reads a list of sample base IDs from a specified text file,
#  finds the corresponding FASTQ files in a source directory, and merges
#  them into a single, combined file in a specified output directory.
#
#  Expected Folder Structure:
#  .
#  ├── merge_replicates.sh  (This script)
#  ├── merge_ids.txt        (File with IDs to merge, one per line)
#  ├── all/                 (Folder containing the original .fastq.gz files)
#  │   ├── 12345.fastq.gz
#  │   └── 12345b.fastq.gz
#  └── merged/               (Folder where merged files will be created)
#
#  How to Use:
#  1. Ensure your files are arranged in the structure described above.
#  2. EDIT the configuration variables below if your filenames or folders differ.
#  3. Open a terminal and run the script: ./merge_replicates.sh
#
# =============================================================================


# --- Step 1: CONFIGURE YOUR PATHS AND FILENAMES ---

# File containing the list of base IDs (one ID per line).
ID_LIST_FILE="merge_ids.txt"

# Directory containing the original FASTQ files to be merged.
SOURCE_DIR="all"

# Directory where the merged FASTQ files will be saved.
OUTPUT_DIR="merged"

# The file extension of your FASTQ files.
FILE_EXTENSION="_msat.fastq"


# --- Step 2: RUN THE SCRIPT (No edits needed below this line) ---
echo "Starting FASTQ file merge process..."
echo "  ID List File:   '${ID_LIST_FILE}'"
echo "  Source Directory: '${SOURCE_DIR}/'"
echo "  Output Directory: '${OUTPUT_DIR}/'"
echo "-------------------------------------"

# Check if the input file exists and is readable
if [ ! -f "$ID_LIST_FILE" ]; then
  echo "Error: Input file '${ID_LIST_FILE}' not found."
  echo "Please create this file or edit the script to point to the correct file path."
  exit 1
fi

# Check if the source directory exists
if [ ! -d "$SOURCE_DIR" ]; then
  echo "Error: Source directory '${SOURCE_DIR}' not found."
  exit 1
fi

# Create the output directory if it doesn't already exist
mkdir -p "$OUTPUT_DIR"

# Loop through each line (ID) in the specified file
while IFS= read -r id || [[ -n "$id" ]]; do
  
  # Skip any empty lines in the file
  if [ -z "$id" ]; then
    continue
  fi
  
  # Define the input file pattern and the name for the final merged file
  input_files_pattern="${SOURCE_DIR}/${id}*${FILE_EXTENSION}"
  output_file="${OUTPUT_DIR}/${id}_merged${FILE_EXTENSION}"
  
  echo "Processing sample base ID: ${id}"
  
  # Check if any files match the pattern before attempting to merge
  # `shopt -s nullglob` prevents errors if no files match
  shopt -s nullglob
  files_to_merge=(${input_files_pattern})
  shopt -u nullglob
  
  if [ ${#files_to_merge[@]} -eq 0 ]; then
      echo "-> WARNING: No files found matching '${input_files_pattern}'. Skipping."
      echo ""
      continue
  fi

  # Use the `cat` command to concatenate all files matching the pattern
  echo "-> Merging ${#files_to_merge[@]} file(s) into '${output_file}'"
  
  # Execute the merge command
  cat "${files_to_merge[@]}" > "${output_file}"
  
  echo "-> Successfully merged files for ${id}."
  echo "" # Add a blank line for cleaner output

done < "$ID_LIST_FILE"

echo "-------------------------------------"
echo "All merging tasks are complete."

