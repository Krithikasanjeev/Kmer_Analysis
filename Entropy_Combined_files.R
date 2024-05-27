# Load necessary libraries
library(ggplot2)    # For data visualization (not used in this code but loaded)
library(dplyr)      # For data manipulation
library(Biostrings) # For handling biological sequences

# Define a function to calculate entropy
calculate_entropy <- function(probabilities) {
  return(-(sum(probabilities * log2(probabilities))))
}

# Define a function to process a single directory
process_directory <- function(dir) {
  # Initialize an empty data frame to store combined k-mer counts
  combined_data <- data.frame(Pattern = character(), Count = integer(), stringsAsFactors = FALSE)
  # Get the list of files within the current subdirectory
  files <- list.files(dir, full.names = TRUE)
  # Loop through each file in the subdirectory
  for (file in files) {
    tryCatch({
      # Read the current file into a data frame
      current_data <- read.table(file, header = FALSE, col.names = c("Pattern", "Count"))
      
      # Combine the current data with the combined data frame
      combined_data <- rbind(combined_data, current_data)
    }, error = function(e) {
      message(sprintf("Error processing file %s: %s", file, e))
    })
  }
  
  # Convert the Count column to numeric
  combined_data$Count <- as.numeric(combined_data$Count)
  combined_data$Count <- as.numeric(as.character(combined_data$Count))
  
  # Filter out rows with NA values in the Count column
  combined_data <- filter(combined_data, !is.na(Count))
  
  # Filter out rows with invalid k-mer patterns (containing characters other than A, T, G, C)
  combined_data <- filter(combined_data, !grepl("[^ATGC]", Pattern))
  
  # Aggregate the counts by k-mer pattern
  summed_data <- aggregate(Count ~ Pattern, data = combined_data, sum)
  
  # Calculate the probabilities of each k-mer
  summed_data$probabilities <- summed_data$Count / sum(summed_data$Count)
  
  # Calculate the entropy for the current subdirectory
  entropy <- calculate_entropy(summed_data$probabilities)
  
  # Define the output file name for the merged data
  output_file_name <- file.path("output", paste0(basename(dir), "_merged.csv"))
  
  # Write the summed data to the output CSV file
  write.csv(summed_data, file = output_file_name, row.names = FALSE, quote = FALSE)
  
  return(entropy)
}

# Set the main directory (configurable)
main_directory <- Sys.getenv("KMER_MAIN_DIRECTORY", "C:/Users/Krithika/Desktop/IBAB/kmer/bascillus/kmer_nht")

# Create an output directory if it doesn't exist
if (!dir.exists("output")) {
  dir.create("output")
}

# Get the list of subdirectories within the main directory
dirs <- list.dirs(main_directory, recursive = FALSE)

# Initialize a data frame to store entropy values for each directory
entropy_data <- data.frame(Directory = character(), Entropy = numeric(), stringsAsFactors = FALSE)

# Loop through each subdirectory
for (dir in dirs) {
  # Process the directory and calculate entropy
  entropy <- process_directory(dir)
  
  # Print the directory name and calculated entropy to the console
  print(paste(basename(dir), entropy))
  
  # Append the entropy value to the entropy data frame
  entropy_data <- rbind(entropy_data, data.frame(Directory = basename(dir), Entropy = entropy))
}

# Write the entropy data frame to a CSV file
write.csv(entropy_data, file = "non_halotolerant_entropy.csv", row.names = FALSE)
