library(dplyr)
library(Biostrings)

# Function to calculate super entropy for sequences in a directory
calculate_super_entropy <- function(main_dir, kmer_lengths) {
  # Loop through each k-mer length
  for (kmer_length in kmer_lengths) {
    # Define the output directory for the current k-mer length
    output_dir <- file.path(getwd(), paste(basename(main_dir), kmer_length, "non_overlapping", sep = "_"))
    
    # Create the output directory if it does not exist
    dir.create(output_dir, showWarnings = FALSE)
    
    # Initialize a data frame to store super entropy results for each file
    super_entropy_results <- data.frame(file = character(), super_entropy = numeric(), stringsAsFactors = FALSE)
    
    # List all files in the main directory
    files <- list.files(main_dir, full.names = TRUE)
    
    # Loop through each file in the directory
    for (file in files) {
      # Read the sequence from the file
      sequence <- readDNAStringSet(file, format = "fasta")
      
      # Initialize a data frame to store entropy results for all k-mers
      all_results <- data.frame(kmer = character(), entropy = numeric(), stringsAsFactors = FALSE)
      
      # Loop through the sequence in steps of k-mer length to get non-overlapping k-mers
      for (i in seq(1, nchar(sequence), by = kmer_length)) {
        # Extract the k-mer of length kmer_length
        kmer <- substr(sequence, i, i + kmer_length - 1)
        
        # Skip if k-mer is shorter than the desired length
        if (nchar(kmer) < kmer_length) { next }
        
        # Count the occurrences of each nucleotide in the k-mer
        sequence_counts <- table(factor(strsplit(kmer, "")[[1]], levels = c("A", "T", "C", "G")))
        
        # Calculate probabilities for each nucleotide
        probabilities <- sequence_counts / sum(sequence_counts)
        
        # Calculate the entropy for the k-mer
        entropy <- -(sum(probabilities * log2(probabilities)))
        
        # Create a data frame with the k-mer and its entropy
        entropy_result <- data.frame(pattern = kmer, entropy = entropy)
        
        # Combine the entropy result with all results
        all_results <- rbind(all_results, entropy_result)
      }
      
      # Initialize a data frame to store entropy probabilities
      prob_results <- data.frame()
      
      # Loop through each k-mer in all results to calculate entropy probabilities
      for (i in 1:nrow(all_results)) {
        pattern <- all_results$pattern[i]
        entropy <- all_results$entropy[i]
        entropy_prob <- all_results$entropy[i] / sum(all_results$entropy)
        entropy_prob_results <- data.frame(pattern = pattern, entropy = entropy, entropy_prob = entropy_prob)
        prob_results <- rbind(prob_results, entropy_prob_results)
      }
      
      # Define the output file name for the current k-mer length
      output_file <- paste(tools::file_path_sans_ext(basename(file)), kmer_length, "mer_no_entropy_results.txt", sep = "_")
      
      # Write the entropy probability results to the output file
      write.table(prob_results, file.path(output_dir, output_file), sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Print a message indicating that the entropy results have been written to the file
      print(paste("Entropy results for", basename(file), "have been written to", output_file))
      
      # Filter out rows where entropy probability is zero
      filtered_data <- prob_results %>%
        filter(entropy_prob != 0)
      
      # Calculate the super entropy for the file
      super_entropy <- -sum(filtered_data$entropy_prob * log2(filtered_data$entropy_prob))
      
      # Create a data frame with the file name and super entropy
      ind_super_entropy <- data.frame(file = sub("_modified.fasta(\\d+)_results.txt", "", basename(file)), super_entropy = super_entropy)
      
      # Print the super entropy (for debugging purposes)
      print(super_entropy)
      
      # Combine the super entropy result with the super entropy results
      super_entropy_results <- rbind(super_entropy_results, ind_super_entropy)
    }
    
    # Write the super entropy results to a CSV file
    write.csv(super_entropy_results, paste(basename(main_dir), kmer_length, "super_entropy_non_halotolerant.csv", sep = "_"))
  }
}

# Set the working directory
main_dir <- "C:/Users/Krithika/Desktop/IBAB/kmer/bascillus/non_halotolerant/"

# Define the range of k-mer lengths to analyze
kmer_lengths <- 2:15

# Call the function to calculate super entropy
calculate_super_entropy(main_dir, kmer_lengths)
