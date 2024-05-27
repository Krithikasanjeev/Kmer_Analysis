# Load necessary libraries 
library(Biostrings)
library(dplyr)

# Set working directory and folder path for fasta files
setwd("C:/Users/Krithika/Desktop/IBAB/kmer/ecoli_kmers/AMK/")
folder_path <- "C:/Users/Krithika/Desktop/IBAB/e.coli_dataset/AMK"

# Define k-mer lengths to analyze
kmer_lengths <- 2:15

# Get list of fasta files
fasta_files <- list.files(folder_path, pattern = "\\.fa$|\\.fasta$", full.names = TRUE)
file_names <- tools::file_path_sans_ext(basename(fasta_files))

# Initialize data frame for entropy values
entropy_df <- data.frame(matrix(NA, nrow = length(fasta_files), ncol = length(kmer_lengths)))
colnames(entropy_df) <- paste0("kmer_", kmer_lengths)
rownames(entropy_df) <- file_names

# Function to calculate entropy
calculate_entropy <- function(counts) {
  probabilities <- counts / sum(counts)
  return(-sum(probabilities * log2(probabilities)))
}

# Loop through each k-mer length and fasta file
for (kmer_length in kmer_lengths) {
  # Create directory for current k-mer length
  kmer_dir <- paste0("amk_", kmer_length, "_mer")
  if (!file.exists(kmer_dir)) dir.create(kmer_dir)
  
  for (fasta_file in fasta_files) {
    sequence <- readDNAStringSet(fasta_file, format = "fasta")
    sequence <- gsub("N", "", as.character(sequence))
    
    # Generate k-mer counts
    kmers <- sapply(1:(nchar(sequence) - kmer_length + 1), 
                    function(i) substr(sequence, i, i + kmer_length - 1))
    kmer_counts <- table(kmers)
    
    # Filter out invalid k-mers
    kmer_counts <- kmer_counts[grepl("^[ATGC]+$", names(kmer_counts))]
    
    # Write k-mer counts to file
    output_file <- file.path(kmer_dir, paste0(tools::file_path_sans_ext(basename(fasta_file)), "_", kmer_length, "mer.txt"))
    write.table(as.data.frame(kmer_counts), file = output_file, col.names = c("pattern", "count"), row.names = FALSE)
    
    # Calculate and store entropy
    entropy <- calculate_entropy(kmer_counts)
    entropy_df[file_names == tools::file_path_sans_ext(basename(fasta_file)), paste0("kmer_", kmer_length)] <- entropy
    
    # Print result
    cat("Filename:", basename(fasta_file), "Kmer:", kmer_length, "Entropy:", entropy, "\n")
  }
}

# Write the entropy data frame to a CSV file
write.csv(entropy_df, paste(basename(folder_path, "csv", sep = ".")))
