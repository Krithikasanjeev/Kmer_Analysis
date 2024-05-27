calculate_emd <- function(folder_path) {
  # Initialize a data frame to store calculated Earth Mover's Distance (EMD) sums
  calculated_sums <- data.frame(kmer = character(), EMD = numeric())
    
  # List the files in the directory 
  files <- list.files(folder_path, full.names = TRUE)
  
  # Loop through different files
  for (file in files) {
    # Read the combined data file for the current k-mer length
    file_data <- read.csv(file)
    
    # Extract probability distributions from the file
    dist1 <- file_data$probabilities_positive
    dist2 <- file_data$probabilities_negative
    
    # Initialize the initial distance variable
    di <- 0
    calculated_values <- data.frame(fi = numeric(), mi = numeric(), di = numeric(), mod_val = numeric())
    
    # Loop through each row in the file
    for (j in 1:nrow(file_data)) {
      # Calculate the normalized probabilities for positive and negative distributions
      fi <- dist1[j] / sum(dist1)
      mi <- dist2[j] / sum(dist2)
      
      # Calculate the difference between positive and negative distributions
      value <- fi + di - mi
      
      # Calculate the absolute value of the difference
      mod_val <- abs(value)
      
      # Append the calculated values to the data frame
      calculated_values <- rbind(calculated_values, data.frame(fi = fi, mi = mi, di = di, mod_val = mod_val))
      
      # Update the distance variable for the next iteration
      di <- value
    }
    
    # Remove the first row which was initialized with numeric() in calculated_values
    calculated_values <- calculated_values[-1, ]
    
    # Write the calculated values to a CSV file
    write.csv(calculated_values, paste0(tools::file_path_sans_ext(basename(file)), "_Wasserstein.csv"), row.names = FALSE)
    
    # Calculate the sum of absolute modified values
    EMD <- sum(calculated_values$mod_val)
    
    # Print the file name and its calculated EMD
    print(paste(basename(file), EMD))
    
    # Append the calculated EMD sum to the data frame
    calculated_sums <- rbind(calculated_sums, data.frame(kmer = basename(file), EMD = EMD))
  }
  
  # Write the calculated EMD sums to a CSV file
  write.csv(calculated_sums, file.path(folder_path, "EMD.csv"), row.names = FALSE)
}

# Set the working directory
folder_path <- "C:/Users/Krithika/Desktop/IBAB/kmer/bascillus/combined_data/"

# Call the function to calculate EMD
calculate_emd(folder_path)
