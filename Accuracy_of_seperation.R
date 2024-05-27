library(dplyr)
library(ggplot2)

# Function to generate file paths
generate_file_path <- function(base_dir, prefix, kmer_length, suffix) {
  file_name <- paste(prefix, kmer_length, suffix, sep = "_")
  return(file.path(base_dir, file_name))
}

# Function to calculate the AOS for a given k-mer length
calculate_aos <- function(file_pos, file_neg, kmer_length) {
  # Read positive and negative datasets
  pos_file <- read.csv(file_pos)
  neg_file <- read.csv(file_neg)
  
  # Calculate density for positive and negative datasets
  density_neg <- density(neg_file$super_entropy)
  density_pos <- density(pos_file$super_entropy)
  
  # Generate a sequence of x values spanning the range of both distributions
  x_values <- seq(min(min(density_neg$x), min(density_pos$x)),
                  max(max(density_neg$x), max(density_pos$x)), by = 0.00001)
  
  # Interpolate density values for both distributions at x_values
  density_values_neg <- approx(density_neg$x, density_neg$y, xout = x_values)$y
  density_values_pos <- approx(density_pos$x, density_pos$y, xout = x_values)$y
  
  # Calculate the maximum density at each x value
  max_densities <- pmax(density_values_neg, density_values_pos, na.rm = TRUE)
  
  # Create a data frame with x_values and max_densities
  max_density_df <- data.frame(x_values = x_values, max_densities = max_densities)
  
  # Filter out rows where max_densities is not equal to 0
  max_density_df <- max_density_df %>%
    filter(max_densities != 0)
  
  # Plot the maximum density curve
  plot <- ggplot(max_density_df, aes(x = x_values, y = max_densities)) +
    geom_line() +
    labs(title = paste("Maximum Density across X-axis for kmer length", kmer_length), 
         x = "X-axis", y = "Maximum Density") +
    theme_minimal()
  
  # Print the plot
  print(plot)
  
  # Calculate the integral of the maximum density curve
  integral_value <- sum(diff(max_density_df$x_values) *
                          (max_density_df$max_densities[-1] +
                             max_density_df$max_densities[-length(max_density_df$max_densities)]) / 2)
  
  # Calculate the AOS
  aos <- integral_value / 2
  
  # Print AOS for debugging
  print(paste("AOS for kmer length", kmer_length, ":", aos))
  
  return(aos)
}

# Main function to execute the AOS calculation
main <- function(kmer_lengths, pos_dir, neg_dir, output_file) {
  # Initialize an empty data frame to store AOS values
  aos_df <- data.frame(kmer_length = numeric(), aos = numeric())
  
  # Iterate over each k-mer length
  for (kmer_length in kmer_lengths) {
    # Generate file paths for positive and negative datasets
    file_pos <- generate_file_path(pos_dir, "crabtree_positive_stripped", kmer_length, "super_entropy_non_halotolerant.csv")
    file_neg <- generate_file_path(neg_dir, "crabtree_negative_stripped", kmer_length, "super_entropy_non_halotolerant.csv")
    
    # Calculate AOS for the given k-mer length
    aos <- calculate_aos(file_pos, file_neg, kmer_length)
    
    # Add kmer length and AOS to aos_df
    aos_df <- rbind(aos_df, data.frame(kmer_length = kmer_length, aos = aos))
  }
  
  # Write AOS values to a CSV file
  write.csv(aos_df, file = output_file, row.names = FALSE)
}

# Run the main function
main(
  kmer_lengths = 2:15,
  pos_dir = "C:/Users/Krithika/Desktop/IBAB/kmer/crabtree/kmer_cp",
  neg_dir = "C:/Users/Krithika/Desktop/IBAB/kmer/crabtree/kmer_cn",
  output_file = "aos_values(super_information)1.csv"
)
