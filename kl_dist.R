library(dplyr) 

# Function to calculate Kullback-Leibler (KL) divergence and Jensen-Shannon Distance (JS) for given positive and negative files
calculate_KL_and_JS <- function(positive_file, negative_file) {
  # Read data from positive and negative files
  positive_data <- read.csv(positive_file, header = TRUE)
  negative_data <- read.csv(negative_file, header = TRUE)
  
  # Get unique patterns from both files
  all_patterns <- union(positive_data$Pattern, negative_data$Pattern)
  
  # Initialize a data frame to store results
  result_row <- data.frame(kmer = numeric(), kl_n_p = numeric(), kl_p_n = numeric(), js_dist = numeric())
  
  # Loop through each pattern and calculate probabilities
  for (pattern in all_patterns) {
    # Find the pattern in the positive file
    positive_prob <- ifelse(pattern %in% positive_data$Pattern, positive_data$probabilities[positive_data$Pattern == pattern], 0)
    
    # Find the pattern in the negative file
    negative_prob <- ifelse(pattern %in% negative_data$Pattern, negative_data$probabilities[negative_data$Pattern == pattern], 0)
    
    # Append the results with probabilities for positive and negative files
    results <- data.frame(Pattern = pattern, probabilities_positive = positive_prob, probabilities_negative = negative_prob, stringsAsFactors = FALSE)
    
    # Filter out rows with zero probabilities
    filtered_data <- results %>% filter(probabilities_positive != 0 | probabilities_negative != 0)
    
    # Calculate Kullback-Leibler (KL) divergence for P/N and N/P
    filtered_data_sum <- filtered_data %>%
      mutate(`prob(P/N)` = probabilities_positive * log2(probabilities_positive / probabilities_negative),
             `prob(N/P)` = probabilities_negative * log2(probabilities_negative / probabilities_positive)) %>%
      filter(!is.nan(`prob(N/P)`) & !is.infinite(`prob(N/P)`) & !is.nan(`prob(P/N)`) & !is.infinite(`prob(P/N)`))
    
    # Calculate Jensen-Shannon Distance (JS)
    js_dist <- (0.5 * sum(filtered_data_sum$`prob(P/N)`)) + (0.5 * sum(filtered_data_sum$`prob(N/P)`))
    
    # Append the results for the current pattern
    result_row <- rbind(result_row, data.frame(kmer = pattern, kl_n_p = sum(filtered_data_sum$`prob(N/P)`), kl_p_n = sum(filtered_data_sum$`prob(P/N)`), js_dist = js_dist))
  }
  
  return(result_row)
}

# Set the working directory
setwd("C:/Users/Krithika/Desktop/IBAB/kmer")

# Define the file paths for positive and negative files
positive_file <- "cp_mer_merged.csv"
negative_file <- "cn_mer_merged.csv"
