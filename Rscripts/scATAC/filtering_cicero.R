library(tidyverse)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check that we have two arguments (for input and output files)
if (length(args) != 3) {
  stop("Please provide both an input and output file name.")
}

input_file <- args[1]
threshold <- args[2]
output_file <- args[3]

# Read the file
data <- read_csv(input_file, col_names = c("id", "Peak1", "Peak2", "value"), 
                 skip = 1)

# Process the data
output <- data %>% 
  select(-id) %>% # Remove the first column
  mutate(generate_name = str_c(Peak1, Peak2, sep = "__")) %>% 
  separate(Peak1, into = c("chr1", "start1", "end1"), sep = "_") %>% # Separate values for Peak1
  separate(Peak2, into = c("chr2", "start2", "end2"), sep = "_") %>% 
  dplyr::mutate(strand1 = ".") %>% 
  dplyr::mutate(strand2 = ".") %>% 
  dplyr::filter(value > as.numeric(threshold)) %>%  # Keep only the top 20% values
  # Drop the quartile columns (optional)
  dplyr::select(chr1, start1, end1, chr2, start2, end2, generate_name, value, strand1, strand2)

# Write the output
write_tsv(output, output_file, col_names = FALSE)
