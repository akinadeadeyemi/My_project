######################## Correlation coefficient of age and DNA methylation of skin samples from mammals #######################
####. Author:Adeyemi Akinade
####  Place : Clemson University, SC
####  lab   : Gopalan's lab
####  Project:  Comparative Epigenetics of Aging
####  Date   : 22 February, 2025


## loading libraries
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)

# selecting the blood and skin samples
# blood_skin_metadata <- sample_metadata %>% filter(tissue == "Blood"| tissue == "Skin")
## remove the NAs in the metadata
# skin_metadata1 <- skin_metadata %>% filter(!age_years == "NA") 

### select samples that are present in the merged betas
# skin_metadata_betas1 <- merged_betas[, colnames(merged_betas) %in% skin_metadata1$geo_accession]

### select the samples that are present in both beta and metadata
# final_skin_metadata <- skin_metadata1[skin_metadata1$geo_accession %in% colnames(skin_metadata_betas1), ]

# write.csv(final_skin_metadata, "final_skin_metadata.csv", row.names = FALSE)


# filtering by the number of species greater than or equal to 14
# skin_species_count <- final_skin_metadata %>% group_by(species) %>%       
#  tally() %>% filter(n >= 50) %>% ungroup
# View(skin_species_count)
# unique (skin_species_count$species)

# Extract the species names that meet the condition
# filtered_species <- skin_species_count$species

# Subset the metadata to include only those species with blood and skin count  50
# final_skin_samples_50 <- final_skin_metadata %>%
#  filter(species %in% filtered_species)

# write out the final blood samples
# write.csv(final_skin_samples_50, "final_skin_samples_50.csv", row.names = FALSE)


####### merging sample IDs in both beta and blood_skin_samples to generate  ########
# skin_metadata_betas50 <- skin_metadata_betas1[, colnames(skin_metadata_betas1) %in% final_skin_samples_50$geo_accession]
# View(skin_metadata_betas50)

# write.csv(skin_metadata_betas50, "skin_metadata_betas50.csv", row.names = FALSE)

#### Ensure metadata is filtered for matching IDs
# final_skin_samples_50$age_years <- as.numeric(as.character(final_skin_samples_50$age_years))




### Reading and loading files 
cat(" file loading ...")
final_skin_samples_50 <- read_csv("~/Phd_data/final_skin_samples_50.csv")
skin_metadata_betas50 <- read_csv("~/Phd_data/skin_metadata_betas50.csv")


cat("Started working on correlation ...")

correlation_test <- function(data, metadata, n_samples, reps = 100) {
  # Get unique species from the metadata
  unique_species <- unique(metadata$species)
  
  # Initialize a list to store results for each species
  species_results <- list()
  
  for (species in unique_species) {
    cat("Processing species:", species, "\n")
    
    # Subset metadata for the current species
    species_metadata <- metadata[metadata$species == species, ]
    
    # Check if there are enough samples for the species
    n_available <- nrow(species_metadata)
    if (n_available < n_samples) {
      warning(paste("Only", n_available, "samples available for species", species, 
                    "but requested", n_samples, "samples. Using all available samples."))
      n_samples <- n_available
    }
    
    # Initialize matrices for storing correlations and p-values
    n_features <- nrow(data)
    correlation_matrix <- matrix(NA, nrow = n_features, ncol = reps,
                                 dimnames = list(rownames(data), paste0("Rep_", 1:reps)))
    p_value_matrix <- matrix(NA, nrow = n_features, ncol = reps,
                             dimnames = list(rownames(data), paste0("Rep_", 1:reps)))
    conf_int_lower <- matrix(NA, nrow = n_features, ncol = reps,
                             dimnames = list(rownames(data), paste0("Rep_", 1:reps)))
    conf_int_upper <- matrix(NA, nrow = n_features, ncol = reps,
                             dimnames = list(rownames(data), paste0("Rep_", 1:reps)))
    
    for (rep in 1:reps) {
      cat("Running repetition:", rep, "\n")
      
      # Randomly sample n_samples individuals for the current species
      sampled_metadata <- species_metadata[sample(1:nrow(species_metadata), size = n_samples, replace = FALSE), ]
      
      # Extract sample IDs and corresponding age data
      sampled_ids <- sampled_metadata$geo_accession
      age_data <- as.numeric(as.character(sampled_metadata$age_years))  # Ensure numeric
      
      # Subset the data file for the sampled IDs
      matched_indices <- match(sampled_ids, colnames(data))
      if (any(is.na(matched_indices))) {
        warning("Some sampled IDs do not match column names in the data.")
      }
      sampled_data <- data[, matched_indices, drop = FALSE]
      
      # Calculate correlation and p-values for each row (feature)
      for (i in seq_len(n_features)) {
        x <- as.numeric(sampled_data[i, ])
        y <- age_data
        
        # Perform correlation test
        cor_test_result <- cor.test(x, y, use = "complete.obs")
        correlation_matrix[i, rep] <- cor_test_result$estimate
        p_value_matrix[i, rep] <- cor_test_result$p.value
        conf_int_lower[i, rep] <- cor_test_result$conf.int[1]  # Lower bound
        conf_int_upper[i, rep] <- cor_test_result$conf.int[2]  # Upper bound
      }
    }
    
    # Store results for the current species
    species_results[[species]] <- list(
      correlations = correlation_matrix,
      p_values = p_value_matrix,
      conf_int_lower = conf_int_lower,
      conf_int_upper = conf_int_upper
    )
  }
  
  # Return the results as a list
  return(species_results)
}

# Set seed for reproducibility
set.seed(123)

# Run the updated function
skin_sample_50_corr_test <- correlation_test(
  data = skin_metadata_betas50,  # The beta values of the CpG sites
  metadata = final_skin_samples_50,  # The metadata data frame
  n_samples = 50,  # Number of individuals to sample
  reps = 100  # Number of bootstrap repetitions
)

cat("Saving the result to the folder")

saveRDS(skin_sample_50_corr_test, file = "~/Phd_data/outputs/skin_sample_50_correlation.rds")

cat("Results have been saved to the current working directory.\n")

